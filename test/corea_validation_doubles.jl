using InferCell
using StaticArrays: SVector

import InferCell: states, parameters, dynamics, contributions, contributed_states,
                  coupling, inputs, module_id, membrane_protein_states,
                  extracellular_states, reduction_notes, formalism, inference_mode

# Executable doubles for spec §11 phase 14a: the assembled Core A′ with one
# thing broken, so each balance check can be shown to fail naming its quantity.
#
# They live here rather than in the test file because two callers read them:
# `test/test_corea_validation.jl` asserts on them, and
# `dev/scripts/corea_validation.jl` records the full-scale numbers. Written
# twice, the recorded mutation and the asserted one could silently differ.
#
# Each mutation targets one balance check. Where the chemistry of the mutation
# also breaks phosphate, the table says so:
#
# | Mutation | Check | Moiety it must break |
# |---|---|---|
# | GAPD draws two NAD⁺ | 3 | redox |
# | LDH makes two lactate | 2 | carbon |
# | recycling mutants (phase 8's `MutatedRecycling`) | 4, 4b | adenylate, guanylate, phosphate |
# | GLCpts1 creates phospho-HPr | 5 | the ptsH carrier (and phosphate, which it also creates) |
# | GAPD's counters charge ten A for ten U | 4, across the boundary | adenylate |
# | translation's debits unclamped | 1 | GTP goes negative |

"""
    glycolysis_mutation(id, species, coefficient) -> CentralGlycolysis

The assembled model's glycolysis (`enzymes = :translated`) with one reaction's
coefficient on one species changed. Phase 6's `mutate_stoichiometry`, applied
to the translated mode.
"""
function glycolysis_mutation(id::Symbol, species::Symbol, coefficient::Int)
    bump(ps) = [s === species ? (s => coefficient) : (s => n) for (s, n) in ps]
    rs = [r.id === id ?
          merge(r, (substrates = bump(r.substrates), products = bump(r.products))) : r
          for r in GLYCOLYTIC_REACTIONS]
    return CentralGlycolysis(enzymes = :translated, reactions = rs)
end

"""
    MeteredCarrierLeak(inner = MeteredPtsTransport())

The metered PTS module with GLCpts1 written as if it *created* phospho-HPr
rather than transferring the phosphate from phospho-EI. This is phase 7's
`PtsCarrierLeak`, on the module the assembled checks compose. The ptsH pair is
no longer conserved, and the other three carriers are.
"""
struct MeteredCarrierLeak{P} <: InferCell.AbstractSubModel
    inner::P
end
MeteredCarrierLeak() = MeteredCarrierLeak(MeteredPtsTransport())
for f in (:states, :parameters, :inputs, :contributed_states, :coupling,
          :membrane_protein_states, :extracellular_states, :reduction_notes,
          :formalism, :inference_mode, :module_id)
    @eval $f(m::MeteredCarrierLeak) = $f(m.inner)
end
contributions(u, p, t, m::MeteredCarrierLeak, ui) = contributions(u, p, t, m.inner, ui)
function dynamics(u, p, t, m::MeteredCarrierLeak, ui)
    du = dynamics(u, p, t, m.inner, ui)
    own = InferCell._pts_own(u)
    _, v1, _, _, _, _ = InferCell._pts_rates(InferCell._live_scalars(m.inner.inner, p), own, ui)
    return Base.setindex(du, du[5] + v1, 5)   # M_ptsh_P_c created, not transferred
end

"""
    shifted_genes(locus; A = +1, U = -1) -> Vector{TranscriptionGene}

The transcription extract with one gene's base counts moved by `A` and `U`,
length unchanged. Given to `CoreATranscription` alone, its counters charge `A`
extra ATP and as many fewer UTP per transcript than the sequence the closure
reads from the extract. The phosphate total is unchanged, since both are one
phosphate per base, so only adenylate should fail. The mutation is in the
boundary accounting, not in any ODE stoichiometry.
"""
function shifted_genes(locus::Symbol; A::Int = 1, U::Int = -1)
    A + U == 0 || throw(ArgumentError("A + U must be 0 so the length is unchanged"))
    return [g.locus === locus ?
            InferCell.TranscriptionGene(g.locus, g.reaction, g.length,
                                        (A = g.counts.A + A, C = g.counts.C,
                                         G = g.counts.G, U = g.counts.U + U),
                                        g.first_two, g.ptn_count, g.mean_mrna) : g
            for g in read_transcription_genes()]
end

"""
    validation_models(; glycolysis, pts, recycling, charging, transcription, translation)

The assembled composition with the carbon meters, and any one module swapped.
Every default is exactly `corea_models(metered = true)`.
"""
function validation_models(; glycolysis = CentralGlycolysis(enzymes = :translated),
                             pts = MeteredPtsTransport(),
                             recycling = NucleotideRecycling(enzymes = :translated),
                             charging = TrnaCharging(),
                             transcription = CoreATranscription(),
                             translation = CoreATranslation())
    return InferCell.AbstractSubModel[glycolysis, pts, recycling, charging,
                                      transcription, CoreATranscriptDecay(), translation]
end

"The kinase-removed recycling module of check 4's second configuration."
kinase_removed() = NucleotideRecycling(enzymes = :translated,
                                       reactions = Tuple(r for r in RECYCLING_REACTIONS
                                                         if r !== :R_ADK1))
