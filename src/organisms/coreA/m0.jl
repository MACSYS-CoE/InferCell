"""
M0, the two-gene reference system (spec/phases/16-recovery.md D16.4, task 16a.8).

Core A′ cut to two genes over 600 s. Transcription, decay and translation run
for GAPD (`JCVISYN3A_0607`, the best-determined promoter) and ptsG
(`JCVISYN3A_0779`, the headline, Core A′'s only membrane protein). Everything
else in the model is unchanged.
- The other twelve enzyme counts the metabolic modules read are held at their
  proteomics copy numbers by [`HeldEnzymeCounts`](@ref), so the catalytic channel
  runs exactly as it does in Core A′, with one live count.
- The other three phosphotransferase carriers are neither translated nor
  credited, so their totals are held.

Both are labelled M0 reductions in `reduction_report`.
"""

"M0's two genes: GAPD and ptsG."
const M0_GENES = [:JCVISYN3A_0607, :JCVISYN3A_0779]

"M0's horizon, 600 s (spec §4 D10)."
const M0_HORIZON_S = 600.0

"""
    M0_TARGETS

M0's free parameters apart from σ, which is drawn separately: D11's six less
`S_PGI`, whose gene M0 cuts (spec/phases/16-recovery.md D16.4).
"""
const M0_TARGETS = [promoter_param(:JCVISYN3A_0607), promoter_param(:JCVISYN3A_0779),
                    :krnadeg, :kcatF_R_ENO, :kcatF_R_FBA]

"""
    HeldEnzymeCounts(loci)

A jump module with no reactions, owning `P_<locus>` for each locus at its
proteomics copy number. It stands in for translation of the genes M0 cuts, so the
metabolic modules' catalytic edges read a constant count. The hold is an M0
reduction, and `reduction_notes` says so.
"""
struct HeldEnzymeCounts <: AbstractSubModel
    genes::Vector{TranscriptionGene}
    params::Vector{InferParameter}
end

function HeldEnzymeCounts(loci::AbstractVector{Symbol})
    all_genes = Dict(g.locus => g for g in read_transcription_genes())
    missing_loci = [l for l in loci if !haskey(all_genes, l)]
    isempty(missing_loci) || throw(ArgumentError("unknown loci $missing_loci"))
    genes = [all_genes[l] for l in loci]
    params = [InferParameter(float(g.ptn_count), Poisson(g.ptn_count), true,
                             Symbol(protein_state(g.locus), "0"), :HeldEnzymeCounts,
                             :initial_condition) for g in genes]
    return HeldEnzymeCounts(genes, params)
end

states(m::HeldEnzymeCounts) = [protein_state(g.locus) for g in m.genes]
parameters(m::HeldEnzymeCounts) = m.params
formalism(::HeldEnzymeCounts) = :jump
inference_mode(::HeldEnzymeCounts) = :simulation
reactions(::HeldEnzymeCounts) = Reaction[]
reduction_notes(m::HeldEnzymeCounts) = [
    "M0 holds $(length(m.genes)) enzyme counts at their proteomics copy numbers " *
    "($(join(string.(g.locus for g in m.genes), ", "))): their genes are cut, so " *
    "nothing translates them. **Ours, M0 only** (spec/phases/16-recovery.md D16.4)",
    "M0 translates no phosphotransferase carrier but ptsG, so the ptsI, ptsH and crr " *
    "totals are held at their initial values. **Ours, M0 only**"]

# The enzyme counts the metabolic modules read through catalytic edges.
_catalytic_species(models) = unique(e.species for m in models for e in coupling(m)
                                    if e isa CatalyticEdge)

"""
    m0_models(; smoothing = nothing) -> Vector{AbstractSubModel}

M0 as composed (spec/phases/16-recovery.md D16.4). It is [`d11_models`](@ref)'s
composition with the three stochastic modules cut to [`M0_GENES`](@ref), and a
[`HeldEnzymeCounts`](@ref) for every catalytic count those genes do not make.
The ENO and FBA forward and reverse constants are free, as in `d11_models`.
"""
function m0_models(; smoothing::Union{Nothing, Real} = nothing)
    drain = smoothing === nothing ? (clip = :clamped_deficit_carried, smoothing = nothing) :
                                    (clip = :smoothed, smoothing = Float64(smoothing))
    all_genes = read_transcription_genes()
    genes = [g for g in all_genes if g.locus in M0_GENES]
    length(genes) == length(M0_GENES) || error("M0's genes are not all in the extract")
    ode = AbstractSubModel[
        CentralGlycolysis(enzymes = :translated,
                          free = [:kcatF_R_ENO, :kcatR_R_ENO, :kcatF_R_FBA, :kcatR_R_FBA]),
        PtsTransport(),
        NucleotideRecycling(enzymes = :translated),
        TrnaCharging(),
    ]
    made = Set(protein_state(l) for l in M0_GENES)
    held = [Symbol(String(s)[3:end]) for s in _catalytic_species(ode) if !(s in made)]
    return AbstractSubModel[
        ode...,
        CoreATranscription(; genes, drain...),
        CoreATranscriptDecay(; genes, drain...),
        CoreATranslation(; genes, drain...),
        HeldEnzymeCounts(held),
    ]
end

"""
    build_m0(; tspan = (0.0, M0_HORIZON_S), smoothing = nothing, kwargs...) -> HandshakeDriver

Build M0 in completeness mode, as [`build_corea`](@ref) builds Core A′.
"""
build_m0(; tspan = (0.0, M0_HORIZON_S), smoothing::Union{Nothing, Real} = nothing,
         kwargs...) =
    build_problem(m0_models(; smoothing); tspan, complete = true, kwargs...)

"""
    M0_PANEL

The metabolites an M0 dataset observes: the 15.7 panel, chosen by 15.6's
influence audit (spec §12 2026-09-29), so M0 and Core A′ are observed alike.
"""
const M0_PANEL = [:M_g6p_c, :M_f6p_c, :M_lac__L_c, :M_g3p_c, :M_dhap_c, :M_trna_chg_c,
                  :M_pi_c, :M_3pg_c, :M_2pg_c, :M_gtp_c, :M_fdp_c, :M_atp_c,
                  :M_lac__L_e, :M_amp_c, :M_gmp_c, :M_gdp_c, :M_adp_c]

"The metabolite noise scale's prior, `Modality`'s default (spec/phases/16-recovery.md §3)."
const SIGMA_MET_PRIOR = LogNormal(log(0.2), 1.0)

"""
    m0_dataset(replicate_seed; ncells = 50, horizon = M0_HORIZON_S) -> NamedTuple

One M0 replicate, everything drawn from one recorded seed
(spec/phases/16-recovery.md task 16a.8, V7):
- a truth over [`M0_TARGETS`](@ref), drawn from the prior Haldane-consistently;
- σ from [`SIGMA_MET_PRIOR`](@ref);
- `ncells` cell seeds and a noise seed.

The cells are the published model, clamped under fractional carry. The
[`M0_PANEL`](@ref) metabolites are observed lognormally at σ, floored at one
particle, and M0's two transcripts exactly. The t = 0 save, the fixed initial
condition, is dropped. Returns [`generate_dataset`](@ref)'s record without t = 0,
with `sigma`, `replicate_seed` and `noise` added.
"""
function m0_dataset(replicate_seed::Integer; ncells::Integer = 50,
                    horizon::Real = M0_HORIZON_S)
    rng = Xoshiro(replicate_seed)
    models = m0_models()
    truth = draw_truth(rng, models; names = M0_TARGETS, purpose = :calibration)
    σ = rand(rng, SIGMA_MET_PRIOR)
    cells = Int.(rand(rng, UInt32, ncells))
    allunique(cells) || error("replicate $replicate_seed drew a repeated cell seed")
    noise_seed = Int(rand(rng, UInt32))
    noise = NoiseModel(Modality(:metabolite, :lognormal, M0_PANEL; floor = 1.0,
                                prior = SIGMA_MET_PRIOR),
                       Modality(:transcript, :poisson,
                                [transcript_state(g) for g in M0_GENES]))
    ds = generate_dataset(truth, cells; models, horizon, noise, noise_seed,
                          scales = Dict(:sigma_metabolite => σ), names = M0_TARGETS)
    keep = findall(>(0), ds.times)
    return (; ds..., times = ds.times[keep], latent = ds.latent[:, :, keep],
            transcripts = ds.transcripts[:, :, keep], observed = ds.observed[:, :, keep],
            sigma = σ, replicate_seed = Int(replicate_seed), noise = noise)
end

export M0_GENES, M0_HORIZON_S, M0_TARGETS, M0_PANEL, SIGMA_MET_PRIOR, HeldEnzymeCounts,
       m0_models, build_m0, m0_dataset
