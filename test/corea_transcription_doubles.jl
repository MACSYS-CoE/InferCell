# Doubles for spec §11 phase 10 (transcription).
#
# `_rate` and `_jic` come from the phase-1 and phase-2 doubles files, which
# `runtests.jl` includes first; `caught` comes from `corea_test_models.jl`.

import InferCell: states, parameters, reactions, formalism, inputs,
                  written_states, coupling, module_id, inference_mode,
                  reduction_notes, rebuilt_params, rate_constants

"""
    ToyTranscriptDecay(genes; scale = 1.0)

First-order decay of every transcript `CoreATranscription` owns, at the
published per-gene constant `(18/452)·88 / n_g`.

**Why this exists, and why it is not scaffolding.** Phase 10's calibration
check compares each gene's time-averaged transcript count against its measured
mean. Transcription alone has no consumer, so the counts only accumulate and
every gene's average is an order of magnitude high — the check would be
vacuous rather than merely weak. This is the same shape as phase 8's charging
drain double and phase 9's translation-demand double: the acceptance criterion
needs a consumer, and the consumer is a double until the phase that owns it
lands. Recorded as a §12 amendment dated 2026-09-10.

**Phase 12 supersedes it** for composed runs and keeps it only for standalone
ones.

It also exercises the amendment of 2026-09-04 from the consumer's side: a
transcript is not a registry species, so no edge can name it, and the write is
gated by `written_states` alone. `scale` is overridable so a test can drive
decay faster than the published rate without rebuilding the gene list.
"""
struct ToyTranscriptDecay <: AbstractSubModel
    genes::Vector{TranscriptionGene}
    params::Vector{InferParameter}
    scale::Float64
    written::Vector{Symbol}
end

function ToyTranscriptDecay(genes = read_transcription_genes(); scale = 1.0,
                            written = nothing)
    genes = collect(TranscriptionGene, genes)
    params = InferParameter[_jic(0, :mRNA_decayed0, :ToyTranscriptDecay)]
    for g in genes
        k = scale * transcript_decay_constant(g)
        push!(params, _rate(k, Symbol("k_deg_", g.locus), :ToyTranscriptDecay))
    end
    names = [transcript_state(g.locus) for g in genes]
    return ToyTranscriptDecay(genes, params, Float64(scale),
                              collect(Symbol, written === nothing ? names : written))
end

states(::ToyTranscriptDecay) = [:mRNA_decayed]
parameters(m::ToyTranscriptDecay) = m.params
formalism(::ToyTranscriptDecay) = :jump
inference_mode(::ToyTranscriptDecay) = :simulation
inputs(m::ToyTranscriptDecay) = [transcript_state(g.locus) for g in m.genes]
written_states(m::ToyTranscriptDecay) = m.written

function reactions(m::ToyTranscriptDecay)
    return [Reaction(
                # `_jic` builds a *fixed* parameter, which does not reach the
                # composed parameter vector, so the i-th free parameter is the
                # i-th decay constant.
                (u, p, t, w) -> p[i] * w[i],
                (u, w) -> (w[i] -= 1; u[1] += 1))
            for i in eachindex(m.genes)]
end
