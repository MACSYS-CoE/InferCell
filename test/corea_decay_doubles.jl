# Doubles for spec §11 phase 12 (transcript decay).
#
# `_jic` comes from `jump_test_models.jl`, which `runtests.jl` includes first.

import InferCell: states, parameters, reactions, formalism, inputs,
                  written_states, coupling, module_id, inference_mode

"""
    TranscriptSource(genes; n0 = 5)

Owns the seventeen transcripts at a fixed initial count and has no reactions.

**Why this exists.** Task 12.4 wants the recyclable monomers to reach the pools
nucleotide recycling closes, and the only proof of that is a hybrid build in
which the hook executes decay's AMP and GMP credits. `CoreATranscription` cannot
be in that build, because its CTP and UTP counters need an owner the registry
forbids (§12, 2026-09-10 E; task 13.9). This stands in for it as the
transcripts' owner, so decay has something to decrement and nothing else moves.
"""
struct TranscriptSource <: AbstractSubModel
    genes::Vector{TranscriptionGene}
    params::Vector{InferParameter}
end

function TranscriptSource(genes = read_transcription_genes(); n0 = 5)
    genes = collect(TranscriptionGene, genes)
    params = InferParameter[_jic(n0, Symbol(transcript_state(g.locus), "0"),
                                 :TranscriptSource) for g in genes]
    return TranscriptSource(genes, params)
end

states(m::TranscriptSource) = [transcript_state(g.locus) for g in m.genes]
parameters(m::TranscriptSource) = m.params
formalism(::TranscriptSource) = :jump
inference_mode(::TranscriptSource) = :simulation
reactions(::TranscriptSource) = Reaction[]
