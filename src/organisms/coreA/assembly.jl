"""
The assembled Core A′: all seven modules, every counter wired (spec §11 phase 13b).

Phase 11's full compositions stripped transcription's and decay's counters,
because their CTP and UTP debits could not compose until phase 13a gave a
chemostat an ownerless path. This is the composition with nothing stripped, and
the one place it is written down, so the tests and the Slurm drivers cannot
drift apart on what "the assembled model" means.
"""

"""
    COREA_CYCLE_S

One cell cycle, 6,300 s — the horizon every full-cycle check runs over. Spec §3
forbids shortening it: 144 s of the adenylate pool looked fine before the dead
end was found.
"""
const COREA_CYCLE_S = CHARGING_CYCLE_S

"""
    corea_models() -> Vector{AbstractSubModel}

The seven Core A′ modules as the assembled model composes them. Glycolysis and
recycling read their enzyme concentrations from translated protein counts
(`enzymes = :translated`, spec §11 phase 11a); transcription, decay and
translation carry every deferred counter they declare.
"""
corea_models() = AbstractSubModel[
    CentralGlycolysis(enzymes = :translated),
    PtsTransport(),
    NucleotideRecycling(enzymes = :translated),
    TrnaCharging(),
    CoreATranscription(),
    CoreATranscriptDecay(),
    CoreATranslation(),
]

"""
    build_corea(; tspan = (0.0, COREA_CYCLE_S), kwargs...) -> HandshakeDriver

Build the assembled Core A′ over one full cycle, in completeness mode: the build
fails if any registry dynamic state is unowned or any moiety is stranded (spec
§11 task 13.1). The remaining keywords are the hybrid driver's
([`build_problem`](@ref)); the solver and tolerances default to the pinned
Rodas5P at `abstol = 1e-10`, `reltol = 1e-8` of spec §3.
"""
build_corea(; tspan = (0.0, COREA_CYCLE_S), kwargs...) =
    build_problem(corea_models(); tspan, complete = true, kwargs...)

export COREA_CYCLE_S, corea_models, build_corea
