"""
Observables and synthetic data for Core A′ (spec §11 phase 15).

What the likelihood may condition on, and where its ground truth comes from.
Both rules are in spec §4 D8, and both are enforced here by a check, not left
to the prose:

- **No observable closes the promoter-proxy loop.** Promoter strength is protein
  copy number over 180 (§4 D5). So a protein count in the observable set, or a
  freed protein initial condition, lets the likelihood re-estimate the quantity
  the promoter prior was built from.
- **A coverage or calibration truth is drawn from the prior, never set to the
  published values.** A truth at the nominal values agrees with a prior centred
  on the proxy by construction. A run at nominal truth is a smoke test and
  must be labelled as one.
"""

# ---------------------------------------------------------------------------
# 15.1: the circularity guard
# ---------------------------------------------------------------------------

"""
    protein_count_states(models) -> Vector{Symbol}

Every state in `models` that is a protein count: translation's per-gene protein
states and its cytosolic ptsG, and the eight phosphotransferase carrier states.
The carriers count here too. Their initial conditions are proteomics counts
times proteomics fractions (§4 D7), so observing them closes the same loop as
observing a translated protein.
"""
function protein_count_states(models::AbstractVector{<:AbstractSubModel})
    out = Symbol[]
    for m in models
        if m isa CoreATranslation
            append!(out, [protein_state(g.locus) for g in m.genes])
            push!(out, TL_PTSG_CYTO)
        elseif m isa PtsTransport
            append!(out, species_in_group(:pts))
        end
    end
    return out
end

"""
    assert_no_circularity(models, observables; free = Symbol[]) -> Nothing

Spec §11 task 15.1 and §4 D8. Throws if `observables` includes a protein count
(see [`protein_count_states`](@ref)), or if a protein count's initial condition
is free, either in `models` or in the extra parameter names `free`. The error
names the seventeen promoter parameters the loop would double-count. Metabolite
concentrations, transcript counts and the cell volume pass.
"""
function assert_no_circularity(models::AbstractVector{<:AbstractSubModel},
                               observables::AbstractVector{Symbol};
                               free::AbstractVector{Symbol} = Symbol[])
    proteins = Set(protein_count_states(models))
    ics = Set(Symbol(s, "0") for s in proteins)
    freed = vcat(collect(free),
                 [p.name for m in models for p in parameters(m)
                  if !p.fixed && p.role === :initial_condition])
    bad_obs = [o for o in observables if o in proteins]
    bad_ic = unique([n for n in freed if n in ics])
    isempty(bad_obs) && isempty(bad_ic) && return nothing
    promoters = [promoter_param(g.locus) for g in read_transcription_genes()]
    what = String[]
    isempty(bad_obs) || push!(what, "observes protein counts $(bad_obs)")
    isempty(bad_ic) || push!(what, "frees protein initial conditions $(bad_ic)")
    throw(ArgumentError(
        "The data model $(join(what, " and ")). Promoter strength is protein copy " *
        "number over $PROMOTER_DIVISOR (spec §4 D5), so this double-counts the " *
        "$(length(promoters)) promoter parameters $(promoters) (spec §4 D8, " *
        "task 15.1)"))
end

# ---------------------------------------------------------------------------
# 15.2: the truth-drawing rule
# ---------------------------------------------------------------------------

"The purposes a truth set can serve. Only `:smoke` may sit at the nominal values."
const TRUTH_PURPOSES = (:coverage, :calibration, :recovery, :smoke)

"""
    draw_truth(rng, models; names = D11_TARGETS, purpose = :coverage) -> NamedTuple

Spec §11 task 15.2 and §4 D8: one ground truth for a synthetic dataset, drawn
from the prior by [`draw_parameters`](@ref) under the Haldane-consistent rule of
task 14c.1. Each reversible reaction a draw touches keeps the equilibrium
constant its nominal Mode values imply (§12 2026-09-28 A). Returns
`(values, purpose, label)`, where `values` is ready for
[`set_parameters!`](@ref) and `label` is the T2 row the rule adds.

A truth for `purpose = :smoke` is not drawn. Use [`nominal_truth`](@ref).
"""
function draw_truth(rng::AbstractRNG, models::AbstractVector{<:AbstractSubModel};
                    names::AbstractVector{Symbol} = D11_TARGETS,
                    purpose::Symbol = :coverage)
    purpose in TRUTH_PURPOSES || throw(ArgumentError(
        "purpose must be one of $TRUTH_PURPOSES, not :$purpose"))
    purpose === :smoke && throw(ArgumentError(
        "a smoke test runs at the nominal values: use nominal_truth"))
    values = draw_parameters(rng, models, names; reverse = :derived)
    truth = (values = values, purpose = purpose, label = truth_label())
    check_truth(models, truth; names)
    return truth
end

"""
    nominal_truth(models; names = D11_TARGETS, purpose = :smoke) -> NamedTuple

A truth at the published nominal values. Spec §4 D8 permits this only for a
smoke test, so any other `purpose` throws.
"""
function nominal_truth(models::AbstractVector{<:AbstractSubModel};
                       names::AbstractVector{Symbol} = D11_TARGETS,
                       purpose::Symbol = :smoke)
    purpose === :smoke || throw(ArgumentError(
        "a truth at the nominal published values is permitted only for a smoke " *
        "test (spec §4 D8); a :$purpose run must draw its truth with draw_truth"))
    nominal = nominal_parameter_values(models)
    return (values = [n => nominal[n] for n in names], purpose = :smoke,
            label = truth_label())
end

"""
    check_truth(models, truth; names = D11_TARGETS) -> Nothing

Refuse a truth set that is not a smoke test and that puts any target at its
nominal published value, and any truth that breaks a reaction's equilibrium
constant (via [`assert_haldane`](@ref)). A value drawn from a continuous prior
lands exactly on the nominal with probability zero, so an exact match means the
truth was set, not drawn.
"""
function check_truth(models::AbstractVector{<:AbstractSubModel}, truth::NamedTuple;
                     names::AbstractVector{Symbol} = D11_TARGETS)
    truth.purpose in TRUTH_PURPOSES || throw(ArgumentError(
        "purpose must be one of $TRUTH_PURPOSES, not :$(truth.purpose)"))
    assert_haldane(models, truth.values)
    truth.purpose === :smoke && return nothing
    nominal = nominal_parameter_values(models)
    given = Dict(truth.values)
    at = [n for n in names if haskey(given, n) && given[n] == nominal[n]]
    isempty(at) || throw(ArgumentError(
        "A :$(truth.purpose) run's truth sets $(at) to the nominal published " *
        "values. The truth must be drawn from the prior (spec §4 D8): a truth at " *
        "the proxy agrees with a prior centred on it by construction. Label the " *
        "run a smoke test, or draw the truth with draw_truth"))
    missing = [n for n in names if !haskey(given, n)]
    isempty(missing) || throw(ArgumentError("the truth gives no value for $(missing)"))
    return nothing
end

"""
    truth_label() -> ReductionLabel

Spec §6 T2's row for the truth-drawing rule: truths keep the equilibrium
constants the nominal Mode values imply, rather than upstream's balanced
equilibrium-constant rows, which disagree with the Modes by more than 1% on 9 of
15 reactions (§12 2026-09-28 A). The rule is ours.
"""
truth_label() = ReductionLabel(
    :model_note, :truth_draws,
    "synthetic truths are drawn from the prior with each reverse constant derived " *
    "through the equilibrium constant the nominal Mode values imply, not " *
    "upstream's balanced equilibrium-constant rows, which disagree with the Modes " *
    "by more than 1% on 9 of 15 reactions (spec §12 2026-09-28 A). The rule is ours")

# ---------------------------------------------------------------------------
# 15.3: the candidate observables at the published cadence
# ---------------------------------------------------------------------------

"""
    reaction_fluxes(u, p, t, m, u_inputs) -> Vector{Pair{Symbol, Float64}}

One ODE module's net reaction rates, in mM/s, named by reaction. Together the
four Core A′ ODE modules give all 22: ten glycolytic, five recycling, the five
phosphotransferase steps plus lactate export, and charging. Spec §4 D8 keeps
fluxes as a reported diagnostic and out of the likelihood.
"""
reaction_fluxes(u, p, t, m::CentralGlycolysis, u_inputs) =
    [r.id for r in GLYCOLYTIC_REACTIONS] .=> collect(reaction_rates(u, p, t, m, u_inputs))
reaction_fluxes(u, p, t, m::NucleotideRecycling, u_inputs) =
    collect(RECYCLING_REACTIONS) .=> collect(recycling_fluxes(u, p, t, m, u_inputs))
reaction_fluxes(u, p, t, m::PtsTransport, u_inputs) =
    PTS_FLUX_NAMES .=> collect(_pts_rates(_live_scalars(m, p), u, u_inputs))
reaction_fluxes(u, p, t, m::TrnaCharging, u_inputs) =
    [:R_trna_chg => charging_flux(u, p, t, m, u_inputs)]

"The phosphotransferase module's six rates, in `_pts_rates` order."
const PTS_FLUX_NAMES = [:R_GLCpts0, :R_GLCpts1, :R_GLCpts2, :R_GLCpts3, :R_GLCpts4,
                        :R_LACt]

# Each ODE module's slices of the block, as `_build_ode_problem` wires them: the
# contexts carry their input maps only once the coupling is resolved.
function _ode_contexts(ode_models)
    contexts = _build_contexts(ode_models)
    _resolve_coupling(ode_models, contexts; validate = false)
    return contexts
end

# Every ODE module's fluxes at the driver's live state. The ODE block is
# `_build_ode_problem(ode_models)`, so `_build_contexts(ode_models)` gives each
# module's slices of it.
function _driver_fluxes(d::HandshakeDriver, ode_models, contexts)
    u, p, t = d.ode.u, d.ode.p, d.ode.t
    out = Pair{Symbol, Float64}[]
    for (m, ctx) in zip(ode_models, contexts)
        ui = u[Int[ctx.input_map[s] for s in inputs(m)]]
        append!(out, reaction_fluxes(u[ctx.state_idxs], p[ctx.param_idxs], t, m, ui))
    end
    return out
end

"""
    emit_observables!(d, models; every = 60, horizon = d.t_end) -> NamedTuple

Spec §11 task 15.3: run the driver to `horizon`, and at t = 0 and every
`every` handshakes (the published 60 s cadence at the 1 s handshake) record the
candidate observables:

- `counts`, species × times: every ODE state in particles, with its carried
  remainder, so both rounding policies report the amount they hold, and every
  jump species as its count. The deferred cost counters are bookkeeping, not
  species, and are left out.
- `fluxes`, reactions × times, in mM/s, from [`reaction_fluxes`](@ref).
- `volume`, in litres.

Returns `(times, species, counts, reactions, fluxes, volume)`.
"""
function emit_observables!(d::HandshakeDriver, models::AbstractVector{<:AbstractSubModel};
                           every::Integer = 60, horizon::Real = d.t_end)
    every >= 1 || throw(ArgumentError("every must be at least one handshake"))
    ode_models = AbstractSubModel[m for m in models if formalism(m) === :ode]
    jump_models = AbstractSubModel[m for m in models if formalism(m) === :jump]
    contexts = _ode_contexts(ode_models)
    ode_names = _block_names(models, :ode)
    jump_names = _block_names(models, :jump)
    counters = Set(e.counter for m in models for e in coupling(m)
                   if e isa DeferredCounterEdge)
    keep = [j for (j, s) in enumerate(jump_names) if !(s in counters)]
    species = vcat(ode_names, jump_names[keep])
    n_steps = round(Int, (horizon - d.ode.t) / d.interval)
    n_steps >= 0 || throw(ArgumentError("horizon $horizon is before t = $(d.ode.t)"))

    times = Float64[]
    counts = Vector{Vector{Float64}}()
    fluxes = Vector{Vector{Float64}}()
    volume = Float64[]
    reactions = Symbol[]
    function record!()
        push!(times, d.ode.t)
        push!(counts, vcat(collect(d.ode.u) .* d.factor .+ d.rounding.remainders,
                           Float64.(d.jump.u[keep])))
        f = _driver_fluxes(d, ode_models, contexts)
        isempty(reactions) && append!(reactions, first.(f))
        push!(fluxes, last.(f))
        push!(volume, d.volume_litres)
    end
    record!()
    for k in 1:n_steps
        handshake_step!(d)
        (k % every == 0 || k == n_steps) && record!()
    end
    return (times = times, species = species, counts = reduce(hcat, counts),
            reactions = reactions, fluxes = reduce(hcat, fluxes), volume = volume)
end

"""
    emit_ensemble(build, seeds; every = 60, horizon = nothing) -> NamedTuple

Spec §11 task 15.3 over replicate cells. `build()` returns `(models, driver)`
for one fresh cell. Each is run by [`emit_observables!`](@ref) under
`Random.seed!(seed)`, seeded before the build and again before the run, as the
phase-14 drivers do, so a seed fixes a cell.

Returns `(times, seeds, species, counts, reactions, fluxes, volume, genes,
transcripts)`. `counts` is replicates × species × times and `fluxes` replicates
× reactions × times. `transcripts` is an `Int` array of replicates × genes ×
times: per-cell counts, not a replicate mean, which is what a count likelihood
needs (§4 D8, §10 R14).
"""
function emit_ensemble(build, seeds::AbstractVector{<:Integer}; every::Integer = 60,
                       horizon::Union{Nothing, Real} = nothing)
    runs = map(seeds) do s
        Random.seed!(s)
        models, d = build()
        Random.seed!(s)
        emit_observables!(d, models; every, horizon = horizon === nothing ? d.t_end : horizon)
    end
    r1 = first(runs)
    all(r -> r.species == r1.species && r.times == r1.times && r.reactions == r1.reactions,
        runs) || error("replicates disagree on their species, reactions or save times")
    stack3(f) = permutedims(cat((f(r) for r in runs)...; dims = 3), (3, 1, 2))
    counts = stack3(r -> r.counts)
    genes = [g.locus for g in read_transcription_genes()]
    rows = state_rows([transcript_state(g) for g in genes], r1.species)
    tx = counts[:, rows, :]
    all(isinteger, tx) || error("a transcript count is not an integer")
    return (times = r1.times, seeds = collect(seeds), species = r1.species,
            counts = counts, reactions = r1.reactions, fluxes = stack3(r -> r.fluxes),
            volume = permutedims(reduce(hcat, (r.volume for r in runs))),
            genes = genes, transcripts = round.(Int, tx))
end

export protein_count_states, assert_no_circularity
export TRUTH_PURPOSES, draw_truth, nominal_truth, check_truth, truth_label
export reaction_fluxes, PTS_FLUX_NAMES, emit_observables!, emit_ensemble
