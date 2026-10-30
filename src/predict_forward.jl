# The forward descriptions for `predict` (KIM-31): one stochastic hybrid cell per
# seed, or one mean-field solve, each read for the predicted quantity Q.
#
# - `:hybrid` runs the reference composition, the handshake driver, one cell per
#   seed. Each seed fixes a cell: the global generator is reseeded before the
#   build and again before the run, as the phase-14 drivers do.
# - `:mean_field` runs route (ii), one ODE per parameter set. It has no cells and
#   no seed, and no single-cell distribution.
#
# Q is read from the run's metered lactate export (`MeteredPtsTransport`), so
# both arms of any comparison are metered (metering moves a trajectory at the
# integrator's tolerance, validation.jl). The mean field's Q is corrected, when
# asked, by δ from a hybrid ensemble at one reference point
# (dev/notes/2026-10-29-cut-predictive-distribution.md §5.4).

"""
    metered(models) -> Vector{AbstractSubModel}

`models` with its `PtsTransport` swapped for a [`MeteredPtsTransport`](@ref), so
that a run records the glucose taken up and the lactate exported. A composition
already metered is returned as it is. One with no PTS module is refused, since
nothing in it exports lactate.
"""
function metered(models::AbstractVector{<:AbstractSubModel})
    any(m -> m isa PtsTransport || m isa MeteredPtsTransport, models) || throw(ArgumentError(
        "This composition has no PtsTransport, so it has no lactate export to meter"))
    return AbstractSubModel[m isa PtsTransport ? MeteredPtsTransport(m) : m for m in models]
end

"""
    CellRun

One forward run, as [`predict`](@ref) keeps it: one hybrid cell, or one
mean-field solve.

Fields:
- `seed`: the cell's seed, or 0 for the mean field, which has none;
- `q`: the predicted quantity, `NaN` when the run failed;
- `failure`: `:none`, or why the run did not finish. It is `:overdraw` when the
  mean field overdrew a pool (its form of starving), `:not_advanced` when the
  hybrid's metabolic block stopped, and `:retcode` when a mean-field solve ended
  without success, such as at its step cap;
- `message`: the failure's text, or empty;
- `drains`: the drains the hybrid made, or 0 for the mean field;
- `clips`: [`clip_summary`](@ref)'s rows, one per counter that clipped. The mean
  field drains unclamped, so it never clips;
- `minima`: each ODE state's smallest particle count over the run, the meters
  aside. It is read at every handshake for the hybrid, and only at the saves for
  the mean field;
- `accounts`: [`carbon_accounts`](@ref) for the run, or `nothing` when it failed
  or is unmetered.
"""
struct CellRun
    seed::Int
    q::Float64
    failure::Symbol
    message::String
    drains::Int
    clips::Vector{NamedTuple}
    minima::Dict{Symbol, Float64}
    accounts::Union{Nothing, NamedTuple}
end

_failed_cell(seed, failure, message) =
    CellRun(seed, NaN, failure, message, 0, NamedTuple[], Dict{Symbol, Float64}(), nothing)

"The quantities `predict` names: `carbon_accounts`' integrated export and exported fraction."
const PREDICT_QUANTITIES = (:lactate_out, :exported)

function _check_quantity(q)
    q isa Symbol && !(q in PREDICT_QUANTITIES) && throw(ArgumentError(
        "quantity :$q is not one of $PREDICT_QUANTITIES; pass a function of the run instead"))
    return q
end

_quantity(q::Symbol, run) = Float64(getproperty(carbon_accounts(run), q))
_quantity(f, run) = Float64(f(run))

_quantity_label(q::Symbol) = q
_quantity_label(f) = Symbol(string(f))

const _METERS = (GLC_UPTAKE_METER, LAC_EXPORT_METER)

function _cell(seed, run::ValidationRun, quantity, drains, clips)
    minima = Dict{Symbol, Float64}(r.species => r.min_particles
                                   for r in particle_floor(run) if !(r.species in _METERS))
    accounts = LAC_EXPORT_METER in run.ode_names ? carbon_accounts(run) : nothing
    return CellRun(seed, _quantity(quantity, run), :none, "", drains,
                   collect(NamedTuple, clips), minima, accounts)
end

# What a forward description needs, built once per prediction.
struct _ForwardSetup
    forward::Symbol
    models::Vector{AbstractSubModel}
    horizon::Float64
    route::Union{Nothing, MeanFieldRoute}
    saveat::Vector{Float64}
end

function _forward_setup(models, forward::Symbol; horizon::Real, every::Real = 60.0,
                        complete::Bool = true)
    forward in (:hybrid, :mean_field) || throw(ArgumentError(
        "forward must be :hybrid or :mean_field, not :$forward"))
    horizon > 0 || throw(ArgumentError("horizon must be positive, got $horizon"))
    ms = collect(AbstractSubModel, models)
    h = Float64(horizon)
    forward === :hybrid && return _ForwardSetup(forward, ms, h, nothing, Float64[])
    rt = build_mean_field_problem(ms; tspan = (0.0, h), complete)
    saveat = sort(unique(vcat(collect(0.0:Float64(every):h), h)))
    return _ForwardSetup(forward, ms, h, rt, saveat)
end

# One hybrid cell, fixed by its seed. A metabolic block that stops is a failed
# cell, kept with its reason; anything else is an error in the call.
function _hybrid_cell(f::_ForwardSetup, values, seed::Integer, quantity)
    ms = f.models
    try
        Random.seed!(seed)
        d = build_problem(ms; tspan = (0.0, f.horizon), complete = true)
        set_parameters!(d, ms, values)
        Random.seed!(seed)
        n = round(Int, f.horizon / d.interval)
        n * d.interval ≈ f.horizon || throw(ArgumentError(
            "horizon $(f.horizon) is not a whole number of $(d.interval) s handshakes"))
        rec = ClipRecord(d)
        run = validation_run!(ms, d, n; moieties = Moiety[], clips = rec)
        return _cell(Int(seed), run, quantity, rec.drains, clip_summary(rec))
    catch e
        e isa MetabolicBlockStalled || rethrow()
        return _failed_cell(Int(seed), :not_advanced, sprint(showerror, e))
    end
end

# Route (ii)'s solution as a `ValidationRun` at its saves, so that Q and the
# minima are read the way the hybrid's are. The membrane count is not a species.
function _route_run(rt::MeanFieldRoute, sol)
    keep = [i for i in eachindex(rt.state_names) if i != rt.n_idx]
    factor = [_route_geometry(rt, rt.n_idx == 0 ? 0.0 : u[rt.n_idx]).factor for u in sol.u]
    ode = [collect(Float64, u[keep]) for u in sol.u]
    samples = NamedTuple{(:t, :u, :p), Tuple{Float64, Vector{Float64}, Vector{Float64}}}[]
    return ValidationRun(rt.state_names[keep], Symbol[], Moiety[], collect(Float64, sol.t),
                         ode, factor, Int[], Dict{Symbol, Vector{Float64}}(), 0.0, samples,
                         maximum(factor), 1e-10, 1e-8, :mean_field, 0)
end

function _mean_field_cell(f::_ForwardSetup, values, quantity)
    rt = f.route
    p = mean_field_parameters(rt, values)
    sol = try
        solve_mean_field(rt; p, saveat = f.saveat)
    catch e
        e isa PoolOverdraw || rethrow()
        return _failed_cell(0, :overdraw, sprint(showerror, e))
    end
    successful_retcode(sol) || return _failed_cell(0, :retcode,
                                                   "the solve ended with retcode $(sol.retcode)")
    return _cell(0, _route_run(rt, sol), quantity, 0, NamedTuple[])
end

# One forward run per seed for the hybrid, or one for the mean field.
function _run_cells(f::_ForwardSetup, values, seeds, quantity; threads::Bool = false)
    f.forward === :hybrid || return [_mean_field_cell(f, values, quantity)]
    out = Vector{CellRun}(undef, length(seeds))
    _each_cause(length(seeds), threads) do i
        out[i] = _hybrid_cell(f, values, seeds[i], quantity)
    end
    return out
end

# Every run of a matrix of parameter sets at once, the pairs with `run[m, k]`
# false left empty: one job per hybrid cell, or one per mean-field solve, so the
# threads share the work however M, K and n_cells fall. Each hybrid cell reseeds
# the generator of the task it runs on (task-local) before its build and its
# run, and the mean field is deterministic, so a threaded run is the serial one.
function _run_matrix(f::_ForwardSetup, values::AbstractMatrix, seeds::AbstractArray{Int, 3},
                     quantity, run::AbstractMatrix{Bool}; threads::Bool)
    M, K = size(values)
    out = [CellRun[] for _ in 1:M, _ in 1:K]
    jobs = Tuple{Int, Int, Int}[]
    for m in 1:M, k in 1:K
        run[m, k] || continue
        n = f.forward === :hybrid ? size(seeds, 1) : 1
        out[m, k] = Vector{CellRun}(undef, n)
        append!(jobs, ((m, k, c) for c in 1:n))
    end
    _each_cause(length(jobs), threads) do j
        m, k, c = jobs[j]
        out[m, k][c] = f.forward === :hybrid ?
                       _hybrid_cell(f, values[m, k], seeds[c, k, m], quantity) :
                       _mean_field_cell(f, values[m, k], quantity)
    end
    return out
end

# `_each`, with an error inside a task rethrown as itself rather than wrapped in
# a `CompositeException` of `TaskFailedException`s, so a refusal keeps its type
# whether or not the work ran threaded.
function _each_cause(f, n::Integer, threads::Bool)
    try
        _each(f, n, threads)
    catch e
        c = _task_cause(e)
        c === e ? rethrow() : throw(c)
    end
    return nothing
end
_task_cause(e::CompositeException) = _task_cause(first(e.exceptions))
_task_cause(e::TaskFailedException) = _task_cause(e.task.exception)
_task_cause(e) = e

# The cells reseed the global generator; a caller's stream is put back after.
function _restoring_global_rng(f)
    saved = copy(Random.default_rng())
    try
        return f()
    finally
        copy!(Random.default_rng(), saved)
    end
end

"""
    forward_cells(models, values, seeds; forward = :hybrid, horizon,
                  quantity = :lactate_out, threads) -> Vector{CellRun}

Run `models` at the parameter writes `values` (`name => value` pairs, written as
given). Unlike [`intervention_cells`](@ref) and [`delta_correction`](@ref), this
derives no reverse constants. Pass `values` through [`haldane_values`](@ref)
first to keep every equilibrium constant. The composition is [`metered`](@ref)
first.

- `forward = :hybrid` runs one cell per seed.
- `forward = :mean_field` runs route (ii) once, and `seeds` is ignored.

`quantity` is `:lactate_out` (the lactate exported over the run, in particles),
`:exported` (the fraction of the lactate formed that left the cell), or a function
of the run's [`ValidationRun`](@ref). The global generator is restored afterwards.
With `threads`, the hybrid's cells run threaded, with the same result.
"""
function forward_cells(models::AbstractVector{<:AbstractSubModel}, values,
                       seeds::AbstractVector{<:Integer}; forward::Symbol = :hybrid,
                       horizon::Real, quantity = :lactate_out,
                       threads::Bool = Threads.nthreads() > 1)
    _check_quantity(quantity)
    f = _forward_setup(metered(models), forward; horizon)
    return _restoring_global_rng(() -> _run_cells(f, values, seeds, quantity; threads))
end

"""
    intervention_cells(models, values, a, seeds; forward = :hybrid, horizon,
                       quantity = :lactate_out, threads) -> (; baseline, intervened, delta)

Two arms on the same seeds: `values`, and `values` under intervention `a`
([`intervene`](@ref)). Each hybrid cell is paired with itself. Both arms derive
their reverse constants from what they run ([`haldane_values`](@ref)), which
[`forward_cells`](@ref) leaves to its caller. `delta` is the paired effect Q_intervened −
Q_baseline for each cell, `NaN` where either arm failed. `threads` is as for
`forward_cells`.
"""
function intervention_cells(models::AbstractVector{<:AbstractSubModel}, values,
                            a::Intervention, seeds::AbstractVector{<:Integer};
                            forward::Symbol = :hybrid, horizon::Real,
                            quantity = :lactate_out,
                            threads::Bool = Threads.nthreads() > 1)
    _check_quantity(quantity)
    f = _forward_setup(metered(models), forward; horizon)
    vals = haldane_values(f.models, Pair{Symbol, Float64}[n => Float64(v) for (n, v) in values])
    base, int = _restoring_global_rng() do
        _run_cells(f, vals, seeds, quantity; threads),
        _run_cells(f, intervene(a, f.models, vals), seeds, quantity; threads)
    end
    return (; baseline = base, intervened = int, delta = _paired_delta(base, int))
end

_paired_delta(base, int) = [b.failure === :none && i.failure === :none ? i.q - b.q : NaN
                            for (b, i) in zip(base, int)]

"""
    DeltaCorrection

The mean field's correction for Q, computed once at a reference point and kept so
it can be reused ([`delta_correction`](@ref)): δ = μ̂_hybrid(θ_ref) − m(θ_ref),
added to the mean field's Q at every draw. It assumes δ is constant over the
posterior's support, which the crossed design tests on M0
(dev/notes/2026-10-29-cut-predictive-distribution.md §5.4).

Fields:
- `values`: θ_ref's parameter writes;
- `horizon` and `quantity`: what Q is, `quantity` as its label;
- `quantity_def`: Q itself, as it was passed (a symbol, or the function of the
  run), which a prediction reusing this correction must match by `isequal`.
  Labels alone do not tell two closures apart;
- `modules`: the composition's module ids;
- `fingerprint`: a content hash of the composition, over each module's type and
  every field, so two compositions of the same modules with different settings
  differ. A prediction reusing this correction must match it, and the ids;
- `seeds` and `hybrid`: the ensemble's cells;
- `mean_field`: the mean field's run at θ_ref;
- `hybrid_mean`, `hybrid_se` and `n_ok`: the mean of Q over the cells that
  finished, its standard error, and how many finished;
- `delta`: `hybrid_mean − mean_field.q`;
- `intervention`: the intervention the correction also covers, or `nothing`.
  The intervened arm runs on the same seeds, and `hybrid_intervened`,
  `mean_field_intervened` and `delta_intervened` are its counterparts, so the
  corrected mean-field effect is `m(g_a θ) − m(θ) + delta_intervened − delta`.
"""
struct DeltaCorrection
    values::Vector{Pair{Symbol, Float64}}
    horizon::Float64
    quantity::Symbol
    quantity_def::Any
    modules::Vector{Symbol}
    fingerprint::UInt
    seeds::Vector{Int}
    hybrid::Vector{CellRun}
    mean_field::CellRun
    hybrid_mean::Float64
    hybrid_se::Float64
    n_ok::Int
    delta::Float64
    intervention::Union{Nothing, Intervention}
    hybrid_intervened::Vector{CellRun}
    mean_field_intervened::Union{Nothing, CellRun}
    delta_intervened::Float64
end

"""
    delta_correction(models, values, seeds; horizon, quantity = :lactate_out,
                     intervention = nothing, threads) -> DeltaCorrection

The hybrid ensemble at the reference writes `values`, one cell per seed, against
the mean field there. Reverse constants are derived from `values`
([`haldane_values`](@ref)). With an `intervention`, the intervened arm runs too,
on the same seeds ([`intervene`](@ref)). A cell that fails is kept but left out of
the mean. If fewer than two cells finish, or the mean field fails at θ_ref, δ is
undefined and the call throws. With `threads`, the ensemble's cells run threaded,
with the same result.
"""
function delta_correction(models::AbstractVector{<:AbstractSubModel}, values,
                          seeds::AbstractVector{<:Integer}; horizon::Real,
                          quantity = :lactate_out,
                          intervention::Union{Nothing, Intervention} = nothing,
                          threads::Bool = Threads.nthreads() > 1)
    _check_quantity(quantity)
    length(seeds) >= 2 && allunique(seeds) || throw(ArgumentError(
        "the ensemble needs at least two distinct seeds"))
    ms = metered(models)
    vals = haldane_values(ms, Pair{Symbol, Float64}[n => Float64(v) for (n, v) in values])
    hy = _forward_setup(ms, :hybrid; horizon)
    mf = _forward_setup(ms, :mean_field; horizon)
    arm(v) = (_run_cells(hy, v, seeds, quantity; threads),
              only(_run_cells(mf, v, Int[], quantity)))
    (cells, m), int = _restoring_global_rng() do
        arm(vals), intervention === nothing ? nothing : arm(intervene(intervention, ms, vals))
    end
    hm, se, n, d = _delta_arm(cells, m)
    cells_i, m_i, d_i = int === nothing ? (CellRun[], nothing, NaN) :
                        (int[1], int[2], _delta_arm(int[1], int[2])[4])
    return DeltaCorrection(vals, Float64(horizon), _quantity_label(quantity), quantity,
                           [module_id(x) for x in ms], _fingerprint(ms), collect(Int, seeds),
                           cells, m, hm, se, n, d, intervention, cells_i, m_i, d_i)
end

# A content hash of a composition: each module's type and every field,
# recursively. Two compositions built alike agree in any session, and two of the
# same modules with different settings differ. Symbols and types are hashed by
# name, since their own hashes are by address.
_fingerprint(models) = foldl((h, m) -> _fingerprint(m, h), models; init = hash(length(models)))
_fingerprint(x::Union{Bool, Integer, AbstractFloat, AbstractString, Char, Nothing, Missing},
             h::UInt) = hash(x, hash(string(typeof(x)), h))
_fingerprint(x::Symbol, h::UInt) = hash(String(x), hash("s", h))
_fingerprint(x::Union{Type, Module}, h::UInt) = hash(string(x), hash("t", h))
_fingerprint(x::Union{AbstractArray, Tuple}, h::UInt) =
    foldl((g, y) -> _fingerprint(y, g), x; init = hash(size(x), hash(string(typeof(x)), h)))
_fingerprint(x::AbstractDict, h::UInt) =
    foldl((g, kv) -> _fingerprint(kv[2], _fingerprint(kv[1], g)),
          sort!(collect(x); by = kv -> string(kv[1])); init = hash(length(x), hash("d", h)))
function _fingerprint(x, h::UInt)
    h = hash(string(typeof(x)), h)
    for f in fieldnames(typeof(x))
        h = isdefined(x, f) ? _fingerprint(getfield(x, f), hash(String(f), h)) : hash("#undef", h)
    end
    return h
end

function _delta_arm(cells, m)
    m.failure === :none || error("The mean field fails at θ_ref ($(m.failure)), so δ is undefined")
    # The cells that count are those `predict`'s μ̂ counts: a finite Q.
    hm, v, n = _moments([c.q for c in cells])
    n >= 2 || error("Fewer than two hybrid cells finished at θ_ref, so δ is undefined")
    return hm, sqrt(v / n), n, hm - m.q
end

export metered, CellRun, PREDICT_QUANTITIES, forward_cells, intervention_cells,
       DeltaCorrection, delta_correction
