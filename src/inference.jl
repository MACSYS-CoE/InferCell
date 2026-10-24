@model function _infercell_model(data_obs, data_times, prob, priors, n_ode,
                                  solver, sensealg, rows)
    n = length(priors)
    theta = Vector{Real}(undef, n)
    for i in 1:n
        theta[i] ~ priors[i]
    end

    # Comprehension preserves Dual element type for AD compatibility
    p_ode = [theta[i] for i in 1:n_ode]

    sol = solve(remake(prob, p=p_ode), solver;
                saveat=data_times, sensealg=sensealg)

    if sol.retcode !== ReturnCode.Success
        Turing.@addlogprob! -Inf
        return
    end

    sigma = theta[n_ode + 1]
    for i in eachindex(data_times)
        data_obs[:, i] ~ MvNormal(sol.u[i][rows], sigma)
    end
end

# The model under a per-modality noise model (spec §11 task 15.4): the model's
# free parameters, then one log-scale σ per :lognormal modality. The models' own
# :observation parameters belong to the single-σ form and are not sampled.
@model function _infercell_model_noise(data, prob, priors, n_ode, solver, sensealg,
                                       rows, noise, idx)
    n = length(priors)
    theta = Vector{Real}(undef, n)
    for i in 1:n
        theta[i] ~ priors[i]
    end
    p_ode = [theta[i] for i in 1:n_ode]
    sol = solve(remake(prob, p=p_ode), solver; saveat=data.times, sensealg=sensealg)
    if sol.retcode !== ReturnCode.Success
        Turing.@addlogprob! -Inf
        return
    end
    scales = [theta[i] for i in (n_ode + 1):n]
    Turing.@addlogprob! observation_loglik(noise, data, i -> sol.u[i][rows], scales, idx)
end

"""
    build_turing_model(models, data, prob; solver, sensealg, priors_override=nothing,
                       noise=nothing)

Construct the Turing.jl probabilistic model used by [`infer`](@ref) for the
differentiable (NUTS) path. `priors_override::Dict{Symbol,Distribution}`, when
supplied, replaces the model-declared prior for each named parameter — used by
the boundary protocol to feed conditioned priors back into NUTS.

Returns `(turing_model, param_names::Vector{Symbol})` where `param_names`
orders the sampled parameters (model free params first, then observation
free params).

`noise`, a [`NoiseModel`](@ref), replaces the single additive σ with one noise
model per modality (spec §11 task 15.4). The observation parameters are then
the `:lognormal` modalities' scales, `sigma_<name>`, and the models' own
`:observation` parameters are not sampled.
"""
function build_turing_model(models::Vector{<:AbstractSubModel}, data::ObservedData, prob;
                             solver=Tsit5(), sensealg=ForwardDiffSensitivity(),
                             priors_override::Union{Nothing, Dict{Symbol, <:Distribution}}=nothing,
                             noise::Union{Nothing, NoiseModel}=nothing)
    all_model_free = unique_params(reduce(vcat, model_free_params.(parameters.(models))))
    if noise !== nothing
        model_priors = [priors_override !== nothing && haskey(priors_override, p.name) ?
                        priors_override[p.name] : p.prior for p in all_model_free]
        scales = noise_scales(noise)
        rows = state_rows(data.species, reduce(vcat, states.(models)))
        idx = _modality_rows(noise, data)   # checked once, at build time
        turing_model = _infercell_model_noise(
            data, prob, vcat(model_priors, last.(scales)), length(all_model_free),
            solver, sensealg, rows, noise, idx)
        return turing_model, vcat([p.name for p in all_model_free], first.(scales))
    end
    all_obs_free = unique_params(reduce(vcat, obs_free_params.(parameters.(models))))
    all_free = vcat(all_model_free, all_obs_free)

    priors = [
        priors_override !== nothing && haskey(priors_override, p.name) ?
            priors_override[p.name] : p.prior
        for p in all_free
    ]
    param_names = [p.name for p in all_free]
    n_ode = length(all_model_free)

    rows = state_rows(data.species, reduce(vcat, states.(models)))
    turing_model = _infercell_model(
        data.observations, data.times, prob,
        priors, n_ode, solver, sensealg, rows
    )

    return turing_model, param_names
end

function _rename_chain(chain, param_names::Vector{Symbol})
    mapping = Dict("theta[$i]" => string(name) for (i, name) in enumerate(param_names))
    return replacenames(chain, mapping...)
end

"""
    infer(models, data; sampler=NUTS(), n_samples=1000, kwargs...)
    infer(model,  data; kwargs...)

Dispatch entry point for parameter inference. Selects NUTS for sub-models with
`inference_mode = :differentiable` and ABC-SMC for `:simulation`. Returns an
`MCMCChains.Chains` object for NUTS or an [`ABCPosterior`](@ref) for ABC-SMC.

The mode is the **composition's**, not its first module's (spec §11 task 13.3):
every module must declare the same one, and a hybrid ODE/jump composition is
refused, because neither backend can run the handshake driver it builds.

`priors_override::Dict{Symbol,Distribution}` is forwarded to
[`build_turing_model`](@ref) on the NUTS path and is the hook used by the
boundary protocol's iterative loop.

The route is [`inference_route`](@ref)'s, so `route`, if given, must be the one this
data takes. A hybrid composition's routes take a cell and bulk record instead, by
the method `infer(models, ds; route)`.
"""
function infer(models::Vector{<:AbstractSubModel}, data::ObservedData;
               sampler=NUTS(), n_samples=1000,
               solver=Tsit5(), sensealg=ForwardDiffSensitivity(),
               prob=nothing, tspan=nothing,
               priors_override::Union{Nothing, Dict{Symbol, <:Distribution}}=nothing,
               noise::Union{Nothing, NoiseModel}=nothing,
               route::Union{Nothing, Symbol}=nothing,
               kwargs...)
    r = inference_route(models, data; route)
    if r === :nuts
        return _infer_nuts(models, data; sampler=sampler, n_samples=n_samples,
                           solver=solver, sensealg=sensealg, prob=prob, tspan=tspan,
                           priors_override=priors_override, noise=noise)
    elseif r === :abc
        return _infer_abc(models, data; n_particles=n_samples, tspan=tspan, kwargs...)
    else
        error("infer has no runner for route :$r with ObservedData")
    end
end

infer(model::AbstractSubModel, data::ObservedData; kwargs...) =
    infer([model], data; kwargs...)

# Which backend a composition takes, read from all of it. Reading the first
# module alone let composition order choose the path: a hybrid listing an ODE
# module first went to NUTS and one listing a jump module first went to ABC,
# and neither can run the `HandshakeDriver` a hybrid builds.
function _composition_inference_mode(models::Vector{<:AbstractSubModel})
    isempty(models) && throw(ArgumentError("infer needs at least one sub-model"))
    if _determine_formalism(models) === :mixed
        ids(f) = join((string(module_id(m)) for m in models if formalism(m) === f), ", ")
        throw(ArgumentError(
            "infer has no path for a hybrid composition: $(ids(:ode)) are :ode " *
            "and $(ids(:jump)) are :jump. With ObservedData it builds a " *
            "HandshakeDriver, which NUTS (an ODEProblem through Turing) and ABC-SMC " *
            "(a JumpProblem re-solved with SSAStepper) both cannot run. A hybrid's " *
            "routes take per-cell transcripts and a bulk record, as m0_dataset " *
            "returns: infer(models, ds; route = :exact) or route = :mixed"))
    end
    modes = unique(inference_mode.(models))
    length(modes) == 1 || throw(ArgumentError(
        "The sub-models in this composition declare different inference modes " *
        "— $(join(("$(module_id(m)) :$(inference_mode(m))" for m in models), ", ")) " *
        "— so no one backend is right for all of them. Declare one mode"))
    return only(modes)
end

# The fields a cell and bulk record needs, as `m0_dataset` returns, and the two
# routes that take one (spec/phases/approximate-method.md §12 2026-10-24, "B6 as
# planned").
const CELL_DATA_FIELDS = (:times, :genes, :seeds, :replicate_seed, :transcripts, :bulk)
const CELL_ROUTES = (:exact, :mixed)

"""
    inference_route(models, data; route = nothing) -> Symbol

The route [`infer`](@ref) takes for a composition and its data, one of the three of
spec/phases/approximate-method.md §1, or ABC:
- `ObservedData` and an all-ODE composition: `:nuts`; an all-jump one: `:abc`. A
  hybrid is refused.
- A cell and bulk record (a `NamedTuple` as [`m0_dataset`](@ref) returns) and a
  hybrid composition with M0's gene expression, one `CoreATranscription` and one
  `CoreATranscriptDecay`: `route`, which must be named, `:exact` (the three-block Gibbs
  sampler, the reference) or `:mixed` (the mixed likelihood, an approximation).

Anything else throws an `ArgumentError` that names what is wrong.
"""
function inference_route(models::Vector{<:AbstractSubModel}, data::ObservedData;
                         route::Union{Nothing, Symbol} = nothing)
    mode = _composition_inference_mode(models)
    r = mode === :differentiable ? :nuts : mode === :simulation ? :abc :
        error("Unknown inference mode: $mode")
    route === nothing || route === r || throw(ArgumentError(
        "route :$route does not apply here: this composition with ObservedData takes :$r"))
    return r
end

function inference_route(models::Vector{<:AbstractSubModel}, ds::NamedTuple;
                         route::Union{Nothing, Symbol} = nothing)
    isempty(models) && throw(ArgumentError("infer needs at least one sub-model"))
    lacks = [f for f in CELL_DATA_FIELDS if !haskey(ds, f)]
    isempty(lacks) || throw(ArgumentError(
        "a cell and bulk record needs the fields $(join(CELL_DATA_FIELDS, ", ")), as " *
        "m0_dataset returns; this one lacks $(join(lacks, ", "))"))
    f = _determine_formalism(models)
    f === :mixed || throw(ArgumentError(
        "cell and bulk data have routes only for a hybrid composition, jump gene " *
        "expression coupled to ODE metabolism; this one is all :$f. Its route takes " *
        "ObservedData"))
    nt = count(m -> m isa CoreATranscription, models)
    nd = count(m -> m isa CoreATranscriptDecay, models)
    nt == 1 && nd == 1 || throw(ArgumentError(
        "the cell-data routes are built for M0's gene expression: one CoreATranscription " *
        "and one CoreATranscriptDecay, as m0_models() holds. This hybrid has $nt and $nd"))
    route === nothing && throw(ArgumentError(
        "a hybrid composition with cell and bulk data has two routes, and one must be " *
        "named: route = :exact (the three-block Gibbs sampler, the reference) or " *
        "route = :mixed (the mixed likelihood, an approximation)"))
    route in CELL_ROUTES || throw(ArgumentError(
        "unknown route :$route for cell and bulk data; the routes are " *
        join((":$r" for r in CELL_ROUTES), " and ")))
    return route
end

"""
    infer(models, ds; route, rng = Random.default_rng(), kwargs...) -> Chains

Inference for a hybrid composition from a cell and bulk record `ds`, as
[`m0_dataset`](@ref) returns, by the named `route` ([`inference_route`](@ref)). It
runs the route end to end in one process:
- `:exact`: [`gibbs_setup`](@ref), a start drawn from the prior as V7's,
  [`gibbs_init`](@ref) and `n_sweeps` of [`gibbs_sweep!`](@ref), the first `n_burn`
  dropped. Keywords: `n_sweeps`, `n_burn = 0`, `N = 20`. Each cell's driver is
  built by `gibbs_setup`'s default, `build_m0`, so `models` must be M0's
  composition: a different one gives drivers and priors from two compositions.
- `:mixed`: B5's pipeline. [`reference_theta`](@ref), the reference ensemble at it
  ([`reference_ensemble`](@ref), [`reference_seeds`](@ref)), [`delta_and_S`](@ref),
  then NUTS on [`mixed_model`](@ref) from θ_ref. Keywords: `n_adapt`, `n_draws`,
  `n_ref = 800`, `target_accept = 0.8`, `max_depth = 10`, `map_maxiters = 200`.

Each call runs one chain, on its own stream forked from `rng`, and leaves the global
generator as the fork left it, so consecutive calls differ; join several with
`chainscat`. Both routes return `MCMCChains.Chains` over
[`mixed_parameter_names`](@ref)' coordinates, `log_<target>` for each of
`M0_TARGETS` and then `log_sigma_b`, with the route's own state in the chain's
`info` (`start`, or `theta_ref`, `delta` and `S`) and `info.route` set.
"""
function infer(models::Vector{<:AbstractSubModel}, ds::NamedTuple;
               route::Union{Nothing, Symbol} = nothing,
               rng::AbstractRNG = Random.default_rng(), kwargs...)
    r = inference_route(models, ds; route)
    # The route draws from its own stream, forked here. Its setup reseeds the global
    # generator once per cell (`gibbs_setup`, `emit_ensemble`), so the global state is
    # put back as it was after the fork: drawing after the setup would give every call
    # the same draws, and leaving it reseeded would make the next call repeat this one.
    stream = Xoshiro(rand(rng, UInt64))
    saved = copy(Random.default_rng())
    try
        return r === :exact ? _infer_exact(models, ds; rng = stream, kwargs...) :
               _infer_mixed(models, ds; rng = stream, kwargs...)
    finally
        copy!(Random.default_rng(), saved)
    end
end

function _infer_exact(models, ds; rng::AbstractRNG,
                      n_sweeps::Integer, n_burn::Integer = 0, N::Integer = 20)
    0 <= n_burn < n_sweeps || throw(ArgumentError(
        "n_burn must be in 0:$(n_sweeps - 1), got $n_burn"))
    g = gibbs_setup(ds; models, N, horizon = ds.times[end])
    order = vcat(g.cme, g.forwards)
    sort(order) == sort(M0_TARGETS) || throw(ArgumentError(
        "the Gibbs sampler's free parameters $(order) are not M0_TARGETS"))
    perm = [findfirst(==(n), order) for n in M0_TARGETS]
    v = Dict(draw_truth(rng, models; names = M0_TARGETS, purpose = :calibration).values)
    start = (cme = [v[n] for n in g.cme], u = [log(v[f]) for f in g.forwards],
             σ = rand(rng, g.sigma_prior))
    s = gibbs_init(rng, g; start...)
    for _ in 1:n_sweeps
        gibbs_sweep!(s, g)
    end
    draws = reduce(vcat, (permutedims(vcat(vcat(log.(t.cme), t.u)[perm], log(t.σ)))
                          for t in s.trace[(n_burn + 1):end]))
    return Chains(reshape(draws, size(draws)..., 1), _log_names(M0_TARGETS);
                  info = (route = :exact, start = start))
end

function _infer_mixed(models, ds; rng::AbstractRNG,
                      n_adapt::Integer, n_draws::Integer, n_ref::Integer = 800,
                      target_accept::Real = 0.8, max_depth::Integer = 10,
                      map_maxiters::Integer = 200)
    s0 = mixed_setup(ds; models)
    ref = reference_theta(s0; maxiters = map_maxiters)
    r = _mixed_solve(s0, exp.(ref.x[1:length(s0.names)]))
    r === nothing && error("route (ii) cannot run at θ_ref ($(ref.method)), so δ and S cannot be formed")
    cells = reference_ensemble(mixed_theta_values(s0, ref.x), reference_seeds(ds, n_ref);
                               models, horizon = ds.times[end])
    d = delta_and_S(cells, mean_field_log_bulk(s0, r.sol); C = size(ds.transcripts, 1))
    s = with_reference(s0, d.delta, d.S)
    chain = sample(rng, mixed_model(s), NUTS(n_adapt, target_accept; max_depth), n_draws;
                   initial_params = ref.x, progress = false)
    chain = replacenames(chain, ("x[$i]" => string(n)
                                 for (i, n) in enumerate(mixed_parameter_names(s)))...)
    return setinfo(chain, merge(chain.info, (route = :mixed, theta_ref = ref, delta = d.delta,
                                             S = d.S)))
end

function _infer_nuts(models, data; sampler=NUTS(), n_samples=1000,
                     solver=Tsit5(), sensealg=ForwardDiffSensitivity(),
                     prob=nothing, tspan=nothing,
                     priors_override::Union{Nothing, Dict{Symbol, <:Distribution}}=nothing,
                     noise::Union{Nothing, NoiseModel}=nothing)
    if prob === nothing
        t = tspan === nothing ? (data.times[1], data.times[end]) : tspan
        prob = build_problem(models; tspan=t)
    end

    turing_model, param_names = build_turing_model(
        models, data, prob; solver=solver, sensealg=sensealg,
        priors_override=priors_override, noise=noise)

    chain = sample(turing_model, sampler, n_samples)

    return _rename_chain(chain, param_names)
end

function _infer_abc(models, data;
                    n_particles=1000, n_populations=10, alpha=0.5,
                    n_replicates=100, tspan=nothing,
                    verbose=false, rng=Random.default_rng(),
                    priors=nothing, param_names=nothing)
    t = tspan === nothing ? (data.times[1], data.times[end]) : tspan

    if priors === nothing
        all_free = unique_params(model_free_params(reduce(vcat, parameters.(models))))
        priors = [p.prior for p in all_free]
        param_names = [p.name for p in all_free]
    end
    names = reduce(vcat, states.(models))
    state_rows(data.species, names)   # throws on a species that is not a state

    # The spread is kept whenever the data carries it (spec §11 task 15.9).
    # Data observed from one noisy ODE solution has none, and then the summary
    # is the per-time mean alone.
    spread = data.spread !== nothing
    observed_stats = spread ?
        vec(vcat(data.observations, data.spread)) : vec(data.observations)
    base_prob = build_problem(models; tspan=t)

    function simulate(theta)
        prob = remake(base_prob, p=theta)
        trajectories = [solve(prob, SSAStepper(); saveat=data.times) for _ in 1:n_replicates]
        return compute_summary_stats(trajectories, data.species; times=data.times,
                                     spread=spread, state_names=names)
    end

    return abc_smc(simulate, observed_stats, priors, param_names;
                   n_particles=n_particles, n_populations=n_populations,
                   alpha=alpha, verbose=verbose, rng=rng)
end

"""
    observe(sol, times, model; sigma=0.1, rng=...) -> ObservedData
    observe(sol, times, models; sigma=0.1, rng=...)
    observe(trajectories, times, model)
    observe(trajectories, times, models)

Synthesise an [`ObservedData`](@ref) record from a forward simulation. The two
`sol`-flavoured methods sample a single ODE solution at `times` and add
i.i.d. Gaussian noise with standard deviation `sigma`. The two
`trajectories`-flavoured methods take the mean across SSA replicates, and keep
the replicates' standard deviation as the record's `spread`, so the ABC summary
carries both (spec §11 task 15.9).
"""
function observe(sol, times, model::AbstractSubModel;
                 sigma=0.1, rng=Random.default_rng())
    pred = Array(sol(times))
    noise = sigma .* randn(rng, size(pred))
    return ObservedData(collect(Float64, times), pred .+ noise, states(model))
end

function observe(sol, times, models::Vector{<:AbstractSubModel};
                 sigma=0.1, rng=Random.default_rng())
    pred = Array(sol(times))
    noise = sigma .* randn(rng, size(pred))
    species = reduce(vcat, states.(models))
    return ObservedData(collect(Float64, times), pred .+ noise, species)
end

"""
    observe(sol, times, models, noise::NoiseModel; scales, rng=...) -> ObservedData

Synthesise data under a per-modality noise model (spec §11 task 15.4): each
`:lognormal` species is `x·exp(σ·ε)` with its modality's σ from `scales`, a
`Dict` keyed by [`scale_name`](@ref), and each `:poisson` species is a Poisson
count with mean `x`. Only the noise model's species are observed.
"""
function observe(sol, times, models::Vector{<:AbstractSubModel}, noise::NoiseModel;
                 scales::AbstractDict, rng=Random.default_rng())
    names = reduce(vcat, states.(models))
    species = reduce(vcat, [m.species for m in noise.modalities])
    rows = state_rows(species, names)
    y = zeros(length(species), length(times))
    for (j, t) in enumerate(times)
        x = sol(t)[rows]
        k = 0
        for m in noise.modalities, s in m.species
            k += 1
            if m.kind === :lognormal
                y[k, j] = max(x[k], m.floor) * exp(scales[scale_name(m)] * randn(rng))
            else
                y[k, j] = rand(rng, Poisson(max(x[k], m.floor)))
            end
        end
    end
    return ObservedData(times, y, species)
end

observe(trajectories::Vector, times, model::AbstractSubModel) =
    observe(trajectories, times, [model])

function observe(trajectories::Vector, times, models::Vector{<:AbstractSubModel})
    species = reduce(vcat, states.(models))
    stats = reshape(compute_summary_stats(trajectories, species; times = times),
                    2 * length(species), length(times))
    n = length(species)
    # One trajectory has no spread to keep; a zero would pull ABC toward
    # deterministic parameter regions.
    spread = length(trajectories) > 1 ? stats[(n + 1):end, :] : nothing
    return ObservedData(times, stats[1:n, :], species, spread)
end

"""
    posterior_predictive(model_or_models, chain_or_posterior; n_samples=100, tspan, saveat=nothing)

Draw forward simulations from the posterior. Accepts either an `MCMCChains`
chain (NUTS path) or an [`ABCPosterior`](@ref) (ABC-SMC path; particles are
resampled by their weights). `tspan` is required; `saveat` defaults to 100
points across `tspan`.

Returns a [`PosteriorPredictive`](@ref).
"""
function posterior_predictive(model::AbstractSubModel, chain;
                               n_samples=100, solver=Tsit5(),
                               tspan=nothing, saveat=nothing)
    tspan === nothing && error("tspan must be provided for posterior_predictive")

    param_names = [string(p.name) for p in model_free_params(parameters(model))]
    n_chain = size(chain, 1)
    n_draw = min(n_samples, n_chain)
    idxs = rand(1:n_chain, n_draw)

    extract_params(idx) = [chain[name].data[idx] for name in param_names]

    return _run_posterior_predictive([model], tspan, saveat, idxs, extract_params, solver)
end

function posterior_predictive(models::Vector{<:AbstractSubModel}, chain;
                               n_samples=100, solver=Tsit5(),
                               tspan=nothing, saveat=nothing)
    tspan === nothing && error("tspan must be provided for posterior_predictive")

    all_free = unique_params(reduce(vcat, model_free_params.(parameters.(models))))
    param_names = [string(p.name) for p in all_free]
    n_chain = size(chain, 1)
    n_draw = min(n_samples, n_chain)
    idxs = rand(1:n_chain, n_draw)

    extract_params(idx) = [chain[name].data[idx] for name in param_names]

    return _run_posterior_predictive(models, tspan, saveat, idxs, extract_params, solver)
end

function posterior_predictive(model::AbstractSubModel, result::ABCPosterior;
                               n_samples=100, tspan=nothing, saveat=nothing)
    tspan === nothing && error("tspan must be provided for posterior_predictive")

    n_available = size(result.particles, 2)
    n_draw = min(n_samples, n_available)
    idxs = [_weighted_sample(result.weights, Random.default_rng()) for _ in 1:n_draw]

    stepper = formalism(model) == :jump ? SSAStepper() : Tsit5()
    extract_params(idx) = result.particles[:, idx]

    return _run_posterior_predictive([model], tspan, saveat, idxs, extract_params, stepper)
end

function _run_posterior_predictive(models::Vector{<:AbstractSubModel}, tspan, saveat,
                                   idxs, extract_params, solver)
    prob = build_problem(models; tspan=tspan)
    save_times = saveat === nothing ? range(tspan[1], tspan[2]; length=100) : saveat

    solutions = Vector{Any}(undef, length(idxs))
    for (j, idx) in enumerate(idxs)
        solutions[j] = solve(remake(prob, p=extract_params(idx)), solver; saveat=save_times)
    end

    species = reduce(vcat, states.(models))
    return PosteriorPredictive(solutions, collect(Float64, save_times), species)
end

"""
    check_identifiability(model_or_models, prob, times; solver=Tsit5())

Local structural-identifiability diagnostic: builds the forward-map Jacobian
∂(sol(times)) / ∂p via ForwardDiff and returns its rank alongside the parameter
and observation counts. `full_rank == true` means every parameter is locally
identifiable from the chosen observation pattern; a deficient rank flags
sloppy / unidentifiable directions that need richer data.
"""
function check_identifiability(model::AbstractSubModel, prob, times;
                                solver=Tsit5())
    all_free = model_free_params(parameters(model))
    p0 = [p.value for p in all_free]

    function forward_map(p)
        sol = solve(remake(prob, p=p), solver; saveat=times)
        return vec(Array(sol))
    end

    J = ForwardDiff.jacobian(forward_map, p0)
    r = rank(J)

    return (rank=r, n_params=length(p0), n_obs=length(states(model)) * length(times),
            full_rank=r >= length(p0), jacobian=J)
end

function check_identifiability(models::Vector{<:AbstractSubModel}, prob, times;
                                solver=Tsit5())
    all_free = unique_params(reduce(vcat, model_free_params.(parameters.(models))))
    p0 = [p.value for p in all_free]
    n_states = sum(length(states(m)) for m in models)

    function forward_map(p)
        sol = solve(remake(prob, p=p), solver; saveat=times)
        return vec(Array(sol))
    end

    J = ForwardDiff.jacobian(forward_map, p0)
    r = rank(J)

    return (rank=r, n_params=length(p0), n_obs=n_states * length(times),
            full_rank=r >= length(p0), jacobian=J)
end
