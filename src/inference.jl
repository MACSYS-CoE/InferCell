@model function _infercell_model(data_obs, data_times, prob, priors, n_ode,
                                  solver, sensealg)
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
        data_obs[:, i] ~ MvNormal(sol[:, i], sigma)
    end
end

function build_turing_model(models::Vector{<:AbstractSubModel}, data::ObservedData, prob;
                             solver=Tsit5(), sensealg=ForwardDiffSensitivity(),
                             priors_override::Union{Nothing, Dict{Symbol, <:Distribution}}=nothing)
    all_model_free = unique_params(reduce(vcat, model_free_params.(parameters.(models))))
    all_obs_free = unique_params(reduce(vcat, obs_free_params.(parameters.(models))))
    all_free = vcat(all_model_free, all_obs_free)

    priors = [
        priors_override !== nothing && haskey(priors_override, p.name) ?
            priors_override[p.name] : p.prior
        for p in all_free
    ]
    param_names = [p.name for p in all_free]
    n_ode = length(all_model_free)

    turing_model = _infercell_model(
        data.observations, data.times, prob,
        priors, n_ode, solver, sensealg
    )

    return turing_model, param_names
end

function _rename_chain(chain, param_names::Vector{Symbol})
    mapping = Dict("theta[$i]" => string(name) for (i, name) in enumerate(param_names))
    return replacenames(chain, mapping...)
end

function infer(models::Vector{<:AbstractSubModel}, data::ObservedData;
               sampler=NUTS(), n_samples=1000,
               solver=Tsit5(), sensealg=ForwardDiffSensitivity(),
               prob=nothing, tspan=nothing,
               priors_override::Union{Nothing, Dict{Symbol, <:Distribution}}=nothing,
               kwargs...)
    mode = inference_mode(models[1])
    if mode == :differentiable
        return _infer_nuts(models, data; sampler=sampler, n_samples=n_samples,
                           solver=solver, sensealg=sensealg, prob=prob, tspan=tspan,
                           priors_override=priors_override)
    elseif mode == :simulation
        return _infer_abc(models, data; n_particles=n_samples, tspan=tspan, kwargs...)
    else
        error("Unknown inference mode: $mode")
    end
end

infer(model::AbstractSubModel, data::ObservedData; kwargs...) =
    infer([model], data; kwargs...)

function _infer_nuts(models, data; sampler=NUTS(), n_samples=1000,
                     solver=Tsit5(), sensealg=ForwardDiffSensitivity(),
                     prob=nothing, tspan=nothing,
                     priors_override::Union{Nothing, Dict{Symbol, <:Distribution}}=nothing)
    if prob === nothing
        t = tspan === nothing ? (data.times[1], data.times[end]) : tspan
        prob = build_problem(models; tspan=t)
    end

    turing_model, param_names = build_turing_model(
        models, data, prob; solver=solver, sensealg=sensealg,
        priors_override=priors_override)

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
    species = reduce(vcat, states.(models))

    observed_stats = vec(data.observations)
    base_prob = build_problem(models; tspan=t)

    function simulate(theta)
        prob = remake(base_prob, p=theta)
        trajectories = [solve(prob, SSAStepper(); saveat=data.times) for _ in 1:n_replicates]
        return compute_summary_stats(trajectories, species; times=data.times)
    end

    return abc_smc(simulate, observed_stats, priors, param_names;
                   n_particles=n_particles, n_populations=n_populations,
                   alpha=alpha, verbose=verbose, rng=rng)
end

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

function observe(trajectories::Vector, times, model::AbstractSubModel)
    # For stochastic models: take the mean across replicate trajectories
    species = states(model)
    n_species = length(species)
    n_times = length(times)
    obs = zeros(n_species, n_times)
    for sol in trajectories
        for (j, t) in enumerate(times)
            obs[:, j] .+= sol(t)
        end
    end
    obs ./= length(trajectories)
    return ObservedData(collect(Float64, times), obs, species)
end

function observe(trajectories::Vector, times, models::Vector{<:AbstractSubModel})
    species = reduce(vcat, states.(models))
    n_species = length(species)
    n_times = length(times)
    obs = zeros(n_species, n_times)
    for sol in trajectories
        for (j, t) in enumerate(times)
            obs[:, j] .+= sol(t)
        end
    end
    obs ./= length(trajectories)
    return ObservedData(collect(Float64, times), obs, species)
end

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
