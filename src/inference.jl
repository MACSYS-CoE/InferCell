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
                             solver=Tsit5(), sensealg=ForwardDiffSensitivity())
    all_ode_free = InferParameter[]
    all_obs_free = InferParameter[]
    for m in models
        append!(all_ode_free, ode_free_params(parameters(m)))
        append!(all_obs_free, obs_free_params(parameters(m)))
    end
    all_free = vcat(all_ode_free, all_obs_free)

    priors = [p.prior for p in all_free]
    param_names = [p.name for p in all_free]
    n_ode = length(all_ode_free)

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
               prob=nothing, tspan=nothing)
    if prob === nothing
        t = tspan === nothing ? (data.times[1], data.times[end]) : tspan
        prob = build_problem(models; tspan=t)
    end

    turing_model, param_names = build_turing_model(
        models, data, prob; solver=solver, sensealg=sensealg)

    chain = sample(turing_model, sampler, n_samples)

    return _rename_chain(chain, param_names)
end

infer(model::AbstractSubModel, data::ObservedData; kwargs...) =
    infer([model], data; kwargs...)

function observe(sol, times, model::AbstractSubModel;
                 sigma=0.1, rng=Random.default_rng())
    pred = Array(sol(times))
    noise = sigma .* randn(rng, size(pred))
    return ObservedData(collect(Float64, times), pred .+ noise, states(model))
end

function posterior_predictive(model::AbstractSubModel, chain;
                               n_samples=100, solver=Tsit5(),
                               tspan=nothing, saveat=nothing)
    tspan === nothing && error("tspan must be provided for posterior_predictive")

    prob = build_problem([model]; tspan=tspan)
    param_names = [string(p.name) for p in ode_free_params(parameters(model))]
    save_times = saveat === nothing ? range(tspan[1], tspan[2]; length=100) : saveat

    n_chain = size(chain, 1)
    n_draw = min(n_samples, n_chain)
    idxs = rand(1:n_chain, n_draw)

    solutions = Vector{Any}(undef, n_draw)
    for (j, idx) in enumerate(idxs)
        p_draw = [chain[name].data[idx] for name in param_names]
        solutions[j] = solve(remake(prob, p=p_draw), solver; saveat=save_times)
    end

    return PosteriorPredictive(solutions, collect(Float64, save_times), states(model))
end

function check_identifiability(model::AbstractSubModel, prob, times;
                                solver=Tsit5())
    all_free = ode_free_params(parameters(model))
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
