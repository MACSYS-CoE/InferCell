function build_problem(models::Vector{<:AbstractSubModel}; tspan=(0.0, 100.0))
    for m in models
        validate_formalism(formalism(m))
    end

    collective = _determine_formalism(models)

    if collective == :ode
        return _build_ode_problem(models; tspan=tspan)
    elseif collective == :jump
        return _build_jump_problem(models; tspan=tspan)
    else
        error("Mixed formalism composition not yet supported")
    end
end

function _determine_formalism(models::Vector{<:AbstractSubModel})
    formalisms = unique(formalism.(models))
    length(formalisms) == 1 && return formalisms[1]
    return :mixed
end

function _build_ode_problem(models::Vector{<:AbstractSubModel}; tspan=(0.0, 100.0))
    contexts = _build_contexts(models)
    _resolve_coupling(models, contexts)

    u0 = _build_u0(models)
    p0 = _build_p0(models)
    rhs = _build_rhs(models, contexts)

    return ODEProblem{false}(rhs, u0, tspan, p0)
end

function _build_jump_problem(models::Vector{<:AbstractSubModel}; tspan=(0.0, 100.0))
    contexts = _build_contexts(models)
    _resolve_coupling(models, contexts)

    u0 = _build_u0_integer(models)
    p0 = _build_p0(models)

    dprob = DiscreteProblem(u0, tspan, p0)
    jumps = JumpSet(; constant_jumps=reduce(vcat, reactions.(models)))

    return JumpProblem(dprob, Direct(), jumps)
end

build_problem(model::AbstractSubModel; kwargs...) = build_problem([model]; kwargs...)

function _build_contexts(models::Vector{<:AbstractSubModel})
    contexts = SubModelContext[]
    state_offset = 0
    param_offset = 0

    for m in models
        n_states = length(states(m))
        n_params = length(ode_free_params(parameters(m)))

        state_idxs = (state_offset + 1):(state_offset + n_states)
        param_idxs = (param_offset + 1):(param_offset + n_params)

        push!(contexts, SubModelContext(state_idxs, param_idxs, Dict{Symbol, Int}()))
        state_offset += n_states
        param_offset += n_params
    end

    return contexts
end

function _resolve_coupling(models::Vector{<:AbstractSubModel},
                           contexts::Vector{SubModelContext})
    state_owners = Dict{Symbol, Int}()
    for (i, m) in enumerate(models)
        for (j, s) in enumerate(states(m))
            global_idx = contexts[i].state_idxs[j]
            haskey(state_owners, s) && error(
                "State :$s is owned by multiple sub-models")
            state_owners[s] = global_idx
        end
    end

    for (i, m) in enumerate(models)
        for inp in inputs(m)
            haskey(state_owners, inp) || error(
                "Input :$inp declared by $(typeof(m)) is not owned by any sub-model")
            contexts[i].input_map[inp] = state_owners[inp]
        end
    end
end

function _build_u0(models::Vector{<:AbstractSubModel})
    u0_vals = Float64[]
    for m in models
        ics = ic_params(parameters(m))
        for s in states(m)
            ic = findfirst(p -> p.name == Symbol(s, "0") || p.name == s, ics)
            push!(u0_vals, ic !== nothing ? ics[ic].value : 0.0)
        end
    end
    return SVector{length(u0_vals)}(u0_vals)
end

function _build_u0_integer(models::Vector{<:AbstractSubModel})
    u0_vals = Int[]
    for m in models
        ics = ic_params(parameters(m))
        for s in states(m)
            ic = findfirst(p -> p.name == Symbol(s, "0") || p.name == s, ics)
            push!(u0_vals, ic !== nothing ? round(Int, ics[ic].value) : 0)
        end
    end
    return u0_vals
end

function _build_p0(models::Vector{<:AbstractSubModel})
    reduce(vcat, [p.value for p in ode_free_params(parameters(m))] for m in models;
           init=Float64[])
end

function _build_rhs(models::Vector{<:AbstractSubModel},
                    contexts::Vector{SubModelContext})
    function rhs(u, p, t)
        du_parts = map(models, contexts) do m, ctx
            u_local = u[ctx.state_idxs]
            p_local = p[ctx.param_idxs]
            dynamics(u_local, p_local, t, m)
        end
        return vcat(du_parts...)
    end
    return rhs
end
