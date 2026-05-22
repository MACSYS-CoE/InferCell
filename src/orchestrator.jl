"""
    build_problem(models; tspan=(0.0, 100.0))
    build_problem(model;  tspan=(0.0, 100.0))

Compose one or more [`AbstractSubModel`](@ref)s into a SciML problem object.
Returns an `ODEProblem` when every sub-model has `formalism = :ode`, or a
`JumpProblem` when every sub-model has `formalism = :jump`; mixed-formalism
composition is not yet supported.

Shared parameters across sub-models are deduplicated by name and must agree
on `value`, `prior`, and `fixed`. Cross-block coupling is resolved via each
sub-model's [`inputs`](@ref) — input symbols must be `states` of some other
sub-model in the composition.
"""
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
    _validate_shared_params(models)
    contexts = _build_contexts(models)
    _resolve_coupling(models, contexts)

    u0 = _build_u0(models)
    p0 = _build_p0(models)
    rhs = _build_rhs(models, contexts)

    return ODEProblem{false}(rhs, u0, tspan, p0)
end

function _validate_shared_params(models::Vector{<:AbstractSubModel})
    seen = Dict{Symbol, InferParameter}()
    for m in models
        for p in model_free_params(parameters(m))
            if haskey(seen, p.name)
                existing = seen[p.name]
                if p.value != existing.value || typeof(p.prior) != typeof(existing.prior) || p.fixed != existing.fixed
                    error("Shared parameter :$(p.name) has inconsistent definitions across modules :$(existing.module_id) and :$(p.module_id)")
                end
            else
                seen[p.name] = p
            end
        end
    end
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

    # Build global parameter map, deduplicating by name
    global_param_map = Dict{Symbol, Int}()
    global_param_count = 0
    model_param_indices = Vector{Vector{Int}}()

    for m in models
        mfp = model_free_params(parameters(m))
        indices = Int[]
        for p in mfp
            if haskey(global_param_map, p.name)
                push!(indices, global_param_map[p.name])
            else
                global_param_count += 1
                global_param_map[p.name] = global_param_count
                push!(indices, global_param_count)
            end
        end
        push!(model_param_indices, indices)
    end

    for (i, m) in enumerate(models)
        n_states = length(states(m))
        state_idxs = (state_offset + 1):(state_offset + n_states)
        push!(contexts, SubModelContext(state_idxs, model_param_indices[i], Dict{Symbol, Int}()))
        state_offset += n_states
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

function _collect_ic_values(models::Vector{<:AbstractSubModel})
    vals = Float64[]
    for m in models
        ics = ic_params(parameters(m))
        for s in states(m)
            ic = findfirst(p -> p.name == Symbol(s, "0") || p.name == s, ics)
            push!(vals, ic !== nothing ? ics[ic].value : 0.0)
        end
    end
    return vals
end

_build_u0(models::Vector{<:AbstractSubModel}) =
    (v = _collect_ic_values(models); SVector{length(v)}(v))

_build_u0_integer(models::Vector{<:AbstractSubModel}) =
    round.(Int, _collect_ic_values(models))

function _build_p0(models::Vector{<:AbstractSubModel})
    seen = Set{Symbol}()
    vals = Float64[]
    for m in models
        for p in model_free_params(parameters(m))
            if !(p.name in seen)
                push!(seen, p.name)
                push!(vals, p.value)
            end
        end
    end
    return vals
end

function _build_rhs(models::Vector{<:AbstractSubModel},
                    contexts::Vector{SubModelContext})
    model_inputs = [inputs(m) for m in models]
    function rhs(u, p, t)
        du_parts = map(models, contexts, model_inputs) do m, ctx, inp_syms
            u_local = u[ctx.state_idxs]
            p_local = p[ctx.param_idxs]
            u_inputs = [u[ctx.input_map[s]] for s in inp_syms]
            dynamics(u_local, p_local, t, m, u_inputs)
        end
        return vcat(du_parts...)
    end
    return rhs
end
