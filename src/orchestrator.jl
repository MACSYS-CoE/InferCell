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

    # Equal values from different source files pass every check above — the
    # values agree, so nothing fires, and deduplication silently keeps whichever
    # module was seen first. That is the cross-file trap in its most easily
    # missed form, so it is reported rather than resolved by first-come.
    all_params = InferParameter[p for m in models for p in parameters(m)]
    for (name, files) in provenance_conflicts(all_params)
        # Once per parameter per session: an infer + posterior-predictive run
        # rebuilds the problem many times, and the conflict does not change.
        @warn "Parameter :$name is imported from more than one source file " *
              "($(join(files, " and "))). Deduplication keeps the first-seen " *
              "definition; declare which file governs rather than letting " *
              "composition order decide." maxlog = 1 _id = Symbol(:provenance_conflict_, name)
    end
    return nothing
end

function _build_jump_problem(models::Vector{<:AbstractSubModel}; tspan=(0.0, 100.0))
    _validate_shared_params(models)
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
    # Validate the Core A′ coupling contract before wiring anything up. This is
    # a no-op for sub-models outside Core A′: they name no registry species and
    # declare no typed edges, so nothing here fires for them.
    resolve_coupling(models)

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
            if !haskey(state_owners, inp)
                # A chemostat can never be owned — the resolver forbids
                # integrating one — so inputs() cannot deliver it. Wiring the
                # registry's held value into dynamics is coupling execution,
                # which wave 2 owns; until then the value travels as a fixed
                # parameter beside a declared ClampedEdge.
                if is_registered(inp) && is_chemostatted(inp)
                    error(
                        "Input :$inp declared by $(typeof(m)) is chemostatted by " *
                        "the Core A′ registry, so no sub-model integrates it and " *
                        "inputs() cannot deliver it. Remove :$inp from inputs(); " *
                        "declare a ClampedEdge with its held_value and carry the " *
                        "value as a fixed parameter instead")
                end
                error(
                    "Input :$inp declared by $(typeof(m)) is not owned by any sub-model")
            end
            contexts[i].input_map[inp] = state_owners[inp]
        end
    end

    # The contribution channel: each declared target resolves to the owner's
    # global index. The resolver has already rejected names outside the
    # registry and chemostatted species, which need no composition to judge;
    # ownership does, and standalone resolution must keep reporting an unowned
    # target rather than failing on it, so that check lives here.
    for (i, m) in enumerate(models)
        targets = contributed_states(m)
        isempty(targets) && continue
        formalism(m) === :jump && error(
            "Module $(module_id(m)) has formalism :jump but declares contributions " *
            "to $(join(string.(":", targets), ", ")). A jump process cannot add " *
            "continuous derivative terms; a jump module's writes to a peer's state " *
            "go through its reactions (spec §11 phase 2)")
        for s in targets
            haskey(state_owners, s) || error(
                "Module $(module_id(m)) contributes to :$s, which no sub-model in " *
                "this composition owns. A contribution needs an owner to add to; " *
                "compose the module that integrates :$s, or drop the contribution")
            push!(contexts[i].contrib_idxs, state_owners[s])
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
    # Which modules contribute is known at build time, so a composition with
    # no contributors runs the pre-phase-1 computation unchanged.
    contributors = [i for (i, ctx) in enumerate(contexts) if !isempty(ctx.contrib_idxs)]
    function rhs(u, p, t)
        du_parts = map(models, contexts, model_inputs) do m, ctx, inp_syms
            u_local = u[ctx.state_idxs]
            p_local = p[ctx.param_idxs]
            u_inputs = [u[ctx.input_map[s]] for s in inp_syms]
            dynamics(u_local, p_local, t, m, u_inputs)
        end
        du = vcat(du_parts...)
        isempty(contributors) && return du
        for i in contributors
            m, ctx = models[i], contexts[i]
            u_local = u[ctx.state_idxs]
            p_local = p[ctx.param_idxs]
            u_inputs = [u[ctx.input_map[s]] for s in model_inputs[i]]
            c = contributions(u_local, p_local, t, m, u_inputs)
            du = _accumulate(du, ctx.contrib_idxs, c)
        end
        return du
    end
    return rhs
end

# Add each contribution to the owner's derivative. `du` is an SVector, so the
# result is a new SVector; the eltype is promoted once so a Float64 derivative
# can receive a dual-number contribution under automatic differentiation.
function _accumulate(du::SVector{N}, idxs, c) where {N}
    length(idxs) == length(c) || error(
        "contributions() returned $(length(c)) term$(length(c) == 1 ? "" : "s") " *
        "but contributed_states() names $(length(idxs)) species")
    out = SVector{N, promote_type(eltype(du), eltype(c))}(du)
    for (j, i) in enumerate(idxs)
        out = setindex(out, out[i] + c[j], i)
    end
    return out
end
