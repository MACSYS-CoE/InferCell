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

    # Each module's reactions are written in local coordinates; the slot maps
    # them onto the composed vectors (spec §11 phase 2). This is the jump-side
    # twin of `_build_rhs`, without the tuple machinery: the SSA re-evaluates
    # every propensity per event through JumpProcesses' own dispatch, so there
    # is no static composed function to keep type-stable.
    jumps = ConstantRateJump[]
    for (m, ctx) in zip(models, contexts)
        slot = _jump_slot(m, ctx)
        for r in reactions(m)
            push!(jumps, _global_jump(r, slot))
        end
    end

    dprob = DiscreteProblem(u0, tspan, p0)
    return JumpProblem(dprob, Direct(), JumpSet(; constant_jumps=jumps))
end

# Where a jump module's states, parameters and inputs sit in the composed
# vectors, and which of its inputs it has declared it writes.
struct _JumpSlot{S, P, I, W, N}
    sidx::S      # global indices of this module's states
    pidx::P      # global indices of its free parameters
    iidx::I      # global indices of its inputs, in inputs(m) order
    writable::W  # per input, whether written_states(m) declares the write
    inames::N    # the input names, in the same order, for the refusal message
    mod::Symbol  # module_id(m), likewise
end

function _jump_slot(m::AbstractSubModel, ctx::SubModelContext)
    ins = inputs(m)
    writes = written_states(m)
    return _JumpSlot(_svec(ctx.state_idxs), _svec(ctx.param_idxs),
                     _svec(Int[ctx.input_map[s] for s in ins]),
                     SVector{length(ins), Bool}(Tuple(s in writes for s in ins)),
                     SVector{length(ins), Symbol}(Tuple(ins)),
                     module_id(m))
end

"""
    PeerView

What a jump module's `Reaction` receives as `u_inputs`: a live view of the
composed state vector at the module's declared inputs, in [`inputs`](@ref)
order. Reading is always allowed. Writing is allowed only where
[`written_states`](@ref) declares it; any other write throws an `ArgumentError`
naming the module and the state, so a peer write that no declaration covers
fails at its first firing rather than landing silently (spec §11 task 2.5, as
amended 2026-09-04).
"""
struct PeerView{T, U <: AbstractVector{T}, S <: _JumpSlot} <: AbstractVector{T}
    u::U
    slot::S
end

Base.size(v::PeerView) = (length(v.slot.iidx),)
Base.IndexStyle(::Type{<:PeerView}) = IndexLinear()
Base.getindex(v::PeerView, i::Int) = v.u[v.slot.iidx[i]]
function Base.setindex!(v::PeerView, x, i::Int)
    sl = v.slot
    sl.writable[i] || throw(ArgumentError(
        "Module $(sl.mod) writes :$(sl.inames[i]), which it does not own and does " *
        "not list in written_states(). A jump module may read any state in " *
        "inputs() but may write only those it declares; add :$(sl.inames[i]) to " *
        "written_states() or make the affect leave it alone"))
    v.u[sl.iidx[i]] = x
    return v
end

# The rate sees the local slice, the local parameters and the inputs as views
# of the global vectors, so `remake(prob; p = θ)` keeps working: nothing here
# captures `p`. The affect mutates the same views, so a write lands on the
# owner's global slot and nowhere else.
function _global_jump(r::Reaction, sl::_JumpSlot)
    rate = (u, p, t) -> r.rate(view(u, sl.sidx), view(p, sl.pidx), t, PeerView(u, sl))
    affect! = integrator -> (r.affect!(view(integrator.u, sl.sidx), PeerView(integrator.u, sl)); nothing)
    return ConstantRateJump(rate, affect!)
end

# Anything else — a `ConstantRateJump`, the pre-phase-2 contract — is refused
# by name rather than composed: its closures index the global vectors, which
# is the aliasing bug this path exists to prevent (spec §2, G3).
_global_jump(r, sl::_JumpSlot) = throw(ArgumentError(
    "Module $(sl.mod) returned a $(nameof(typeof(r))) from reactions(). " *
    "Its rate and affect closures index the composed state and parameter " *
    "vectors at this module's local positions, so a second jump module in " *
    "the composition would read and write the first module's slice. Return " *
    "Reaction(rate, affect!) with rate(u, p, t, u_inputs) and " *
    "affect!(u, u_inputs) written in local coordinates instead"))

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
    # registry, chemostatted species and contributions from a :jump module,
    # none of which needs a composition to judge; ownership does, and
    # standalone resolution must keep reporting an unowned target rather than
    # failing on it, so that check lives here.
    for (i, m) in enumerate(models)
        targets = contributed_states(m)
        isempty(targets) && continue
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

# The composed right-hand side.
#
# State-container decision, re-taken at 32 states (spec §11 task 1.6; numbers in
# dev/scripts/bench_rhs_containers_result.md): the problem stays out-of-place
# over an `SVector`. What makes that allocation-free is that every index the
# closure uses is a static vector fixed at build time — a `UnitRange` slice of an
# `SVector` returns a heap `Vector`, a `Vector{Int}` index allocates, and a
# splatted `vcat` over a `Vector` of parts is not inferable — and that the
# modules are held as a `Tuple`, so `map` and `reduce` unroll. Base stops
# unrolling tuple `map` at 32 elements (`Base.Any32` matches any tuple of 32 or
# more), hence the cap of 31; Core A′ is seven modules. The framework adds no
# allocation; the composed function is allocation-free iff every `dynamics` and
# `contributions` is.
struct _Slot{S, P, I, C}
    sidx::S      # global indices of this module's states
    pidx::P      # global indices of its free parameters
    iidx::I      # global indices of its inputs, in inputs(m) order
    cidx::C      # global indices of its contributed states, in contributed_states(m) order
end

_svec(v) = SVector{length(v), Int}(Tuple(v))

function _build_rhs(models::Vector{<:AbstractSubModel},
                    contexts::Vector{SubModelContext})
    length(models) <= 31 || error(
        "build_problem composes at most 31 sub-models: Base.map over a tuple of " *
        "32 or more is not unrolled (Base.Any32), so the right-hand side would " *
        "allocate and lose type stability. $(length(models)) were given")
    slots = Tuple(_Slot(_svec(ctx.state_idxs), _svec(ctx.param_idxs),
                        _svec(Int[ctx.input_map[s] for s in inputs(m)]),
                        _svec(ctx.contrib_idxs))
                  for (m, ctx) in zip(models, contexts))
    # Build-time bookkeeping ends here; _make_rhs holds only the runtime closure.
    return _make_rhs(Tuple(models), slots)
end

function _make_rhs(models::Tuple, slots::Tuple)
    function rhs(u, p, t)
        parts = map(models, slots) do m, sl
            u_local = u[sl.sidx]
            p_local = p[sl.pidx]
            u_inputs = u[sl.iidx]
            (dynamics(u_local, p_local, t, m, u_inputs),
             contributions(u_local, p_local, t, m, u_inputs))
        end
        du = reduce(vcat, map(first, parts))
        return _fold_contributions(du, slots, map(last, parts))
    end
    return rhs
end

@inline _fold_contributions(du, ::Tuple{}, ::Tuple{}) = du
@inline _fold_contributions(du, slots::Tuple, cs::Tuple) =
    _fold_contributions(_accumulate(du, first(slots).cidx, first(cs)),
                        Base.tail(slots), Base.tail(cs))

# A module that contributes nothing leaves `du` untouched — same object, same
# type — which is what keeps a non-contributing composition byte-identical to
# the pre-phase-1 path.
@inline _accumulate(du::SVector, ::SVector{0, Int}, ::SVector{0}) = du

# Add each contribution to the owner's derivative. The eltype is promoted once
# so a Float64 derivative can receive a dual-number contribution under
# automatic differentiation; `setindex` would otherwise refuse the conversion.
@inline function _accumulate(du::SVector{N}, idxs::SVector{C, Int}, c::SVector{C}) where {N, C}
    out = SVector{N, promote_type(eltype(du), eltype(c))}(du)
    for j in 1:C
        i = idxs[j]
        out = setindex(out, out[i] + c[j], i)
    end
    return out
end

# Anything else — the wrong number of terms, or a non-static container — is a
# protocol error; say so rather than failing inside StaticArrays.
_accumulate(::SVector, ::SVector{C, Int}, c::AbstractVector) where {C} = error(
    "contributions() returned a $(typeof(c)) with $(length(c)) " *
    "term$(length(c) == 1 ? "" : "s") but contributed_states() names $C species; " *
    "return a static vector of length $C (e.g. SA[...])")
