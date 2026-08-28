"""
    AbstractSubModel

Abstract supertype for every InferCell block (TX/TL, metabolism, stochastic gene
expression, …). Concrete sub-types implement the protocol functions
[`states`](@ref), [`parameters`](@ref), [`dynamics`](@ref) and/or
[`reactions`](@ref), and optionally override [`inputs`](@ref),
[`formalism`](@ref), and [`inference_mode`](@ref).
"""
abstract type AbstractSubModel end

"""
    states(m::AbstractSubModel) -> Vector{Symbol}

Ordered state-variable names for sub-model `m`. The orchestrator concatenates
state vectors across sub-models in this order.
"""
function states end

"""
    parameters(m::AbstractSubModel) -> Vector{InferParameter}

All parameters of `m` — rates, initial conditions, and observation parameters
combined. The inference layer filters by role and `fixed` to decide which are
sampled.
"""
function parameters end

"""
    dynamics(u, p, t, m::AbstractSubModel[, u_inputs])

Out-of-place right-hand side for ODE sub-models. Returns `du` for the local
state slice `u` given local parameters `p`, time `t`, and (optionally) the
coupling inputs `u_inputs` resolved by the orchestrator.
"""
function dynamics end

"""
    inputs(m::AbstractSubModel) -> Vector{Symbol}

State names that `m` reads from *other* sub-models (its coupling inputs).
Defaults to empty — sub-models with no cross-block coupling need not override.
"""
inputs(::AbstractSubModel) = Symbol[]

"""
    coupling(m::AbstractSubModel) -> Vector{CouplingEdge}

Typed coupling declarations for `m`: not only *which* state crosses a boundary
but *how*, as one of the seven kinds in [`EDGE_KINDS`](@ref). Defaults to empty,
so a sub-model with no cross-block coupling — and every sub-model written before
the Core A′ interface contract — needs no modification.

[`inputs`](@ref) keeps its own meaning and its own default. A sub-model may
declare both; where it declares any coupling edge, [`resolve_coupling`](@ref)
checks the two against each other in both directions — every registry-species
input has an inbound edge, and every inbound mass or currency edge on a state
the module does not integrate is listed in `inputs` — so a typed declaration
that drifts out of step with `inputs` is caught rather than silently
disagreeing. Inputs naming non-registry state (the legacy blocks' mRNA and
protein) are outside the typed contract and pass through unchecked.
"""
coupling(::AbstractSubModel) = CouplingEdge[]

"""
    module_id(m::AbstractSubModel) -> Symbol

Identity of a sub-model, used to name it in coupling errors and to match the
`peer` of a [`CouplingEdge`](@ref). Defaults to the type name.

Each Core A′ module is a distinct type, so the default is usually right. A
sub-model override matters when one type is instantiated more than once in a
composition, where the type name would name both instances identically.
"""
module_id(m::AbstractSubModel) = nameof(typeof(m))

"""
    reduction_notes(m::AbstractSubModel) -> Vector{String}

Simplifications this sub-model introduces that the published model does not
make. Defaults to empty.

Core A′ carries several treatments that are ours rather than the source model's
— the lumped tRNA charging step most of all — and each is a place where a result
could depend on our choice rather than on the published model. Registering them
here is what lets [`reduction_declarations`](@ref) enumerate them, so the label
can reach any report that depends on one.
"""
reduction_notes(::AbstractSubModel) = String[]

"""
    formalism(m::AbstractSubModel) -> Symbol

One of `:ode`, `:sde`, or `:jump`. Controls how the orchestrator assembles the
problem object. Defaults to `:ode`.
"""
formalism(::AbstractSubModel) = :ode

"""
    inference_mode(m::AbstractSubModel) -> Symbol

`:differentiable` (NUTS via Turing) or `:simulation` (ABC-SMC). Controls which
backend [`infer`](@ref) dispatches to. Defaults to `:differentiable`.
"""
inference_mode(::AbstractSubModel) = :differentiable

"""
    reactions(m::AbstractSubModel) -> Vector{ConstantRateJump}

Jump-process reactions for SSA sub-models. Required for `formalism = :jump`;
defaults to an error to catch incomplete sub-model definitions.
"""
reactions(m::AbstractSubModel) = error("reactions() not implemented for $(typeof(m))")

# Fallback: models with no inputs ignore u_inputs
dynamics(u, p, t, m::AbstractSubModel, u_inputs) = dynamics(u, p, t, m)

"""
    SubModelContext(state_idxs, param_idxs, input_map)

Per-sub-model bookkeeping built by the orchestrator: where this sub-model's
states live in the global state vector, which global parameter indices it reads,
and how its declared input symbols map onto global state indices.
"""
struct SubModelContext
    state_idxs::UnitRange{Int}
    param_idxs::Vector{Int}
    input_map::Dict{Symbol, Int}
end

const VALID_FORMALISMS = (:ode, :sde, :jump)

function validate_formalism(f::Symbol)
    f in VALID_FORMALISMS || throw(ArgumentError(
        "Invalid formalism :$f. Must be one of $VALID_FORMALISMS"))
end
