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
    contributed_states(m::AbstractSubModel) -> Vector{Symbol}

Registry states integrated by *other* sub-models to which `m` adds derivative
terms — the outbound counterpart of [`inputs`](@ref). Defaults to empty.

Where `inputs` lets a module read a state it does not own, this lets it write
into one: the orchestrator adds each term returned by [`contributions`](@ref)
to the owner's derivative at the owner's global index, so a shared pool
receives terms from every module that produces into or draws from it. The
names must be Core A′ registry species and may not be chemostatted, which
[`resolve_coupling`](@ref) checks, and must be owned by some module in the
composition, which the orchestrator checks at build time; each rejection names
the species and this module.

For an ODE module the resolver holds this list to the coupling edges in both
directions: every mass or currency edge, inbound or outbound, on a dynamic
registry species `m` does not own must appear here, and every entry here must
have such an edge — a module listing contributions with no edges is checked,
not skipped. A consumer of a foreign pool contributes a negative term; a
producer a positive one. A `:jump` module may not list contributions at all;
its writes to a peer's state go through its reactions.
"""
contributed_states(::AbstractSubModel) = Symbol[]

"""
    contributions(u, p, t, m::AbstractSubModel[, u_inputs])

Signed derivative terms `m` adds to the states named by
[`contributed_states`](@ref), in that order, as a static vector of the same
length. Receives exactly what [`dynamics`](@ref) receives: the local state
slice, local parameters, time and resolved inputs. Defaults to an empty static
vector, so a sub-model that contributes nothing needs no method and takes the
same code path through the right-hand side as before this channel existed.

This is the phase-1 mechanism of `spec/spec.md` (§11, task 1.1). It was chosen
over a second return value from `dynamics`, which would change the return
arity of every existing model, and over an in-place global write buffer, which
would force the composed problem in-place and every `dynamics` to a mutating
form. Being a separate function with an empty default, it adds nothing to a
module that does not use it.
"""
contributions(u, p, t, ::AbstractSubModel) = SVector{0, Float64}()
# Fallback, as for `dynamics`: a module that defines only the four-argument
# form is still reached from the orchestrator's five-argument call.
contributions(u, p, t, m::AbstractSubModel, u_inputs) = contributions(u, p, t, m)

"""
    SubModelContext(state_idxs, param_idxs, input_map[, contrib_idxs])

Per-sub-model bookkeeping built by the orchestrator: where this sub-model's
states live in the global state vector, which global parameter indices it reads,
how its declared input symbols map onto global state indices, and the global
index of each state in [`contributed_states`](@ref), in that order.
"""
struct SubModelContext
    state_idxs::UnitRange{Int}
    param_idxs::Vector{Int}
    input_map::Dict{Symbol, Int}
    contrib_idxs::Vector{Int}
end

SubModelContext(state_idxs, param_idxs, input_map) =
    SubModelContext(state_idxs, param_idxs, input_map, Int[])

const VALID_FORMALISMS = (:ode, :sde, :jump)

function validate_formalism(f::Symbol)
    f in VALID_FORMALISMS || throw(ArgumentError(
        "Invalid formalism :$f. Must be one of $VALID_FORMALISMS"))
end
