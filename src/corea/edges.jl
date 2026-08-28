"""
Coupling edge kinds.

Seven ways state crosses a module boundary in Core A′, one per entry in the
legend of `dev/notes/figures/reduced-syn3a-coupling/fig1r_state_graph_reduced.pdf`.
The figure is the source of truth for this list: if the figure and these types
disagree, one of them is a bug.

The kinds are closed. An eighth way for state to cross an interface is a change
to the figure and to this file together, not a new subtype declared elsewhere —
[`resolve_coupling`](@ref) dispatches on these types and fails on anything else
rather than falling through to a generic method.

Two of the fields carry decisions the scoping note leaves open. A
[`DeferredCounterEdge`](@ref) declares how it handles a debit larger than the
pool, and a [`RateConstantEdge`](@ref) declares how often it refreshes. Both
default to what the published model does, so departing from the published model
is something an author has to write down.
"""

const EDGE_DIRECTIONS = (:in, :out)

"""
Clipping policies for a deferred counter.

- `:clamped_deficit_carried` — the published model: the pool floors at zero and
  the shortfall is carried to the next step. A `max(0, ·)` on the interface, and
  therefore non-differentiable.
- `:unclamped` — the pool is allowed to go negative.
- `:smoothed` — the discontinuity is replaced by a differentiable approximation.

The published clamped drain is the default; either other policy is a deviation
from the published model, and labelled as one.
"""
const CLIP_POLICIES = (:clamped_deficit_carried, :unclamped, :smoothed)

"""
Refresh cadences for a rate-constant edge.

- `:piecewise_constant` — the published model: rate constants are recomputed on
  a fixed interval and held between refreshes. The production driver rebuilds
  the stochastic block every 60 s.
- `:continuous` — propensities track the upstream pools continuously. Arguably a
  better model, definitely a different one, so it is labelled as a deviation.
"""
const RATE_CADENCES = (:piecewise_constant, :continuous)

"""
Whether a simplification is the published model's or this reduction's.

`:ours` is the label that has to survive into any result that depends on the
simplification — the chemostatted CTP, UTP and amino-acid pools, and the lumped
tRNA charging step, are ours rather than the published model's.
"""
const EDGE_ORIGINS = (:published, :ours)

"""
    CouplingEdge

Supertype for the seven coupling-edge kinds. Every edge names a registry
species, a direction (`:in` where the declaring module reads or consumes,
`:out` where it writes or produces), and optionally a peer module.

A `peer` of `nothing` means "whichever module owns this species", which is the
usual case inside the ODE block; naming a peer explicitly is what lets the
resolver report a counterpart that is missing from the composition.
"""
abstract type CouplingEdge end

_require(kind::Symbol, field::Symbol, value) =
    value === nothing && throw(ArgumentError(
        "$kind requires the field `$field`, which was not given. " *
        "Every $kind must name it; see src/corea/edges.jl"))

function _check_common(kind::Symbol, species, direction)
    _require(kind, :species, species)
    _require(kind, :direction, direction)
    _check_vocab(kind, :direction, direction, EDGE_DIRECTIONS)
    return nothing
end

"""
    MassEdge(; species, direction, peer=nothing)

Shared state, continuous. Gradients cross inside the ODE block. The plainest
kind: two modules integrating a common species.
"""
struct MassEdge <: CouplingEdge
    species::Symbol
    direction::Symbol
    peer::Union{Symbol, Nothing}

    # Validation lives in the inner constructor so that positional
    # construction cannot bypass it. Same pattern on all seven kinds.
    function MassEdge(species, direction, peer)
        _check_common(:MassEdge, species, direction)
        return new(species, direction, peer)
    end
end

MassEdge(; species=nothing, direction=nothing, peer=nothing) =
    MassEdge(species, direction, peer)

"""
    CurrencyEdge(; species, direction, peer=nothing, pool=species)

Routed via a shared pool node rather than directly module to module — how ATP
and GTP traffic is drawn in fig 1. `pool` names the registry species acting as
the routing node, which for the energy currencies is the species itself.
"""
struct CurrencyEdge <: CouplingEdge
    species::Symbol
    direction::Symbol
    peer::Union{Symbol, Nothing}
    pool::Symbol

    function CurrencyEdge(species, direction, peer, pool)
        _check_common(:CurrencyEdge, species, direction)
        _require(:CurrencyEdge, :pool, pool)
        return new(species, direction, peer, pool)
    end
end

CurrencyEdge(; species=nothing, direction=nothing, peer=nothing, pool=nothing) =
    CurrencyEdge(species, direction, peer, pool === nothing ? species : pool)

"""
    DeferredCounterEdge(; species, direction, counter, peer=nothing,
                          clip=:clamped_deficit_carried, smoothing=nothing)

The stochastic block accrues a cost in `counter`; the hook debits it against
`species` one step later. This is the kind that carries the interface's
`max(0, ·)`: see [`CLIP_POLICIES`](@ref), and [`obstructs_gradients`](@ref) for
what a clamped policy costs a gradient-based sampler.

`smoothing` is the width of the differentiable approximation and is required
exactly when `clip = :smoothed` — the parameter controlling the smoothing is
exposed rather than hidden, so the deviation is fully specified where it is
declared. The other policies take no smoothing and reject one.
"""
struct DeferredCounterEdge <: CouplingEdge
    species::Symbol
    direction::Symbol
    peer::Union{Symbol, Nothing}
    counter::Symbol
    clip::Symbol
    smoothing::Union{Float64, Nothing}

    function DeferredCounterEdge(species, direction, peer, counter, clip, smoothing)
        _check_common(:DeferredCounterEdge, species, direction)
        _require(:DeferredCounterEdge, :counter, counter)
        _check_vocab(:DeferredCounterEdge, :clip, clip, CLIP_POLICIES)
        if clip === :smoothed
            smoothing === nothing && throw(ArgumentError(
                "DeferredCounterEdge with clip = :smoothed requires the field " *
                "`smoothing` — the width of the differentiable approximation is " *
                "part of the deviation, not an implementation detail"))
            smoothing > 0 || throw(ArgumentError(
                "DeferredCounterEdge field `smoothing` must be positive, got $smoothing"))
            return new(species, direction, peer, counter, clip, Float64(smoothing))
        end
        smoothing === nothing || throw(ArgumentError(
            "DeferredCounterEdge with clip = :$clip takes no `smoothing`; " *
            "carrying one would imply an approximation it does not make"))
        return new(species, direction, peer, counter, clip, nothing)
    end
end

DeferredCounterEdge(; species=nothing, direction=nothing, peer=nothing,
                    counter=nothing, clip=:clamped_deficit_carried,
                    smoothing=nothing) =
    DeferredCounterEdge(species, direction, peer, counter, clip, smoothing)

"""
    CatalyticEdge(; species, direction, param_slot, peer=nothing)

Counts enter a rate law as parameters. **No mass flows**: a catalytic edge
contributes to no balance, and [`mass_contribution`](@ref) rejects it rather
than counting it as zero. `param_slot` names the rate-law parameter the count
fills.
"""
struct CatalyticEdge <: CouplingEdge
    species::Symbol
    direction::Symbol
    peer::Union{Symbol, Nothing}
    param_slot::Symbol

    function CatalyticEdge(species, direction, peer, param_slot)
        _check_common(:CatalyticEdge, species, direction)
        _require(:CatalyticEdge, :param_slot, param_slot)
        return new(species, direction, peer, param_slot)
    end
end

CatalyticEdge(; species=nothing, direction=nothing, peer=nothing,
              param_slot=nothing) =
    CatalyticEdge(species, direction, peer, param_slot)

"""
    RateConstantEdge(; species, direction, peer=nothing,
                       cadence=:piecewise_constant, interval=60.0)

Pools re-enter the stochastic block as recomputed rate constants. The only
ODE→stochastic channel in the model, and so the edge that decides whether the
coupling is bidirectional at all.

No mass moves — the pool is read, not consumed — so `direction` follows the
information: the module whose rate constants are rebuilt from the pool declares
`:in`, and the module owning the pool declares `:out` where it names the
channel from its own side.

Defaults to the published model's 60 s rebuild: an unspecified `interval` is
60 s for a `:piecewise_constant` cadence and absent for a `:continuous` one,
which takes no interval and is marked as a deviation. An interval passed
explicitly is always validated, so `:continuous` alongside one is an error
rather than a silent discard.
"""
struct RateConstantEdge <: CouplingEdge
    species::Symbol
    direction::Symbol
    peer::Union{Symbol, Nothing}
    cadence::Symbol
    interval::Union{Float64, Nothing}

    function RateConstantEdge(species, direction, peer, cadence, interval)
        _check_common(:RateConstantEdge, species, direction)
        _check_vocab(:RateConstantEdge, :cadence, cadence, RATE_CADENCES)
        if cadence === :piecewise_constant
            interval === nothing && throw(ArgumentError(
                "RateConstantEdge requires the field `interval` when cadence is " *
                ":piecewise_constant — the published model rebuilds every 60 s"))
            interval > 0 || throw(ArgumentError(
                "RateConstantEdge field `interval` must be positive, got $interval"))
            return new(species, direction, peer, cadence, Float64(interval))
        end
        interval === nothing || throw(ArgumentError(
            "RateConstantEdge with cadence = :continuous takes no `interval`; " *
            "carrying one would imply a cadence it does not have"))
        return new(species, direction, peer, cadence, nothing)
    end
end

# The interval default adapts to the cadence; an explicit value — including an
# explicit `nothing` — reaches the inner constructor's validation untouched.
RateConstantEdge(; species=nothing, direction=nothing, peer=nothing,
                 cadence=:piecewise_constant,
                 interval=(cadence === :continuous ? nothing : 60.0)) =
    RateConstantEdge(species, direction, peer, cadence, interval)

"""
    VolumeEdge(; species, direction, peer=nothing)

Counts set surface area, which sets volume, which rescales every concentration.
The second coupling channel from the stochastic block into the deterministic
one, and the cheapest.
"""
struct VolumeEdge <: CouplingEdge
    species::Symbol
    direction::Symbol
    peer::Union{Symbol, Nothing}

    function VolumeEdge(species, direction, peer)
        _check_common(:VolumeEdge, species, direction)
        return new(species, direction, peer)
    end
end

VolumeEdge(; species=nothing, direction=nothing, peer=nothing) =
    VolumeEdge(species, direction, peer)

"""
    ClampedEdge(; species, direction, origin, peer=nothing, held_value=nothing)

A real dependence replaced by a constant. `origin` is required and has no
default on purpose: whether a clamp is the published model's or ours is the
thing this kind exists to record, so an author has to state it.
"""
struct ClampedEdge <: CouplingEdge
    species::Symbol
    direction::Symbol
    peer::Union{Symbol, Nothing}
    held_value::Union{Float64, Nothing}
    origin::Symbol

    function ClampedEdge(species, direction, peer, held_value, origin)
        _check_common(:ClampedEdge, species, direction)
        _require(:ClampedEdge, :origin, origin)
        _check_vocab(:ClampedEdge, :origin, origin, EDGE_ORIGINS)
        return new(species, direction, peer,
                   held_value === nothing ? nothing : Float64(held_value),
                   origin)
    end
end

ClampedEdge(; species=nothing, direction=nothing, peer=nothing,
            held_value=nothing, origin=nothing) =
    ClampedEdge(species, direction, peer, held_value, origin)

"""
    EDGE_KINDS

The seven legend names, in the figure's order.
"""
const EDGE_KINDS = (:mass, :currency, :deferred_counter, :catalytic,
                    :rate_constant, :volume, :clamped)

"""
    edge_kind(e::CouplingEdge) -> Symbol

The legend name of an edge's kind: `:mass`, `:currency`, `:deferred_counter`,
`:catalytic`, `:rate_constant`, `:volume` or `:clamped`.
"""
edge_kind(::MassEdge) = :mass
edge_kind(::CurrencyEdge) = :currency
edge_kind(::DeferredCounterEdge) = :deferred_counter
edge_kind(::CatalyticEdge) = :catalytic
edge_kind(::RateConstantEdge) = :rate_constant
edge_kind(::VolumeEdge) = :volume
edge_kind(::ClampedEdge) = :clamped
edge_kind(e::CouplingEdge) = throw(ArgumentError(
    "$(typeof(e)) is not one of the seven Core A′ edge kinds. " *
    "The kinds are closed: $(EDGE_KINDS)"))

"""
    is_consumer(e::CouplingEdge) -> Bool
    is_producer(e::CouplingEdge) -> Bool

Whether the declaring module draws the species down or supplies it. Direction
`:in` consumes, `:out` produces. Only the mass-carrying kinds count for the
dead-end check; see [`carries_mass`](@ref).
"""
is_consumer(e::CouplingEdge) = e.direction === :in
is_producer(e::CouplingEdge) = e.direction === :out

"""
    carries_mass(e::CouplingEdge) -> Bool

Whether the edge moves matter. True for mass, currency and deferred-counter
edges; false for the four that carry information rather than material —
catalytic, rate-constant, volume and clamped.

This is what the dead-end and conservation checks filter on: an edge that
carries no mass cannot create or strand a moiety.
"""
carries_mass(::MassEdge) = true
carries_mass(::CurrencyEdge) = true
carries_mass(::DeferredCounterEdge) = true
carries_mass(::CatalyticEdge) = false
carries_mass(::RateConstantEdge) = false
carries_mass(::VolumeEdge) = false
carries_mass(::ClampedEdge) = false

"""
    mass_contribution(e::CouplingEdge) -> Int

The sign an edge contributes to a mass or moiety balance: `-1` where it
consumes, `+1` where it produces.

Throws for a [`CatalyticEdge`](@ref). Counts entering a rate law move no mass,
so including one in a conservation check is a category error, and returning a
silent zero would let that error pass unnoticed.
"""
function mass_contribution(e::CouplingEdge)
    carries_mass(e) || throw(ArgumentError(
        "A $(edge_kind(e)) edge on :$(e.species) carries no mass, so it cannot " *
        "appear in a conservation check. Counts entering a rate law as " *
        "parameters move no matter; filter with carries_mass first"))
    return is_producer(e) ? 1 : -1
end

"""
    obstructs_gradients(e::CouplingEdge) -> Bool

Whether the edge is non-differentiable at the boundary. True only for a
deferred counter under the published clamped policy — the `max(0, ·)` that will
obstruct NUTS if the ODE block is sampled with a gradient-based method.
"""
obstructs_gradients(e::DeferredCounterEdge) = e.clip === :clamped_deficit_carried
obstructs_gradients(::CouplingEdge) = false

"""
    deviates_from_published(e::CouplingEdge) -> Bool

Whether the edge departs from the published model. True for a smoothed or
unclamped deferred counter, a continuous rate-constant edge, and a clamp this
reduction introduced — exactly the edges [`deviation_reason`](@ref) has a
reason for.

This is the predicate behind the "what is ours" enumeration: a deviation that
reaches a result unlabelled is the failure mode it exists to prevent.
"""
deviates_from_published(e::CouplingEdge) = deviation_reason(e) !== nothing

"""
    deviation_reason(e::CouplingEdge) -> Union{String, Nothing}

A one-line description of how an edge departs from the published model, or
`nothing` where it does not. Written for a report a human reads.
"""
function deviation_reason(e::DeferredCounterEdge)
    e.clip === :smoothed &&
        return "deferred counter on :$(e.species) uses a smoothed clip of width " *
               "$(e.smoothing), replacing the published model's max(0, ·) with a " *
               "differentiable approximation"
    e.clip === :unclamped &&
        return "deferred counter on :$(e.species) is unclamped — the pool may go " *
               "negative, where the published model floors it at zero and " *
               "carries the deficit forward"
    return nothing
end

function deviation_reason(e::RateConstantEdge)
    e.cadence === :continuous || return nothing
    return "rate-constant edge on :$(e.species) is continuous, where the " *
           "published model recomputes on a fixed interval and holds between refreshes"
end

function deviation_reason(e::ClampedEdge)
    e.origin === :ours || return nothing
    return "clamp on :$(e.species) is this reduction's simplification, not the " *
           "published model's"
end

deviation_reason(::CouplingEdge) = nothing

"""
    deviation_category(e::CouplingEdge) -> Symbol

The [`ReductionLabel`](@ref) category a deviating edge files under. Defined only
for the kinds that can deviate — asking for the category of an edge whose
[`deviation_reason`](@ref) is `nothing` is a caller error, and a new deviating
kind must add its method here rather than fall through to a wrong default.
"""
deviation_category(e::DeferredCounterEdge) =
    e.clip === :smoothed ? :smoothed_counter : :unclamped_counter
deviation_category(::RateConstantEdge) = :continuous_rebuild
deviation_category(::ClampedEdge) = :clamp

export CouplingEdge, MassEdge, CurrencyEdge, DeferredCounterEdge, CatalyticEdge,
       RateConstantEdge, VolumeEdge, ClampedEdge, EDGE_KINDS,
       edge_kind, is_consumer, is_producer, carries_mass, mass_contribution,
       obstructs_gradients, deviates_from_published, deviation_reason,
       deviation_category
