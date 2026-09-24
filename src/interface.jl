"""
    AbstractSubModel

Abstract supertype for every InferCell block (TX/TL, metabolism, stochastic gene
expression, …). Concrete sub-types implement the protocol functions
[`states`](@ref), [`parameters`](@ref), [`dynamics`](@ref) and/or
[`reactions`](@ref), and optionally override [`inputs`](@ref),
[`written_states`](@ref), [`formalism`](@ref), and [`inference_mode`](@ref).
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
    written_states(m::AbstractSubModel) -> Vector{Symbol}

States of *other* sub-models that a `:jump` module's reaction affects modify —
the jump-side twin of [`contributed_states`](@ref). Defaults to empty. Every
name must also appear in [`inputs`](@ref), since the write goes through the same
`u_inputs` view the read does; the orchestrator refuses a write to an input not
listed here, naming the module and the state, so removing the declaration makes
the write fail rather than silently land. Where the state is a Core A′ registry
species [`resolve_coupling`](@ref) additionally holds the declaration to a mass
or currency edge in either direction — `:in` where the module consumes the
pool, `:out` where it produces into it, and an outbound edge alone satisfies
the inputs contract for that species, since the pool is an input only as the
write channel. A non-registry state, such as a transcript, is gated by the
declaration alone (spec §12, amendment of 2026-09-04).
"""
written_states(::AbstractSubModel) = Symbol[]

"""
    coupling(m::AbstractSubModel) -> Vector{CouplingEdge}

Typed coupling declarations for `m`: not only *which* state crosses a boundary
but *how*, as one of the seven kinds in [`EDGE_KINDS`](@ref). Defaults to empty,
so a sub-model with no cross-block coupling — and every sub-model written before
the Core A′ interface contract — needs no modification.

[`inputs`](@ref) keeps its own meaning and its own default. A sub-model may
declare both; where it declares any coupling edge, [`resolve_coupling`](@ref)
checks the two against each other in both directions — every registry-species
input has an inbound edge (or, for a jump module producing into the pool, an
outbound mass or currency edge plus a [`written_states`](@ref) entry), and
every inbound mass or currency edge on a state the module does not integrate is
listed in `inputs` — so a typed declaration that drifts out of step with
`inputs` is caught rather than silently disagreeing. Inputs naming non-registry state (the legacy blocks' mRNA and
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
    reduction_notes(m::AbstractSubModel) -> Vector

Simplifications this sub-model introduces that the published model does not
make. Defaults to empty.

Each entry is a `String`, reported under `:model_note`, or a `category =>
description` pair naming one of [`REDUCTION_CATEGORIES`](@ref) — which is how
the lumped charging step is reported as `:lumping` and its deterministic
integration as `:formalism`, two separate declarations (spec §11 task 13.6).

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
    Reaction(rate, affect!)

One jump reaction of a `:jump` sub-model, written in the module's **own**
coordinates. `rate(u, p, t, u_inputs)` returns the propensity; `affect!(u,
u_inputs)` applies the state change by mutating `u`. Both receive the module's
local state slice `u`, its local parameters `p` and its resolved coupling inputs
`u_inputs`, in [`inputs`](@ref) order — exactly what [`dynamics`](@ref)
receives — never the composed global vectors.

The orchestrator wraps each `Reaction` in a `ConstantRateJump` over views of the
global state and parameter vectors, so a module cannot address a slot outside
its own slice: `u[1]` is this module's first state whether it is composed first
or seventh. A peer's state is read through `u_inputs` and may be written through
it only where [`written_states`](@ref) declares the write.
"""
struct Reaction{R, A}
    rate::R
    affect!::A
end

"""
    reactions(m::AbstractSubModel) -> Vector{Reaction}

Jump-process reactions for SSA sub-models, each a [`Reaction`](@ref) in local
coordinates. Required for `formalism = :jump`; defaults to an error to catch
incomplete sub-model definitions.

Before spec §11 phase 2 this returned `ConstantRateJump`s whose closures indexed
the global vectors at local positions, so a second jump module in a composition
wrote the first module's slice. `build_problem` now rejects that form by name.
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
    rebuilt_params(m::AbstractSubModel) -> Vector{Symbol}

Names of `m`'s **own free parameters** that the 60 s rebuild recomputes from
live ODE pools. Defaults to empty.

This is the declaration side of the rate-constant channel — the only
ODE→stochastic channel in the model, and so the one that decides whether the
coupling is bidirectional at all (spec §11 phase 4). It is the rate-constant
twin of [`contributed_states`](@ref): the names say *which* slots the hook
fills, and [`rate_constants`](@ref) says with what.

The values live in the composed **parameter vector**, not on the sub-model
struct, so `remake(prob; p = θ)` and [`model_free_params`](@ref) still see
them; a constant held on the struct would be invisible to both, and so could
never appear in a posterior, a provenance table or an identifiability Jacobian.
Each name must therefore be one of this module's own free parameters, and must
not collide with a free parameter of any other module in the composition —
parameters deduplicate by name into one slot, so a rebuild into a shared name
would drive the other module's rate law too. Both are refused at build time,
by name.

A module that declares names here must also declare at least one inbound
[`RateConstantEdge`](@ref), and the converse: an inbound rate-constant edge on
a module that rebuilds nothing is a channel declared and never executed.

**A rebuilt parameter is a derived quantity, and the inference layer does not
know that yet.** `build_turing_model` and the ABC path both sample every entry
of [`model_free_params`](@ref), a rebuilt slot included — and the hook
overwrites the drawn value at the first refresh, so it influences the
trajectory only over the first interval. At seventeen transcription constants
that is seventeen posterior dimensions costing sampler effort and coming back
shaped like their priors, which must not be read as an identifiability result.
Excluding derived names from the sampled set belongs to the inference phases
(spec §11 phases 15 to 17), which are the ones that will first compose a
hybrid driver with `infer`; until then this is a documented consequence of the
mechanism, not a defect in it.
"""
rebuilt_params(::AbstractSubModel) = Symbol[]

"""
    rate_constants(p, t, m::AbstractSubModel, pools)

The values the rebuild writes into the slots [`rebuilt_params`](@ref) names, in
that order, as a static vector of the same length. Defaults to an empty static
vector, so a module that rebuilds nothing needs no method.

`p` is the module's local parameter slice and `t` the rebuild time. `pools` are
the concentrations, **in mM**, of the species this module's inbound
[`RateConstantEdge`](@ref)s name, in declaration order. They arrive as an
argument rather than through [`inputs`](@ref) because a jump module may not
name an ODE-block state in `inputs()` — the two blocks hold separate state
vectors, and `inputs()` is resolved within a block.

Note the slot a rebuild fills is read from `p` on the next call, so a rebuild
law written as a function of its own previous value compounds. Write it as a
function of `pools` and of the module's *other* parameters.
"""
rate_constants(p, t, ::AbstractSubModel, pools) = SVector{0, Float64}()

"""
    membrane_protein_states(m::AbstractSubModel) -> Vector{Symbol}

States of `m` whose **counts** set the cell's membrane surface area, and so its
radius, its volume, and every count-to-concentration conversion. Defaults to
empty.

This is the declaration side of the volume channel (spec §11 phase 5), the one
coupling channel that runs on neither the 1 s nor the 60 s clock. The published
model computes `CellSA_Prot = 28.0 nm² × Σ(membrane protein counts)` every hook
and derives radius and volume from it
(`dev/notes/well-stirred-minimal-cell.md:178-196`); the driver does the same
over whatever a composition flags here.

The flag is **data on the module** rather than knowledge in the growth code, so
that a composition is *swept* for its membrane proteins rather than the driver
knowing which sub-model to ask. Every name must therefore be one of the
declaring module's own [`states`](@ref), and no two modules may flag the same
species — it would be counted into the area twice. Both are refused by name.

Either block may declare it. A `:jump` module's state is already a count; an
`:ode` module's is a concentration, and its count is `u × factor` at the
conversion factor in force at that hook, which is the published semantics —
the particle map is master at the hook instant. The count is used as a
continuous value rather than rounded, so the geometry injects no drift of its
own into check 0's round trip.

A module that declares names here must also declare an outbound
[`VolumeEdge`](@ref) on each of them, and the converse: an *outbound* volume
edge with nothing flagged is a channel declared and never executed. Both are
refused at build time. The converse holds for the outbound direction only —
an inbound edge fills a rate-law slot rather than the surface area, flags
nothing, and is paired against a `param_slot` instead.

**The chain runs only under the hybrid handshake driver**, since it converts
counts to concentrations and back and only a mixed composition has both sides.
A homogeneous `build_problem` resolves this declaration and executes nothing
from it — as it does for the deferred-counter and rate-constant edges, which
are equally cross-block — so a module smoke-tested
alone is not evidence that its volume channel works. The same is true of an
inbound edge: standalone, its slot keeps the module's declared value, which
should therefore be the registry's initial radius so that a standalone run is
the published rate law at the published geometry rather than a placeholder.
A catalytic edge is the exception: a build with no jump block refuses it,
because its slot must be free and would be sampled with nothing filling it
(spec §11 phase 11a).
Note that the driver writes `radius_from_area_nm(502831) = 200.03505` nm and
not a round 200.0, so the two agree to the four significant figures the source
states the radius in and differ beyond them.
"""
membrane_protein_states(::AbstractSubModel) = Symbol[]

"""
    membrane_protein_states(models::Vector{<:AbstractSubModel}) -> Vector{Symbol}

Every membrane-protein state in a composition, in composition order, found by
sweeping the models rather than by knowing which one declares them.

This is the query spec §11 tasks 5.1 and 7.6 are written against: a module that
did not declare the flag can find every contributing state without naming any
other module's type. Refuses, by name, a flag on a state the declaring module
does not own, and the same species flagged twice.
"""
function membrane_protein_states(models::Vector{<:AbstractSubModel})
    flagged = Symbol[]
    declarers = Dict{Symbol, Symbol}()
    for m in models
        id = module_id(m)
        own = states(m)
        for s in membrane_protein_states(m)
            s in own || throw(ArgumentError(
                "Module $id flags :$s as a membrane protein, but does not own " *
                "it — states($id) is $own. The flag is data on the module that " *
                "integrates the state, so that a composition can be swept for " *
                "its membrane proteins; flag it on the module that owns it"))
            haskey(declarers, s) && throw(ArgumentError(
                "Modules $(declarers[s]) and $id both flag :$s as a membrane " *
                "protein, so its count would enter the surface area twice. One " *
                "declaration per species"))
            declarers[s] = id
            push!(flagged, s)
        end
    end
    return flagged
end

"""
    extracellular_states(m::AbstractSubModel) -> Vector{Symbol}

States of `m` whose concentration is referred to a volume other than the
cell's, so that cell growth does not dilute them. Defaults to empty.

The growth chain rescales every ODE concentration when the volume changes,
holding the particle count fixed — which is right for everything inside the
cell and wrong for anything outside it. Core A′ integrates external lactate as
a dynamic state whose millimolar is a *medium* concentration, scaled by a
medium-to-cell volume ratio (spec §4 D6); diluting it by the cell's volume
ratio would destroy lactate that has already left, and carbon balance (check 2)
is the check that would fail. There is no extracellular marker in the registry
— external lactate is an ordinary dynamic species there — so the exemption is
declared by the module that integrates it.

Only an `:ode` module may declare it: a `:jump` module's states are counts, and
a count is not referred to a volume at all. Every name must be one of the
declaring module's own [`states`](@ref), and a state may not be both exempt and
flagged by [`membrane_protein_states`](@ref) — a membrane protein's count is
read *from* its concentration, so exempting it would make the area depend on
the volume it sets. All three are refused at build time, by name.
"""
extracellular_states(::AbstractSubModel) = Symbol[]

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
