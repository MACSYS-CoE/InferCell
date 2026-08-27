"""
Coupling resolution.

Takes a set of sub-models and turns their declared [`CouplingEdge`](@ref)s into
a validated [`CouplingGraph`](@ref), or fails with an error that names the
modules and species involved.

The pass is standalone on purpose. A module author working alone in a worktree
can call [`resolve_coupling`](@ref) on a single module and get the checks, and a
report can be produced without building an `ODEProblem` or running a simulation.

It is also deliberately lenient about sub-models that are not part of Core A′.
A module whose states are not registry species — the TX/TL and metabolism blocks
that predate this contract — declares no coupling edges, so resolution is a
no-op for it and composition behaves exactly as it did before.
"""

"""
    ResolvedEdge

One declared edge after resolution: which module declared it, which species it
touches, where that species sits in the registry, and which module it couples
to. `state_index` is always a registry position — an edge naming a species
outside the registry fails resolution before a `ResolvedEdge` exists.
"""
struct ResolvedEdge
    edge::CouplingEdge
    kind::Symbol
    species::Symbol
    declared_by::Symbol
    peer::Union{Symbol, Nothing}
    state_index::Int
end

"""
    DeadEnd

A mass-carrying species whose declared flow stops: consumed with no declared
producer (`missing_role = :producer`) — the structural form of the failure that
exhausted the adenylate pool after 2.3% of a cell cycle, and of the GMP that
transcription buries and nothing returns — or produced with no declared
consumer (`missing_role = :consumer`), mass accumulating with nothing drawing
it down.

`moiety` is the species' registry group, which is what makes stranded AMP and
stranded GMP the same class of error rather than two unrelated ones. `modules`
names the declarers on the side of the flow that exists.
"""
struct DeadEnd
    species::Symbol
    moiety::Symbol
    missing_role::Symbol   # :producer or :consumer
    modules::Vector{Symbol}
end

"""
    CouplingGraph

The resolved boundary of a composition.

- `edges` — every declared edge, resolved.
- `unowned_states` — registry dynamic states no module integrates. A report, not
  an error: a single module under test legitimately leaves most of them unowned.
- `dead_ends` — mass-carrying species consumed with no declared producer, or
  produced with no declared consumer.
- `chemostat_exemptions` — species whose costs were exempted from the dead-end
  check because the registry chemostats them. Recorded so that a later change
  making a pool live reinstates the check rather than inheriting the exemption
  silently.
- `gradient_obstructions` — edges that are non-differentiable at the boundary.
- `deviations` — edges that depart from the published model.
"""
struct CouplingGraph
    edges::Vector{ResolvedEdge}
    unowned_states::Vector{Symbol}
    dead_ends::Vector{DeadEnd}
    chemostat_exemptions::Vector{Symbol}
    gradient_obstructions::Vector{ResolvedEdge}
    deviations::Vector{ResolvedEdge}
end

"""
    resolve_coupling(models) -> CouplingGraph
    resolve_coupling(model)  -> CouplingGraph

Resolve and validate the declared coupling of a composition.

Throws when the boundary is inconsistent: an edge naming an unregistered
species or pool, an edge whose named peer is absent, two modules declaring the
same species and direction with different kinds, a module integrating a
chemostat, a state integrated twice, two clamps holding one species at
different values, a clamp holding a chemostat away from the registry's value,
or an `inputs` list that has drifted from the module's own inbound edges in
either direction.

Returns a graph reporting what is merely incomplete rather than wrong: unowned
states and dead ends — including a declared cost whose paying state is absent —
which a partial composition is expected to have.
"""
function resolve_coupling(models::Vector{<:AbstractSubModel})
    _check_state_ownership(models)
    _check_inputs_consistency(models)

    resolved = _resolve_edges(models)
    _check_kind_agreement(resolved)
    _check_clamp_agreement(resolved)

    owned = _owned_states(models)
    unowned = [s for s in dynamic_species() if !haskey(owned, s)]

    dead_ends, exemptions = _find_dead_ends(resolved)

    obstructions = [r for r in resolved if obstructs_gradients(r.edge)]
    deviations = [r for r in resolved if deviates_from_published(r.edge)]

    return CouplingGraph(resolved, unowned, dead_ends, exemptions,
                         obstructions, deviations)
end

resolve_coupling(model::AbstractSubModel) = resolve_coupling([model])

function _owned_states(models::Vector{<:AbstractSubModel})
    owned = Dict{Symbol, Symbol}()   # species => declaring module
    for m in models
        for s in states(m)
            is_registered(s) && (owned[s] = module_id(m))
        end
    end
    return owned
end

# A module may not integrate a species the registry holds constant, and no two
# modules may integrate the same one. Both are checked only for registry
# species, so a sub-model outside Core A′ is unaffected.
function _check_state_ownership(models::Vector{<:AbstractSubModel})
    seen = Dict{Symbol, Symbol}()
    for m in models
        name = module_id(m)
        for s in states(m)
            is_registered(s) || continue

            if is_chemostatted(s)
                throw(ArgumentError(
                    "Module $name declares dynamics for :$s, which the Core A′ " *
                    "registry holds at a fixed concentration. A chemostatted " *
                    "species is not integrated; remove it from states($name) or " *
                    "change its treatment in the registry"))
            end

            if haskey(seen, s)
                throw(ArgumentError(
                    "Dynamic state :$s is integrated by two modules, $(seen[s]) " *
                    "and $name. Exactly one module owns each registry state"))
            end
            seen[s] = name
        end
    end
    return nothing
end

# Where a module declares typed coupling at all, its untyped `inputs` must not
# drift from its inbound edges — in either direction. Modules that declare no
# coupling — every sub-model predating this contract — are skipped entirely.
function _check_inputs_consistency(models::Vector{<:AbstractSubModel})
    for m in models
        edges = coupling(m)
        isempty(edges) && continue

        inbound = Set(e.species for e in edges if is_consumer(e))
        for s in inputs(m)
            # Only registry species are governed by the typed contract. A hybrid
            # module may also read legacy non-registry state (mRNA, protein)
            # through inputs(), which no edge kind can name, so demanding an
            # edge for those would make typed coupling and the legacy channel
            # mutually exclusive.
            is_registered(s) || continue
            s in inbound || throw(ArgumentError(
                "Module $(module_id(m)) lists :$s in inputs() but declares no " *
                "inbound coupling edge for it. A typed declaration that has " *
                "drifted from inputs() is the disagreement this check exists to " *
                "catch; add the edge or drop the input"))
        end

        # The converse, for the one kind whose inbound side is wired through
        # `inputs`: an inbound MassEdge on a dynamic species the module does not
        # itself integrate reaches its dynamics only via inputs(), so leaving it
        # out means the declared coupling silently never arrives.
        own = Set(states(m))
        declared = Set(inputs(m))
        for e in edges
            e isa MassEdge || continue
            is_consumer(e) || continue
            is_registered(e.species) && is_dynamic(e.species) || continue
            e.species in own && continue
            e.species in declared || throw(ArgumentError(
                "Module $(module_id(m)) declares an inbound mass edge on " *
                ":$(e.species) but does not list it in inputs(). Only inputs() " *
                "wires a state into dynamics, so the declared coupling would " *
                "silently never arrive; add :$(e.species) to inputs() or drop " *
                "the edge"))
        end
    end
    return nothing
end

function _resolve_edges(models::Vector{<:AbstractSubModel})
    present = Set(module_id(m) for m in models)
    resolved = ResolvedEdge[]

    for m in models
        name = module_id(m)
        for e in coupling(m)
            kind = edge_kind(e)   # throws for a subtype outside the seven

            is_registered(e.species) || throw(ArgumentError(
                "Module $name declares a $kind edge on :$(e.species), which is " *
                "not a Core A′ registry species. See src/corea/registry.jl for " *
                "the $(length(COREA_SPECIES)) registered species"))

            if e.peer !== nothing && !(e.peer in present)
                throw(ArgumentError(
                    "Module $name declares a $kind edge on :$(e.species) to peer " *
                    "$(e.peer), which is not in this composition. The composition " *
                    "holds $(join(sort(collect(present)), ", ")). An edge with an " *
                    "unnamed peer resolves against whichever module owns the " *
                    "species, which is what a single-module composition under " *
                    "test should use"))
            end

            if e isa CurrencyEdge && !is_registered(e.pool)
                throw(ArgumentError(
                    "Module $name declares a currency edge on :$(e.species) " *
                    "routed via pool :$(e.pool), which is not a Core A′ registry " *
                    "species. The pool is the routing node; a typo here would " *
                    "route the currency nowhere"))
            end

            # A clamp on a chemostatted species must hold it at the registry's
            # value: two modules clamping one medium concentration differently
            # would each pass alone and silently disagree at the join.
            if e isa ClampedEdge && e.held_value !== nothing &&
               is_chemostatted(e.species)
                registry_value = held_value(e.species)
                if registry_value !== nothing && e.held_value != registry_value
                    throw(ArgumentError(
                        "Module $name clamps :$(e.species) at $(e.held_value) mM, " *
                        "but the registry chemostats it at $registry_value mM. " *
                        "One held value per species; change the edge or the " *
                        "registry, not one of two copies"))
                end
            end

            push!(resolved, ResolvedEdge(e, kind, e.species, name, e.peer,
                                         species_index(e.species)))
        end
    end
    return resolved
end

# Two modules may not describe the same crossing differently. Mass and a
# deferred counter are not two spellings of one thing: one is continuous shared
# state, the other is a cost accrued in one block and debited in another a step
# later.
function _check_kind_agreement(resolved::Vector{ResolvedEdge})
    seen = Dict{Tuple{Symbol, Symbol}, ResolvedEdge}()
    for r in resolved
        key = (r.species, r.edge.direction)
        prior = get(seen, key, nothing)
        if prior === nothing
            seen[key] = r
        elseif prior.kind !== r.kind
            who = prior.declared_by === r.declared_by ?
                "Module $(prior.declared_by) describes" :
                "Modules $(prior.declared_by) and $(r.declared_by) disagree about"
            throw(ArgumentError(
                "$who how :$(r.species) crosses the boundary in direction " *
                ":$(r.edge.direction) in two ways. $(prior.declared_by) declares it " *
                ":$(prior.kind) — $(_kind_gloss(prior.kind)); $(r.declared_by) " *
                "declares it :$(r.kind) — $(_kind_gloss(r.kind)). These are " *
                "different semantics, not two spellings of one"))
        end
    end
    return nothing
end

# The registry check in _resolve_edges covers a clamp against a recorded
# chemostat value, but two modules can clamp a species the registry records no
# value for. Their held values must still agree: two clamps that each pass
# alone and disagree at the join would let composition order decide a held
# concentration.
function _check_clamp_agreement(resolved::Vector{ResolvedEdge})
    seen = Dict{Symbol, ResolvedEdge}()
    for r in resolved
        r.edge isa ClampedEdge || continue
        r.edge.held_value === nothing && continue
        prior = get(seen, r.species, nothing)
        if prior === nothing
            seen[r.species] = r
        elseif prior.edge.held_value != r.edge.held_value
            who = prior.declared_by === r.declared_by ?
                "Module $(prior.declared_by) clamps" :
                "Modules $(prior.declared_by) and $(r.declared_by) clamp"
            throw(ArgumentError(
                "$who :$(r.species) at two different values, " *
                "$(prior.edge.held_value) mM and $(r.edge.held_value) mM. One " *
                "held value per species; composition order must not decide " *
                "which one wins"))
        end
    end
    return nothing
end

function _kind_gloss(kind::Symbol)
    kind === :mass && return "shared state, continuous, gradients cross inside the ODE block"
    kind === :currency && return "routed via a shared pool node"
    kind === :deferred_counter && return "a cost accrued in one block and debited a step later"
    kind === :catalytic && return "counts entering a rate law as parameters, no mass flowing"
    kind === :rate_constant && return "pools re-entering the stochastic block as rate constants"
    kind === :volume && return "counts setting surface area, hence volume, hence every concentration"
    kind === :clamped && return "a real dependence replaced by a constant"
    return "unknown kind"
end

# Every declared cost needs a declared producer that returns it — and the
# mirror image, a produced species needs a consumer drawing it down. Both are
# reports rather than errors: a partial composition legitimately lacks the
# counterpart, and a module author validating a single module alone imports
# species whose owners and producers are absent by construction. The ownership
# gap itself is reported through `unowned_states`.
function _find_dead_ends(resolved::Vector{ResolvedEdge})
    mass_edges = [r for r in resolved if carries_mass(r.edge)]

    producers = Dict{Symbol, Vector{Symbol}}()
    consumers = Dict{Symbol, Vector{Symbol}}()
    for r in mass_edges
        bucket = is_producer(r.edge) ? producers : consumers
        push!(get!(bucket, r.species, Symbol[]), r.declared_by)
    end

    exemptions = Symbol[]
    dead_ends = DeadEnd[]

    for (species, consuming) in consumers
        if is_chemostatted(species)
            # The chemostat is the return path. Record the exemption so that a
            # later change making the pool live reinstates the check.
            push!(exemptions, species)
            continue
        end

        if !haskey(producers, species)
            push!(dead_ends,
                  DeadEnd(species, species_group(species), :producer, sort(consuming)))
        end
    end

    for (species, producing) in producers
        is_chemostatted(species) && continue   # the chemostat absorbs it
        if !haskey(consumers, species)
            push!(dead_ends,
                  DeadEnd(species, species_group(species), :consumer, sort(producing)))
        end
    end

    return sort!(dead_ends; by = d -> (d.species, d.missing_role)),
           sort!(unique!(exemptions))
end

"""
    dead_end_report(graph::CouplingGraph) -> String

The dead ends of a resolved graph, written for a human. Names the conserved
moiety of each stranded species, because that is what identifies AMP with no
route back and stranded GMP as one class of error rather than two.
"""
function dead_end_report(graph::CouplingGraph)
    isempty(graph.dead_ends) &&
        return "No dead ends: every consumed species has a declared producer, " *
               "and every produced species a declared consumer."
    lines = ["$(length(graph.dead_ends)) species whose declared flow stops:"]
    for d in graph.dead_ends
        push!(lines, d.missing_role === :producer ?
              "  :$(d.species) ($(d.moiety) moiety) consumed by " *
              "$(join(d.modules, ", ")) with nothing returning it" :
              "  :$(d.species) ($(d.moiety) moiety) produced by " *
              "$(join(d.modules, ", ")) with nothing drawing it down")
    end
    return join(lines, "\n")
end

"""
    gradient_report(graph::CouplingGraph) -> String

Edges that obstruct gradients, for a composition about to be sampled with a
gradient-based method. A diagnosable warning rather than a silent success, so
that proceeding with a clamped interface is a deliberate choice.
"""
function gradient_report(graph::CouplingGraph)
    isempty(graph.gradient_obstructions) &&
        return "No gradient obstructions: no clamped deferred counter on the boundary."
    lines = ["$(length(graph.gradient_obstructions)) edge(s) non-differentiable at the boundary:"]
    for r in graph.gradient_obstructions
        push!(lines, "  :$(r.species) in $(r.declared_by) clips at zero and carries " *
                     "the deficit forward — a max(0, ·) that will obstruct a " *
                     "gradient-based sampler")
    end
    return join(lines, "\n")
end

"""
    check_gradient_safety(models) -> Vector{ResolvedEdge}

Report the non-differentiable boundary edges of a composition whose sub-models
declare `inference_mode = :differentiable`. Returns the obstructing edges and
emits a warning naming them; it does not throw, because a clamped interface is
a legitimate model, just not a differentiable one.
"""
function check_gradient_safety(models::Vector{<:AbstractSubModel})
    graph = resolve_coupling(models)
    isempty(graph.gradient_obstructions) && return ResolvedEdge[]

    any(m -> inference_mode(m) === :differentiable, models) || return graph.gradient_obstructions

    @warn """
    A composition with a differentiable sub-model contains a non-differentiable boundary edge.
    $(gradient_report(graph))
    Declare clip = :smoothed (with its smoothing width) on the edge to sample it with gradients, \
    and label the result as departing from the published model.
    """
    return graph.gradient_obstructions
end

export ResolvedEdge, DeadEnd, CouplingGraph, resolve_coupling,
       dead_end_report, gradient_report, check_gradient_safety
