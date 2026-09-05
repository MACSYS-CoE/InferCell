"""
The 1 s handshake: executing the coupling the edge kinds declare.

Spec §11 phases 3 and 4. Phase 1 made mass and currency edges execute inside
the ODE block; phase 3 added the catalytic and deferred-counter channels and
phase 4 the rate-constant one, so five of the seven [`CouplingEdge`](@ref)
kinds now run and two — volume and clamped — remain declare-and-validate only. It
advances a jump block and an ODE block by a split-operator exchange, mirroring
the published model's `hookSimulation` (`dev/notes/well-stirred-minimal-cell.md`):
every `delt = 1.0` s, counts become initial conditions, the metabolic block
integrates one second, and the accrued costs are debited back.

**Mechanism, chosen in task 3.1: an outer loop over two stepped integrators.**
The two rejected alternatives, and why, are recorded in the pull request. The
short version is that a `JumpProblem` over an `ODEProblem` — the SciML-native
hybrid — makes propensities track the ODE pools *continuously*, which is a
different model from the published one and is not representable alongside the
60 s piecewise-constant rebuild phase 4 adds. An outer loop is the only one of
the three under which "the propensities are bitwise unchanged between
refreshes" is even a statement that can be true.

The two blocks keep **separate state vectors** — the ODE block in mM over
`Float64`, the jump block in particles over `Int` — which is what makes the
conversion below a real boundary rather than a units comment.
"""

# ---------------------------------------------------------------------------
# Counts and concentrations (spec §11 task 3.3)
# ---------------------------------------------------------------------------

"""
Avogadro's constant, mol⁻¹. Exact by the 2019 SI definition, so this is a
definition rather than a measurement and carries no uncertainty.
"""
const AVOGADRO = 6.02214076e23

"""
The published initial cell radius, in nm.

`dev/notes/reduced-syn3a-scoping.md` records the derived figure this must
reproduce: at 200 nm, 1 mM is 20,180 particles. Phase 5 makes the radius grow;
until then it is the one fixed number the conversion closes over, and it is
passed as an argument rather than baked in so that phase 5 changes a call site
and not an equation.
"""
const COREA_INITIAL_RADIUS_NM = 200.0

"""
    cell_volume_litres(radius_nm) -> Float64

Volume of a sphere of the given radius, in litres. The published model's growth
law is surface-area accounting on a sphere, so a radius is the only geometric
quantity the conversion needs.
"""
cell_volume_litres(radius_nm) = (4 / 3) * π * (radius_nm * 1e-9)^3 * 1e3

"""
    particles_per_mM(volume_litres) -> Float64

How many particles one millimolar comes to at this volume. Derived from
[`AVOGADRO`](@ref) rather than hard-coded, so that the 20,180 of the scoping
note is a *reproduced* number and not a transcribed one — a test asserts the
derivation returns it at the published radius.
"""
particles_per_mM(volume_litres) = 1e-3 * volume_litres * AVOGADRO

"""
    corea_particles_per_mM(radius_nm = COREA_INITIAL_RADIUS_NM) -> Float64

The conversion factor at a given radius. 20,180 particles per mM at 200 nm.
"""
corea_particles_per_mM(radius_nm = COREA_INITIAL_RADIUS_NM) =
    particles_per_mM(cell_volume_litres(radius_nm))

"""
    counts_to_mM(n, factor) -> Float64

Particles to millimolar. Exact division: no policy applies in this direction,
because the target is continuous.
"""
counts_to_mM(n, factor) = n / factor

"""
Rounding policies for writing a continuous quantity back as whole particles.

Spec §3, check 0. Across the ~6,300 handshakes of a full cycle the three
policies have three different *signatures*, and the check is on the signature
rather than on any single magnitude:

- `:fractional_carry` — the remainder is carried to the next handshake, so the
  quantisation injected at each hook **telescopes**: the emitted counts track
  the exact running total and what is left is one carry, bounded by half a
  particle and not growing with the handshake count. The default, and the only
  policy under which the tolerance principle of §3 tests the integrator rather
  than the rounding. Note the telescoping is a statement about the injected
  perturbations, not about the residual of the trajectory they are injected
  into: a nonlinear right-hand side does not preserve it, which is why task
  3.4 asserts exactness on a frozen pool.
- `:stochastic` — round up with probability equal to the fraction. Unbiased, so
  the residual accumulates as the **square root** of the handshake count, about
  0.002 mM over a cycle.
- `:deterministic` — round to nearest, discarding the remainder. Biased, so the
  residual accumulates **linearly**: up to 0.5 particle × 6,300 = 3,150
  particles, about 0.156 mM per species, roughly 140× the integrator bound.
  **Rejected by the spec.** It is implemented only so that check 0 can
  demonstrate the signature that rejects it; selecting it is a labelled
  departure from this project's own policy.
"""
const ROUNDING_POLICIES = (:fractional_carry, :stochastic, :deterministic)

"""
    RoundingState(policy; nspecies)

The carried remainders, one per species, and the policy that governs them.

The remainder vector is live state, not bookkeeping: under `:fractional_carry`
it is exactly what makes the round trip lossless, and it is the reason the
policy is a field on the driver rather than a call-site argument. The other two
policies carry no state and leave it at zero, which is itself worth being able
to assert.
"""
struct RoundingState
    policy::Symbol
    remainders::Vector{Float64}

    function RoundingState(policy::Symbol, remainders::Vector{Float64})
        _check_vocab(:RoundingState, :policy, policy, ROUNDING_POLICIES)
        return new(policy, remainders)
    end
end

RoundingState(policy::Symbol = :fractional_carry; nspecies::Integer = 0) =
    RoundingState(policy, zeros(Float64, nspecies))

"""
    round_to_counts!(state, i, x) -> Int

Emit whole particles for the continuous quantity `x` at slot `i`, under the
state's policy, updating its carried remainder.

Fractional carry rounds the *carry-adjusted* value to nearest rather than
flooring it. Flooring is the textbook formulation and is wrong here: a quantity
that is already an exact integer arrives one ulp low after a division and a
multiplication, so `floor` would emit `n - 1` and leave a carry of nearly one,
making the counts oscillate around the right answer instead of sitting on it.
Rounding to nearest keeps the carry in [-0.5, 0.5] and returns `n` exactly,
which is what task 3.4 asserts over 600 handshakes.
"""
function round_to_counts!(state::RoundingState, i::Integer, x::Real)
    if state.policy === :fractional_carry
        total = x + state.remainders[i]
        n = round(Int, total)
        state.remainders[i] = total - n
        return n
    elseif state.policy === :stochastic
        base = floor(x)
        n = Int(base) + (rand() < (x - base) ? 1 : 0)
        return n
    else  # :deterministic — biased, and rejected; kept so check 0 can show why
        return round(Int, x)
    end
end

export AVOGADRO, COREA_INITIAL_RADIUS_NM, cell_volume_litres, particles_per_mM,
       corea_particles_per_mM, counts_to_mM, ROUNDING_POLICIES, RoundingState,
       round_to_counts!

# ---------------------------------------------------------------------------
# Lowered exchange records (spec §11 tasks 3.5 and 3.6)
#
# The resolved graph is a description; these are the executable form of it.
# They are built once, at build time, with every index already resolved, so
# that the per-handshake loop touches no symbol lookup and no edge object. The
# wave-0 review left this as an advisory ("lower ResolvedEdges into concrete
# typed callback state before any hot path iterates them"); at 6,300 handshakes
# per trajectory it is what keeps task 3.8's measurement honest.
# ---------------------------------------------------------------------------

"""
    CatalyticExchange

One [`CatalyticEdge`](@ref) made executable: a protein count in the jump block
fills a named rate-law parameter of an ODE module, converted to mM. No mass
moves — that is the edge kind's defining property — so this writes a parameter
and never a state.
"""
struct CatalyticExchange
    count_idx::Int        # the protein's slot in the jump block's state vector
    param_idx::Int        # the slot it fills in the ODE block's parameter vector
    species::Symbol
    param_slot::Symbol
    declared_by::Symbol
end

"""
    DeferredDebit

One [`DeferredCounterEdge`](@ref) made executable. The stochastic block accrues
a cost in whole particles; the hook debits it against an ODE pool one step
later and, under the published `:clamped_deficit_carried` policy, floors the
pool at zero and carries the shortfall to the next handshake.

`deficit` is the carried shortfall, in particles, and is the interface's
`max(0, ·)` in its stateful form — the quantity check 7's census counts and the
one [`obstructs_gradients`](@ref) is about.
"""
mutable struct DeferredDebit
    counter_idx::Int      # the accrual counter's slot in the jump state vector
    pool_idx::Int         # the debited pool's slot in the ODE state vector
    sign::Int             # -1 where the module consumes the pool, +1 where it produces
    clip::Symbol
    smoothing::Union{Float64, Nothing}
    deficit::Float64
    species::Symbol
    counter::Symbol
    declared_by::Symbol
end

# One counter may feed several pools — the published charged-tRNA transfer is
# exactly that shape, debiting the charged pool and crediting the uncharged one
# from a single accrual (spec §3). So the counter is read and cleared *once* per
# handshake, before any pool is touched, and every debit on it sees the same
# accrued value. Reading it inside the per-pool loop would let the first debit
# zero it and hand every later one a silent zero.
struct _CounterRead
    counter_idx::Int
    debit_idxs::Vector{Int}   # positions in `driver.debits` fed by this counter
end

"""
    RateConstantRebuild

One module's [`RateConstantEdge`](@ref) declarations made executable: the 60 s
rebuild, and the only ODE→stochastic channel in the model.

Every `steps` handshakes the hook reads the pools the module's inbound
rate-constant edges name, calls its [`rate_constants`](@ref), and writes the
results into the jump block's **parameter vector** at the slots
[`rebuilt_params`](@ref) names. Between refreshes nothing touches them, which
is what makes the coupling the published piecewise-constant one rather than an
accidental continuous one.

The record is per module rather than per edge: one pool feeds many constants —
the nucleotide pools set all seventeen transcription rate constants — so the
edges name the inputs and the module names the outputs.
"""
mutable struct RateConstantRebuild
    model::AbstractSubModel
    p_idxs::Vector{Int}       # global jump-parameter slots of the module's own free params
    fill_idxs::Vector{Int}    # the slots the rebuild fills, in rebuilt_params order
    pool_idxs::Vector{Int}    # global ODE state indices of the pools, in edge order
    steps::Int                # handshakes between refreshes
    interval::Float64
    names::Vector{Symbol}     # rebuilt_params(m), for reports and refusals
    species::Vector{Symbol}   # the pools, in the same order as `pool_idxs`
    declared_by::Symbol
    n_refreshes::Int
end

# ---------------------------------------------------------------------------
# The driver
# ---------------------------------------------------------------------------

"""
    HandshakeDriver

A mixed ODE/jump composition, advanced by a split-operator exchange.

Holds two live integrators and the lowered exchange records between them.
`interval` is the exchange period — 1.0 s, the published model's `delt` — and
`drain_interval` how often the deferred counters are debited, which defaults to
`interval` and is a labelled reduction when it is coarser (spec §4 D10). The
60 s rate-constant rebuild runs on each [`RateConstantRebuild`](@ref)'s own
declared schedule. The volume chain (phase 5) hooks into the same loop and is
deliberately absent here.

`n_handshakes`, `n_drains` and `n_clipped` are check 7's census: how many
exchanges ran, how many of them applied a debit, and at how many of *those* any
deferred counter carried a deficit. K5 is scored from the last two, and the
denominator is the drain count rather than the handshake count — see
[`clipping_census`](@ref).
"""
mutable struct HandshakeDriver{OI, JI}
    ode::OI
    jump::JI
    interval::Float64
    drain_interval::Float64               # how often the deferred counters are debited
    steps_per_drain::Int                  # that interval, in handshakes
    factor::Float64                       # particles per mM, at fixed volume
    rounding::RoundingState
    catalytic::Vector{CatalyticExchange}
    debits::Vector{DeferredDebit}
    counters::Vector{_CounterRead}        # debits grouped by the counter feeding them
    rebuilds::Vector{RateConstantRebuild} # the 60 s rate-constant channel
    t_end::Float64                        # the declared end of both blocks' tspan
    n_handshakes::Int
    n_drains::Int
    n_clipped::Int
end

"""
    clipping_census(driver) -> NamedTuple

Check 7 on this driver: the number of handshakes run, the number of them that
applied a debit, the number of *those* at which any deferred counter carried a
deficit, the fraction, the carried deficits and the accrual still sitting in
the counters. Spec §3 check 7 requires **zero at published parameters**; more
than 5% of prior draws clipping puts the non-smooth drain in the operating
regime and fires K5.

**`fraction` is per drain, not per handshake.** Only a drain can clip, so at a
`drain_interval` coarser than the exchange a per-handshake denominator would
dilute the fraction by exactly `steps_per_drain` — at 60 s, a run in which
*every* debit clipped would report 1.7% and pass K5's 5% gate. `pending` exists
for the same reason: between drains the accrual sits in the jump block's
counters rather than in any `deficit`, so without it the ledger reads closed
while a period of cost is unpaid.

Note that a coarse drain changes what this census is *about*, in both
directions: each debit is `steps_per_drain` times larger against the same pool,
so it clips more readily, while there are that many fewer of them. Check 7 and
K5 are scored at the published 1 s drain (spec §12, 2026-09-05).
"""
clipping_census(d::HandshakeDriver) = (handshakes = d.n_handshakes,
                                       drains = d.n_drains,
                                       clipped = d.n_clipped,
                                       fraction = d.n_drains == 0 ? 0.0 :
                                                  d.n_clipped / d.n_drains,
                                       deficits = [b.deficit for b in d.debits],
                                       pending = [float(d.jump.u[c.counter_idx])
                                                  for c in d.counters])

"""
    rebuild_census(driver) -> Vector{NamedTuple}

How many times each rate-constant rebuild has fired, and on what schedule. One
row per rebuilding module: its id, the parameters it rebuilds, the pools it
reads, its declared interval and the refresh count.
"""
rebuild_census(d::HandshakeDriver) =
    [(module_id = r.declared_by, params = copy(r.names), pools = copy(r.species),
      interval = r.interval, refreshes = r.n_refreshes) for r in d.rebuilds]

export CatalyticExchange, DeferredDebit, RateConstantRebuild, HandshakeDriver,
       clipping_census, rebuild_census

# ---------------------------------------------------------------------------
# Building a hybrid composition (spec §11 task 3.2)
# ---------------------------------------------------------------------------

# A pool of Float64s as a static vector. `_svec` in the orchestrator does the
# same for indices; this is its floating-point twin, and it exists because the
# same construction appears at four sites, three of them in the per-handshake
# path.
_fvec(v) = SVector{length(v), Float64}(Tuple(v))

# A hook fires only at a handshake — the one instant at which the two blocks
# agree on a state — so a period between handshakes is unreachable, and rounding
# one silently would run a cadence nobody declared. Written once because the
# drain interval and the rebuild interval are held to the same rule, and a
# silent divergence between two copies of it is exactly the failure the rule
# exists to prevent.
function _handshake_steps(period, interval, subject)
    n = round(Int, period / interval)
    (n >= 1 && isapprox(n * interval, period; rtol = 1e-12)) || throw(ArgumentError(
        "$subject of $period s is not a whole number of $(interval) s handshakes. " *
        "It is applied at a handshake — the only instant at which the two blocks " *
        "agree on a state — so a period between handshakes is unreachable"))
    return n
end

# Where a species sits in one block's composed state vector, or `nothing`.
# Blocks are built by the existing homogeneous paths, so the layout is
# `_build_contexts`': models in order, each contributing its `states` in order.
function _block_state_index(models::Vector{<:AbstractSubModel}, s::Symbol)
    i = 0
    for m in models, st in states(m)
        i += 1
        st === s && return i
    end
    return nothing
end

# The ODE block's parameter layout, taken from the builder's own contexts
# rather than re-derived. `_build_p0` and `_build_contexts` already agree on the
# first-seen-by-name dedup; a third implementation of that rule here would be one
# more thing to keep in step, and the failure if it drifted would be a silent
# write to the wrong rate-law slot.

# Lower every edge that crosses the formalism boundary into an executable
# record, with every index resolved at build time. An edge whose two ends are in
# one block is that block's own business — phase 1's contribution channel or
# phase 2's peer writes — and is deliberately not the driver's.
function _lower_exchanges(models, ode_models, jump_models, ode_contexts)
    catalytic = CatalyticExchange[]
    debits = DeferredDebit[]

    # A rate-law slot must belong to exactly one ODE module. Two modules
    # declaring a free parameter of the same name share one slot in the composed
    # vector (`_build_p0` dedups by name), so a catalytic write into it would
    # silently drive the other module's rate law too.
    # Scanned over *both* blocks, for the reason `_lower_rebuilds` scans both:
    # within the ODE block a shared name is one deduplicated slot, so the count
    # would drive the other module's rate law; across the boundary the two
    # blocks hold separate parameter vectors, so the write would land in one and
    # not the other and a single named parameter would carry two values —
    # invisible to `_validate_shared_params`, which compares declared values and
    # finds them equal.
    slot_owners = Dict{Symbol, Vector{Symbol}}()
    for m in vcat(ode_models, jump_models), q in model_free_params(parameters(m))
        push!(get!(slot_owners, q.name, Symbol[]), module_id(m))
    end

    for (i, m) in enumerate(ode_models)
        id = module_id(m)
        own = model_free_params(parameters(m))
        for e in coupling(m)
            e isa CatalyticEdge || continue
            # A catalytic channel is lowered from the ODE side only. The jump
            # module owning the count may name the same channel from its own
            # side, but only the consuming module has a `param_slot` — the slot
            # is a position in *its* rate law — so the peer declaration is the
            # same channel described twice and is not lowered again.
            count_idx = _block_state_index(jump_models, e.species)
            count_idx === nothing && throw(ArgumentError(
                "Module $id declares a CatalyticEdge on :$(e.species), but no " *
                "jump module in this composition owns that state. A catalytic " *
                "edge across the boundary reads a count from the stochastic " *
                "block; compose the module that owns :$(e.species), or drop " *
                "the edge"))
            j = findfirst(q -> q.name === e.param_slot, own)
            j === nothing && throw(ArgumentError(
                "Module $id declares a CatalyticEdge filling the parameter " *
                "slot :$(e.param_slot), which is not one of its own free " *
                "parameters. The slot names the rate-law parameter the count " *
                "fills, so it must exist on the module whose rate law it " *
                "belongs to; its free parameters are $([q.name for q in own])"))
            holders = slot_owners[e.param_slot]
            length(holders) == 1 || throw(ArgumentError(
                "Module $id fills the parameter slot :$(e.param_slot) from a " *
                "catalytic edge, but $(join(holders, " and ")) each declare a " *
                "free parameter of that name. Within the metabolic block those " *
                "are one deduplicated slot, so the count would drive every one " *
                "of those rate laws; across the boundary they are two slots in " *
                "two parameter vectors, so the write would land in one and leave " *
                "the other. Give the slot a name unique to $id"))
            push!(catalytic, CatalyticExchange(count_idx, ode_contexts[i].param_idxs[j],
                                               e.species, e.param_slot, id))
        end
    end

    # Deferred counters are lowered from the side that owns the counter — the
    # block that accrues the cost. That side's `direction` is the one that means
    # something: `:in` where it draws the pool down, `:out` where it produces
    # into it. The pool's owner may name the same channel from its own side with
    # the opposite direction; taking the sign from whichever declaration came
    # first would let composition order decide whether a cost is debited or
    # credited.
    seen = Set{Tuple{Symbol, Symbol}}()
    for m in models
        id = module_id(m)
        for e in coupling(m)
            e isa DeferredCounterEdge || continue
            # Not a declaration this module can speak for: it does not own the
            # counter, so it is either the pool owner's mirror of a channel the
            # accruing side lowers, or an orphan. The loop below decides which.
            _block_state_index([m], e.counter) === nothing && continue
            counter_idx = _block_state_index(jump_models, e.counter)
            counter_idx === nothing && throw(ArgumentError(
                "Module $id accrues into :$(e.counter), which it owns itself but " *
                "which is not a jump-block state. A deferred counter is the " *
                "stochastic block's accrual slot, cleared by the hook every " *
                "handshake; an ODE module's own state cannot be one"))
            is_registered(e.counter) && throw(ArgumentError(
                "Module $id accrues into :$(e.counter), which is a Core A′ " *
                "registry species. A deferred counter is an accrual slot the " *
                "hook clears every handshake, not a modelled pool; clearing a " *
                "registry species would destroy it once per second"))
            pool_idx = _block_state_index(ode_models, e.species)
            pool_idx === nothing && throw(ArgumentError(
                "Module $id declares a DeferredCounterEdge debiting " *
                ":$(e.species), but no ODE module in this composition " *
                "integrates that pool. The hook debits a continuous pool; " *
                "compose the module that owns :$(e.species)"))
            # One channel is one debit however many times it is declared.
            key = (e.counter, e.species)
            key in seen && continue
            push!(seen, key)
            push!(debits, DeferredDebit(counter_idx,
                                        pool_idx, mass_contribution(e), e.clip,
                                        e.smoothing, 0.0, e.species, e.counter, id))
        end
    end

    # Every declared channel must have been lowered by the side that owns its
    # counter. Two ways it might not have been, and both are silent without this:
    # no module owns the counter at all, or the only module that declared the
    # channel was the pool's owner, whose declaration this driver does not act
    # on. A cost that is declared and never debited is the failure the deferred
    # counter exists to prevent, and it leaves a clean census behind it.
    for m in models, e in coupling(m)
        e isa DeferredCounterEdge || continue
        (e.counter, e.species) in seen && continue
        owner = _block_state_index(jump_models, e.counter) === nothing ? nothing :
                first(md for md in jump_models
                      if _block_state_index([md], e.counter) !== nothing)
        owner === nothing && throw(ArgumentError(
            "Module $(module_id(m)) declares a DeferredCounterEdge accruing " *
            "into :$(e.counter), but no jump module in this composition owns " *
            "that state. The counter is the stochastic block's accrual slot; " *
            "it must be a state of a :jump module"))
        throw(ArgumentError(
            "Module $(module_id(m)) declares the deferred channel " *
            ":$(e.counter) → :$(e.species), but $(module_id(owner)), which owns " *
            "the counter, declares no matching edge. The channel is lowered " *
            "from the accruing side, so this one would never fire: the counter " *
            "would grow without bound, the pool would never be debited, and the " *
            "clipping census would report nothing wrong. Declare it on " *
            "$(module_id(owner)) as well"))
    end

    # The same trap on the catalytic side: an edge declared only by the jump
    # module that owns the count names no rate law, so nothing is written and
    # the ODE module silently keeps its nominal enzyme concentration.
    lowered_catalytic = Set{Symbol}(c.species for c in catalytic)
    for m in jump_models, e in coupling(m)
        e isa CatalyticEdge || continue
        e.species in lowered_catalytic && continue
        throw(ArgumentError(
            "Module $(module_id(m)) declares a CatalyticEdge on :$(e.species), " *
            "which it owns, but no ODE module in this composition declares the " *
            "matching channel. `param_slot` names a position in the *consuming* " *
            "module's rate law, so the channel is lowered from the ODE side; " *
            "declared only from here it would never write anything"))
    end

    # One accrual can be split among several consumers — they share what the
    # counter holds — but two producers would each be credited the whole of it,
    # which creates matter. Spec §3's only shape is one of each; refuse the rest
    # rather than leave a silent doubling for a later phase to meet.
    for c in unique(b.counter for b in debits)
        producers = [b for b in debits if b.counter === c && b.sign > 0]
        length(producers) <= 1 || throw(ArgumentError(
            "Counter :$c credits more than one pool — " *
            "$(join((string(":", b.species) for b in producers), ", ")) — and each " *
            "would receive the whole accrual, creating matter from one cost. " *
            "Split the accrual across separate counters, one per credited pool"))
    end

    # Group the debits by counter, so the hook reads and clears each counter
    # once and every pool it feeds sees the same accrued value.
    groups = _CounterRead[]
    for (k, b) in enumerate(debits)
        g = findfirst(c -> c.counter_idx == b.counter_idx, groups)
        g === nothing ? push!(groups, _CounterRead(b.counter_idx, [k])) :
                        push!(groups[g].debit_idxs, k)
    end
    return catalytic, debits, groups
end

# Lower the rate-constant channel (spec §11 phase 4). One record per jump module
# that rebuilds anything, with every index resolved at build time.
#
# The channel is lowered from the module whose constants are rebuilt — the side
# that declares `:in`, by the direction convention wave 0 pinned. Unlike a
# deferred counter, no mass moves, so the pool's owner has nothing to execute
# and is not required to mirror the declaration; an `:out` edge with no `:in`
# counterpart is refused below, because that one *is* a channel declared and
# never run.
function _lower_rebuilds(ode_models, jump_models, jump_contexts, interval)
    rebuilds = RateConstantRebuild[]

    # A rebuilt slot must belong to exactly one jump module, for the reason a
    # catalytic `param_slot` must: `_build_p0` deduplicates free parameters by
    # name into one slot, so rebuilding a shared name would drive the other
    # module's propensities too. With seventeen genes and a handful of shared
    # gene-expression globals, reusing a name is the obvious mistake.
    # Scanned over *both* blocks. Within the jump block a shared name is one
    # slot, so the rebuild would drive the other module's propensities. Across
    # the boundary the two blocks hold separate parameter vectors, so the hook
    # would rewrite only the jump block's copy and the two would silently carry
    # different values for one named parameter — which `_validate_shared_params`
    # cannot catch, because it compares *declared* values and they agree.
    slot_owners = Dict{Symbol, Vector{Symbol}}()
    for m in vcat(jump_models, ode_models), q in model_free_params(parameters(m))
        push!(get!(slot_owners, q.name, Symbol[]), module_id(m))
    end

    for m in ode_models
        isempty(rebuilt_params(m)) || throw(ArgumentError(
            "Module $(module_id(m)) is an :ode sub-model and declares " *
            "rebuilt_params. The rate-constant channel rebuilds the *stochastic* " *
            "block's constants from live pools; an ODE module reads a pool it " *
            "does not own directly, through inputs() and a mass or currency edge"))
    end

    for (i, m) in enumerate(jump_models)
        id = module_id(m)
        names = rebuilt_params(m)
        edges = [e for e in coupling(m) if e isa RateConstantEdge && is_consumer(e)]

        # The direction convention wave 0 pinned: the module whose constants are
        # rebuilt declares `:in`, the pool's owner `:out`. So an outbound edge on
        # a jump module is the convention read backwards, and without this it is
        # the *quiet* mistake — the inbound form throws for a missing
        # `rebuilt_params`, while the outbound one falls through the filter above
        # and composes with the stochastic block silently keeping its nominal
        # constants for the whole trajectory.
        for e in coupling(m)
            e isa RateConstantEdge && is_producer(e) || continue
            throw(ArgumentError(
                "Module $id is a :jump module and declares an *outbound* " *
                "RateConstantEdge on :$(e.species). The direction follows the " *
                "information: the module whose rate constants are rebuilt from " *
                "the pool declares :in, and the pool's owner declares :out. As " *
                "written this channel is lowered by nobody and $id would keep " *
                "its nominal rate constants for the whole trajectory; declare " *
                "direction = :in"))
        end

        if isempty(names) && isempty(edges)
            continue
        end
        isempty(edges) && throw(ArgumentError(
            "Module $id declares rebuilt_params $(names) but no inbound " *
            "RateConstantEdge. The rebuild reads the pools its edges name, so " *
            "with none declared there is nothing to rebuild from; declare the " *
            "edge, or drop the rebuilt parameters"))
        isempty(names) && throw(ArgumentError(
            "Module $id declares an inbound RateConstantEdge on " *
            "$(join((string(":", e.species) for e in edges), ", ")) but no " *
            "rebuilt_params, so the channel would be declared and never " *
            "executed: the pool would be read by nothing and the module would " *
            "keep its nominal rate constants for the whole trajectory. Name the " *
            "parameters the rebuild fills, or drop the edge"))

        # Task 4.4. The mechanism chosen in task 3.1 — an outer loop over two
        # stepped integrators — holds the constants between refreshes by
        # construction, which is precisely what makes task 4.3's assertion
        # possible. A genuinely continuous cadence needs propensities that track
        # the pools between handshakes, which is the `JumpProblem` over
        # `ODEProblem` shape task 3.1 rejected. Refuse it by name rather than
        # silently discard the cadence field, or run it at 1 s and call that
        # continuous.
        for e in edges
            e.cadence === :continuous && throw(ArgumentError(
                "Module $id declares a RateConstantEdge on :$(e.species) with " *
                "cadence = :continuous, which the handshake driver cannot " *
                "execute. It advances the two blocks by an outer split-operator " *
                "loop (spec §11 task 3.1), so a rate constant is by construction " *
                "held between refreshes; propensities tracking a pool " *
                "continuously would need a JumpProblem over an ODEProblem, which " *
                "that task rejected. Declare cadence = :piecewise_constant with " *
                "a shorter interval instead — that is what a cadence-sensitivity " *
                "comparison varies (spec §10 R11) — or drop the edge"))
        end

        # Every rebuilt name must be one of this module's own free parameters:
        # the slot is a position in *its* propensities, exactly as a catalytic
        # `param_slot` is a position in the consuming module's rate law.
        own = model_free_params(parameters(m))
        ctx = jump_contexts[i]
        fill_idxs = Int[]
        for n in names
            j = findfirst(q -> q.name === n, own)
            j === nothing && throw(ArgumentError(
                "Module $id names :$n in rebuilt_params, which is not one of its " *
                "own free parameters. The rebuild fills a slot in this module's " *
                "own parameter vector; its free parameters are " *
                "$([q.name for q in own])"))
            holders = slot_owners[n]
            length(holders) == 1 || throw(ArgumentError(
                "Module $id rebuilds :$n, but $(join(holders, " and ")) each " *
                "declare a free parameter of that name. Within the stochastic " *
                "block those are one deduplicated slot, so the rebuild would " *
                "drive every one of those modules' propensities; across the " *
                "boundary they are two slots in two parameter vectors, so the " *
                "hook would rewrite one and leave the other, and one named " *
                "parameter would carry two values. Give the parameter a name " *
                "unique to $id"))
            push!(fill_idxs, ctx.param_idxs[j])
        end

        # The refresh interval is the edge's, not the driver's, and the rebuild
        # can only fire at a handshake — that is the one instant at which the two
        # blocks agree on a state. So an interval that is not a whole number of
        # handshakes is unreachable, and rounding it silently would run a cadence
        # nobody declared.
        declared = unique(e.interval for e in edges)
        length(declared) == 1 || throw(ArgumentError(
            "Module $id declares rate-constant edges with different intervals " *
            "($(join(declared, ", "))). One module rebuilds its constants on one " *
            "schedule; split the pools across modules, or declare one interval"))
        iv = first(declared)
        steps = _handshake_steps(iv, interval, "Module $id's rate-constant interval")

        pool_idxs = Int[]
        for e in edges
            k = _block_state_index(ode_models, e.species)
            k === nothing && throw(ArgumentError(
                "Module $id declares an inbound RateConstantEdge on " *
                ":$(e.species), but no ODE module in this composition integrates " *
                "that pool. The rebuild reads a live concentration from the " *
                "metabolic block; compose the module that owns :$(e.species), or " *
                "drop the edge"))
            push!(pool_idxs, k)
        end

        push!(rebuilds, RateConstantRebuild(m, collect(Int, ctx.param_idxs), fill_idxs,
                                            pool_idxs, steps, Float64(iv), collect(names),
                                            [e.species for e in edges], id, 0))
    end

    # The mirror trap: the pool's owner names the channel, nothing consumes it.
    # Nothing would be rebuilt and nothing would say so.
    consumed = Set{Symbol}(s for r in rebuilds for s in r.species)
    for m in ode_models, e in coupling(m)
        e isa RateConstantEdge || continue
        is_producer(e) || continue
        e.species in consumed && continue
        throw(ArgumentError(
            "Module $(module_id(m)) declares an outbound RateConstantEdge on " *
            ":$(e.species), but no jump module in this composition rebuilds any " *
            "constant from it. The channel is lowered from the side whose " *
            "constants are rebuilt, so this declaration would never fire and the " *
            "stochastic block would keep its nominal rate constants. Compose the " *
            "module that reads :$(e.species), or drop the edge"))
    end

    return rebuilds
end

"""
    _build_hybrid_problem(models; tspan, interval, rounding, radius_nm,
                          ode_solver, abstol, reltol)

The mixed-formalism branch of [`build_problem`](@ref), which is how callers
reach it. Returns a [`HandshakeDriver`](@ref) rather than a
SciML problem, because the two blocks are advanced separately and there is no
single problem object that means what the published coupling means.

Keyword arguments beyond `tspan` are the driver's declared policy: the exchange
`interval` (1 s, the published `delt`), the `rounding` policy (spec check 0),
the cell `radius_nm` behind the count↔concentration conversion (fixed until
phase 5), and the ODE solver with its pinned tolerances (spec §3, the tolerance
principle).
"""
function _build_hybrid_problem(models::Vector{<:AbstractSubModel};
                               tspan = (0.0, 100.0),
                               interval = 1.0,
                               drain_interval = nothing,
                               rounding = :fractional_carry,
                               radius_nm = COREA_INITIAL_RADIUS_NM,
                               ode_solver = Rodas5P(),
                               abstol = 1e-10, reltol = 1e-8)
    interval > 0 || throw(ArgumentError(
        "The handshake interval must be positive, got $interval"))
    drain = drain_interval === nothing ? Float64(interval) : Float64(drain_interval)
    drain > 0 || throw(ArgumentError(
        "The drain interval must be positive, got $drain"))
    steps_per_drain = _handshake_steps(drain, interval, "The drain interval")
    _validate_shared_params(models)
    # The contract is validated over the whole composition, not per block: a
    # boundary crossing is by definition not visible from one side.
    resolve_coupling(models)

    ode_models = AbstractSubModel[m for m in models if formalism(m) === :ode]
    jump_models = AbstractSubModel[m for m in models if formalism(m) === :jump]
    isempty(ode_models) && error(
        "A hybrid composition needs at least one :ode sub-model; this one has none")
    isempty(jump_models) && error(
        "A hybrid composition needs at least one :jump sub-model; this one has none")

    # `_check_state_ownership` already spans both blocks — it iterates every
    # model in the composition regardless of formalism — but it skips any name
    # the registry does not know. So a *non-registry* name owned on both sides
    # passes it, and the transcripts phase 10 owns are non-registry by task
    # 10.2. A name owned in both blocks is two different states wearing one
    # name, and every cross-block reference to it — a catalytic count, a debited
    # pool — would silently resolve to whichever block happened to be asked.
    # That is the aliasing bug phase 2 fixed inside a block, arriving between
    # them.
    ode_owned = Dict{Symbol, Symbol}()
    for m in ode_models, s in states(m)
        ode_owned[s] = module_id(m)
    end
    jump_owned = Dict{Symbol, Symbol}()
    for m in jump_models, s in states(m)
        jump_owned[s] = module_id(m)
        haskey(ode_owned, s) && error(
            "State :$s is owned in both blocks of this hybrid composition — by " *
            "$(ode_owned[s]) on the ODE side and $(module_id(m)) on the jump " *
            "side. The blocks hold separate state vectors, so one name for two " *
            "states makes every crossing that mentions :$s ambiguous. Rename " *
            "one, or compose only the module that should own it")
    end

    # `inputs()` and `written_states()` are wired inside a block, by that block's
    # own builder, from that block's own state vector. Naming a state the other
    # block owns therefore cannot work — and left to the block builders it fails
    # as "not owned by any sub-model", which is true of that block and useless as
    # a diagnosis. Say what is actually wrong, and what the boundary offers
    # instead.
    for (block, foreign, other) in ((ode_models, jump_owned, "jump"),
                                    (jump_models, ode_owned, "ODE"))
        for m in block
            for s in inputs(m)
                haskey(foreign, s) && error(
                    "Module $(module_id(m)) lists :$s in inputs(), but :$s is " *
                    "owned by $(foreign[s]) in the $other block. inputs() is " *
                    "resolved within a block, so it cannot reach across the " *
                    "boundary. State crosses it through a declared edge that " *
                    "the handshake executes: a CatalyticEdge to read a count " *
                    "into a rate-law parameter, or a DeferredCounterEdge to " *
                    "debit a pool")
            end
            for s in written_states(m)
                haskey(foreign, s) && error(
                    "Module $(module_id(m)) lists :$s in written_states(), but " *
                    ":$s is owned by $(foreign[s]) in the $other block. A jump " *
                    "module's affect writes its own block's vector and cannot " *
                    "reach the other; debit an ODE pool through a " *
                    "DeferredCounterEdge, which the hook applies between steps")
            end
        end
    end

    factor = corea_particles_per_mM(radius_nm)
    catalytic, debits, counters = _lower_exchanges(models, ode_models, jump_models,
                                                   _build_contexts(ode_models))
    rebuilds = _lower_rebuilds(ode_models, jump_models, _build_contexts(jump_models),
                               Float64(interval))

    # The contract was resolved above, over the whole composition. The block
    # builders must not resolve it again on their own module subset: from one
    # side a boundary crossing looks like an edge naming an absent peer.
    ode_prob = _build_ode_problem(ode_models; tspan = tspan, validate = false)
    jump_prob = _build_jump_problem(jump_models; tspan = tspan, validate = false)

    ode_integ = init(ode_prob, ode_solver; abstol = abstol, reltol = reltol,
                     save_everystep = false)
    jump_integ = init(jump_prob, SSAStepper())

    return HandshakeDriver(ode_integ, jump_integ, Float64(interval), drain,
                           steps_per_drain, factor,
                           RoundingState(rounding; nspecies = length(ode_prob.u0)),
                           catalytic, debits, counters, rebuilds,
                           Float64(tspan[2]), 0, 0, 0)
end

# ---------------------------------------------------------------------------
# The exchange itself (spec §11 tasks 3.4, 3.5, 3.6)
# ---------------------------------------------------------------------------

# The ODE block is built out-of-place over an `SVector` (the state-container
# decision of task 1.6), so a debit reassigns the state rather than mutating
# it, and the integrator is told its state moved under it.
function _set_ode_state!(integ, i::Int, v)
    u = integ.u
    if u isa SVector
        integ.u = setindex(u, v, i)
    else
        u[i] = v
    end
    u_modified!(integ, true)
    return nothing
end

# max(0, x), smoothed. Used only under `clip = :smoothed`, which is a labelled
# departure from the published model; `smoothing` is its width and the edge
# refuses to be constructed without one.
function _soft_excess(x, width)
    z = x / width
    # `exp(z)` overflows for large z. At the cut the two branches differ by
    # `width * log1p(exp(-30))`, about fourteen ulps of the result — negligible
    # against a quantity that is about to be rounded to whole particles, but not
    # exact, so it is stated rather than claimed away.
    return z > 30 ? x : width * log1p(exp(z))
end

"""
    handshake_step!(driver) -> driver

One exchange, in the published model's order: counts become rate-law
parameters, the metabolic block integrates one interval, the accrued costs are
debited against the pools, and then the stochastic block advances one interval.

The debit runs *after* the ODE step and before the jump step, which is what
makes it deferred: the cost accrued during one stochastic interval is paid at
the next hook, exactly as `hookSimulation` does it.
"""
function handshake_step!(d::HandshakeDriver)
    step = d.n_handshakes + 1

    # 1. The enzyme-concentration channel: counts fill rate-law parameter slots.
    for c in d.catalytic
        d.ode.p[c.param_idx] = counts_to_mM(d.jump.u[c.count_idx], d.factor)
    end

    # 2. Integrate the metabolic block one interval.
    target = d.ode.t + d.interval
    step!(d.ode, d.interval, true)
    _assert_ode_advanced(d, target)

    # 3. Debit the deferred counters against their pools.
    #
    # Each counter is read and cleared once, before any pool it feeds is
    # touched, so that a counter feeding several pools hands the same accrued
    # value to all of them. The published charged-tRNA transfer is exactly that
    # shape — one accrual debiting the charged pool and crediting the uncharged
    # one (spec §3) — and reading the counter inside the pool loop would give
    # the second pool a silent zero.
    #
    # `drain_interval` decides how often. At the default of one handshake this
    # is the published 1 s drain; coarser, the counters accrue across several
    # handshakes and the pool sees one aggregate debit, which is the reduction
    # spec §4 D10 measures the cost of and `driver_declarations` labels.
    clipped = false
    drained = step % d.steps_per_drain == 0
    if drained
        for g in d.counters
            accrued_now = float(d.jump.u[g.counter_idx])
            d.jump.u[g.counter_idx] = 0

            # Consumers first, and remember what they actually paid. A credit on the
            # same accrual must match the debit, not the demand: if the drawn pool
            # clips, only what left it may arrive anywhere else. Crediting the raw
            # accrual would create matter out of a shortfall — and on the charged-tRNA
            # transfer this shape exists for, it would do so every time the charged
            # pool ran dry.
            paid_total = 0.0
            consumers = false
            for k in g.debit_idxs
                b = d.debits[k]
                b.sign < 0 || continue
                consumers = true
                accrued = accrued_now + b.deficit
                pool = d.ode.u[b.pool_idx] * d.factor     # the pool, in particles

                if b.clip === :clamped_deficit_carried
                    paid = min(accrued, max(pool, 0.0))
                elseif b.clip === :unclamped
                    paid = accrued
                else                                      # :smoothed
                    paid = accrued - _soft_excess(accrued - pool, b.smoothing)
                end
                b.deficit = accrued - paid
                paid_total += paid

                # Whether this handshake clipped, judged against the cost rather
                # than against zero. Under `:smoothed` the softplus leaves a
                # strictly positive residue at every step however deep the pool, so
                # a `deficit > 0` test would report every handshake as clipping —
                # and `:smoothed` is exactly the policy adopted to escape the
                # clamped counter's gradient obstruction, so the census that scores
                # K5 would be uninformative precisely where it is needed.
                b.deficit > _CLIP_ATOL * max(accrued, 1.0) && (clipped = true)
                _write_pool!(d, b.pool_idx, pool - paid)
            end

            # Then the producers, crediting what was taken rather than what was
            # asked for. With no consumer on this counter there is nothing to match,
            # and the accrual is a pure production.
            credit = consumers ? paid_total : accrued_now
            for k in g.debit_idxs
                b = d.debits[k]
                b.sign > 0 || continue
                b.deficit = 0.0
                _write_pool!(d, b.pool_idx, d.ode.u[b.pool_idx] * d.factor + credit)
            end
        end
    end

    # 3b. The 60 s rebuild: live pools become the stochastic block's rate
    # constants, on each module's own declared schedule (spec §11 phase 4).
    #
    # It runs after the debit, so the pools it reads are the ones the hook just
    # left behind — the published order. The values go into the jump block's
    # parameter vector rather than onto the sub-model, so `remake(prob; p = θ)`
    # and `model_free_params` still see them.
    rebuilt = false
    for r in d.rebuilds
        step % r.steps == 0 || continue
        pools = _fvec([d.ode.u[k] for k in r.pool_idxs])
        vals = rate_constants(view(d.jump.p, r.p_idxs), d.ode.t, r.model, pools)
        length(vals) == length(r.fill_idxs) || error(
            "rate_constants() for $(r.declared_by) returned $(length(vals)) " *
            "value$(length(vals) == 1 ? "" : "s") but rebuilt_params() names " *
            "$(length(r.fill_idxs)) parameter$(length(r.fill_idxs) == 1 ? "" : "s") " *
            "($(r.names)); return a static vector of the same length")
        for (k, v) in zip(r.fill_idxs, vals)
            d.jump.p[k] = v
        end
        r.n_refreshes += 1
        rebuilt = true
    end

    # 4. Advance the stochastic block one interval.
    #
    # The hook cleared the counters, and `Direct` caches propensities and their
    # total between events. Nothing in the phase-3 toy's rate laws reads a
    # counter, so the cache would happen to stay correct — but that is a
    # property of the toy, not of the mechanism, and a module whose propensity
    # reads a state the hook touches would be silently wrong. A rebuild is not
    # a matter of luck at all: it rewrites the parameters every propensity is
    # computed from. Rebuild the aggregation whenever either happened.
    if rebuilt || (drained && !isempty(d.debits))
        reset_aggregated_jumps!(d.jump)
    end
    step!(d.jump, d.interval, true)

    d.n_handshakes += 1
    drained && (d.n_drains += 1)
    clipped && (d.n_clipped += 1)
    return d
end

# A shortfall this far below the cost is the smoothing's own residue rather
# than a pool that ran dry. Loose enough to stay right at a smoothing width
# comparable to the pool, tight enough that missing a real one-particle
# shortfall would need an accrual above a million particles.
const _CLIP_ATOL = 1e-6

# Back to a concentration through whole particles, under the driver's stated
# policy. This is the round trip check 0 is about: a pool is a continuous
# quantity between hooks and a particle count across them.
function _write_pool!(d::HandshakeDriver, idx::Int, particles)
    n = round_to_counts!(d.rounding, idx, particles)
    _set_ode_state!(d.ode, idx, counts_to_mM(n, d.factor))
    return nothing
end

# `step!` returns normally when the solve has failed: it breaks out of its loop
# on a non-success retcode, leaving `t` where it was. Without this the driver
# would keep looping over a frozen — or NaN — state, incrementing the handshake
# count and reporting a clean census. The clamped debit pins a pool at exactly
# zero every time it clips, which is enough to put a rate law with that pool in
# a denominator into `Unstable`, so this is a reachable state and not a
# defensive flourish.
function _assert_ode_advanced(d::HandshakeDriver, target)
    ok = d.ode.sol.retcode === ReturnCode.Default || d.ode.sol.retcode === ReturnCode.Success
    (ok && isapprox(d.ode.t, target; atol = 1e-9, rtol = 1e-12)) || error(
        "The metabolic block failed to advance to t = $target at handshake " *
        "$(d.n_handshakes + 1): it stopped at t = $(d.ode.t) with retcode " *
        "$(d.ode.sol.retcode). The trajectory from here would be a frozen " *
        "state advanced by a counter, so the run stops instead")
    return nothing
end

"""
    run_handshake!(driver, n_steps) -> NamedTuple

Run `n_steps` exchanges, recording the two blocks' states after each.

Returns `(; t, ode, jump, jump_p, census, rebuilds)` — the handshake times, the
ODE state in mM and the jump state in particles at each of them, the jump
block's parameter vector at each of them, check 7's census and
[`rebuild_census`](@ref). The record is per handshake rather than per solver
step because the handshake is the only instant at which the two blocks agree on
a state.

`jump_p` is what makes the rate-constant channel observable from the outside:
a refresh count is read off it by counting the handshakes at which a rebuilt
slot changed, rather than asserted from the schedule that produced it.
"""
function run_handshake!(d::HandshakeDriver, n_steps::Integer)
    (d.n_handshakes + n_steps) % d.steps_per_drain == 0 || error(
        "$n_steps handshakes from handshake $(d.n_handshakes) would end at " *
        "handshake $(d.n_handshakes + n_steps), which is not a multiple of the " *
        "$(d.steps_per_drain)-handshake drain interval. The run would stop with up " *
        "to $(d.steps_per_drain - 1) exchanges of accrued cost still in the " *
        "counters and never debited, so the recorded pools would be high by an " *
        "amount nothing in the census names. Run a whole number of drain " *
        "intervals, or step with handshake_step! if a mid-period state is what " *
        "you want")
    finish = d.ode.t + n_steps * d.interval
    finish <= d.t_end + 1e-9 || error(
        "$n_steps handshakes of $(d.interval) s from t = $(d.ode.t) would reach " *
        "t = $finish, past the tspan both blocks were built with, which ends at " *
        "$(d.t_end). Build the driver with a tspan that covers the run")
    t = Vector{Float64}(undef, n_steps)
    ode = Vector{Vector{Float64}}(undef, n_steps)
    jump = Vector{Vector{Int}}(undef, n_steps)
    jump_p = Vector{Vector{Float64}}(undef, n_steps)
    for i in 1:n_steps
        handshake_step!(d)
        t[i] = d.ode.t
        ode[i] = collect(Float64, d.ode.u)
        jump[i] = collect(Int, d.jump.u)
        jump_p[i] = collect(Float64, d.jump.p)
    end
    return (t = t, ode = ode, jump = jump, jump_p = jump_p,
            census = clipping_census(d), rebuilds = rebuild_census(d))
end

# ---------------------------------------------------------------------------
# The channel's gain (spec §11 task 4.5)
# ---------------------------------------------------------------------------

"""
    rate_constant_elasticity(driver; rel = 1e-4) -> Vector{NamedTuple}

The gain of the rate-constant channel: the elasticity

```
ε = d ln k / d ln pool
```

of every rebuilt rate constant with respect to every pool its module reads, at
the driver's **current** pool concentrations, by a central difference in log
space.

One row per (module, parameter, pool) triple, carrying the pool concentration
it was measured at. This is the diagnostic behind §6 F1's reverse-channel rows:
`dev/notes/reduced-syn3a-scoping.md` measures the real transcription channel at
0.044–0.051, and phase 10 compares against that.

An elasticity is a local derivative, so it is a function of the pools and not a
constant of the model. Spec §10 R8 requires the channel gains to be reported as
functions of the quantities we assert rather than as point values, which is why
the pool concentration travels with the number.
"""
function rate_constant_elasticity(d::HandshakeDriver; rel = 1e-4)
    rows = NamedTuple[]
    for r in d.rebuilds
        pools = [d.ode.u[k] for k in r.pool_idxs]
        p_local = collect(Float64, view(d.jump.p, r.p_idxs))
        k_at(v) = rate_constants(p_local, d.ode.t, r.model, _fvec(v))
        base = collect(Float64, k_at(pools))
        for (j, species) in enumerate(r.species)
            x = pools[j]
            x > 0 || throw(ArgumentError(
                "The elasticity of $(r.declared_by)'s rate constants to :$species " *
                "is a log derivative, and :$species is at $x mM. Measure it at a " *
                "positive pool concentration"))
            up, down = copy(pools), copy(pools)
            up[j] = x * (1 + rel)
            down[j] = x * (1 - rel)
            ku, kd = k_at(up), k_at(down)
            for (i, name) in enumerate(r.names)
                base[i] > 0 || throw(ArgumentError(
                    "The elasticity of :$name is a log derivative and :$name " *
                    "rebuilds to $(base[i]) at the current pools. Measure it " *
                    "where the rate constant is positive"))
                eps = (log(ku[i]) - log(kd[i])) / (log(up[j]) - log(down[j]))
                push!(rows, (module_id = r.declared_by, param = name,
                             pool = species, concentration = x,
                             rate_constant = base[i], elasticity = eps))
            end
        end
    end
    return rows
end

# ---------------------------------------------------------------------------
# Driver-level departures from the published model (spec §6 T2)
# ---------------------------------------------------------------------------

"""
    driver_declarations(driver) -> Vector{ReductionLabel}

Departures from the published model carried by the **driver's policy** rather
than by any module's declarations.

[`reduction_declarations`](@ref) enumerates what the sub-models declare, and
cannot see this: the drain granularity and the count-to-concentration rounding
policy are fields on the driver. Both are rows in spec §6 T2, and both are
departures a result could depend on unremarked — which is the failure the
labelling exists to prevent. Call it beside `reduction_declarations`, not
instead of it.
"""
function driver_declarations(d::HandshakeDriver)
    labels = ReductionLabel[]
    if d.drain_interval > d.interval
        push!(labels, ReductionLabel(
            :coarse_drain, :drain_interval,
            "the deferred counters are debited every $(d.drain_interval) s rather " *
            "than at every $(d.interval) s handshake, so the pools see one " *
            "aggregate drain where the published model applies " *
            "$(d.steps_per_drain) — a labelled reduction whose cost spec §4 D10 " *
            "requires to be measured, not assumed"))
    end
    if d.rounding.policy !== :fractional_carry
        push!(labels, ReductionLabel(
            :rounding_policy, d.rounding.policy,
            "counts are written back under the :$(d.rounding.policy) rounding " *
            "policy rather than fractional carry, which is the only policy whose " *
            "round trip is exact (spec §3, check 0)" *
            (d.rounding.policy === :deterministic ?
             " — and deterministic rounding is rejected by the spec, its residual " *
             "accumulating linearly in the handshake count" : "")))
    end
    return labels
end

"""
    reduction_declarations(models, driver) -> Vector{ReductionLabel}
    reduction_report(models, driver) -> String

The composition's declared departures together with the driver's own policy.

[`reduction_declarations`](@ref) enumerates what the sub-models declare and
structurally cannot see a field on the driver; [`driver_declarations`](@ref)
covers the rest. Spec §6 F11's enumeration panel and §6 T2 are generated from
these two-argument forms, because a result produced under a coarse drain or a
non-carry rounding policy that carried only the one-argument enumeration would
be exactly the unlabelled departure §6 T2 exists to prevent.
"""
reduction_declarations(models::Vector{<:AbstractSubModel}, d::HandshakeDriver) =
    vcat(reduction_declarations(models), driver_declarations(d))

reduction_report(models::Vector{<:AbstractSubModel}, d::HandshakeDriver) =
    _reduction_report(reduction_declarations(models, d); driver_seen = true)

export handshake_step!, run_handshake!, rate_constant_elasticity,
       driver_declarations
