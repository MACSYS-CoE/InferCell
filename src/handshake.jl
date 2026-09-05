"""
The 1 s handshake: executing the coupling the edge kinds declare.

Spec §11 phase 3. Phase 1 made mass and currency edges execute inside the ODE
block; this file adds the catalytic and deferred-counter channels, so four of
the seven [`CouplingEdge`](@ref) kinds now run and three — rate constant,
volume and clamped — remain declare-and-validate only. It
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

# ---------------------------------------------------------------------------
# The driver
# ---------------------------------------------------------------------------

"""
    HandshakeDriver

A mixed ODE/jump composition, advanced by a split-operator exchange.

Holds two live integrators and the lowered exchange records between them.
`interval` is the exchange period — 1.0 s, the published model's `delt`. The
60 s rate-constant rebuild (phase 4) and the volume chain (phase 5) hook into
the same loop and are deliberately absent here.

`n_handshakes` and `n_clipped` are check 7's census: how many exchanges ran,
and at how many of them any deferred counter carried a deficit. K5 is scored
from the second number.
"""
mutable struct HandshakeDriver{OI, JI}
    ode::OI
    jump::JI
    interval::Float64
    factor::Float64                       # particles per mM, at fixed volume
    rounding::RoundingState
    catalytic::Vector{CatalyticExchange}
    debits::Vector{DeferredDebit}
    counters::Vector{_CounterRead}        # debits grouped by the counter feeding them
    t_end::Float64                        # the declared end of both blocks' tspan
    n_handshakes::Int
    n_clipped::Int
end

"""
    clipping_census(driver) -> NamedTuple

Check 7 on this driver: the number of handshakes run, the number at which any
deferred counter carried a deficit, and the fraction. Spec §3 check 7 requires
**zero at published parameters**; more than 5% of prior draws clipping puts the
non-smooth drain in the operating regime and fires K5.
"""
clipping_census(d::HandshakeDriver) = (handshakes = d.n_handshakes,
                                       clipped = d.n_clipped,
                                       fraction = d.n_handshakes == 0 ? 0.0 :
                                                  d.n_clipped / d.n_handshakes,
                                       deficits = [b.deficit for b in d.debits])

export CatalyticExchange, DeferredDebit, HandshakeDriver, clipping_census

# ---------------------------------------------------------------------------
# Building a hybrid composition (spec §11 task 3.2)
# ---------------------------------------------------------------------------

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
    slot_owners = Dict{Symbol, Vector{Symbol}}()
    for m in ode_models, q in model_free_params(parameters(m))
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
                "catalytic edge, but $(join(holders, " and ")) all declare a " *
                "free parameter of that name. Parameters are deduplicated by " *
                "name into one slot, so the count would drive every one of " *
                "those rate laws. Give the slot a name unique to $id"))
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
            push!(debits, DeferredDebit(_block_state_index(jump_models, e.counter),
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
                               rounding = :fractional_carry,
                               radius_nm = COREA_INITIAL_RADIUS_NM,
                               ode_solver = Rodas5P(),
                               abstol = 1e-10, reltol = 1e-8)
    interval > 0 || throw(ArgumentError(
        "The handshake interval must be positive, got $interval"))
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

    # The contract was resolved above, over the whole composition. The block
    # builders must not resolve it again on their own module subset: from one
    # side a boundary crossing looks like an edge naming an absent peer.
    ode_prob = _build_ode_problem(ode_models; tspan = tspan, validate = false)
    jump_prob = _build_jump_problem(jump_models; tspan = tspan, validate = false)

    ode_integ = init(ode_prob, ode_solver; abstol = abstol, reltol = reltol,
                     save_everystep = false)
    jump_integ = init(jump_prob, SSAStepper())

    return HandshakeDriver(ode_integ, jump_integ, Float64(interval), factor,
                           RoundingState(rounding; nspecies = length(ode_prob.u0)),
                           catalytic, debits, counters, Float64(tspan[2]), 0, 0)
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
    clipped = false
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

    # 4. Advance the stochastic block one interval.
    #
    # The hook cleared the counters, and `Direct` caches propensities and their
    # total between events. Nothing in the phase-3 toy's rate laws reads a
    # counter, so the cache would happen to stay correct — but that is a
    # property of the toy, not of the mechanism, and a module whose propensity
    # reads a state the hook touches would be silently wrong. Rebuild the
    # aggregation instead of relying on it.
    isempty(d.debits) || reset_aggregated_jumps!(d.jump)
    step!(d.jump, d.interval, true)

    d.n_handshakes += 1
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

Returns `(; t, ode, jump, census)` — the handshake times, the ODE state in mM
and the jump state in particles at each of them, and check 7's census. The
record is per handshake rather than per solver step because the handshake is
the only instant at which the two blocks agree on a state.
"""
function run_handshake!(d::HandshakeDriver, n_steps::Integer)
    finish = d.ode.t + n_steps * d.interval
    finish <= d.t_end + 1e-9 || error(
        "$n_steps handshakes of $(d.interval) s from t = $(d.ode.t) would reach " *
        "t = $finish, past the tspan both blocks were built with, which ends at " *
        "$(d.t_end). Build the driver with a tspan that covers the run")
    t = Vector{Float64}(undef, n_steps)
    ode = Vector{Vector{Float64}}(undef, n_steps)
    jump = Vector{Vector{Int}}(undef, n_steps)
    for i in 1:n_steps
        handshake_step!(d)
        t[i] = d.ode.t
        ode[i] = collect(Float64, d.ode.u)
        jump[i] = collect(Int, d.jump.u)
    end
    return (t = t, ode = ode, jump = jump, census = clipping_census(d))
end

export handshake_step!, run_handshake!
