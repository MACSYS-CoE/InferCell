using InferCell
using Distributions: LogNormal, Normal
using StaticArrays

import InferCell: states, parameters, dynamics, reactions, formalism, inputs,
                  coupling, module_id, inference_mode, reduction_notes

# Executable doubles for the 1 s handshake (spec §11 phase 3).
#
# Purpose-built rather than assembled from the shipped models, for the reason
# phases 1 and 2 built their own: tasks 3.5 and 3.6 need a real
# `DeferredCounterEdge` and a real `CatalyticEdge`, and no shipped model
# declares any typed coupling at all.
#
# Both doubles name **registry** species where an edge touches them, because
# `_resolve_edges` refuses an edge on a species outside the Core A′ registry.
# The protein is `M_ptsg_c` — Core A′'s only membrane protein, and the species
# `test/test_edges.jl` already uses as its catalytic example — and the debited
# pool is `M_atp_c`. `_rate` and `_jic` come from the phase-1 and phase-2
# doubles files, which `runtests.jl` includes first.

# `_jic` under a name that reads right in an ODE double: an initial condition
# is the same `InferParameter` either way, and it is the composed builder that
# rounds one block's to whole particles and leaves the other's in mM.
_toy_ic(x0, name) = _jic(x0, name, :ToyPool)

"""
    ToyExpression(; k_tx = 2.0, gamma_m = 0.5, k_tl = 1.0, cost = 10,
                    mrna0 = 0, protein0 = 0, edges = <deferred counter on ATP>)

The stochastic half of the toy: one gene, transcribed, decayed and translated.
Owns its transcript, the protein `M_ptsg_c`, and the accrual counter
`:atp_cost`. Translation charges `cost` particles of ATP into that counter,
which the hook debits against the ODE block's pool one handshake later.

`edges` is overridable so a test can build a *mis*-declared variant — a counter
no module owns, a pool no module integrates — and witness the refusal.
"""
struct ToyExpression <: AbstractSubModel
    params::Vector{InferParameter}
    cost::Int
    edges::Vector{CouplingEdge}
end

function ToyExpression(; k_tx = 2.0, gamma_m = 0.5, k_tl = 1.0, cost = 10,
                       mrna0 = 0, protein0 = 0,
                       edges = CouplingEdge[DeferredCounterEdge(species = :M_atp_c,
                                                                direction = :in,
                                                                counter = :atp_cost)])
    # `edges` is overridable so a test can drive the two labelled clip policies,
    # a producing counter, and one counter feeding two pools — the shape spec §3
    # gives the charged-tRNA transfer.
    params = [_rate(k_tx, :k_tx_toy, :ToyExpression),
              _rate(gamma_m, :gamma_m_toy, :ToyExpression),
              _rate(k_tl, :k_tl_toy, :ToyExpression),
              _jic(mrna0, :toy_mrna0, :ToyExpression),
              _jic(protein0, :M_ptsg_c0, :ToyExpression),
              _jic(0, :atp_cost0, :ToyExpression)]
    return ToyExpression(params, Int(cost), collect(CouplingEdge, edges))
end

states(::ToyExpression) = [:toy_mrna, :M_ptsg_c, :atp_cost]
parameters(m::ToyExpression) = m.params
formalism(::ToyExpression) = :jump
inference_mode(::ToyExpression) = :simulation
coupling(m::ToyExpression) = m.edges
function reactions(m::ToyExpression)
    cost = m.cost
    return [
        # ∅ -> mRNA
        Reaction((u, p, t, _) -> p[1], (u, _) -> (u[1] += 1)),
        # mRNA -> ∅
        Reaction((u, p, t, _) -> p[2] * u[1], (u, _) -> (u[1] -= 1)),
        # mRNA -> mRNA + protein, charging `cost` particles of ATP to the counter
        Reaction((u, p, t, _) -> p[3] * u[1], (u, _) -> (u[2] += 1; u[3] += cost)),
    ]
end

"""
    ToyPool(; kcat = 30.0, km = 1.0, enzyme = 1.0, atp0 = 3.6529, adp0 = 0.2178,
              edges = <catalytic edge on M_ptsg_c>)

The metabolic half of the toy: one irreversible reaction turning ATP into ADP,
Michaelis–Menten in ATP and first order in the enzyme concentration.

`enzyme` is the `:enzyme_conc` parameter the catalytic edge fills, so its
declared value is only the value before the first handshake overwrites it.
`kcat = 0.0` gives the **zero-derivative** variant task 3.4 needs, which is why
the rate constant is a keyword rather than a constant: with no derivative, any
drift over 600 handshakes is the exchange's and nothing else's.

Initial conditions are the registry's own for both species.
"""
struct ToyPool <: AbstractSubModel
    params::Vector{InferParameter}
    edges::Vector{CouplingEdge}
end

function ToyPool(; kcat = 30.0, km = 1.0, enzyme = 1.0,
                 atp0 = 3.6529, adp0 = 0.2178,
                 edges = CouplingEdge[CatalyticEdge(species = :M_ptsg_c,
                                                    direction = :in,
                                                    param_slot = :enzyme_conc)])
    params = [_rate(kcat, :kcat_toy, :ToyPool),
              _rate(km, :km_toy, :ToyPool),
              # The enzyme concentration must be a *free* parameter: only
              # those reach the composed parameter vector, and the catalytic
              # edge needs a slot in it to fill. Its declared value is nominal
              # and overridable, which is what task 6.4 asks of the real
              # modules too.
              _rate(enzyme, :enzyme_conc, :ToyPool),
              _toy_ic(atp0, :M_atp_c0),
              _toy_ic(adp0, :M_adp_c0)]
    return ToyPool(params, collect(CouplingEdge, edges))
end

states(::ToyPool) = [:M_atp_c, :M_adp_c]
parameters(m::ToyPool) = m.params
formalism(::ToyPool) = :ode
coupling(m::ToyPool) = m.edges
function dynamics(u, p, t, ::ToyPool)
    kcat, km, enzyme = p[1], p[2], p[3]
    atp = u[1]
    v = kcat * enzyme * atp / (km + atp)
    return SA[-v, v]
end

"""
    toy_flux(prob_or_p, u) -> Float64

The toy pool's reaction rate at a given parameter vector and state — the
quantity task 3.6 asserts scales with the protein count. Written once here so
the test asserts the model's own rate law rather than a copy of it.
"""
toy_flux(p, u) = p[1] * p[3] * u[1] / (p[2] + u[1])

"""
    ToyPool2()

A second ODE pool module that also declares a free parameter named
`:enzyme_conc`. It exists to witness one refusal: parameters deduplicate by
name into a single slot, so a catalytic edge filling `:enzyme_conc` would drive
this module's rate law as well as its own. With seventeen genes and several
enzyme-catalysed modules, reusing the obvious name is the obvious mistake.
"""
struct ToyPool2 <: AbstractSubModel
    params::Vector{InferParameter}
end
ToyPool2() = ToyPool2([_rate(2.0, :kcat_toy2, :ToyPool2),
                       _rate(1.0, :enzyme_conc, :ToyPool2),
                       _jic(1.0, :M_gtp_c0, :ToyPool2)])
states(::ToyPool2) = [:M_gtp_c]
parameters(m::ToyPool2) = m.params
formalism(::ToyPool2) = :ode
dynamics(u, p, t, ::ToyPool2) = SA[-p[1] * p[2] * u[1]]

# --- Phase 4: the 60 s rebuild ---------------------------------------------
#
# A second stochastic double, whose transcription rate constant is *rebuilt*
# from live ODE pools rather than declared once. It owns the same three states
# as `ToyExpression` — the protein `M_ptsg_c` so `ToyPool`'s catalytic edge
# still has a count to read, and the accrual counter `:atp_cost` so the
# deferred debit still runs — which is what lets one composition exercise all
# three executed channels at once, the shape the assembled model has. The two
# doubles are never composed together, so the shared state names cannot
# collide.
#
# Written as its own type rather than as a keyword on `ToyExpression` so that
# every phase-3 test keeps running exactly the model it was written against.

import InferCell: rebuilt_params, rate_constants

"""
    ToyRebuiltExpression(; k_tx = 2.0, k_tx_max = 6.0, km_tx = 1.0,
                           gamma_m = 0.5, k_tl = 1.0, cost = 10,
                           mrna0 = 0, protein0 = 0,
                           pools = [:M_atp_c], interval = 60.0,
                           cadence = :piecewise_constant,
                           rebuilt = [:k_tx_rb],
                           edges = <rate-constant + deferred counter>)

`ToyExpression` with its rate constants on the 60 s rebuild. The `i`-th name in
`rebuilt` takes

```
k_i = (k_tx_max / i) * Π_j pool_j / (km_tx + pool_j)
```

A Michaelis-shaped law in each pool, chosen because its elasticity has a closed
form — `d ln k_i / d ln pool_j = km_tx / (km_tx + pool_j)`, the same for every
`i` — so task 4.5's diagnostic is checked against an analytic value rather than
merely reported. `km_tx` therefore tunes the channel gain, which is what lets a
test drive it to the 0.044–0.051 range the real transcription channel is
measured at. The `1/i` scaling is arbitrary and exists only so that two rebuilt
constants take different values.

Rebuilding more than one name matters for a reason beyond exercising the
plural path. With `rebuilt = [:k_tx_rb, :k_tl_rb]` the translation propensity
`k_tl_rb · mRNA` becomes a *state-dependent* function of a rebuilt constant,
which is what makes "the propensities change only at interval boundaries"
assertable at a held state as something more than a restatement of the
parameter slot. Transcription is zeroth order — in Core A′ as here — so its
propensity simply *is* its rate constant.

Every keyword that names a declaration — `pools`, `interval`, `cadence`,
`rebuilt`, `edges` — is overridable so a test can build a mis-declared variant
and witness the refusal.
"""
struct ToyRebuiltExpression <: AbstractSubModel
    params::Vector{InferParameter}
    cost::Int
    edges::Vector{CouplingEdge}
    rebuilt::Vector{Symbol}
end

function ToyRebuiltExpression(; k_tx = 2.0, k_tx_max = 6.0, km_tx = 1.0,
                              gamma_m = 0.5, k_tl = 1.0, cost = 10,
                              mrna0 = 0, protein0 = 0,
                              pools = [:M_atp_c], interval = 60.0,
                              cadence = :piecewise_constant,
                              rebuilt = [:k_tx_rb],
                              edges = nothing)
    params = [_rate(k_tx, :k_tx_rb, :ToyRebuiltExpression),
              _rate(k_tx_max, :k_tx_max_rb, :ToyRebuiltExpression),
              _rate(km_tx, :km_tx_rb, :ToyRebuiltExpression),
              _rate(gamma_m, :gamma_m_rb, :ToyRebuiltExpression),
              _rate(k_tl, :k_tl_rb, :ToyRebuiltExpression),
              _jic(mrna0, :toy_mrna0, :ToyRebuiltExpression),
              _jic(protein0, :M_ptsg_c0, :ToyRebuiltExpression),
              _jic(0, :atp_cost0, :ToyRebuiltExpression)]
    default = CouplingEdge[DeferredCounterEdge(species = :M_atp_c,
                                               direction = :in,
                                               counter = :atp_cost)]
    for s in pools
        push!(default, cadence === :continuous ?
              RateConstantEdge(species = s, direction = :in, cadence = :continuous) :
              RateConstantEdge(species = s, direction = :in, cadence = cadence,
                               interval = interval))
    end
    return ToyRebuiltExpression(params, Int(cost),
                                collect(CouplingEdge, edges === nothing ? default : edges),
                                collect(Symbol, rebuilt))
end

states(::ToyRebuiltExpression) = [:toy_mrna, :M_ptsg_c, :atp_cost]
parameters(m::ToyRebuiltExpression) = m.params
formalism(::ToyRebuiltExpression) = :jump
inference_mode(::ToyRebuiltExpression) = :simulation
coupling(m::ToyRebuiltExpression) = m.edges
rebuilt_params(m::ToyRebuiltExpression) = m.rebuilt

# The rebuild reads only the pools and the module's *other* parameters. Reading
# a slot it is about to fill would compound that slot's own previous value.
function rate_constants(p, t, m::ToyRebuiltExpression, pools)
    kmax, km = p[2], p[3]
    base = kmax * prod(x / (km + x) for x in pools)
    n = length(m.rebuilt)
    return SVector{n, Float64}(ntuple(i -> base / i, n))
end

"""
    toy_rebuilt_k(m, pools[, i]) -> Float64

The rebuild law for the `i`-th rebuilt name, written once so a test asserts the
model's own arithmetic rather than a copy of it.
"""
toy_rebuilt_k(m::ToyRebuiltExpression, pools, i = 1) =
    m.params[2].value * prod(x / (m.params[3].value + x) for x in pools) / i

"""
    toy_rebuilt_elasticity(m, pool) -> Float64

The closed-form elasticity of the law above with respect to one pool,
`km / (km + pool)`, which task 4.5's finite-difference diagnostic is checked
against.
"""
toy_rebuilt_elasticity(m::ToyRebuiltExpression, pool) =
    m.params[3].value / (m.params[3].value + pool)

function reactions(m::ToyRebuiltExpression)
    cost = m.cost
    return [
        # ∅ -> mRNA, at the rebuilt rate constant
        Reaction((u, p, t, _) -> p[1], (u, _) -> (u[1] += 1)),
        # mRNA -> ∅
        Reaction((u, p, t, _) -> p[4] * u[1], (u, _) -> (u[1] -= 1)),
        # mRNA -> mRNA + protein, charging `cost` particles of ATP to the
        # counter. First order in the transcript, as translation is in Core A′,
        # so with :k_tl_rb rebuilt this propensity depends on both a rebuilt
        # constant and the state.
        Reaction((u, p, t, _) -> p[5] * u[1], (u, _) -> (u[2] += 1; u[3] += cost)),
    ]
end

"""
    toy_slow_pool(; kwargs...)

`ToyPool` with a slow enzyme and a large pool: ATP declines at every handshake
and ends around two thirds full, so an upstream pool that genuinely *moves* is
available without the drain ever clipping. Phase 3's long-horizon tests freeze
the pool instead (`kcat = 0.0`), which is the opposite of what phase 4 needs.
Defined here rather than in each caller so the test suite and
`dev/scripts/rebuild_channel.jl` cannot drift apart on what "slow" means.
"""
toy_slow_pool(; kwargs...) = ToyPool(; kcat = 0.02, atp0 = 20.0, kwargs...)

"""
    ToyPoolRebuildPeer(; species = :M_gtp_c)

An ODE double that declares the *outbound* half of a rate-constant channel and
owns `:M_gtp_c`. It exists to witness one refusal: a pool's owner naming a
channel no jump module consumes, which would be declared and never executed.
"""
struct ToyPoolRebuildPeer <: AbstractSubModel
    params::Vector{InferParameter}
    edges::Vector{CouplingEdge}
end

ToyPoolRebuildPeer(; species = :M_gtp_c) =
    ToyPoolRebuildPeer([_rate(0.1, :k_peer_rb, :ToyPoolRebuildPeer),
                        _jic(1.0, :M_gtp_c0, :ToyPoolRebuildPeer)],
                       CouplingEdge[RateConstantEdge(species = species,
                                                     direction = :out)])

states(::ToyPoolRebuildPeer) = [:M_gtp_c]
parameters(m::ToyPoolRebuildPeer) = m.params
formalism(::ToyPoolRebuildPeer) = :ode
coupling(m::ToyPoolRebuildPeer) = m.edges
dynamics(u, p, t, ::ToyPoolRebuildPeer) = SA[-p[1] * u[1]]

"""
    ToyPoolRebuilder()

An `:ode` double that wrongly declares `rebuilt_params`. The rate-constant
channel rebuilds the *stochastic* block's constants; an ODE module reads a pool
directly.
"""
struct ToyPoolRebuilder <: AbstractSubModel
    params::Vector{InferParameter}
end
ToyPoolRebuilder() = ToyPoolRebuilder([_rate(1.0, :k_odereb, :ToyPoolRebuilder),
                                       _jic(1.0, :M_gtp_c0, :ToyPoolRebuilder)])
states(::ToyPoolRebuilder) = [:M_gtp_c]
parameters(m::ToyPoolRebuilder) = m.params
formalism(::ToyPoolRebuilder) = :ode
rebuilt_params(::ToyPoolRebuilder) = [:k_odereb]
dynamics(u, p, t, ::ToyPoolRebuilder) = SA[-p[1] * u[1]]

"""
    ToyRebuildClash()

A second jump double that declares a free parameter named `:k_tx_rb`. It exists
to witness one refusal: free parameters deduplicate by name into a single slot,
so a rebuild into a shared name would drive this module's propensities too —
the trap `param_slot` has on the catalytic side, arriving on the rebuild side.
"""
struct ToyRebuildClash <: AbstractSubModel
    params::Vector{InferParameter}
end
# The declared value matches `ToyRebuiltExpression`'s, so `_validate_shared_params`
# passes and the refusal under test is the one this double exists for rather
# than the shared-parameter disagreement check firing first.
ToyRebuildClash() = ToyRebuildClash([_rate(2.0, :k_tx_rb, :ToyRebuildClash),
                                     _jic(0, :toy_clash0, :ToyRebuildClash)])
states(::ToyRebuildClash) = [:toy_clash]
parameters(m::ToyRebuildClash) = m.params
formalism(::ToyRebuildClash) = :jump
reactions(::ToyRebuildClash) = [Reaction((u, p, t, _) -> p[1], (u, _) -> (u[1] += 1))]

"""
    ToyOdeNameClash()

An `:ode` double declaring a free parameter named `:k_tx_rb`, with the same
declared value as `ToyRebuiltExpression`'s so `_validate_shared_params` passes.
It witnesses the cross-block half of the rebuild name clash: the two blocks
hold separate parameter vectors, so the hook would rewrite the jump block's
copy and leave this one, and a single named parameter would carry two values.
"""
struct ToyOdeNameClash <: AbstractSubModel
    params::Vector{InferParameter}
end
ToyOdeNameClash() = ToyOdeNameClash([_rate(2.0, :k_tx_rb, :ToyOdeNameClash),
                                     _jic(1.0, :M_gtp_c0, :ToyOdeNameClash)])
states(::ToyOdeNameClash) = [:M_gtp_c]
parameters(m::ToyOdeNameClash) = m.params
formalism(::ToyOdeNameClash) = :ode
dynamics(u, p, t, ::ToyOdeNameClash) = SA[-p[1] * u[1]]

# ---------------------------------------------------------------------------
# Doubles for the volume channel (spec §11 phase 5)
# ---------------------------------------------------------------------------

import InferCell: membrane_protein_states, extracellular_states

# The conversion factor at the published *initial surface area*, which is the
# primitive once a composition grows. Written out here rather than imported so
# that the doubles and their assertions do not both come from the code under
# test.
_toy_factor0() = 1e-3 * ((4 / 3) * π * (sqrt(502831.0 / (4 * π)) * 1e-9)^3 * 1e3) *
                 6.02214076e23

"""
    toy_growth(n, n0; total = 502831.0, footprint = 28.0)

The published growth law, written out longhand: a baseline frozen so that `n0`
membrane proteins come to `total` nm², then `n` of them, then radius, volume and
volume against the initial one.

Independent of `src/handshake.jl` on purpose — it is the closed form the
diagnostic is checked against, as `toy_rebuilt_k` is for the rebuild channel, so
that a test asserts the law and not the implementation restated.
"""
function toy_growth(n, n0; total = 502831.0, footprint = 28.0)
    base = total - footprint * n0
    vol(a) = (4 / 3) * π * (sqrt(a / (4 * π)) * 1e-9)^3 * 1e3
    area = base + footprint * n
    return (baseline_nm2 = base,
            area_nm2 = area,
            radius_nm = sqrt(area / (4 * π)),
            volume_litres = vol(area),
            fractional = vol(area) / vol(total))
end

"""
    ToyGrowingExpression(; k_tx = 0.0, gamma_m = 0.5, k_tl = 1.0, cost = 10,
                           mrna0 = 0, protein0 = 831,
                           membrane = [:M_ptsg_c], edges = <deferred + volume>)

The stochastic half of the growth toy: `ToyExpression` with its membrane protein
flagged and an outbound `VolumeEdge` on it, so the protein count drives surface
area, radius, volume and every conversion.

Its own type rather than a keyword on `ToyExpression`, for the reason
`ToyRebuiltExpression` is: every phase-3 and phase-4 test keeps running exactly
the model it was written against.

`k_tx` defaults to **zero**, which is what makes the growth arithmetic testable:
with no transcription nothing fires, the count stays where a test puts it, and
the area is a function of that count alone rather than of a random path.
`membrane`, `extracellular` and `edges` are overridable so a test can build a
mis-declared variant — a flag with no edge, an edge with nothing flagged, a flag
on a state this module does not own, a jump module claiming an extracellular
state — and witness the refusal.
"""
struct ToyGrowingExpression <: AbstractSubModel
    params::Vector{InferParameter}
    cost::Int
    edges::Vector{CouplingEdge}
    membrane::Vector{Symbol}
    extracellular::Vector{Symbol}
end

function ToyGrowingExpression(; k_tx = 0.0, gamma_m = 0.5, k_tl = 1.0, cost = 10,
                              mrna0 = 0, protein0 = 831,
                              membrane = [:M_ptsg_c],
                              extracellular = Symbol[],
                              edges = CouplingEdge[
                                  DeferredCounterEdge(species = :M_atp_c,
                                                      direction = :in,
                                                      counter = :atp_cost),
                                  VolumeEdge(species = :M_ptsg_c,
                                             direction = :out)])
    params = [_rate(k_tx, :k_tx_gr, :ToyGrowingExpression),
              _rate(gamma_m, :gamma_m_gr, :ToyGrowingExpression),
              _rate(k_tl, :k_tl_gr, :ToyGrowingExpression),
              _jic(mrna0, :toy_mrna0, :ToyGrowingExpression),
              _jic(protein0, :M_ptsg_c0, :ToyGrowingExpression),
              _jic(0, :atp_cost0, :ToyGrowingExpression)]
    return ToyGrowingExpression(params, Int(cost), collect(CouplingEdge, edges),
                                collect(Symbol, membrane),
                                collect(Symbol, extracellular))
end

states(::ToyGrowingExpression) = [:toy_mrna, :M_ptsg_c, :atp_cost]
parameters(m::ToyGrowingExpression) = m.params
formalism(::ToyGrowingExpression) = :jump
inference_mode(::ToyGrowingExpression) = :simulation
coupling(m::ToyGrowingExpression) = m.edges
membrane_protein_states(m::ToyGrowingExpression) = m.membrane
# A :jump module may not declare this; the keyword exists only so a test can
# witness the refusal.
extracellular_states(m::ToyGrowingExpression) = m.extracellular
function reactions(m::ToyGrowingExpression)
    cost = m.cost
    return [
        Reaction((u, p, t, _) -> p[1], (u, _) -> (u[1] += 1)),
        Reaction((u, p, t, _) -> p[2] * u[1], (u, _) -> (u[1] -= 1)),
        Reaction((u, p, t, _) -> p[3] * u[1], (u, _) -> (u[2] += 1; u[3] += cost)),
    ]
end

"""
    ToyMembranePool(; kcat = 0.0, km = 1.0, enzyme = 1.0, atp0 = 3.6529,
                      adp0 = 0.2178, ptsi_particles = 353, lac_e0 = 0.5,
                      membrane = [:M_ptsi_c], extracellular = [:M_lac__L_e],
                      edges = <catalytic + volume>)

The metabolic half of the growth toy, and the only double that exercises the
*ODE* side of both new declarations: it owns a membrane protein as a
concentration, and a state whose millimolar is referred to the medium rather
than to the cell.

`ptsi_particles` is given in particles and converted at the initial factor, so a
test can say what the membrane count is rather than what concentration would
produce it. `kcat` defaults to zero so that a change in an ODE state is the
dilution and nothing else.
"""
struct ToyMembranePool <: AbstractSubModel
    params::Vector{InferParameter}
    edges::Vector{CouplingEdge}
    membrane::Vector{Symbol}
    extracellular::Vector{Symbol}
end

function ToyMembranePool(; kcat = 0.0, km = 1.0, enzyme = 1.0,
                         atp0 = 3.6529, adp0 = 0.2178,
                         ptsi_particles = 353, lac_e0 = 0.5,
                         membrane = [:M_ptsi_c],
                         extracellular = [:M_lac__L_e],
                         edges = CouplingEdge[
                             CatalyticEdge(species = :M_ptsg_c,
                                           direction = :in,
                                           param_slot = :enzyme_conc),
                             VolumeEdge(species = :M_ptsi_c,
                                        direction = :out)])
    params = [_rate(kcat, :kcat_toy, :ToyMembranePool),
              _rate(km, :km_toy, :ToyMembranePool),
              _rate(enzyme, :enzyme_conc, :ToyMembranePool),
              _toy_ic(atp0, :M_atp_c0),
              _toy_ic(adp0, :M_adp_c0),
              _toy_ic(ptsi_particles / _toy_factor0(), :M_ptsi_c0),
              _toy_ic(lac_e0, :M_lac__L_e0)]
    return ToyMembranePool(params, collect(CouplingEdge, edges),
                           collect(Symbol, membrane),
                           collect(Symbol, extracellular))
end

states(::ToyMembranePool) = [:M_atp_c, :M_adp_c, :M_ptsi_c, :M_lac__L_e]
parameters(m::ToyMembranePool) = m.params
formalism(::ToyMembranePool) = :ode
coupling(m::ToyMembranePool) = m.edges
membrane_protein_states(m::ToyMembranePool) = m.membrane
extracellular_states(m::ToyMembranePool) = m.extracellular
function dynamics(u, p, t, ::ToyMembranePool)
    kcat, km, enzyme = p[1], p[2], p[3]
    atp = u[1]
    v = kcat * enzyme * atp / (km + atp)
    return SA[-v, v, 0.0, 0.0]
end
