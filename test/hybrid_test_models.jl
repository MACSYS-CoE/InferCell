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
