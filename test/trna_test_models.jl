using InferCell
using StaticArrays: SA
using Distributions: Normal

import InferCell: states, parameters, dynamics, contributions, contributed_states,
                  coupling, inputs, module_id, reduction_notes

# Executable doubles for spec §11 phase 9.
#
# The charging step owns the tRNA pair, and standalone nothing consumes charged
# tRNA: the pool saturates, the flux goes to zero, and every conservation check
# passes trivially. So its acceptance test is a composition with a stand-in for
# phase 11's translation, and this file holds it. Phase 8's composition and
# helpers (`nucleotide_test_models.jl`) are included first and reused.

"""
The published demand, recomputed here from the residue count rather than read
from the module under test, so a test comparing the two compares two
derivations and not one number with itself (spec §12, 2026-09-11 and
2026-09-23). Deliberately not phase 8's `RECYCLING_DRAIN_*` either.
"""
const TL_DEMAND_PER_S = 3_484_518 / 6300
const TL_DEMAND_MM_PER_S = TL_DEMAND_PER_S / corea_particles_per_mM()

"""
    TranslationDemand(; k_tl)

Phase 11's translation reduced to what phase 9 needs from it: charged tRNA
consumed and uncharged tRNA returned, one each per residue.

**First order in the charged pool, not a constant drain, and deliberately.**
Against a constant drain the steady charging flux equals the drain's rate by
mass balance whatever `k_chg` is, so the check would test the double. A
first-order consumer makes the steady split, and so the flux, depend on both
modules. It is also what translation becomes once its rate constant reads the
charged pool (task 11.3). `k_tl` defaults to the published demand divided by
the nominal charged pool, so at the nominal split it consumes exactly
553.10 residues per second.

It owns one non-registry state, the millimolar of residues translated, so the
delivered flux is read off the trajectory as an exact integral rather than by
quadrature over save points.
"""
struct TranslationDemand <: AbstractSubModel
    params::Vector{InferParameter}
    k_tl::Float64
end
function TranslationDemand(; k_tl = TL_DEMAND_MM_PER_S /
                                   (CHARGING_POOL_DEFAULTS.charged_fraction *
                                    CHARGING_POOL_DEFAULTS.pool_mM))
    ic = InferParameter(0.0, Normal(0.0, 0.1), true, :tl_residues_mM0,
                        :TranslationDemand, :initial_condition)
    return TranslationDemand([ic], k_tl)
end
states(::TranslationDemand) = [:tl_residues_mM]
parameters(m::TranslationDemand) = m.params
module_id(::TranslationDemand) = :TranslationDemand
inputs(::TranslationDemand) = [:M_trna_chg_c]
contributed_states(::TranslationDemand) = [:M_trna_c, :M_trna_chg_c]
coupling(::TranslationDemand) = CouplingEdge[
    MassEdge(species = :M_trna_chg_c, direction = :in),
    MassEdge(species = :M_trna_c, direction = :out),
]

@inline _tl_rate(m::TranslationDemand, chg) = m.k_tl * chg
dynamics(u, p, t, m::TranslationDemand, u_inputs) = SA[_tl_rate(m, u_inputs[1])]
function contributions(u, p, t, m::TranslationDemand, u_inputs)
    v = _tl_rate(m, u_inputs[1])
    return SA[v, -v]          # uncharged returned, charged consumed
end

"""
    ChargingCounter(inner)

The real charging module with a cumulative counter of residues charged, so the
flux it delivered over an interval is read off the trajectory exactly. It wraps
rather than modifies `TrnaCharging`, so the module under test is not changed to
be testable. The counter reads the same flux function the module integrates.
"""
struct ChargingCounter{M <: AbstractSubModel} <: AbstractSubModel
    inner::M
end
states(m::ChargingCounter) = [states(m.inner); :chg_residues_mM]
parameters(m::ChargingCounter) = [parameters(m.inner);
    InferParameter(0.0, Normal(0.0, 0.1), true, :chg_residues_mM0,
                   module_id(m.inner), :initial_condition)]
module_id(m::ChargingCounter) = module_id(m.inner)
inputs(m::ChargingCounter) = inputs(m.inner)
contributed_states(m::ChargingCounter) = contributed_states(m.inner)
coupling(m::ChargingCounter) = coupling(m.inner)
reduction_notes(m::ChargingCounter) = reduction_notes(m.inner)
function dynamics(u, p, t, m::ChargingCounter, u_inputs)
    d = dynamics(SA[u[1], u[2]], p, t, m.inner, u_inputs)
    return SA[d[1], d[2], d[2]]           # the charged gain, accumulated
end
contributions(u, p, t, m::ChargingCounter, u_inputs) =
    contributions(SA[u[1], u[2]], p, t, m.inner, u_inputs)

"""
    charging_models(; charging, demand, reactions, counter = true)

The four-module composition every phase 9 full-cycle check runs: phase 8's
recycling and glycolytic double, the charging module (wrapped in a counter by
default), and the translation-demand double. `demand = nothing` drops the
consumer, which is the saturating configuration task 9.5 shows is vacuous.
"""
function charging_models(; charging = TrnaCharging(), demand = TranslationDemand(),
                         reactions = RECYCLING_REACTIONS, counter = true)
    ms = AbstractSubModel[NucleotideRecycling(reactions = reactions), HeldGlycolytic(),
                          counter ? ChargingCounter(charging) : charging]
    demand === nothing || push!(ms, demand)
    return ms
end

"Composed state index of `sp` in the composition `ms`, by name."
layout_index(ms, sp) = findfirst(==(sp), reduce(vcat, states.(ms)))

"""
    MutatedCharging(inner)

The real module with one coefficient wrong: each charging event makes a charged
tRNA without consuming an uncharged one, so it creates tRNA rather than
transferring it. Everything else delegates. A conservation check that cannot
fail is not evidence (spec §3), and this is the mutation task 9.6 names.
"""
struct MutatedCharging <: AbstractSubModel
    inner::TrnaCharging
end
states(m::MutatedCharging) = states(m.inner)
parameters(m::MutatedCharging) = parameters(m.inner)
module_id(::MutatedCharging) = :TrnaCharging
inputs(m::MutatedCharging) = inputs(m.inner)
contributed_states(m::MutatedCharging) = contributed_states(m.inner)
coupling(m::MutatedCharging) = coupling(m.inner)
contributions(u, p, t, m::MutatedCharging, u_inputs) =
    contributions(u, p, t, m.inner, u_inputs)
function dynamics(u, p, t, m::MutatedCharging, u_inputs)
    v = charging_flux(u, p, t, m.inner, u_inputs)
    return SA[zero(v), v]                  # the uncharged pool is never drawn
end

"The tRNA pair's sum in the composition `ms`, as a function of the state."
trna_total(ms) = (U = layout_index(ms, :M_trna_c); C = layout_index(ms, :M_trna_chg_c);
                  u -> u[U] + u[C])

"""
    TwoAtpCharging(; pool_mM, charged_fraction)

The alternative lumping the scoping note rejected, built for the one comparison
task 9.9 asks for:

    M_trna_c + 2 ATP  ->  M_trna_chg_c + 2 ADP + 2 Pi
    v = k2 · [M_trna_c] · [M_atp_c]²

Mass action in each unit of stoichiometry, so ATP appears squared. It spends
the same two high-energy phosphates per residue as the published form, but it
closes adenylate and phosphate without AMP or pyrophosphate ever appearing. So
it sends no traffic through the adenylate kinase or the pyrophosphatase, which
is the difference 9.9 measures.

**Matched on flux, not on event count** (spec §4 D14): `k2` is derived from the
same demand and nominal pools as `k_chg`, so both forms deliver 553.10
residues/s at nominal ATP, and what differs between them is only how each rate
law responds away from nominal. The demand is computed here, not read from
`TrnaCharging`.
"""
struct TwoAtpCharging <: AbstractSubModel
    params::Vector{InferParameter}
    k2::Float64
end
function TwoAtpCharging(; pool_mM = CHARGING_POOL_DEFAULTS.pool_mM,
                        charged_fraction = CHARGING_POOL_DEFAULTS.charged_fraction)
    atp = species_entry(:M_atp_c).initial_value
    k2 = TL_DEMAND_MM_PER_S / ((1 - charged_fraction) * pool_mM * atp^2)
    ics = [InferParameter(v, Normal(v, 0.1), true, Symbol(s, "0"), :TrnaCharging,
                          :initial_condition)
           for (s, v) in zip(CHARGING_STATES, ((1 - charged_fraction) * pool_mM,
                                               charged_fraction * pool_mM))]
    return TwoAtpCharging(ics, k2)
end
states(::TwoAtpCharging) = CHARGING_STATES
parameters(m::TwoAtpCharging) = m.params
module_id(::TwoAtpCharging) = :TrnaCharging
inputs(::TwoAtpCharging) = [:M_atp_c]
contributed_states(::TwoAtpCharging) = [:M_atp_c, :M_adp_c, :M_pi_c]
coupling(::TwoAtpCharging) = CouplingEdge[
    CurrencyEdge(species = :M_atp_c, direction = :in),
    CurrencyEdge(species = :M_adp_c, direction = :out),
    CurrencyEdge(species = :M_pi_c, direction = :out),
]
@inline _two_atp_rate(m::TwoAtpCharging, trna, atp) = m.k2 * trna * atp^2
function dynamics(u, p, t, m::TwoAtpCharging, u_inputs)
    v = _two_atp_rate(m, u[1], u_inputs[1])
    return SA[-v, v]
end
function contributions(u, p, t, m::TwoAtpCharging, u_inputs)
    v = _two_atp_rate(m, u[1], u_inputs[1])
    return SA[-2v, 2v, 2v]
end
