using InferCell
using StaticArrays: SA, SVector
using Distributions: LogNormal, Normal

import InferCell: states, parameters, dynamics, contributions, contributed_states,
                  coupling, inputs, module_id

# Executable doubles for spec §11 phase 8.
#
# Nucleotide recycling owns the pools every other module routes energy through,
# so standalone it conserves both moieties trivially: nothing consumes them.
# Its acceptance criterion is therefore a composition, and these two doubles are
# what stand in for the modules that are not built yet — the lumped charging
# step of phase 9, and central glycolysis of phase 6.

const RECYCLING_DRAIN_PER_S = 553.1        # residues charged per second
const RECYCLING_PARTICLES_PER_MM_TEST = 20180
"""
The published charging demand as a concentration flux: 3,484,518 residues over a
6,300 s cycle is 553.1/s, and 553.1 particles per second at 20,180 particles per
millimolar is 0.0274 mM/s. Spec §11 task 8.8 records this as the *demand* phase
9's real module must meet, never as a value phase 9 may calibrate against and
then re-check.
"""
const RECYCLING_DRAIN_MM_PER_S = RECYCLING_DRAIN_PER_S / RECYCLING_PARTICLES_PER_MM_TEST

_drain_ic(name, value) = InferParameter(value, Normal(value, 0.1), true,
                                        Symbol(name, "0"), :ChargingDrain,
                                        :initial_condition)

"""
    ChargingDrain(; rate_mm_per_s, k_half)

The lumped tRNA charging step of spec §11 phase 9, reduced to what phase 8 needs
from it: `ATP -> AMP + PPi` at the published residue demand. It owns one
non-registry state counting the millimolar of ATP it has charged, so a test can
read the delivered flux off the trajectory rather than assume it.

**Saturating rather than strictly constant, and this is a deliberate
departure.** §4 D14 observes that a constant drain drives ATP to zero in finite
time; a configuration with the adenylate kinase removed would then integrate
into negative concentrations and the 6,300 s run could not complete at all. So
the drain carries a Michaelis factor with `k_half` three orders of magnitude
below the nominal ATP pool: at 3.6529 mM it delivers 99.97% of the demand, and
it goes to zero with ATP instead of through it. The threshold spec task 8.6
asserts — ATP below 1% of its initial value — is crossed at 0.0365 mM, where the
factor is still 0.97, so the assertion is made where the drain is effectively
the published constant one. Phase 9's real module is mass-action in both
reactants and needs no such device.
"""
struct ChargingDrain <: AbstractSubModel
    params::Vector{InferParameter}
    rate::Float64
    k_half::Float64
end
function ChargingDrain(; rate_mm_per_s = RECYCLING_DRAIN_MM_PER_S, k_half = 1e-3)
    return ChargingDrain([_drain_ic(:chg_charged_mM, 0.0)], rate_mm_per_s, k_half)
end
states(::ChargingDrain) = [:chg_charged_mM]
parameters(m::ChargingDrain) = m.params
module_id(::ChargingDrain) = :ChargingDrain
inputs(::ChargingDrain) = [:M_atp_c]
contributed_states(::ChargingDrain) = [:M_atp_c, :M_amp_c, :M_ppi_c]
coupling(::ChargingDrain) = CouplingEdge[
    CurrencyEdge(species = :M_atp_c, direction = :in),
    CurrencyEdge(species = :M_amp_c, direction = :out),
    CurrencyEdge(species = :M_ppi_c, direction = :out),
]

@inline _drain_rate(m::ChargingDrain, atp) = m.rate * atp / (m.k_half + atp)

dynamics(u, p, t, m::ChargingDrain, u_inputs) = SA[_drain_rate(m, u_inputs[1])]
function contributions(u, p, t, m::ChargingDrain, u_inputs)
    v = _drain_rate(m, u_inputs[1])
    # One ATP consumed, one AMP and one pyrophosphate produced: three phosphates
    # in, one plus two out, so the transfer closes phosphate internally and
    # needs no correction in the closure check (spec §3, check 4b).
    return SA[-v, v, v]
end

# The four held species at their registry values, then a cumulative counter.
# The counter is not decoration: the phosphate closure check needs the *integral*
# of the phosphorylation flux, and integrating it as a state makes that number
# exact where quadrature over 60 s save points would not be.
const HELD_GLYCOLYTIC_STATES = [:M_13dpg_c, :M_3pg_c, :M_pep_c, :M_pyr_c,
                                :slp_phosphorylated_mM]
const HELD_GLYCOLYTIC_VALUES = [0.0098, 1.1015, 0.0409, 3.3660, 0.0]   # registry.jl

"""
    HeldGlycolytic(; k_slp)

The part of central glycolysis nucleotide recycling cannot run without, and
nothing more. Two jobs:

1. **It owns the four glycolytic species PGK3 and PYK3 read**, held at their
   registry values with a zero derivative. Without an owner the composition does
   not build at all, because `inputs` resolves against integrated states.

2. **It rephosphorylates ADP as `ADP + Pi -> ATP`**, which is the substrate-level
   phosphorylation phase 6's GAPD and PGK will supply. This is not decoration.
   Nothing in the recycling module produces ATP — the adenylate kinase converts
   AMP back to ADP *at the cost of an ATP* — so without a source ATP falls to
   zero in about 130 s whatever else is true, and spec task 8.6's "ATP positive"
   could not hold for any configuration. Worse, the kinase-removed run's
   threshold crossing would then be dominated by ATP draining into ADP rather
   than by adenylate stranding as AMP, which is the mechanism the scoping note's
   144 s actually describes.

**The phosphate comes from the free pool, and getting that wrong is instructive.**
A first version drew it from the held 13DPG pool, which is where the published
PGK reaction takes it from — but 13DPG is held here, so phosphate entered the
composition and nothing removed it: 345 mM of it over a cycle, free phosphate
climbing without bound, and pyrophosphate riding up with it to 51 mM because
PPA's reverse term grows as the square of the phosphate it produces. That is an
artefact of holding 13DPG, not a property of the module. Glycolysis from G3P is
`G3P + Pi + ADP -> 3PG + ATP` once GAPD and PGK are composed, so the lumped
stand-in takes the phosphate from `M_pi_c`. Then the cycle closes exactly: each
charging turnover frees two phosphates through the pyrophosphatase and the two
ATP the kinase arithmetic needs consume exactly two.

`k_slp` defaults to the value that matches the charging demand at the registry's
initial ADP. The Michaelis factor in phosphate is there for the same reason as
the drain's: so a configuration that runs the pool down stops rather than
integrating through zero.
"""
struct HeldGlycolytic <: AbstractSubModel
    params::Vector{InferParameter}
    k_slp::Float64
    k_half::Float64
end
function HeldGlycolytic(; k_slp = RECYCLING_DRAIN_MM_PER_S / 0.2178, k_half = 1e-3)
    params = [InferParameter(v, Normal(v, 0.1), true, Symbol(s, "0"),
                             :HeldGlycolytic, :initial_condition)
              for (s, v) in zip(HELD_GLYCOLYTIC_STATES, HELD_GLYCOLYTIC_VALUES)]
    return HeldGlycolytic(params, k_slp, k_half)
end
states(::HeldGlycolytic) = HELD_GLYCOLYTIC_STATES
parameters(m::HeldGlycolytic) = m.params
module_id(::HeldGlycolytic) = :HeldGlycolytic
inputs(::HeldGlycolytic) = [:M_adp_c, :M_pi_c]
contributed_states(::HeldGlycolytic) = [:M_atp_c, :M_adp_c, :M_pi_c]
coupling(::HeldGlycolytic) = CouplingEdge[
    CurrencyEdge(species = :M_atp_c, direction = :out),
    CurrencyEdge(species = :M_adp_c, direction = :in),
    CurrencyEdge(species = :M_pi_c, direction = :in),
]

# The four species it owns are held: a clamp in all but name, and the reason
# every number this composition produces is the recycling module's rather than a
# stand-in glycolysis's.
@inline _slp_rate(m::HeldGlycolytic, adp, pin) =
    m.k_slp * adp * pin / (m.k_half + pin)

dynamics(u, p, t, m::HeldGlycolytic, u_inputs) =
    SA[0.0, 0.0, 0.0, 0.0, _slp_rate(m, u_inputs[1], u_inputs[2])]
# Declared unconditionally, so `k_slp = 0` is the same composition running at
# zero rate rather than a different one — which is what lets the phosphate
# check's exact and flux-corrected forms differ in one number and nothing else.
function contributions(u, p, t, m::HeldGlycolytic, u_inputs)
    v = _slp_rate(m, u_inputs[1], u_inputs[2])
    return SA[v, -v, -v]
end

"""
    MutatedRecycling(inner, mutation)

The real module with one stoichiometric coefficient wrong, and nothing else
changed: states, parameters, edges, inputs and contributions all delegate. A
conservation test that cannot fail is not evidence (spec §3), so each of the
three mutations breaks exactly one check and leaves the others passing, which is
what says a failure localises.

| `mutation` | What it does | Which check should fail |
|---|---|---|
| `:adk1_adp_coefficient` | the kinase returns one ADP where the reaction makes two | adenylate |
| `:gk1_gdp_created` | guanylate kinase creates a GDP rather than transferring | guanylate |
| `:ppa_phosphate_coefficient` | pyrophosphate hydrolyses to one phosphate, not two | phosphate closure |
"""
struct MutatedRecycling <: AbstractSubModel
    inner::NucleotideRecycling
    mutation::Symbol
end
const RECYCLING_MUTATIONS = (:adk1_adp_coefficient, :gk1_gdp_created,
                             :ppa_phosphate_coefficient)
function MutatedRecycling(mutation::Symbol; kwargs...)
    mutation in RECYCLING_MUTATIONS || throw(ArgumentError(
        "unknown mutation :$mutation; try $(join(string.(":", RECYCLING_MUTATIONS), ", "))"))
    return MutatedRecycling(NucleotideRecycling(; kwargs...), mutation)
end
states(m::MutatedRecycling) = states(m.inner)
parameters(m::MutatedRecycling) = parameters(m.inner)
coupling(m::MutatedRecycling) = coupling(m.inner)
inputs(m::MutatedRecycling) = inputs(m.inner)
contributed_states(m::MutatedRecycling) = contributed_states(m.inner)
module_id(::MutatedRecycling) = :NucleotideRecycling
contributions(u, p, t, m::MutatedRecycling, u_inputs) =
    contributions(u, p, t, m.inner, u_inputs)

function dynamics(u, p, t, m::MutatedRecycling, u_inputs)
    v = recycling_fluxes(u, p, 0.0, m.inner, u_inputs)
    pgk3, pyk3, adk1, gk1, ppa = v
    adp_per_kinase = m.mutation === :adk1_adp_coefficient ? 1 : 2
    gdp_from_gk1 = m.mutation === :gk1_gdp_created ? 2 : 1
    pi_per_ppa = m.mutation === :ppa_phosphate_coefficient ? 1 : 2
    return SA[
        -adk1 - gk1,
        adp_per_kinase * adk1 + gk1,
        -adk1,
        pi_per_ppa * ppa,
        pgk3 + pyk3,
        -pgk3 - pyk3 + gdp_from_gk1 * gk1,
        -gk1,
        -ppa,
    ]
end
