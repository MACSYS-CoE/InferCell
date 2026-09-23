using InferCell
using StaticArrays: SA, SVector
using Distributions: LogNormal, Normal
using OrdinaryDiffEq: solve, Rodas5P, Tsit5

import InferCell: states, parameters, dynamics, contributions, contributed_states,
                  coupling, inputs, module_id, _recycling_ddt, RECYCLING_INPUTS

# Executable doubles for spec §11 phase 8.
#
# Nucleotide recycling owns the pools every other module routes energy through,
# so standalone it conserves both moieties trivially: nothing consumes them.
# Its acceptance criterion is therefore a composition, and these two doubles are
# what stand in for the modules that are not built yet — the lumped charging
# step of phase 9, and central glycolysis of phase 6.

const RECYCLING_DRAIN_PER_S = 553.1        # residues charged per second
"""
The published charging demand as a concentration flux: 3,484,518 residues over a
6,300 s cycle is 553.1/s, and 553.1 particles per second at the registry's
20,180.39 particles per millimolar is 0.0274 mM/s. Spec §11 task 8.8 records
this as the *demand* phase 9's real module must meet, never as a value phase 9
may calibrate against and then re-check.

The conversion is [`corea_particles_per_mM`](@ref), not a transcribed 20,180,
for the reason spec §12's 2026-09-10 entry gives: a rounded copy puts a module
out of step with the conversion the handshake itself performs.
"""
const RECYCLING_DRAIN_MM_PER_S = RECYCLING_DRAIN_PER_S / corea_particles_per_mM()

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

**Superseded for composed runs by phase 9's `TrnaCharging`** (with
`TranslationDemand` in `trna_test_models.jl`), and kept for phase 8's
standalone checks, whose recorded numbers belong to this drain.
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
# The four species recycling reads, in the order `states` requires (increasing
# `species_index`), plus the counter. Derived from the module's own input list
# so the two cannot drift.
const HELD_GLYCOLYTIC_STATES = [sort(RECYCLING_INPUTS, by = species_index);
                                :slp_phosphorylated_mM]

# Read from the registry rather than transcribed, so a setpoint cannot drift
# from the initial condition it is supposed to be.
const HELD_GLYCOLYTIC_VALUES =
    Float64[[species_entry(s).initial_value for s in HELD_GLYCOLYTIC_STATES[1:4]]; 0.0]

"""
    HeldGlycolytic(; k_slp)

The part of central glycolysis nucleotide recycling cannot run without, and
nothing more. Two jobs:

1. **It owns the four glycolytic species PGK3 and PYK3 read**, and holds them
   near their registry values by resupplying what recycling draws down and
   consuming what it makes. Without an owner the composition does not build at
   all, because `inputs` resolves against integrated states.

   **A zero derivative does not hold them, and a first version assumed it did.**
   `NucleotideRecycling` names all four in `contributed_states`, and
   `_accumulate` folds a contribution into the *owner's* derivative — so with
   `dynamics` returning zeros the pools moved at exactly ±v_PGK3 and ±v_PYK3.
   13DPG fell below 1% of its initial value at **t = 0.5 s** and PEP at
   **t = 18.5 s**, and the GTP branch then died of substrate starvation inside
   the first save interval, in a run whose stated purpose is to show that PGK3
   and PYK3 make GTP a live product of glycolysis. Each pool now relaxes to its
   setpoint at `k_gly`, which is the four lumped reactions phase 6 supplies:
   GAPD resupplies 13DPG, PGM and ENO carry 3PG to PEP, and LDH drains
   pyruvate. None of the four is in the tracked phosphate sum, so the resupply
   adds nothing to the closed moiety; the only tracked crossing is still the
   one phosphate per GTP that check 4b's correction subtracts.

   **What now bounds the branch is guanylate, which is the honest limit.**
   Nothing in Core A′ consumes GTP until phase 9's translation, so PGK3 and
   PYK3 run until GDP is spent and then stop. That is a structural statement
   about the composition rather than an artefact of the double.

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
integrating through zero. `k_gly` is fast against the GTP branch, so the pools
sit near their setpoints rather than tracking them exactly; the suite asserts
how near.
"""
struct HeldGlycolytic <: AbstractSubModel
    params::Vector{InferParameter}
    k_slp::Float64
    k_half::Float64
    k_gly::Float64
end
function HeldGlycolytic(; k_slp = RECYCLING_DRAIN_MM_PER_S / species_entry(:M_adp_c).initial_value,
                        k_half = 1e-3, k_gly = 1.0)
    params = [InferParameter(v, Normal(v, 0.1), true, Symbol(s, "0"),
                             :HeldGlycolytic, :initial_condition)
              for (s, v) in zip(HELD_GLYCOLYTIC_STATES, HELD_GLYCOLYTIC_VALUES)]
    return HeldGlycolytic(params, k_slp, k_half, k_gly)
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

@inline _slp_rate(m::HeldGlycolytic, adp, pin) =
    m.k_slp * adp * pin / (m.k_half + pin)

# Relaxation to the registry setpoint, one lumped reaction per pool. This has to
# be a derivative term and not a clamp, because the term it is cancelling is
# `NucleotideRecycling`'s contribution, which arrives in this module's `du`.
function dynamics(u, p, t, m::HeldGlycolytic, u_inputs)
    k, sp = m.k_gly, HELD_GLYCOLYTIC_VALUES
    return SA[k * (sp[1] - u[1]), k * (sp[2] - u[2]),
              k * (sp[3] - u[3]), k * (sp[4] - u[4]),
              _slp_rate(m, u_inputs[1], u_inputs[2])]
end
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
conservation test that cannot fail is not evidence (spec §3), so each mutation
names a moiety and drives that moiety's residual to pool scale.

**Two of the three also break phosphate closure, and that is arithmetic rather
than a leaky test.** The phosphate sum spans both nucleotide groups, so a
coefficient wrong on a phosphorylated species shows up there too: returning one
ADP where two are made loses two phosphates per turnover, and creating a GDP
creates two. Only the named moiety is claimed to localise; the third column
records what actually moves, measured over 600 s.

| `mutation` | What it does | Named moiety | Also breaks |
|---|---|---|---|
| `:adk1_adp_coefficient` | the kinase returns one ADP where the reaction makes two | adenylate, 3.70 mM | phosphate, 7.40 mM |
| `:gk1_gdp_created` | guanylate kinase creates a GDP rather than transferring | guanylate, 9.89e-3 mM | phosphate, 1.98e-2 mM |
| `:ppa_phosphate_coefficient` | pyrophosphate hydrolyses to one phosphate, not two | phosphate, 16.5 mM | — |
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

dynamics(u, p, t, m::MutatedRecycling, u_inputs) = _recycling_ddt(
    recycling_fluxes(u, p, t, m.inner, u_inputs);
    adp_per_adk1 = m.mutation === :adk1_adp_coefficient ? 1 : 2,
    gdp_from_gk1 = m.mutation === :gk1_gdp_created ? 2 : 1,
    pi_per_ppa   = m.mutation === :ppa_phosphate_coefficient ? 1 : 2)


# ----------------------------------------------------------------------
# The composition, and the quantities read off it.
#
# These live here rather than in either caller because BOTH callers read them:
# `test/test_corea_nucleotide_recycling.jl` asserts them and
# `dev/scripts/full_cycle_recycling.jl` records them in the results artefact.
# Written twice they were two texts for one measurement, and a moiety weight or
# a state index drifting between them would have made the recorded number and
# the asserted number different quantities while both still passed.
# ----------------------------------------------------------------------

const CYCLE_S = 6300.0
const SAVE_EVERY = 60.0
const ABSTOL_R, RELTOL_R = 1e-10, 1e-8

"""
    recycling_models(; reactions, mutation) -> Vector{AbstractSubModel}

The three-module composition every full-cycle check runs: the recycling module
(or its mutant), the glycolytic double, and the charging drain.
"""
# `k_slp === nothing` means "whatever `HeldGlycolytic` derives", so the default
# is stated in exactly one place. Task 8.7's exact form passes 0.0 to switch the
# substrate-level phosphorylation off.
recycling_models(; reactions = RECYCLING_REACTIONS, mutation = nothing,
                 k_slp = nothing) =
    AbstractSubModel[
        mutation === nothing ? NucleotideRecycling(reactions = reactions) :
            MutatedRecycling(mutation; reactions = reactions),
        k_slp === nothing ? HeldGlycolytic() : HeldGlycolytic(k_slp = k_slp),
        ChargingDrain(),
    ]

"""
    recycling_solve(ms; horizon, abstol, reltol, saveat, alg)

`alg` is a keyword because the mutation runs do not need a stiff solver, and
letting them use `Tsit5` avoids a second `Rodas5P` specialisation that costs
about fifteen seconds of suite compilation for values that agree to ten
significant figures.
"""
recycling_solve(ms; horizon = CYCLE_S, abstol = ABSTOL_R, reltol = RELTOL_R,
                saveat = SAVE_EVERY, alg = Rodas5P()) =
    solve(build_problem(ms; tspan = (0.0, horizon)), alg;
          abstol = abstol, reltol = reltol, saveat = saveat)

# Composed state indices, resolved by name off the composition itself rather
# than written as integers in two files that must agree.
const RECYCLING_LAYOUT = reduce(vcat, states.(recycling_models()))
_ridx(sp) = findfirst(==(sp), RECYCLING_LAYOUT)

const ATP_I, ADP_I, AMP_I, PI_I, GTP_I, GDP_I, GMP_I, PPI_I =
    _ridx.((:M_atp_c, :M_adp_c, :M_amp_c, :M_pi_c,
            :M_gtp_c, :M_gdp_c, :M_gmp_c, :M_ppi_c))
const DPG_I, PG3_I, PEP_I, PYR_I =
    _ridx.((:M_13dpg_c, :M_3pg_c, :M_pep_c, :M_pyr_c))
const SLP_CUM_I = _ridx(:slp_phosphorylated_mM)
const DRAIN_CUM_I = _ridx(:chg_charged_mM)

adenylate_of(u) = u[ATP_I] + u[ADP_I] + u[AMP_I]
guanylate_of(u) = u[GTP_I] + u[GDP_I] + u[GMP_I]
phosphate_of(u) = u[PI_I] + 3u[ATP_I] + 2u[ADP_I] + u[AMP_I] +
                  3u[GTP_I] + 2u[GDP_I] + u[GMP_I] + 2u[PPI_I]

max_drift(sol, f) = maximum(abs(f(u) - f(sol.u[1])) for u in sol.u)

"""
Phosphate crossing the two inbound mass edges. Standalone nothing consumes GTP,
so the crossing is exactly the GTP gain — a state difference, not a quadrature,
which is what keeps the corrected closure a linear functional of the state.
"""
inbound_phosphate(u, u0) = u[GTP_I] - u0[GTP_I]

corrected_phosphate_drift(sol) =
    maximum(abs((phosphate_of(u) - inbound_phosphate(u, sol.u[1])) - phosphate_of(sol.u[1]))
            for u in sol.u)

"""
    tolerance_bound(sol, coeffs; abstol, reltol)

`tol_C` for one invariant: `Σ|nᵢ| · max(abstol, reltol · maxₜ|xᵢ|)`, spec §3.
"""
tolerance_bound(sol, coeffs; abstol = ABSTOL_R, reltol = RELTOL_R) =
    sum(abs(n) * max(abstol, reltol * maximum(abs(u[i]) for u in sol.u))
        for (i, n) in coeffs)

const ADENYLATE_COEFFS = ((ATP_I, 1), (ADP_I, 1), (AMP_I, 1))

"First save point at which ATP has fallen below `frac` of its initial value."
function crossing_time(sol, frac)
    target = frac * sol.u[1][ATP_I]
    i = findfirst(u -> u[ATP_I] < target, sol.u)
    return i === nothing ? nothing : sol.t[i]
end

"The five reactions with `r` removed, for the two knock-out configurations."
without(r) = Tuple(x for x in RECYCLING_REACTIONS if x !== r)
