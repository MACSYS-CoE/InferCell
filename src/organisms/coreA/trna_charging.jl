# Core A′ lumped tRNA charging (spec §11 phase 9).
#
# One reaction, and the only one in the ODE block that is ours rather than the
# published model's:
#
#     M_trna_c + ATP  ->  M_trna_chg_c + AMP + PPi
#     v = k_chg · [M_trna_c] · [M_atp_c]
#
# It replaces the published model's twenty per-amino-acid chains of five
# reactions each, whose four constants are hard-coded inline with no priors to
# inherit (spec §3). It is integrated deterministically, in the ODE block, per
# §4 D13: the registry marks both tRNA species `:metabolite` regime, and the
# scoping note always placed the step on the ODE side.
#
# This module owns the tRNA pair. It reads ATP and contributes to ATP, AMP and
# pyrophosphate, which nucleotide recycling owns, through phase 1's contribution
# channel. Standalone it has no consumer of charged tRNA, so the pool saturates
# and the flux goes to zero; its acceptance test is therefore a composition with
# a translation-demand double (spec task 9.5).

"""
The two states this module integrates, derived from the registry's `:trna`
group so membership and order are the registry's (`states` must increase in
`species_index`).
"""
const CHARGING_STATES = species_in_group(:trna)

"The rate constant's parameter name. Its only kinetic parameter."
const CHARGING_K_ID = :k_chg

# ----------------------------------------------------------------------
# The demand, and the derivation of k_chg from it (spec §4 D14).
#
# The residue count and the cycle length are written here, in this file, and
# nowhere imported from. Phase 8's drain double carries the same demand as
# `RECYCLING_DRAIN_*`; reading it from there would let this module calibrate
# against the number it is then checked against (spec §12, 2026-09-11).
# ----------------------------------------------------------------------

"Residues charged over one cycle: a full proteome doubling of Core A′'s loci."
const CHARGING_RESIDUES_PER_CYCLE = 3_484_518
"The published cell cycle, in seconds."
const CHARGING_CYCLE_S = 6300.0

"""
    charging_demand_per_s() -> Float64

The published residue demand, 3,484,518 / 6,300 = 553.10 residues per second.
"""
charging_demand_per_s() = CHARGING_RESIDUES_PER_CYCLE / CHARGING_CYCLE_S

"""
    charging_demand_mM_per_s() -> Float64

[`charging_demand_per_s`](@ref) as a concentration flux, at the registry's
volume through [`corea_particles_per_mM`](@ref) rather than a transcribed
20,180 (spec §12, 2026-09-10): 0.027408 mM/s.
"""
charging_demand_mM_per_s() = charging_demand_per_s() / corea_particles_per_mM()

"""
Asserted defaults for the tRNA pool (spec §12, 2026-09-23). Neither comes from
the published model, which has twenty per-amino-acid pools and no lumped one.

- `pool_mM = 0.25` — the total, charged plus uncharged: 5,045 particles, the
  order of bacterial tRNA concentrations scaled to Syn3A's volume.
- `charged_fraction = 0.8` — a typical bacterial charged fraction.
"""
const CHARGING_POOL_DEFAULTS = (pool_mM = 0.25, charged_fraction = 0.8)

"""
    derive_k_chg(; pool_mM, charged_fraction, atp_mM) -> Float64

The rate constant that makes `k_chg·[M_trna_c]·[M_atp_c]` equal the published
demand at the nominal pools, where `[M_trna_c] = (1 − charged_fraction)·pool_mM`:

    k_chg = charging_demand_mM_per_s() / ((1 − charged_fraction) · pool_mM · atp_mM)

At the defaults and the registry's 3.6529 mM ATP this is 0.15006 /(mM·s).
`k_chg` is the **derived** quantity here (D14): change the pool and it follows.
"""
function derive_k_chg(; pool_mM::Real = CHARGING_POOL_DEFAULTS.pool_mM,
                      charged_fraction::Real = CHARGING_POOL_DEFAULTS.charged_fraction,
                      atp_mM::Real = species_entry(:M_atp_c).initial_value)
    pool_mM > 0 || throw(ArgumentError("pool_mM must be positive, got $pool_mM"))
    0 <= charged_fraction < 1 || throw(ArgumentError(
        "charged_fraction must lie in [0, 1), got $charged_fraction; at 1 there " *
        "is no uncharged tRNA to charge and no k_chg meets the demand"))
    return charging_demand_mM_per_s() / ((1 - charged_fraction) * pool_mM * atp_mM)
end

_asserted(identifier) = ParameterSource("asserted: spec §4 D14";
                                              table = "tRNA charging",
                                              identifier = identifier,
                                              informedness = :asserted)

"""
    TrnaCharging(; pool_mM, charged_fraction, k_chg = nothing, free_k_chg = false)

The lumped charging step. The initial uncharged and charged pools are
`(1 − charged_fraction)·pool_mM` and `charged_fraction·pool_mM`; the registry
records no value for either (spec §4 D14).

`k_chg = nothing`, the default, derives it through [`derive_k_chg`](@ref). A
number overrides the derivation, which is for doubles and hand-checked tests,
and the override is recorded in the parameter's provenance rather than hidden.

The pool size, the charged fraction and `k_chg` are all `:asserted` parameters,
so they reach [`reduction_declarations`](@ref)'s asserted-prior enumeration.
Only `k_chg` can be freed, through `free_k_chg = true`, which moves it from the
struct into the parameter vector. The pool quantities set only the initial
conditions, and a freed initial condition is sampled and then ignored (spec
§12, 2026-09-10), so they are always fixed.
"""
struct TrnaCharging <: AbstractSubModel
    params::Vector{InferParameter}
    k_held::Float64
    k_slot::Int
end

function TrnaCharging(; pool_mM::Real = CHARGING_POOL_DEFAULTS.pool_mM,
                      charged_fraction::Real = CHARGING_POOL_DEFAULTS.charged_fraction,
                      k_chg::Union{Real, Nothing} = nothing,
                      free_k_chg::Bool = false)
    derived = derive_k_chg(; pool_mM, charged_fraction)
    k = k_chg === nothing ? derived : Float64(k_chg)
    k > 0 || throw(ArgumentError("k_chg must be positive, got $k"))
    k_note = k_chg === nothing ? "derived from demand and pool (D14)" :
                                 "override of the D14 derivation ($(derived))"
    trna0 = (1 - charged_fraction) * pool_mM
    trna_chg0 = charged_fraction * pool_mM
    params = InferParameter[
        InferParameter(k, LogNormal(log(k), log(2.0)), !free_k_chg,
                       CHARGING_K_ID, :TrnaCharging, :rate,
                       _asserted("k_chg: " * k_note)),
        InferParameter(pool_mM, LogNormal(log(pool_mM), log(2.0)), true,
                       :trna_pool_mM, :TrnaCharging, :rate,
                       _asserted("total tRNA pool")),
        InferParameter(charged_fraction, Normal(charged_fraction, 0.1), true,
                       :trna_charged_fraction, :TrnaCharging, :rate,
                       _asserted("nominal charged fraction")),
        InferParameter(trna0, Normal(trna0, 0.1), true,
                       Symbol(CHARGING_STATES[1], "0"), :TrnaCharging,
                       :initial_condition),
        InferParameter(trna_chg0, Normal(trna_chg0, 0.1), true,
                       Symbol(CHARGING_STATES[2], "0"), :TrnaCharging,
                       :initial_condition),
    ]
    # `k_chg` is this module's only freeable parameter, so if it is free it is
    # the first entry of the module's free list.
    return TrnaCharging(params, k, free_k_chg ? 1 : 0)
end

"""
    charging_derivation(m::TrnaCharging) -> NamedTuple

What `k_chg` was derived from, so a report can show the chain rather than the
number: demand, pool, split, the ATP it was calibrated at, and the result.
"""
function charging_derivation(m::TrnaCharging)
    val(n) = only(p.value for p in m.params if p.name === n)
    return (demand_per_s = charging_demand_per_s(),
            demand_mM_per_s = charging_demand_mM_per_s(),
            pool_mM = val(:trna_pool_mM),
            charged_fraction = val(:trna_charged_fraction),
            atp_mM = species_entry(:M_atp_c).initial_value,
            k_chg = m.k_held)
end

"""
The boundary: three currency edges on pools nucleotide recycling owns. The
charging step is the principal consumer of ATP's turnover (~81%, spec §4 D13)
and a principal producer of AMP and pyrophosphate, so it declares all three. As
a jump module it could not have executed these against ODE states (D13); in the
ODE block they run through phase 1's contribution channel.
"""
const CHARGING_EDGES = CouplingEdge[
    CurrencyEdge(species = :M_atp_c, direction = :in),
    CurrencyEdge(species = :M_amp_c, direction = :out),
    CurrencyEdge(species = :M_ppi_c, direction = :out),
]

states(::TrnaCharging) = CHARGING_STATES
coupling(::TrnaCharging) = CHARGING_EDGES
parameters(m::TrnaCharging) = m.params
inputs(::TrnaCharging) = [:M_atp_c]
contributed_states(::TrnaCharging) = [:M_atp_c, :M_amp_c, :M_ppi_c]
formalism(::TrnaCharging) = :ode

"""
    reduction_notes(m::TrnaCharging)

Two declarations, not one, because two separate things here are ours (spec task
9.4):

1. **The lumping.** One reaction replaces the published model's twenty
   per-amino-acid chains of five reactions each. It keeps the published product
   stoichiometry, AMP plus pyrophosphate, over a `trna + 2 ATP -> trna_chg +
   2 ADP + 2 Pi` form. The latter closes the same moieties by fiat and adds two
   more entries to the list of what is ours, and it would put no traffic through
   the adenylate kinase that the published form sends ~553/s through.
2. **The formalism.** Integrating the step deterministically is ours, even
   though its placement in the ODE block is the scoping note's (§4 D13). The
   published model fires it as ~3.49 million stochastic events per cycle.

Both land under `:lumping`, the only category `reduction_notes` feeds. A
separate `:formalism` category would be a change to `src/labels.jl`, which
§10 R15 freezes against a module branch.
"""
reduction_notes(::TrnaCharging) = [
    "tRNA charging lumped: one reaction M_trna_c + ATP -> M_trna_chg_c + AMP + " *
    "PPi replaces the published 20 per-amino-acid chains of 5 reactions each (100 " *
    "reactions), with one effective tRNA pool whose size and charged fraction are " *
    "asserted by this project. The published AMP + PPi product stoichiometry is " *
    "kept in preference to a 2 ATP -> 2 ADP + 2 Pi form, which would close the " *
    "same moieties by fiat and route no traffic through the adenylate kinase",
    "tRNA charging integrated deterministically in the ODE block: the placement " *
    "is the scoping note's, but the formalism is this project's; the published " *
    "model fires charging as ~3.49 million stochastic events per cycle (spec §4 D13)",
]

# From the parameter vector where free, from the struct where fixed, converted
# to the vector's element type so the right-hand side stays inferable under
# ForwardDiff. The same accessor pattern as `NucleotideRecycling`'s `_k`.
@inline _k_chg(p, m::TrnaCharging) =
    m.k_slot == 0 ? convert(eltype(p), m.k_held) : @inbounds(p[m.k_slot])

"""
    charging_flux(u, p, t, m, u_inputs) -> Real

`k_chg · [M_trna_c] · [M_atp_c]`, in mM/s.
"""
@inline charging_flux(u, p, t, m::TrnaCharging, u_inputs) =
    _k_chg(p, m) * u[1] * u_inputs[1]

# One uncharged tRNA becomes one charged tRNA.
function dynamics(u, p, t, m::TrnaCharging, u_inputs)
    v = charging_flux(u, p, t, m, u_inputs)
    return SA[-v, v]
end

# One ATP consumed; one AMP and one pyrophosphate produced. Three phosphates in,
# one plus two out, so the transfer closes adenylate and phosphate internally
# (spec §3, checks 4 and 4b).
function contributions(u, p, t, m::TrnaCharging, u_inputs)
    v = charging_flux(u, p, t, m, u_inputs)
    return SA[-v, v, v]
end

export TrnaCharging, CHARGING_STATES, CHARGING_EDGES, CHARGING_RESIDUES_PER_CYCLE, CHARGING_CYCLE_S,
       CHARGING_POOL_DEFAULTS, charging_flux, charging_demand_per_s,
       charging_demand_mM_per_s, derive_k_chg, charging_derivation
