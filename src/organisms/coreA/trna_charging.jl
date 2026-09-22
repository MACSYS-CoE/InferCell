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

"""
    TrnaCharging(; k_chg, trna0, trna_chg0, free_k_chg = false)

The lumped charging step. `trna0` and `trna_chg0` are the initial uncharged and
charged pools in mM; the registry records no value for either (spec §4 D14).

`k_chg` is fixed by default (D14) and is then read from the struct rather than
the parameter vector; `free_k_chg = true` puts it in the vector so a sampler
moves it. The initial conditions are never freed: a freed initial condition is
sampled and then ignored (spec §12, 2026-09-10).
"""
struct TrnaCharging <: AbstractSubModel
    params::Vector{InferParameter}
    k_held::Float64
    k_slot::Int
end

function TrnaCharging(; k_chg::Real, trna0::Real, trna_chg0::Real,
                      free_k_chg::Bool = false)
    k_chg > 0 || throw(ArgumentError("k_chg must be positive, got $k_chg"))
    (trna0 >= 0 && trna_chg0 >= 0) || throw(ArgumentError(
        "tRNA initial conditions must be non-negative, got $trna0 and $trna_chg0"))
    params = InferParameter[
        InferParameter(k_chg, LogNormal(log(k_chg), log(2.0)), !free_k_chg,
                       CHARGING_K_ID, :TrnaCharging, :rate),
        InferParameter(trna0, Normal(trna0, 0.1), true,
                       Symbol(CHARGING_STATES[1], "0"), :TrnaCharging,
                       :initial_condition),
        InferParameter(trna_chg0, Normal(trna_chg0, 0.1), true,
                       Symbol(CHARGING_STATES[2], "0"), :TrnaCharging,
                       :initial_condition),
    ]
    # `k_chg` is this module's only kinetic parameter, so if it is free it is
    # the first entry of the module's free list.
    return TrnaCharging(params, Float64(k_chg), free_k_chg ? 1 : 0)
end

states(::TrnaCharging) = CHARGING_STATES
parameters(m::TrnaCharging) = m.params
inputs(::TrnaCharging) = [:M_atp_c]
contributed_states(::TrnaCharging) = [:M_atp_c, :M_amp_c, :M_ppi_c]
formalism(::TrnaCharging) = :ode

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

export TrnaCharging, CHARGING_STATES, charging_flux
