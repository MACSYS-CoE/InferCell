using InferCell
using StaticArrays: SA, SVector
using Distributions: LogNormal, Normal

import InferCell: states, parameters, dynamics, contributions, contributed_states,
                  coupling, inputs, module_id

# Executable doubles for the contribution channel (spec §11 phase 1).
#
# `CoreAStub` implements no dynamics, so it cannot show a contribution being
# *executed*. These four do. They use real registry species names because the
# channel is registry-resolved, but their rate laws are placeholders — linear
# decay and constant or linear sources — chosen so every expected derivative
# can be written down in the test and compared bit for bit.
#
# Between them `EnergyPools` and `CarbonBlock` own all 32 dynamic registry
# species, which is the composition size the allocation assertion is made at.

const ENERGY_SPECIES = [:M_atp_c, :M_adp_c, :M_amp_c, :M_pi_c,
                        :M_gtp_c, :M_gdp_c, :M_gmp_c,
                        :M_nad_c, :M_nadh_c, :M_ppi_c, :M_lac__L_e,
                        :M_trna_c, :M_trna_chg_c]
const CARBON_SPECIES = [:M_g6p_c, :M_f6p_c, :M_fdp_c, :M_dhap_c, :M_g3p_c, :M_13dpg_c,
                        :M_3pg_c, :M_2pg_c, :M_pep_c, :M_pyr_c, :M_lac__L_c,
                        :M_ptsi_c, :M_ptsi_P_c, :M_ptsh_c, :M_ptsh_P_c,
                        :M_crr_c, :M_crr_P_c, :M_ptsg_c, :M_ptsg_P_c]

_ic_params(species, id) = [InferParameter(1.0 + i, Normal(1.0 + i, 0.1), true,
                                          Symbol(s, "0"), id, :initial_condition)
                           for (i, s) in enumerate(species)]

"""
    EnergyPools(; gamma = 0.05)

Owns the thirteen energy, redox, external-lactate and tRNA species. Every state
decays at one rate `gamma_o`. It reads nothing and contributes nothing: it is
the pool other modules write into.
"""
struct EnergyPools <: AbstractSubModel
    params::Vector{InferParameter}
end
function EnergyPools(; gamma = 0.05)
    params = vcat(InferParameter[InferParameter(gamma, LogNormal(log(gamma), 0.5), false, :gamma_o, :EnergyPools, :rate)],
                  _ic_params(ENERGY_SPECIES, :EnergyPools))
    return EnergyPools(params)
end
states(::EnergyPools) = ENERGY_SPECIES
parameters(m::EnergyPools) = m.params
dynamics(u, p, t, ::EnergyPools) = -p[1] .* u

"""
    CarbonBlock(; k_prod, k_cons, gamma, contribs, edges, ins)

Owns the nineteen glycolytic and phosphotransferase species. Every state decays
at `gamma_c`. It produces ATP at `k_prod · [g6p]` through an outbound currency
edge and draws phosphate at `k_cons · [pi]` through an inbound mass edge, so it
exercises both directions of the channel: a positive contribution to a pool it
reads nothing of, and a negative one to a pool it reads through `inputs`.

The declarations are keyword-overridable so the drift and rejection tests can
build a module whose contributions disagree with its edges.
"""
struct CarbonBlock <: AbstractSubModel
    params::Vector{InferParameter}
    contribs::Vector{Symbol}
    edges::Vector{CouplingEdge}
    ins::Vector{Symbol}
end
function CarbonBlock(; k_prod = 0.3, k_cons = 0.2, gamma = 0.01,
                     contribs = [:M_atp_c, :M_pi_c],
                     edges = CouplingEdge[CurrencyEdge(species = :M_atp_c, direction = :out),
                                          MassEdge(species = :M_pi_c, direction = :in)],
                     ins = [:M_pi_c])
    params = vcat(InferParameter[InferParameter(k_prod, LogNormal(log(k_prod), 0.5), false, :k_prod, :CarbonBlock, :rate),
              InferParameter(k_cons, LogNormal(log(k_cons), 0.5), false, :k_cons, :CarbonBlock, :rate),
              InferParameter(gamma, LogNormal(log(gamma), 0.5), false, :gamma_c, :CarbonBlock, :rate)],
              _ic_params(CARBON_SPECIES, :CarbonBlock))
    return CarbonBlock(params, collect(Symbol, contribs), collect(CouplingEdge, edges),
                       collect(Symbol, ins))
end
states(::CarbonBlock) = CARBON_SPECIES
parameters(m::CarbonBlock) = m.params
coupling(m::CarbonBlock) = m.edges
inputs(m::CarbonBlock) = m.ins
contributed_states(m::CarbonBlock) = m.contribs
dynamics(u, p, t, ::CarbonBlock, u_inputs) = -p[3] .* u
# Entry 1: ATP produced from g6p (u[1]). Entry 2: phosphate drawn, read via inputs.
contributions(u, p, t, ::CarbonBlock, u_inputs) = SA[p[1] * u[1], -p[2] * u_inputs[1]]

"""
    InertPool()

Owns GTP with zero dynamics and no free parameters. Whatever it holds at the
end of a run is exactly what was contributed to it, which is the done-when
check of phase 1 in its simplest form. Its derivative is `Float64` whatever the
parameter eltype, so composing it with a differentiated `Source` exercises the
mixed-eltype accumulation path.
"""
struct InertPool <: AbstractSubModel
    params::Vector{InferParameter}
end
InertPool() = InertPool(_ic_params([:M_gtp_c], :InertPool))
states(::InertPool) = [:M_gtp_c]
parameters(m::InertPool) = m.params
dynamics(u, p, t, ::InertPool) = SA[zero(u[1])]

"""
    Source(; k = 0.7)

Owns GDP with zero dynamics and contributes a constant `k` into GTP through an
outbound mass edge. The owner's state must then rise linearly at exactly `k`.
"""
struct Source <: AbstractSubModel
    params::Vector{InferParameter}
end
Source(; k = 0.7) = Source(vcat(InferParameter[InferParameter(k, LogNormal(log(k), 0.5), false, :k_src, :Source, :rate)],
                            _ic_params([:M_gdp_c], :Source)))
states(::Source) = [:M_gdp_c]
parameters(m::Source) = m.params
coupling(::Source) = CouplingEdge[MassEdge(species = :M_gtp_c, direction = :out)]
contributed_states(::Source) = [:M_gtp_c]
dynamics(u, p, t, ::Source) = SA[zero(u[1])]
contributions(u, p, t, ::Source, u_inputs) = SA[p[1]]

"""
    MiscountedSource()

Declares two contributed species and returns one term. Exists only so the
length-mismatch error has a witness.
"""
struct MiscountedSource <: AbstractSubModel
    params::Vector{InferParameter}
end
MiscountedSource() = MiscountedSource(vcat(InferParameter[InferParameter(0.5, LogNormal(log(0.5), 0.5), false, :k_mis, :MiscountedSource, :rate)],
                                       _ic_params([:M_g6p_c], :MiscountedSource)))
states(::MiscountedSource) = [:M_g6p_c]
parameters(m::MiscountedSource) = m.params
coupling(::MiscountedSource) = CouplingEdge[MassEdge(species = :M_gtp_c, direction = :out),
                                            MassEdge(species = :M_gmp_c, direction = :out)]
contributed_states(::MiscountedSource) = [:M_gtp_c, :M_gmp_c]
dynamics(u, p, t, ::MiscountedSource) = SA[zero(u[1])]
contributions(u, p, t, ::MiscountedSource, u_inputs) = SA[p[1]]
