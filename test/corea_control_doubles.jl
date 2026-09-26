using InferCell
using StaticArrays: SA
using Distributions: LogNormal, Normal

import InferCell: states, parameters, dynamics, contributions, contributed_states,
                  coupling, inputs, module_id

# The frozen-expression variant of spec §3 check 8 (§11 task 14.8).
#
# Check 8's summation theorems hold at a steady state, and the assembled model
# has none: expression grows the cell, and translation drains charged tRNA and
# GTP through handshake counters that are not in the ODE right-hand side. With
# expression frozen, nothing consumes either, so the tRNA pool would saturate
# and charging's flux would go to zero, which §3 names as the trap. So the
# variant is the four ODE modules at nominal enzyme counts, plus a
# `FrozenDemand` that stands in for translation as two first-order consumers:
#
#     M_trna_chg_c -> M_trna_c              at k_tl  · [M_trna_chg_c]
#     M_gtp_c      -> M_gdp_c + M_pi_c      at k_gtp · [M_gtp_c]
#
# First order rather than constant, so each is a reaction with a rate the
# steady state sets, and its rate constant is the multiplier the theorem needs.
# At the nominal pools they consume the published demand, 553.10 residues/s,
# and two GTP per residue, as translation's counters charge.
#
# Shared by `test/test_corea_validation_14b.jl` and
# `dev/scripts/corea_validation_14b.jl`, so the tested variant and the recorded
# one cannot differ.

"""
    FrozenDemand(; tl_scale = 1.0, gtp_scale = 1.0)

Translation's two ODE-side demands, frozen: charged tRNA returned uncharged,
and GTP hydrolysed to GDP and phosphate. The two rate constants are free
parameters so a control analysis can scale them. It owns one accumulator,
`frozen_residues_mM`, the residues consumed, which has no steady state and is
held by the analysis.
"""
struct FrozenDemand <: AbstractSubModel
    params::Vector{InferParameter}
end
function FrozenDemand(; tl_scale::Real = 1.0, gtp_scale::Real = 1.0)
    (; pool_mM, charged_fraction) = CHARGING_POOL_DEFAULTS
    demand = charging_demand_mM_per_s()
    k_tl = tl_scale * demand / (charged_fraction * pool_mM)
    k_gtp = gtp_scale * 2 * demand / species_entry(:M_gtp_c).initial_value
    return FrozenDemand([
        InferParameter(0.0, Normal(0.0, 0.1), true, :frozen_residues_mM0,
                       :FrozenDemand, :initial_condition),
        InferParameter(k_tl, LogNormal(log(k_tl), log(2.0)), false, :k_frozen_tl,
                       :FrozenDemand, :rate),
        InferParameter(k_gtp, LogNormal(log(k_gtp), log(2.0)), false, :k_frozen_gtp,
                       :FrozenDemand, :rate),
    ])
end
states(::FrozenDemand) = [:frozen_residues_mM]
parameters(m::FrozenDemand) = m.params
module_id(::FrozenDemand) = :FrozenDemand
inputs(::FrozenDemand) = [:M_trna_chg_c, :M_gtp_c]
contributed_states(::FrozenDemand) = [:M_trna_c, :M_trna_chg_c, :M_gtp_c, :M_gdp_c, :M_pi_c]
coupling(::FrozenDemand) = CouplingEdge[
    CurrencyEdge(species = :M_trna_chg_c, direction = :in),
    CurrencyEdge(species = :M_trna_c, direction = :out),
    CurrencyEdge(species = :M_gtp_c, direction = :in),
    CurrencyEdge(species = :M_gdp_c, direction = :out),
    CurrencyEdge(species = :M_pi_c, direction = :out),
]
@inline _frozen_rates(p, ui) = (p[1] * ui[1], p[2] * ui[2])
function dynamics(u, p, t, m::FrozenDemand, ui)
    vt, _ = _frozen_rates(p, ui)
    return SA[vt]
end
function contributions(u, p, t, m::FrozenDemand, ui)
    vt, vg = _frozen_rates(p, ui)
    return SA[vt, -vt, -vg, vg, vg]
end

"""
    HeldProteins()

Owns the ten glycolytic protein counts `CentralGlycolysis` declares as inputs
at nominal enzymes, each held at its published copy number with zero
derivative. It is phase 6's `HeldEnergyPool` without the currencies, which the
frozen variant integrates. Expression is frozen, so nothing moves them.
"""
struct HeldProteins <: AbstractSubModel
    params::Vector{InferParameter}
end
HeldProteins() = HeldProteins([InferParameter(Float64(r.copies), LogNormal(log(r.copies), 1.0),
                                              true, Symbol(s, "0"), :HeldProteins,
                                              :initial_condition)
                               for (s, r) in zip(default_protein_sources(), GLYCOLYTIC_REACTIONS)])
states(::HeldProteins) = default_protein_sources()
parameters(m::HeldProteins) = m.params
module_id(::HeldProteins) = :HeldProteins
dynamics(u, p, t, ::HeldProteins) = zero(u)

# Every rate-linear parameter of the four ODE modules, freed so the control
# analysis can scale it. Freeing changes no value.
const _GLYC_IDS = Tuple(r.id for r in GLYCOLYTIC_REACTIONS)
const _PTS_STEPS = (:glcpts0, :glcpts1, :glcpts2, :glcpts3, :glcpts4)
_kcats(ids) = Symbol[s for r in ids for s in (Symbol("kcatF_", r), Symbol("kcatR_", r))]

"""
    frozen_models() -> Vector{AbstractSubModel}

The frozen-expression variant: glycolysis and recycling at their nominal
enzyme concentrations, which are the published counts at the registry's cell;
PTS transport; charging; and [`FrozenDemand`](@ref). Every parameter a
multiplier scales is free.
"""
frozen_models() = AbstractSubModel[
    CentralGlycolysis(free = _kcats(_GLYC_IDS)),
    PtsTransport(free = vcat([Symbol(k, "_", s) for s in _PTS_STEPS for k in (:kf, :kr)],
                             [:p_lact2r])),
    NucleotideRecycling(free = _kcats(RECYCLING_REACTIONS)),
    TrnaCharging(free_k_chg = true),
    FrozenDemand(),
    HeldProteins(),
]

"""
    frozen_multipliers(pnames) -> Vector{Pair{Symbol, Vector{Int}}}

One multiplier per reaction of the frozen variant, as slots of the parameter
vector whose names are `pnames`: the 22 ODE reactions and the two demands
(spec §3 check 8, amended 2026-09-26).
"""
function frozen_multipliers(pnames::AbstractVector{Symbol})
    function slot(s)
        i = findfirst(==(s), pnames)
        i === nothing && throw(ArgumentError("no free parameter :$s"))
        return i
    end
    out = Pair{Symbol, Vector{Int}}[]
    for r in (_GLYC_IDS..., RECYCLING_REACTIONS...)
        push!(out, r => [slot(Symbol("kcatF_", r)), slot(Symbol("kcatR_", r))])
    end
    for s in _PTS_STEPS
        push!(out, Symbol("R_GLCpts", last(string(s))) => [slot(Symbol(:kf_, s)), slot(Symbol(:kr_, s))])
    end
    push!(out, :R_LACt => [slot(:p_lact2r)])
    push!(out, :R_charging => [slot(:k_chg)])
    push!(out, :frozen_translation => [slot(:k_frozen_tl)])
    push!(out, :frozen_gtp => [slot(:k_frozen_gtp)])
    return out
end

"""
    frozen_gene_groups() -> Vector{Pair{Symbol, Vector{Symbol}}}

The thirteen genes whose protein multiplies a rate, each with its reactions.
PGK and PYK each catalyse the glycolytic step and its GTP branch. The four PTS
genes are carrier totals and are reported through the conserved totals
instead.
"""
frozen_gene_groups() = vcat(
    [Symbol(string(r)[3:end]) => (r === :R_PGK ? [:R_PGK, :R_PGK3] :
                                  r === :R_PYK ? [:R_PYK, :R_PYK3] : [r])
     for r in _GLYC_IDS],
    [Symbol(string(r)[3:end]) => [r] for r in (:R_ADK1, :R_GK1, :R_PPA)])

"""
The states check 8 holds: external lactate, a boundary species; the demand's
accumulator; and the held protein counts.
"""
const FROZEN_HELD = vcat([:M_lac__L_e, :frozen_residues_mM], default_protein_sources())

"""
    frozen_control_problem(; models = frozen_models()) -> ControlProblem

The frozen variant as a [`ControlProblem`](@ref), built through `build_problem`
so the right-hand side is the composed one the model integrates.
"""
function frozen_control_problem(; models = frozen_models(),
                                multipliers = frozen_multipliers)
    prob = build_problem(models; tspan = (0.0, 1.0))
    pnames = unique(Symbol[q.name for m in models for q in model_free_params(parameters(m))])
    length(pnames) == length(prob.p) || error("parameter layout mismatch")
    snames = Symbol[s for m in models for s in states(m)]
    f(u, p) = prob.f(u, p, 0.0)
    return ControlProblem(f, collect(prob.u0), prob.p, snames;
                          held = FROZEN_HELD, multipliers = multipliers(pnames))
end

