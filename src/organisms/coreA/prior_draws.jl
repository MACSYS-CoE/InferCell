"""
Prior draws that keep each reaction's equilibrium constant (spec §11 task 14c.1).

A reversible modular rate law carries its equilibrium constant implicitly, in the
Haldane relation

    Keq = kcatF / kcatR · Π KmPⱼ^nⱼ / Π KmSᵢ^nᵢ

so drawing the forward constant, the reverse constant and the Michaelis
constants independently moves the reaction's thermodynamics as well as its
kinetics. Phase 14b's broad census did exactly that, moving the equilibrium
constants by ×5 to 40 at 1σ (§12, 2026-09-28). Under the rule here the reverse
constant is derived from the others through the equilibrium constant, so a draw
changes rates and never where a reaction settles.

**The equilibrium constant kept is the one the nominal Mode values imply**, not
upstream's `equilibrium constant` rows, which disagree with the Modes on 7 of the
15 reactions, by up to 1.6e4 on TPI (§12, 2026-09-28, phase 14c planning A).
Keeping the Mode-implied constant is what keeps the nominal draw the published
model.
"""

"""
    HaldaneRelation

One reversible modular reaction's equilibrium constant as a function of its
parameters, by name: `forward` and `reverse` catalytic constants, and each side's
Michaelis constants with the exponent the rate law raises them to (one per unit
of stoichiometry, so a coefficient of two is an exponent of two).
"""
struct HaldaneRelation
    reaction::Symbol
    forward::Symbol
    reverse::Symbol
    substrates::Vector{Pair{Symbol, Int}}
    products::Vector{Pair{Symbol, Int}}
end

"""
    haldane_relations(m) -> Vector{HaldaneRelation}
    haldane_relations(models) -> Vector{HaldaneRelation}

The Haldane relation of every reversible modular reaction a module runs. Empty
for a module with none, which is every module but glycolysis and recycling: the
PTS steps are mass-action pairs with asserted constants, export is diffusion,
and charging is irreversible (spec §3).

Glycolysis's are read off the compiled rate laws — the same parameter indices
`_reaction_rate` evaluates — so they cannot describe a different reaction from
the one that runs. Recycling's rate laws index their constants by position,
and [`RECYCLING_HALDANE`](@ref) restates those positions. The suite ties each
relation to its rate law by setting a reaction's mass-action ratio to the
relation's equilibrium constant and asserting the flux vanishes.
"""
haldane_relations(::AbstractSubModel) = HaldaneRelation[]
haldane_relations(models::AbstractVector{<:AbstractSubModel}) =
    reduce(vcat, (haldane_relations(m) for m in models); init = HaldaneRelation[])

function haldane_relations(m::CentralGlycolysis)
    name(k) = m.params[k].name
    return [HaldaneRelation(Symbol(chopprefix(String(name(r.kf)), "kcatF_")),
                            name(r.kf), name(r.kr),
                            [name(k) => n for (k, n) in zip(r.skm, r.scoef)],
                            [name(k) => n for (k, n) in zip(r.pkm, r.pcoef)])
            for r in m.rates]
end

"""
    RECYCLING_HALDANE

Positions in [`RECYCLING_KINETIC_IDS`](@ref) of each recycling reaction's
forward and reverse constants and of its Michaelis constants with their
exponents, in [`RECYCLING_REACTIONS`](@ref) order. Restates the positions
`recycling_fluxes` reads: ADK1 and PPA each square one product's saturation
ratio, so their products carry exponent 2.
"""
const RECYCLING_HALDANE = (
    (1, 2, (11 => 1, 12 => 1), (13 => 1, 14 => 1)),    # PGK3: 13dpg + gdp ⇌ 3pg + gtp
    (3, 4, (15 => 1, 16 => 1), (17 => 1, 18 => 1)),    # PYK3: gdp + pep ⇌ gtp + pyr
    (5, 6, (19 => 1, 20 => 1), (21 => 2,)),            # ADK1: amp + atp ⇌ 2 adp
    (7, 8, (22 => 1, 23 => 1), (24 => 1, 25 => 1)),    # GK1:  atp + gmp ⇌ adp + gdp
    (9, 10, (26 => 1,), (27 => 2,)),                   # PPA:  ppi ⇌ 2 pi
)

function haldane_relations(m::NucleotideRecycling)
    id(i) = RECYCLING_KINETIC_IDS[i]
    return [HaldaneRelation(r, id(kf), id(kr),
                            [id(i) => n for (i, n) in subs],
                            [id(i) => n for (i, n) in prods])
            for (r, active, (kf, kr, subs, prods)) in
                zip(RECYCLING_REACTIONS, m.active, RECYCLING_HALDANE) if active]
end

"""
    equilibrium_constant(rel, values) -> Float64

`kcatF/kcatR · Π KmP^n / Π KmS^n`, with `values` mapping each parameter name
to its value.
"""
function equilibrium_constant(rel::HaldaneRelation, values::AbstractDict{Symbol})
    k = values[rel.forward] / values[rel.reverse]
    for (km, n) in rel.products
        k *= values[km]^n
    end
    for (km, n) in rel.substrates
        k /= values[km]^n
    end
    return k
end

# The reverse constant that puts `rel` at equilibrium constant `keq`, given the
# forward and Michaelis constants in `values`.
function _derived_reverse(rel::HaldaneRelation, values::AbstractDict{Symbol}, keq)
    k = values[rel.forward] / keq
    for (km, n) in rel.products
        k *= values[km]^n
    end
    for (km, n) in rel.substrates
        k /= values[km]^n
    end
    return k
end

"""
    nominal_parameter_values(models) -> Dict{Symbol, Float64}

Every parameter's stored value, by name, over a composition. A name two modules
share must carry one value, or the composition would run one parameter at two.
"""
function nominal_parameter_values(models::AbstractVector{<:AbstractSubModel})
    v = Dict{Symbol, Float64}()
    for m in models, p in parameters(m)
        if haskey(v, p.name)
            v[p.name] == p.value || throw(ArgumentError(
                "Parameter :$(p.name) is $(v[p.name]) in one module and $(p.value) " *
                "in another; a composition cannot run it at both"))
        else
            v[p.name] = p.value
        end
    end
    return v
end

function _parameter_by_name(models, name::Symbol)
    for m in models, p in parameters(m)
        p.name === name && return p
    end
    throw(ArgumentError("No module in the composition has a parameter :$name"))
end

"""
    draw_parameters(rng, models, names; reverse = :derived) -> Vector{Pair{Symbol, Float64}}

One prior draw of the parameters `names`, as `name => value` pairs ready for
[`set_parameters!`](@ref).

Under `reverse = :derived`, the default and the rule of spec §11 task 14c.1,
every other name is drawn from its prior, and then each reversible reaction a
drawn parameter touches has its reverse constant derived from its forward and
Michaelis constants through its nominal equilibrium constant. A reverse constant
named in `names` is derived rather than drawn. So drawing only a forward
constant scales that reaction's reverse constant by the same factor, and every
reaction keeps its equilibrium constant, which [`assert_haldane`](@ref) checks.

`reverse = :independent` draws every name from its own prior, reverse constants
included. That is phase 14b's rule, which breaks the Haldane relations. It is
kept to reproduce 14b's census and as the mutation 14c.1's check must fail on.
"""
function draw_parameters(rng::AbstractRNG, models::AbstractVector{<:AbstractSubModel},
                         names::AbstractVector{Symbol}; reverse::Symbol = :derived)
    reverse in (:derived, :independent) || throw(ArgumentError(
        "reverse must be :derived or :independent, not :$reverse"))
    allunique(names) || throw(ArgumentError("names must not repeat"))
    if reverse === :independent
        return [n => rand(rng, _parameter_by_name(models, n).prior) for n in names]
    end

    rels = haldane_relations(models)
    reverses = Set(r.reverse for r in rels)
    drawn = [n => rand(rng, _parameter_by_name(models, n).prior)
             for n in names if !(n in reverses)]

    nominal = nominal_parameter_values(models)
    values = merge(nominal, Dict(drawn))
    moved = Set(first.(drawn))
    derived = Pair{Symbol, Float64}[]
    for r in rels
        touched = r.reverse in names || r.forward in moved ||
                  any(km in moved for (km, _) in r.substrates) ||
                  any(km in moved for (km, _) in r.products)
        touched || continue
        push!(derived, r.reverse => _derived_reverse(r, values,
                                                     equilibrium_constant(r, nominal)))
    end
    return vcat(drawn, derived)
end

"""
    assert_haldane(models, draw) -> Vector{NamedTuple}

Check that `draw` keeps every reaction's nominal equilibrium constant, to
roundoff. Returns one row per reaction, `(reaction, keq, relative, bound)`, and
throws naming the first reaction whose relative change exceeds its bound.

The bound is `4n·eps`, where `n` counts the factors in the Haldane quotient
(two catalytic constants plus one per Michaelis-constant power): deriving the
reverse constant and recomputing the quotient each round at most once per
factor. A reverse constant drawn independently moves the constant by the ratio
of two log-normal draws, so a real break misses the bound by many orders.
"""
function assert_haldane(models::AbstractVector{<:AbstractSubModel},
                        draw::AbstractVector{<:Pair{Symbol}})
    nominal = nominal_parameter_values(models)
    values = merge(nominal, Dict(draw))
    rows = NamedTuple[]
    for r in haldane_relations(models)
        k0 = equilibrium_constant(r, nominal)
        k = equilibrium_constant(r, values)
        rel = abs(k / k0 - 1)
        nfac = 2 + sum(last, r.substrates) + sum(last, r.products)
        bound = 4 * nfac * eps()
        push!(rows, (reaction = r.reaction, keq = k0, relative = rel, bound = bound))
        rel <= bound || throw(ArgumentError(
            "Draw moves $(r.reaction)'s equilibrium constant from $k0 to $k, a " *
            "relative change of $rel against a roundoff bound of $bound. Its " *
            "reverse constant :$(r.reverse) was not derived through the Haldane " *
            "relation (spec §11 task 14c.1)"))
    end
    return rows
end

"""
    D11_TARGETS

Spec §4 D11's six inference targets, by parameter name, in the table's order:
the GAPD, PGI and ptsG promoter strengths, the decay constant, and the enolase
and aldolase forward constants.
"""
const D11_TARGETS = [promoter_param(:JCVISYN3A_0607), promoter_param(:JCVISYN3A_0445),
                     promoter_param(:JCVISYN3A_0779), :krnadeg,
                     :kcatF_R_ENO, :kcatF_R_FBA]

"""
    d11_models() -> Vector{AbstractSubModel}

The assembled composition with D11's two ODE targets freed, and the reverse
constants their draws derive. The four stochastic targets are free already.
"""
function d11_models()
    ms = corea_models()
    ms[1] = CentralGlycolysis(enzymes = :translated,
                              free = [:kcatF_R_ENO, :kcatR_R_ENO, :kcatF_R_FBA, :kcatR_R_FBA])
    return ms
end

export HaldaneRelation, RECYCLING_HALDANE, haldane_relations, equilibrium_constant,
       nominal_parameter_values, draw_parameters, assert_haldane, D11_TARGETS,
       d11_models
