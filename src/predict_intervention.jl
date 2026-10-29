# Interventions for `predict` (KIM-32): do(a), applied after the parameter draw.
#
# p(Q | D, do(a)) = ∫ p(Q | g_a(θ)) p(θ | D) dθ, and the reported effect is
# paired: Δ = Q(g_a(θ), ξ) − Q(θ, ξ), both arms on the same cell seeds, so most
# of the intrinsic noise cancels (dev/notes/2026-10-29-cut-predictive-distribution.md §4).
# The assumptions the paper states: the intervention acts only through the
# named parameter; the posterior fitted on unperturbed data carries over; and
# the model stays valid at the dose, which the dose screen (KIM-37) checks.

"An intervention `predict` applies after the parameter draw. See [`Knockdown`](@ref)."
abstract type Intervention end

"""
    Knockdown(; param = :S_JCVISYN3A_0779, dose)

A promoter knockdown: `param` is scaled by `1 − dose` after the parameter draw,
and nothing else changes ([`intervene`](@ref)). The default is the ptsG promoter,
the first case (KIM-32). It is an inferred parameter, and its effect crosses the
gene-expression–metabolism boundary. `dose` is in [0, 1], and the dose screen
(KIM-37) picks it. A dose of 0 changes nothing, so its effect is exactly zero.
"""
struct Knockdown <: Intervention
    param::Symbol
    dose::Float64
    function Knockdown(param::Symbol, dose::Real)
        0 <= dose <= 1 || throw(ArgumentError("dose must be in [0, 1], got $dose"))
        return new(param, Float64(dose))
    end
end
Knockdown(; param::Symbol = promoter_param(:JCVISYN3A_0779), dose::Real) = Knockdown(param, dose)

"""
    intervene(a, models, values) -> Vector{Pair{Symbol, Float64}}

The parameter writes `values` with intervention `a` applied. For a
[`Knockdown`](@ref), the parameter's drawn value is scaled by `1 − dose`. If
`values` does not set it, the scaling starts from its nominal value. Reverse
constants are not re-derived here. `predict` intervenes before deriving them.
"""
function intervene(a::Knockdown, models::AbstractVector{<:AbstractSubModel}, values)
    vals = Pair{Symbol, Float64}[n => Float64(v) for (n, v) in values]
    i = findfirst(p -> first(p) === a.param, vals)
    if i === nothing
        nominal = nominal_parameter_values(models)
        haskey(nominal, a.param) || throw(ArgumentError(
            "No module in the composition has the parameter :$(a.param) to knock down"))
        push!(vals, a.param => nominal[a.param] * (1 - a.dose))
    else
        vals[i] = a.param => last(vals[i]) * (1 - a.dose)
    end
    return vals
end

export Intervention, Knockdown, intervene
