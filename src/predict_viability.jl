# Viability-conditioned completion (KIM-30). The viability rule V is data: its
# thresholds are passed in, and none is fixed here. V is frozen in M4 (KIM-36,
# KIM-39), before any prediction is read
# (dev/notes/2026-10-29-cut-predictive-distribution.md §3.3).
#
# A cell is one of five things, and they are not interchangeable: it failed, it
# starved, it fell below a floor, it clipped, or it ran clean. Clipping is
# routine: 81.5% of 14c.2's D11 prior draws clip at least once, while 4% starve.

"""
    ViabilityRule(; label, allow_clipped, min_viable_fraction, starved_fraction = 0.5,
                  floors = Pair{Symbol, Float64}[])

A declared viability rule V, passed to `predict` as data.

Each cell is put in one category, in this order of precedence ([`cell_status`](@ref)):
1. `:failed`: the run did not finish;
2. `:starved`: some counter clipped at more than `starved_fraction` of the drains.
   The default 0.5 is 14c's census rule. For the mean field, an overdrawn pool
   is its form of starving;
3. `:below_floor`: some state in `floors` (`species => particles`) went below its
   floor at some point in the run;
4. `:clipped`: some counter clipped at least once;
5. `:clean`.

A cell is viable when it is `:clean`, or `:clipped` with `allow_clipped`.

V is a property of θ, not of a cell. θ is viable when at least
`min_viable_fraction` of its cells are. A per-cell V would condition the cells
rather than the parameters, and would change what the population mean means
(note §3.3). The mean field has one run per θ, so the fraction there is 0 or 1.
`label` names the rule in a prediction's provenance.
"""
struct ViabilityRule
    label::String
    allow_clipped::Bool
    min_viable_fraction::Float64
    starved_fraction::Float64
    floors::Vector{Pair{Symbol, Float64}}
    function ViabilityRule(label, allow_clipped, min_viable_fraction, starved_fraction, floors)
        0 < min_viable_fraction <= 1 || throw(ArgumentError(
            "min_viable_fraction must be in (0, 1], got $min_viable_fraction"))
        0 < starved_fraction <= 1 || throw(ArgumentError(
            "starved_fraction must be in (0, 1], got $starved_fraction"))
        fl = Pair{Symbol, Float64}[Symbol(s) => Float64(x) for (s, x) in floors]
        allunique(first.(fl)) || throw(ArgumentError("a floor names a state twice"))
        return new(String(label), allow_clipped, Float64(min_viable_fraction),
                   Float64(starved_fraction), fl)
    end
end
ViabilityRule(; label::AbstractString, allow_clipped::Bool, min_viable_fraction::Real,
              starved_fraction::Real = 0.5, floors = Pair{Symbol, Float64}[]) =
    ViabilityRule(label, allow_clipped, min_viable_fraction, starved_fraction, floors)

"The five cell categories, from worst to best."
const CELL_STATUSES = (:failed, :starved, :below_floor, :clipped, :clean)

"""
    cell_status(cell, V) -> Symbol
    cell_status(cell; starved_fraction = 0.5) -> Symbol

The category of one [`CellRun`](@ref) under the rule `V` ([`ViabilityRule`](@ref)),
one of `CELL_STATUSES`. Without a rule, the cell is read with no floors and 14c's
starvation threshold, which is the unfiltered ensemble's diagnostic.
"""
function cell_status(c::CellRun, V::ViabilityRule)
    c.failure === :overdraw && return :starved
    c.failure === :none || return :failed
    any(r -> r.drains > V.starved_fraction * c.drains, c.clips) && return :starved
    for (s, x) in V.floors
        haskey(c.minima, s) || throw(ArgumentError(
            "the rule $(repr(V.label)) sets a floor on :$s, which the run does not have"))
        c.minima[s] < x && return :below_floor
    end
    return isempty(c.clips) ? :clean : :clipped
end
cell_status(c::CellRun; starved_fraction::Real = 0.5) =
    cell_status(c, ViabilityRule("unfiltered", true, 1.0, starved_fraction, Pair{Symbol, Float64}[]))

_viable(status::Symbol, V::ViabilityRule) =
    status === :clean || (status === :clipped && V.allow_clipped)

"""
    assess_viability(cells, V) -> (; viable, fraction, counts)

V's decision on one parameter set from its cells (one [`CellRun`](@ref) each):
whether it is viable, the fraction of its cells that are, and the count of its
cells in each category, which is the reason a rejected set is kept with.
"""
function assess_viability(cells::AbstractVector{CellRun}, V::ViabilityRule)
    isempty(cells) && throw(ArgumentError("no cells to assess"))
    statuses = [cell_status(c, V) for c in cells]
    fraction = count(s -> _viable(s, V), statuses) / length(cells)
    counts = Dict(s => count(==(s), statuses) for s in CELL_STATUSES)
    return (; viable = fraction >= V.min_viable_fraction, fraction, counts)
end

_kish(w) = sum(w) > 0 ? sum(w)^2 / sum(abs2, w) : 0.0

"""
    viability_weights(v::AbstractMatrix{Bool}) -> NamedTuple

The two conditioned variants of note §3.3 as weights on one crossed matrix, where
`v[m, k]` is V's decision for posterior draw m beside held draw k. M3 does not
pick between them; M4's freeze does.
- `per_theta_I`, variant (b): `w = v / (M′ · Σ_k v)`. Each θ_I row that keeps a
  viable θ_H weighs 1/M′, so θ_I's marginal stays the posterior's, conditional
  on at least one viable θ_H per row. M′ counts those rows.
- `joint`, variant (a): `w ∝ v`. A θ_I row weighs its acceptance rate
  P̂(V | θ_I), so the filter acts on θ_I as data.

Also returned: the overall `acceptance` rate; `row_acceptance`, which is
P̂(V | θ_I^m) per row; `failed_rows`, the θ_I rows with no viable θ_H, which (b)
leaves out and reports; and the Kish effective sample size of each weight set,
`ess_per_theta_I` and `ess_joint`. With no viable cell anywhere, both weight sets
are zero and both ESS are 0.
"""
function viability_weights(v::AbstractMatrix{Bool})
    M, K = size(v)
    rows = vec(sum(v; dims = 2))
    kept = count(>(0), rows)
    wb = zeros(M, K)
    for m in 1:M, k in 1:K
        v[m, k] && (wb[m, k] = 1 / (kept * rows[m]))
    end
    total = sum(v)
    wa = total == 0 ? zeros(M, K) : v ./ total
    return (; per_theta_I = wb, joint = wa, acceptance = total / (M * K),
            row_acceptance = rows ./ K, failed_rows = findall(==(0), rows),
            ess_per_theta_I = _kish(wb), ess_joint = _kish(wa))
end

export ViabilityRule, CELL_STATUSES, cell_status, assess_viability, viability_weights
