# The variance split on the crossed matrix (KIM-66): the variance of the
# population mean μ under one of a prediction's weight sets, and its nested
# splits between the inferred θ_I (rows) and the held θ_H (columns), in both
# orders (dev/notes/2026-10-29-cut-predictive-distribution.md §5.2).
#
# Every moment is a weighted U-statistic: a product of two μ̂ enters only when
# the two come from different draws of what is being averaged over, so the cell
# noise ε, which is independent between pairs, cancels from every product. The
# one diagonal term, E[μ²], has its noise τ = σ̂²/n removed. Without a mask the
# estimates are unbiased. Under V the weights are ratio estimates of the
# conditioned law, biased at O(1/K); dev/scripts/variance_split_check_kim66.jl
# measures that bias on a known split.

# Σ over pairs of entries in distinct rows and distinct columns of z_i z_j, from
# the row and column sums: every pair, less those sharing a row or a column,
# plus those sharing both (each entry with itself), which were taken twice.
_offpairs(z) = sum(z)^2 - sum(abs2, sum(z; dims = 2)) - sum(abs2, sum(z; dims = 1)) +
               sum(abs2, z)
_offrows(z) = sum(z)^2 - sum(abs2, sum(z; dims = 2))       # distinct rows only
_offcols(z) = sum(z)^2 - sum(abs2, sum(z; dims = 1))       # distinct columns only

# E[c²] for the conditional mean c along each row (or column): within each line,
# a product of two distinct entries, weighted by the line's mass. A line needs
# two entries with weight. Returns the estimate and the lines it used.
function _line_second_moment(z, w, dims)
    num = vec(sum(z; dims) .^ 2 .- sum(abs2, z; dims))
    den = vec(sum(w; dims) .^ 2 .- sum(abs2, w; dims))
    mass = vec(sum(w; dims))
    use = den .> 0
    any(use) || return NaN, 0
    return sum(mass[use] .* num[use] ./ den[use]) / sum(mass[use]), count(use)
end

"""
    variance_split(p::Prediction; weights = first(keys(p.weights))) -> NamedTuple

The variance of the population mean μ under the weight set `weights` of the
prediction `p` (`:uniform`, or `:per_theta_I` or `:joint` under V), and its nested
splits between the inferred parameters θ_I (the rows) and the held θ_H (the
columns), in both orders (dev/notes/2026-10-29-cut-predictive-distribution.md
§5.2):

```text
I first:  Var[μ] = Var_I[ E(μ | θ_I) ] + E_I[ Var(μ | θ_I) ]
H first:  Var[μ] = Var_H[ E(μ | θ_H) ] + E_H[ Var(μ | θ_H) ]
```

Each moment is a weighted U-statistic on the matrix of μ̂, centred first on its
weighted mean, so the split does not depend on where μ sits. A product of two μ̂
enters only when they come from distinct draws of what is averaged over, so the
cells' Monte Carlo noise cancels from it. The one diagonal term, E[μ²], has the
noise τ = σ̂²/n removed. The crossed design's shared columns then bias nothing.
Without a mask (`:uniform`) every estimate is unbiased. Under V the weights are
ratio estimates of the conditioned law, biased at O(1/K).
`dev/scripts/variance_split_check_kim66.jl` measures that bias on known splits,
and M4's freeze entry records the result.

Returns:
- `total`: Var[μ];
- `i_first = (; between, within)`: Var_I[E(μ | θ_I)] and E_I[Var(μ | θ_I)];
- `h_first = (; between, within)`: Var_H[E(μ | θ_H)] and E_H[Var(μ | θ_H)];
- `held_share = (; i_first, h_first)`: the held parameters' share of `total` in
  each order: `i_first.within / total`, the held parameters with their
  interaction, and `h_first.between / total`, their first-order effect. The gap
  between the orders, `gap = h_first.within − i_first.between`, is the
  interaction plus any dependence V induced. It is taken from the components
  before they are clipped at 0, so it matches the reported ones only when
  nothing is `flagged`;
- `mc`: the noise removed, Σ w τ;
- `flagged`: which of the five components came out negative. Each is reported
  as 0, as the note asks;
- `mc_undefined`: true when no weighted hybrid pair has two finished cells, as
  at `n_cells = 1`. Then σ̂² does not exist, so `total` and both `within` terms
  are NaN. The two `between` terms do not need it;
- `dropped`: the weighted pairs left out because they had under two finished
  cells while others had two, so that one such pair does not blank the split.
  Leaving them out changes the weights as a missing μ̂ does, and they are
  counted here, not hidden;
- `rows` and `columns`: with both groups drawn, the lines with two weighted
  entries that the conditional second moments used; with one, the lines of the
  drawn group that carry weight. `pairs` is the number of weighted pairs;
- `inferred_varies` and `held_varies`: whether each group is drawn. θ_H is
  drawn only with completion and at least two held draws, and θ_I only at two
  or more rows. A group that is not drawn takes no share: its `between` term is 0.

Pairs with no finished run carry no μ̂, so they leave the weights. For the mean
field, τ = 0.
"""
function variance_split(p::Prediction; weights::Symbol = first(keys(p.weights)))
    haskey(p.weights, weights) || throw(ArgumentError(
        "this prediction has weight sets $(keys(p.weights)), not :$weights"))
    τ = p.forward === :hybrid ? p.sigma2 ./ p.n_ok : zeros(size(p.mu))
    return (; weights, _variance_split(p.mu, τ, p.weights[weights]; held = !isempty(p.held))...)
end

# The split on a matrix of μ̂ with noise τ and weights w0, θ_I down the rows and
# θ_H across the columns; `held` says whether the columns are θ_H draws at all.
function _variance_split(mu::AbstractMatrix, τ::AbstractMatrix, w0::AbstractMatrix;
                         held::Bool)
    M, K = size(mu)
    ok = isfinite.(mu)
    # A pair with under two finished cells has no τ, so its noise cannot be taken
    # out of E[μ²]. If other pairs have one, it leaves the weights and is counted.
    weighted = ok .& (w0 .> 0)
    lacks = weighted .& .!isfinite.(τ)
    dropped = any(lacks) && any(weighted .& .!lacks) ? count(lacks) : 0
    dropped > 0 && (ok = ok .& .!lacks)
    w = w0 .* ok
    W = sum(w)
    W > 0 || throw(ArgumentError("the weight set weighs no pair with a finite μ̂"))
    w = w ./ W
    # Centred on the weighted mean first. Under a mask the U-statistics' pair
    # weights are not `w`'s, so they are not invariant to a shift of μ, and a
    # mean far from zero (Q in particles) would swamp its variance. Without a
    # mask they are invariant, and centring changes nothing.
    y = ifelse.(ok, mu, 0.0)
    y = ifelse.(ok, y .- sum(w .* y), 0.0)
    z = w .* y
    mc_undefined = any(!isfinite(τ[i]) for i in eachindex(w) if w[i] > 0)
    mc = mc_undefined ? NaN : sum(w[i] * τ[i] for i in eachindex(w) if w[i] > 0)
    m2 = sum(w[i] * y[i]^2 for i in eachindex(w) if w[i] > 0) - mc

    rows_vary = M >= 2
    cols_vary = held && K >= 2
    rows_vary || cols_vary || throw(ArgumentError(
        "one row and no held draws: nothing varies, so there is no variance to split"))
    # (E μ)², from products of entries independent in every group that varies.
    sq = rows_vary && cols_vary ? _offpairs(z) / _offpairs(w) :
         rows_vary ? _offrows(z) / _offrows(w) : _offcols(z) / _offcols(w)
    total = m2 - sq
    if rows_vary && cols_vary
        c2, nr = _line_second_moment(z, w, 2)          # along rows: θ_I fixed
        d2, nc = _line_second_moment(z, w, 1)          # along columns: θ_H fixed
        a, c = c2 - sq, d2 - sq
    elseif rows_vary
        a, c, nr, nc = total, 0.0, count(>(0), sum(w; dims = 2)), 0
    else
        a, c, nr, nc = 0.0, total, 0, count(>(0), sum(w; dims = 1))
    end
    b, d = total - a, total - c
    clip(x) = isnan(x) ? x : max(x, 0.0)
    t = clip(total)
    return (; total = t, i_first = (; between = clip(a), within = clip(b)),
            h_first = (; between = clip(c), within = clip(d)),
            held_share = (; i_first = clip(b) / t, h_first = clip(c) / t),
            gap = d - a, mc,
            flagged = (; total = total < 0, i_between = a < 0, i_within = b < 0,
                       h_between = c < 0, h_within = d < 0),
            mc_undefined, dropped, rows = nr, columns = nc, pairs = count(>(0), w),
            inferred_varies = rows_vary, held_varies = cols_vary)
end

export variance_split
