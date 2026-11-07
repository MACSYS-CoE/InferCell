# Bands and the variance split over time, from the trajectories `predict` recorded
# (KIM-86). Nothing here runs a cell. Three layers of uncertainty are kept apart,
# each a quantile over what it spreads over:
# - `rows`, the inferred parameters: each row's mean, the held parameters averaged out;
# - `pairs`, inferred and held: each pair's mean over its cells;
# - `cells`, a new cell: every cell.

# Type 7 generalised by weight: order the values, and interpolate on
# s_k = (k − 1) w_k + (n − 1) Σ_{i<k} w_i at q s_n. Equal weights give s_k ∝ k − 1,
# which is `Statistics.quantile`. A value with zero weight is not there. It is not the
# quantile of the data repeated by weight, and at very unequal weights it interpolates
# across a heavy value. KIM-86's check-PR simulated that on a 40 × 50 matrix at 11–12%
# acceptance: 1–3% of an SD on average and inward, so bands come out narrower, 5–13% at
# the 95th percentile of a single band, and larger as the kept pairs fall.
function _weighted_quantile(x::AbstractVector, w::AbstractVector, q::Real)
    keep = findall(>(0), w)
    isempty(keep) && throw(ArgumentError("no value carries weight"))
    o = sortperm(view(x, keep))
    xs, ws = x[keep][o], w[keep][o]
    n = length(xs)
    n == 1 && return Float64(xs[1])
    s = Vector{Float64}(undef, n)
    cum = 0.0
    for k in 1:n
        s[k] = (k - 1) * ws[k] + (n - 1) * cum
        cum += ws[k]
    end
    target = q * s[n]
    j = clamp(searchsortedlast(s, target), 1, n - 1)
    f = clamp((target - s[j]) / (s[j + 1] - s[j]), 0.0, 1.0)
    return xs[j] + f * (xs[j + 1] - xs[j])
end

# The recorded species as one cells × times matrix per pair, and the times. A cell
# that holds no trajectory, because it failed, is not in its pair's matrix.
function _trajectory_values(p::Prediction, species::Symbol, baseline::Bool)
    held = (c.trajectory for cs in p.cells for c in cs if c.trajectory !== nothing)
    isempty(held) && throw(ArgumentError(
        "no run of this prediction holds a trajectory: call `predict` with `trajectory`"))
    tr = first(held)
    j = findfirst(==(species), tr.species)
    j === nothing && throw(ArgumentError(
        "this prediction recorded $(tr.species), not :$species"))
    vals = map(p.cells) do cs
        rs = [c.trajectory.values for c in cs if c.trajectory !== nothing]
        X = Matrix{Float64}(undef, length(rs), length(tr.t))
        for (r, v) in enumerate(rs)
            X[r, :] .= baseline ? view(v, :, j) .- v[1, j] : view(v, :, j)
        end
        X
    end
    return tr.t, vals
end

function _weight_matrix(p::Prediction, weights)
    if weights isa Symbol
        haskey(p.weights, weights) || throw(ArgumentError(
            "this prediction has weight sets $(keys(p.weights)), not :$weights"))
        return p.weights[weights]
    end
    size(weights) == size(p.cells) || throw(ArgumentError(
        "weights are $(size(weights)), but the matrix is $(size(p.cells))"))
    return weights
end

function _layer(x, w, levels)
    keep = w .> 0
    v = w[keep] ./ sum(w[keep])
    return (; mean = sum(v .* x[keep]), lo = [_weighted_quantile(x, w, (1 - l) / 2) for l in levels],
            hi = [_weighted_quantile(x, w, (1 + l) / 2) for l in levels])
end

_layer_out(rows, levels) =
    (; mean = [r.mean for r in rows],
     lo = Dict(l => [r.lo[i] for r in rows] for (i, l) in enumerate(levels)),
     hi = Dict(l => [r.hi[i] for r in rows] for (i, l) in enumerate(levels)))

"""
    trajectory_bands(p::Prediction, species; weights = first(keys(p.weights)),
                     levels = (0.5, 0.9), baseline = false)
        -> (; t, n_failed, pairs, rows, cells)

Central bands over time for a `species` recorded by [`predict`](@ref)'s `trajectory`,
as empirical quantiles and not mean ± SD, since a quantity like Q is skewed once
cells starve. Each of the three layers is `(; mean, lo, hi)`, with `lo` and `hi`
`Dict`s from each level to a vector over `t`:

- `rows`: over the inferred parameters. The quantiles of each row's mean, with the
  held parameters averaged out;
- `pairs`: over the inferred and held parameters. The quantiles of each pair's mean
  μ̂(t). This layer carries the Monte Carlo noise of the pair's finite cells, about
  σ̂²/n. [`trajectory_split`](@ref) removes that noise and this does not;
- `cells`: a new cell. The quantiles over every cell.

`weights` is a weight set's name, or an M × K matrix, so `ones(M, K)` gives the
unconditioned band under a viability rule. A pair with zero weight is absent. In
`cells`, each cell carries its pair's weight over the pair's finished cells; in
`rows`, each row carries its pairs' total.

`baseline = true` subtracts each cell's own value at t = 0, so the metered export
gives Q(t). A failed cell holds no trajectory: it is skipped at every time, so the
bands at early times are those of the cells that finished, and the weighted pairs'
failed cells are counted in `n_failed`. A pair with no finished cell leaves the
weights. Weighted quantiles are type 7 generalised by weight, and equal `Statistics.quantile`
at equal weights.

On the mean field each pair has one run, so `cells` repeats `pairs` and is no new-cell
band. The bands read the runs as they are, so δ, a scalar at the horizon, is not in them.
"""
function trajectory_bands(p::Prediction, species::Symbol; weights = first(keys(p.weights)),
                          levels = (0.5, 0.9), baseline::Bool = false)
    all(l -> 0 < l < 1, levels) || throw(ArgumentError(
        "levels are fractions of the distribution, between 0 and 1, not $levels"))
    t, vals = _trajectory_values(p, species, baseline)
    w0 = _weight_matrix(p, weights)
    w = [size(vals[i], 1) > 0 ? Float64(w0[i]) : 0.0 for i in CartesianIndices(vals)]
    sum(w) > 0 || throw(ArgumentError("these weights weigh no pair with a finished cell"))
    used = findall(>(0), w)
    n_failed = sum(count(c -> c.trajectory === nothing, p.cells[i])
                   for i in eachindex(w0) if w0[i] > 0)

    cell_w = reduce(vcat, [fill(w[i] / size(vals[i], 1), size(vals[i], 1)) for i in used])
    rowmass = vec(sum(w; dims = 2))
    M, K = size(w)
    out = map(eachindex(t)) do a
        mu = [size(vals[i], 1) > 0 ? mean(view(vals[i], :, a)) : 0.0 for i in CartesianIndices(vals)]
        rowmean = [rowmass[m] > 0 ? sum(w[m, k] * mu[m, k] for k in 1:K) / rowmass[m] : 0.0
                   for m in 1:M]
        pooled = reduce(vcat, [view(vals[i], :, a) for i in used])
        (; pairs = _layer(vec(mu), vec(w), levels), rows = _layer(rowmean, rowmass, levels),
         cells = _layer(pooled, cell_w, levels))
    end
    return (; t, n_failed, pairs = _layer_out([o.pairs for o in out], levels),
            rows = _layer_out([o.rows for o in out], levels),
            cells = _layer_out([o.cells for o in out], levels))
end

"""
    trajectory_split(p::Prediction, species; weights = first(keys(p.weights)),
                     baseline = false)
        -> (; t, held_share = (; i_first, h_first), total, mc, flagged, mc_undefined, dropped)

[`variance_split`](@ref) at each recorded time after t = 0, on the pair means μ̂(t)
with their noise τ(t) = σ̂²(t)/n. Each result is a vector over `t`. With
`baseline = true` and the metered export, the last time is `variance_split`'s result
for Q. `weights` is as for [`trajectory_bands`](@ref).

Early in a run the cells' noise can outweigh the spread, so a component comes out negative
and is reported as 0. A share there can exceed 1, or be `Inf` or `NaN` when `total` is 0.
`flagged` says some component did, and `mc`, `mc_undefined` and `dropped` are
`variance_split`'s, at each time. A share lies between 0 and 1 where `flagged` and
`mc_undefined` are false and `total > 0`, though very sparse weights can still leave it `NaN`.
That does not make it signal: `total` can still be cell noise, so read a share only where
`total` also exceeds `mc`.
"""
function trajectory_split(p::Prediction, species::Symbol; weights = first(keys(p.weights)),
                          baseline::Bool = false)
    t, vals = _trajectory_values(p, species, baseline)
    w = _weight_matrix(p, weights)
    times = findall(>(0), t)
    splits = map(times) do a
        mom = [_moments(view(X, :, a)) for X in vals]
        mu, s2, nok = first.(mom), getindex.(mom, 2), last.(mom)
        τ = p.forward === :hybrid ? s2 ./ nok : zeros(size(mu))
        _variance_split(mu, τ, w; held = !isempty(p.held))
    end
    return (; t = t[times],
            held_share = (; i_first = [s.held_share.i_first for s in splits],
                          h_first = [s.held_share.h_first for s in splits]),
            total = [s.total for s in splits], mc = [s.mc for s in splits],
            flagged = [any(values(s.flagged)) for s in splits],
            mc_undefined = [s.mc_undefined for s in splits],
            dropped = [s.dropped for s in splits])
end

export trajectory_bands, trajectory_split
