# KIM-95, item 3: the scores (spec/prereg/kim-92.md §4), read only from the plain TSVs
# that corea_predict_kim95.jl's `collect` writes. It needs no InferCell and no .jls, so
# it reruns at any commit on the archive alone.
#
# Usage: julia --project dev/scripts/corea_predict_kim95_score.jl <dir>
#        julia --project dev/scripts/corea_predict_kim95_score.jl --selfcheck
# Writes <dir>/score.md and <dir>/score_targets.tsv.
#
# For each target, the predictive is its draws' μ̂ (predictive_<arm>_mu.tsv) and the
# truth y its mean over the truth cells (truth.tsv). Declared before any run:
#   - hit50, hit90: y inside the central 50% / 90% interval, by Statistics.quantile
#     (type 7) at 0.25–0.75 and 0.05–0.95, inclusive;
#   - CRPS: the sample CRPS, E|X − y| − ½ E|X − X′| over the draws (the V-statistic,
#     computed exactly by sorting);
#   - skill: 1 − CRPS_post / CRPS_prior;
#   - per species (34: 17 pools, 17 transcripts), over its forecast saves: the hit
#     fractions, the aggregate skill 1 − ΣCRPS_post / ΣCRPS_prior, and the median of
#     the per-save skills; per family, the hit fraction over all its saves and the
#     median of its species' aggregate skills;
#   - beside each, the truth's SE over the posterior predictive's SD, the rule's check,
#     and each arm's MC share: mean(σ̂²/n) over the draws, over var(μ̂) across them,
#     the part of the predictive variance that is cell noise rather than θ (KIM-69's
#     rule asks ≤ 0.1; it can exceed 1 by chance when Var_θ is near zero).
# A draw whose μ̂ is NaN (no finished cell) is left out of that target and counted.

using Printf
using Random
using Statistics: mean, std, var, median, quantile

function crps(x::AbstractVector{<:Real}, y::Real)
    s = sort(x)
    n = length(s)
    e1 = sum(abs(v - y) for v in s) / n
    e2 = 2 * sum((2i - n - 1) * s[i] for i in 1:n) / n^2     # E|X − X′|
    return e1 - e2 / 2
end

hits(x, y, lo, hi) = quantile(x, lo) <= y <= quantile(x, hi)

# CRPS of N(μ, σ) at y, in closed form.
Φ(z) = 0.5 * erfc_(-z / sqrt(2))
φ(z) = exp(-z^2 / 2) / sqrt(2π)
crps_normal(μ, σ, y) = (z = (y - μ) / σ; σ * (z * (2Φ(z) - 1) + 2φ(z) - 1 / sqrt(π)))
# erfc without SpecialFunctions: Numerical Recipes' erfcc, |relative error| < 1.2e-7.
function erfc_(x)
    z = abs(x)
    t = 1 / (1 + z / 2)
    r = t * exp(-z^2 - 1.26551223 + t * (1.00002368 + t * (0.37409196 + t * (0.09678418 +
        t * (-0.18628806 + t * (0.27886807 + t * (-1.13520398 + t * (1.48851587 +
        t * (-0.82215223 + t * 0.17087277)))))))))
    return x >= 0 ? r : 2 - r
end

function selfcheck()
    rng = Xoshiro(95)
    x = randn(rng, 200_000)
    for y in (0.0, 0.7, -2.5)
        a, b = crps(x, y), crps_normal(0.0, 1.0, y)
        abs(a - b) < 5e-3 || error("CRPS of N(0,1) at $y: sample $a, closed form $b")
        println(@sprintf("selfcheck: CRPS N(0,1) at %+.1f: sample %.5f, closed form %.5f", y, a, b))
    end
    crps([3.0], 1.0) == 2.0 || error("a point mass's CRPS is not |x − y|")
    crps(fill(3.0, 5), 1.0) == 2.0 || error("a repeated point's CRPS is not |x − y|")
    v = collect(1.0:100.0)
    hits(v, 50.0, 0.25, 0.75) && !hits(v, 99.0, 0.05, 0.95) && hits(v, 95.0, 0.05, 0.95) &&
        !hits(v, 10.0, 0.25, 0.75) || error("the interval hits are wrong")
    println("selfcheck: point mass, intervals: ok")
end

function read_wide(path)
    lines = readlines(path)
    head = split(lines[1], '\t')
    X = permutedims(reduce(hcat, [parse.(Float64, split(l, '\t')) for l in lines[2:end]]))
    return String.(head[4:end]), X[:, 4:end]
end

function read_table(path)
    lines = readlines(path)
    head = split(lines[1], '\t')
    rows = [split(l, '\t') for l in lines[2:end]]
    return Dict(String(h) => [String(r[i]) for r in rows] for (i, h) in enumerate(head))
end

f(x) = isfinite(x) ? @sprintf("%.4g", x) : "—"
pct(x) = @sprintf("%.0f%%", 100x)

function score(dir)
    ids, P = read_wide(joinpath(dir, "predictive_posterior_mu.tsv"))
    idp, R = read_wide(joinpath(dir, "predictive_prior_mu.tsv"))
    ids == idp || error("the arms' targets differ")
    _, P2 = read_wide(joinpath(dir, "predictive_posterior_s2.tsv"))
    _, Pn = read_wide(joinpath(dir, "predictive_posterior_n.tsv"))
    _, R2 = read_wide(joinpath(dir, "predictive_prior_s2.tsv"))
    _, Rn = read_wide(joinpath(dir, "predictive_prior_n.tsv"))
    function mcshare(M, S2, N, j)
        ok = [i for i in axes(M, 1) if isfinite(M[i, j]) && isfinite(S2[i, j]) && N[i, j] >= 2]
        return mean(S2[ok, j] ./ N[ok, j]) / var(M[ok, j])
    end
    tr = read_table(joinpath(dir, "truth.tsv"))
    tg = read_table(joinpath(dir, "targets.tsv"))
    tr["target"] == ids == tg["target"] || error("truth.tsv or targets.tsv is not in the predictive's order")
    y = parse.(Float64, tr["mean"])
    se = parse.(Float64, tr["se"])
    rows = map(eachindex(ids)) do j
        xp = filter(isfinite, P[:, j])
        xr = filter(isfinite, R[:, j])
        cp, cr = crps(xp, y[j]), crps(xr, y[j])
        (; id = ids[j], family = tg["family"][j], species = tg["species"][j], t = tg["t"][j],
         y = y[j], se = se[j], n_post = length(xp), n_prior = length(xr),
         nan_post = size(P, 1) - length(xp), nan_prior = size(R, 1) - length(xr),
         post_lo = quantile(xp, 0.05), post_med = median(xp), post_hi = quantile(xp, 0.95),
         prior_lo = quantile(xr, 0.05), prior_med = median(xr), prior_hi = quantile(xr, 0.95),
         post50 = hits(xp, y[j], 0.25, 0.75), post90 = hits(xp, y[j], 0.05, 0.95),
         prior50 = hits(xr, y[j], 0.25, 0.75), prior90 = hits(xr, y[j], 0.05, 0.95),
         crps_post = cp, crps_prior = cr, skill = 1 - cp / cr, se_over_sd = se[j] / std(xp),
         mc_post = mcshare(P, P2, Pn, j), mc_prior = mcshare(R, R2, Rn, j))
    end
    open(joinpath(dir, "score_targets.tsv"), "w") do io
        ks = keys(rows[1])
        println(io, join(ks, '\t'))
        for r in rows
            println(io, join([v isa AbstractFloat ? @sprintf("%.17g", v) : string(v) for v in values(r)], '\t'))
        end
    end
    open(joinpath(dir, "score.md"), "w") do io
        println(io, "# KIM-95 scores\n\nFrom the TSVs beside this file. ", size(P, 1), " posterior and ",
                size(R, 1), " prior draws; one truth, so these are scores, not coverage rates.\n")
        println(io, "## Q and the knockdown\n")
        println(io, "| target | truth ± SE | posterior 5–50–95% | prior 5–50–95% | post 50/90 | prior 50/90 | CRPS post | CRPS prior | skill | SE/SD | MC share post | MC share prior |")
        println(io, "|---|---|---|---|---|---|---|---|---|---|---|---|")
        yn(b) = b ? "in" : "out"
        for r in rows
            r.family == "Q" || continue
            println(io, "| ", r.id, " | ", f(r.y), " ± ", f(r.se), " | ", f(r.post_lo), " – ", f(r.post_med),
                    " – ", f(r.post_hi), " | ", f(r.prior_lo), " – ", f(r.prior_med), " – ", f(r.prior_hi),
                    " | ", yn(r.post50), "/", yn(r.post90), " | ", yn(r.prior50), "/", yn(r.prior90), " | ",
                    f(r.crps_post), " | ", f(r.crps_prior), " | ", f(r.skill), " | ", f(r.se_over_sd), " | ",
                    f(r.mc_post), " | ", f(r.mc_prior), " |")
        end
        fams = ("pool", "transcript")
        println(io, "\n## The later horizon, by family\n")
        println(io, "| family | species | saves | post hit 50 | post hit 90 | prior hit 50 | prior hit 90 | median species skill | SE/SD median | SE/SD max | MC share post (median) | MC share prior (median) |")
        println(io, "|---|---|---|---|---|---|---|---|---|---|---|---|")
        bysp(fam) = unique(r.species for r in rows if r.family == fam)
        agg(rs) = 1 - sum(r.crps_post for r in rs) / sum(r.crps_prior for r in rs)
        for fam in fams
            rs = [r for r in rows if r.family == fam]
            sk = [agg([r for r in rs if r.species == sp]) for sp in bysp(fam)]
            println(io, "| ", fam, " | ", length(bysp(fam)), " | ", length(rs), " | ",
                    pct(mean(r.post50 for r in rs)), " | ", pct(mean(r.post90 for r in rs)), " | ",
                    pct(mean(r.prior50 for r in rs)), " | ", pct(mean(r.prior90 for r in rs)), " | ",
                    f(median(sk)), " | ", f(median(r.se_over_sd for r in rs)), " | ",
                    f(maximum(r.se_over_sd for r in rs)), " | ", f(median(r.mc_post for r in rs)), " | ",
                    f(median(r.mc_prior for r in rs)), " |")
        end
        println(io, "\n## The later horizon, by species\n")
        println(io, "| species | family | saves | post hit 50 | post hit 90 | prior hit 50 | prior hit 90 | skill (aggregate) | skill (median per save) |")
        println(io, "|---|---|---|---|---|---|---|---|---|")
        for fam in fams, sp in bysp(fam)
            rs = [r for r in rows if r.species == sp]
            println(io, "| ", sp, " | ", fam, " | ", length(rs), " | ", pct(mean(r.post50 for r in rs)), " | ",
                    pct(mean(r.post90 for r in rs)), " | ", pct(mean(r.prior50 for r in rs)), " | ",
                    pct(mean(r.prior90 for r in rs)), " | ", f(agg(rs)), " | ",
                    f(median(r.skill for r in rs)), " |")
        end
        nn = sum(r.nan_post + r.nan_prior for r in rows)
        println(io, "\nDraw-target values with no finished cell, left out: ", nn, ".")
    end
    println("score: ", length(ids), " targets to ", joinpath(dir, "score.md"))
end

function main(args)
    selfcheck()
    isempty(args) || args[1] == "--selfcheck" || score(args[1])
end

if abspath(PROGRAM_FILE) == @__FILE__
    main(ARGS)
end
