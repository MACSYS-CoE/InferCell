# KIM-66's done-when: the variance-split estimator on the crossed matrix
# (`variance_split`, src/predict_variance.jl), checked by simulation on known
# splits before M5 uses it (dev/notes/2026-10-29-cut-predictive-distribution.md
# §5.2). No cell runs; the matrices are made up.
#
# Usage: julia --project dev/scripts/variance_split_check_kim66.jl <outdir>
#   (dev/scripts/variance_split_check_kim66.sh submits it.)
#
# The design, declared before the run (M4 plan, 2026-10-30):
# - θ_I = x down the rows and θ_H = y across the columns, both standard normal.
#   μ(x, y) = o + a x + b y + g x y + q (x² − 1), and each μ̂ is the mean of n = 4
#   cells with noise sd 1, so τ = σ̂²/n is estimated, as the hybrid's is.
# - Three shapes of μ: "balanced" (a, b, g, q) = (1, 0.7, 0.5, 0.3), and
#   "held-dominant" (0.3, 1, 0.3, 0), both with o = 0; and "offset", balanced
#   with o = 10⁵. A real Q sits far from zero (lactate out is ~10⁵ particles),
#   and the split must not depend on where. (Added after a smoke pass at R = 30
#   showed the masked estimate moving with the mean; the estimator now centres
#   μ̂ first.)
# - Three masks, each a viability decision v_mk ~ Bernoulli(π(x, y)):
#   none (π = 1, the `uniform` weights); independent (π = 0.7); and dependent on
#   θ (π = logistic(1 + x + y), about 0.7 on average, falling with both groups).
#   Under a mask both conditioned variants are checked, `per_theta_I` and `joint`,
#   with `viability_weights`' weights.
# - Two sizes: the planned (M, K) = (100, 20), and a small (30, 10).
# - R = 2,000 replicates of each.
# - The truth under each law is computed on a 1,601-point grid in each of x and
#   y over [−8, 8] (trapezoid weights): the law is φ(x) φ(y) ρ(x, y), with ρ = 1
#   (uniform), π / π̄(x) (per_theta_I) or π / Π (joint).
#   Under per_theta_I, the estimator's within-row moments use rows with two or more
#   viable columns; the truth is the plain law, so that shortfall is part of the
#   measured bias.
# - **Pass:** for each of the five components (total; I first, between and
#   within; H first, between and within), the mean estimate is within ±10% of the
#   truth, the 10% taken of max(|truth|, 0.05 × total) so a component near zero
#   is judged on the total's scale. The MC standard error of each mean is printed
#   beside it.
# - If H first fails under the θ-dependent mask, the freeze reports I first only
#   (note §5.2's fallback).

using InferCell
using Random
using Statistics: mean, std, var
using Printf
using Dates

const OUT = ARGS[1]
const COMMIT = strip(read(`git rev-parse --short HEAD`, String))
const JOB = get(ENV, "SLURM_JOB_ID", "none")
const R = parse(Int, get(ENV, "R", "2000"))   # R=20 for a smoke pass
const N_CELLS = 4
const SHAPES = (balanced = (o = 0.0, a = 1.0, b = 0.7, g = 0.5, q = 0.3),
                held_dominant = (o = 0.0, a = 0.3, b = 1.0, g = 0.3, q = 0.0),
                offset = (o = 1e5, a = 1.0, b = 0.7, g = 0.5, q = 0.3))
const MASKS = (none = (x, y) -> 1.0, independent = (x, y) -> 0.7,
               dependent = (x, y) -> 1 / (1 + exp(-(1 + x + y))))
const SIZES = ((100, 20), (30, 10))

μf(s, x, y) = s.o + s.a * x + s.b * y + s.g * x * y + s.q * (x^2 - 1)

# The five components under the law φ(x) φ(y) ρ(x, y), on a grid.
function truth(s, π, law)
    g = range(-8, 8; length = 1601)
    h = step(g)
    φ = [exp(-x^2 / 2) / sqrt(2pi) * h for x in g]
    P = [π(x, y) for x in g, y in g]
    ρ = law === :uniform ? ones(size(P)) :
        law === :per_theta_I ? P ./ (P * φ) :           # π / π̄(x), π̄(x) = Σ_y φ(y) π
        P
    L = (φ * φ') .* ρ
    L ./= sum(L)
    m = [μf(s, x, y) for x in g, y in g]
    Em, Em2 = sum(L .* m), sum(L .* m .^ 2)
    rx, ry = vec(sum(L; dims = 2)), vec(sum(L; dims = 1))
    cx = vec(sum(L .* m; dims = 2)) ./ rx
    dy = vec(sum(L .* m; dims = 1)) ./ ry
    total = Em2 - Em^2
    ib = sum(rx .* cx .^ 2) - Em^2
    hb = sum(ry .* dy .^ 2) - Em^2
    return (; total, i_between = ib, i_within = total - ib, h_between = hb,
            h_within = total - hb)
end

function replicate(rng, s, π, M, K)
    x, y = randn(rng, M), randn(rng, K)
    cells = [μf(s, x[m], y[k]) .+ randn(rng, N_CELLS) for m in 1:M, k in 1:K]
    mu, τ = mean.(cells), var.(cells) ./ N_CELLS
    v = [rand(rng) < π(x[m], y[k]) for m in 1:M, k in 1:K]
    return mu, τ, v
end

components(e) = (; total = e.total, i_between = e.i_first.between, i_within = e.i_first.within,
                 h_between = e.h_first.between, h_within = e.h_first.within)

function main()
    mkpath(OUT)
    println("KIM-66 variance-split check at $COMMIT, job $JOB, $(now()), R = $R")
    flush(stdout)
    tsv = open(joinpath(OUT, "variance_split_check_kim66.tsv"), "w")
    println(tsv, "# commit $COMMIT, job $JOB, $(now()), R = $R, n_cells = $N_CELLS")
    println(tsv, join(["shape", "mask", "weights", "M", "K", "component", "truth", "mean",
                       "mc_se", "bias", "tolerance", "pass", "flagged_fraction"], '\t'))
    all_pass = true
    seed = 18_436_000                                  # M4's reserved block
    for (sname, s) in pairs(SHAPES), (mname, π) in pairs(MASKS), (M, K) in SIZES
        laws = mname === :none ? (:uniform,) : (:per_theta_I, :joint)
        seed += 1
        rng = Xoshiro(seed)
        est = Dict(l => NamedTuple[] for l in laws)
        flags = Dict(l => 0 for l in laws)
        for r in 1:R
            mu, τ, v = replicate(rng, s, π, M, K)
            for l in laws
                w = l === :uniform ? fill(1 / (M * K), M, K) : getproperty(viability_weights(v), l)
                e = InferCell._variance_split(mu, τ, w; held = true)
                push!(est[l], components(e))
                flags[l] += any(values(e.flagged))
            end
        end
        for l in laws
            t = truth(s, π, l)
            for c in keys(t)
                vals = [getproperty(e, c) for e in est[l]]
                m, se = mean(vals), std(vals) / sqrt(R)
                tol = 0.1 * max(abs(getproperty(t, c)), 0.05 * t.total)
                bias = m - getproperty(t, c)
                pass = abs(bias) <= tol
                all_pass &= pass
                println(tsv, join([sname, mname, l, M, K, c, @sprintf("%.6g", getproperty(t, c)),
                                   @sprintf("%.6g", m), @sprintf("%.3g", se),
                                   @sprintf("%.4g", bias), @sprintf("%.4g", tol), pass,
                                   @sprintf("%.4f", flags[l] / R)], '\t'))
                @printf("%-14s %-12s %-12s M=%3d K=%2d %-10s truth %8.4f mean %8.4f ± %.4f  bias %+8.4f (tol %.4f) %s\n",
                        sname, mname, l, M, K, c, getproperty(t, c), m, se, bias, tol,
                        pass ? "pass" : "FAIL")
            end
        end
        println("  $(now()) done $sname / $mname / $M × $K")
        flush(stdout)
    end
    close(tsv)
    println(all_pass ? "ALL PASS" : "SOME FAIL")
end

main()
