# KIM-178: what the truths with a large promoter offset share. Exploratory, not gated
# (prereg §6 and Tom's call of 2026-11-28, after gate 2 was read). Reads the run directory
# or its archive and changes nothing.
#
# Usage: julia --project dev/scripts/corea_replicates/corea_replicates_diagnose.jl <out>
#
# Per replicate: where the truth sits in the prior (in prior SDs), the errors of the
# non-promoter coordinates, the size of δ and S, and from (−H)⁻¹ the correlation between
# the promoters' mean and each other coordinate. Then the correlation across the ten of
# each of these with the promoters' mean z. Ten points: a pointer for a hypothesis only.

include(joinpath(@__DIR__, "..", "pilot_common.jl"))
using Printf
using LinearAlgebra
using Statistics: mean, std, cor

const SEEDS = 17_200_002:17_200_011

pilot_setup(ds) = mixed_setup(ds; models = d11_models(), targets = COREA_TARGETS)

function row(out, seed)
    dir = joinpath(out, "seed_$seed")
    ds = load(dir, "dataset"); p = load(dir, "prep"); l = load(dir, "laplace")
    j = deserialize(joinpath(dir, "fit.jls")).result
    s = pilot_setup(ds)
    med, sd = mixed_prior_median(s), mixed_prior_sd(s)
    names = string.(j.names)
    idx(n) = findfirst(==(n), names)
    pr = j.prom
    tz = (j.truth .- med) ./ sd                  # the truth, in prior SDs from the median
    C = inv(Symmetric(-l.H))
    w = zeros(length(j.x)); w[pr] .= 1 / length(pr)
    vbar = w' * C * w                            # the promoters' mean, its posterior variance
    corr_with(k) = (w' * C[:, k]) / sqrt(vbar * C[k, k])
    others = ["log_krnadeg", "log_kcatF_R_ENO", "log_kcatF_R_FBA", "log_sigma_b"]
    return (; seed, m = j.m, m_ref = j.m_ref, nneg = j.nneg,
            tz_prom = mean(tz[pr]), tz_prom_sd = std(tz[pr]),
            tz = Dict(n => tz[idx(n)] for n in others),
            z = Dict(n => j.z[idx(n)] for n in others),
            err_prom = mean(j.x[pr] .- j.truth[pr]),           # in log units
            err_krnadeg = j.x[idx("log_krnadeg")] - j.truth[idx("log_krnadeg")],
            corr = Dict(n => corr_with(idx(n)) for n in others),
            maxdelta = maximum(abs, p.delta), trS = tr(p.S),
            sd_prom = mean(l.lap_sd[pr]), lp_gain = load(dir, "map").fit.logdensity - load(dir, "map").lp_ref)
end

function main(out)
    rows = [row(out, s) for s in SEEDS]
    others = ["log_krnadeg", "log_kcatF_R_ENO", "log_kcatF_R_FBA", "log_sigma_b"]
    short = Dict("log_krnadeg" => "krnadeg", "log_kcatF_R_ENO" => "ENO",
                 "log_kcatF_R_FBA" => "FBA", "log_sigma_b" => "σ_b")
    println("# KIM-178 diagnostic: ", out, "\n")
    println("## The truth in prior SDs, and the MAP's errors\n")
    println("| seed | mean z | promoters' truth (mean, SD) | ", join(["truth " * short[n] for n in others], " | "),
            " | ", join(["z " * short[n] for n in others], " | "),
            " | promoter error (log) | krnadeg error (log) | mean Laplace SD | max abs δ | tr S | Δ log p̃ |")
    println("|", repeat("---|", 5 + 2length(others) + 5))
    for r in rows
        @printf("| %d | %+.3f | %+.2f, %.2f | %s | %s | %+.4f | %+.4f | %.4f | %.3g | %.3g | %.2f |\n",
                r.seed, r.m, r.tz_prom, r.tz_prom_sd,
                join([@sprintf("%+.2f", r.tz[n]) for n in others], " | "),
                join([@sprintf("%+.2f", r.z[n]) for n in others], " | "),
                r.err_prom, r.err_krnadeg, r.sd_prom, r.maxdelta, r.trS, r.lp_gain)
    end
    println("\n## Posterior correlation of the promoters' mean with each coordinate, from (−H)⁻¹\n")
    println("| seed | ", join([short[n] for n in others], " | "), " |")
    println("|", repeat("---|", 1 + length(others)))
    for r in rows
        @printf("| %d | %s |\n", r.seed, join([@sprintf("%+.3f", r.corr[n]) for n in others], " | "))
    end
    m = [r.m for r in rows]
    feats = vcat([("promoters' truth, mean", [r.tz_prom for r in rows]),
                  ("promoters' truth, SD", [r.tz_prom_sd for r in rows])],
                 [("truth " * short[n], [r.tz[n] for r in rows]) for n in others],
                 [("z " * short[n], [r.z[n] for r in rows]) for n in others],
                 [("max abs δ", [r.maxdelta for r in rows]), ("tr S", [r.trS for r in rows]),
                  ("Δ log p̃ (MAP − θ_ref)", [r.lp_gain for r in rows]),
                  ("mean z at θ_ref", [r.m_ref for r in rows])])
    println("\n## Correlation across the ten with the promoters' mean z (exploratory, n = 10)\n")
    println("| feature | r |\n|---|---|")
    for (name, v) in feats
        @printf("| %s | %+.2f |\n", name, cor(v, m))
    end
end

main(abspath(ARGS[1]))
