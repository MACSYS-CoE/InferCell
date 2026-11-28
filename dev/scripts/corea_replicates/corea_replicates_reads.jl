# KIM-178: gate 2 and the result note's numbers, regenerated from the run directory or
# its archive. It reads the replicates' fit.jls, and kim92/kim92.jls when present, and
# changes nothing.
#
# Usage: julia --project dev/scripts/corea_replicates/corea_replicates_reads.jl <out>
#
# Prereg §3: the seeds are taken in order (the ten, then the spares), skipping failed fits,
# until ten fits succeed. Prereg §4: gate 2 on those ten.

using Serialization
using Printf
using Statistics: mean, std

const SEEDS = 17_200_002:17_200_011
const SPARES = 17_200_012:17_200_016
const T_MAX, NEG_MIN, N = 3.0, 9, 10

ld(path) = deserialize(path).result

function main(out)
    fits = Dict{Int, Any}()
    for seed in vcat(SEEDS, SPARES)
        f = joinpath(out, "seed_$seed", "fit.jls")
        isfile(f) && (fits[seed] = ld(f))
    end
    println("# KIM-178 reads: ", out, "\n")
    # "mean z at θ_ref" divides by the fitted MAP's Laplace SD, as the driver's judge does: no
    # Hessian exists at θ_ref on a replicate. Prereg §6 asked for θ_ref's own SD.
    println("| seed | fit | promoters' mean z | below the truth | correlation-aware SE | mean z / SE | mean z at θ_ref (fitted MAP's SD) | reason |")
    println("|---|---|---|---|---|---|---|---|")
    used = Int[]
    for seed in vcat(SEEDS, SPARES)
        haskey(fits, seed) || continue
        j = fits[seed]
        @printf("| %d | %s | %+.3f | %d of %d | %.3f | %+.1f | %+.3f | %s |\n", seed, j.ok ? "ok" : "FAILED",
                j.m, j.nneg, length(j.prom), j.se_corr, j.m / j.se_corr, j.m_ref,
                j.ok ? "" : join(j.reasons, "; "))
    end
    # Prereg §6, reported not gated: the optimisers, −H and the time per stage.
    println("\n| seed | θ_ref by | in s | ensemble s | fitted MAP by | in s | Laplace s | −H eigenvalues | asymmetry | box margin (prior SD) |")
    println("|---|---|---|---|---|---|---|---|---|---|")
    for seed in vcat(SEEDS, SPARES)
        haskey(fits, seed) || continue
        d = joinpath(out, "seed_$seed")
        p, mp, l = ld(joinpath(d, "prep.jls")), ld(joinpath(d, "map.jls")), ld(joinpath(d, "laplace.jls"))
        @printf("| %d | %s | %.0f | %.0f | %s | %.0f | %.0f | %.4g to %.4g | %.2g | %.2f |\n", seed,
                p.ref.method, p.seconds.map, p.seconds.ensemble, mp.fit.method, mp.seconds, l.seconds,
                l.eigenvalues[1], l.eigenvalues[end], l.asym, fits[seed].edge_min)
    end
    # One pass in order (prereg §3): a failed fit is skipped, a missing one stops the read.
    for seed in vcat(SEEDS, SPARES)
        length(used) == N && break
        if !haskey(fits, seed)
            println("\nGate 2 is not read: seed $seed has no fit yet, with $(length(used)) of $N succeeded before it.")
            return
        end
        fits[seed].ok && push!(used, seed)
    end
    if length(used) < N
        println("\nGate 2 is not read: the spares ran out with $(length(used)) of $N fits (prereg §3). A stop, for Tom.")
        return
    end
    m = [fits[s].m for s in used]
    nneg = count(<(0), m)
    mbar, sm = mean(m), std(m)
    t = mbar / (sm / sqrt(N))
    se = [fits[s].se_corr for s in used]
    a, b = nneg >= NEG_MIN, abs(t) > T_MAX
    println("\n## Gate 2 (prereg §4)\n")
    @printf("- Replicates: %s\n", join(used, ", "))
    @printf("- (a) mean z below zero in %d of %d (fails at %d or more): %s\n", nneg, N, NEG_MIN, a ? "FIRES" : "does not fire")
    @printf("- (b) pooled mean %+.3f, SD %.3f, SE %.3f, t = %+.2f on %d df (fails at |t| > %.0f): %s\n",
            mbar, sm, sm / sqrt(N), t, N - 1, T_MAX, b ? "FIRES" : "does not fire")
    println("- **Verdict: the transcript factor ", a || b ? "is MISCALIBRATED (gate 2 fails)" : "passes gate 2", ".**")
    # Not gated (Decision log 2026-11-28): the spread of the ten against their own SEs.
    @printf("\nNot gated: Σ(mean z / SE)² = %.1f on %d df; SD of the ten means / mean SE = %.2f; %d of %d beyond 3 SE.\n",
            sum(abs2, m ./ se), N, sm / mean(se), count(>(3), abs.(m ./ se)), N)
    k = joinpath(out, "kim92", "kim92.jls")
    if isfile(k)
        r = ld(k)
        println("\n## KIM-92's dataset (reported, not pooled)\n")
        @printf("- promoters' mean z: %+.3f by MAP and Laplace, %+.3f by NUTS, %d of %d below the truth\n",
                r.m, mean(r.z_nuts[r.prom]), r.nneg, length(r.prom))
        @printf("- Laplace SD / NUTS SD: %.3f to %.3f at the fitted MAP, %.3f to %.3f at θ_ref\n",
                extrema(r.ratio)..., extrema(r.ratio_ref)...)
        @printf("- fitted log density: %.4f at the MAP, %.4f at θ_ref\n", r.lp_map, r.lp_ref)
        lk = ld(joinpath(out, "kim92", "laplace.jls"))
        @printf("- −H at the fitted MAP: eigenvalues %.4g to %.4g, asymmetry %.2g\n",
                lk.eigenvalues[1], lk.eigenvalues[end], lk.asym)
        @printf("- promoters' mean z at θ_ref by θ_ref's own Laplace SD (KIM-92's metric): %+.3f\n",
                mean((r.z_ref .* lk.lap_sd ./ r.lap_ref)[r.prom]))   # (θ_ref − truth) / θ_ref's SD
    end
end

main(abspath(ARGS[1]))
