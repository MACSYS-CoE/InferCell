# spec/phases/16-recovery.md task 16a.7: is the true path a typical draw from
# the window proposal?
#
# csmc_collapse_16a7 found the full-scale reference at ranks 0.03 to 0.10 on
# M_pi_c, M_3pg_c and M_2pg_c in windows 2 and 3. That could be chance, or the
# proposal could be wrong at full scale. At each window, K proposals start from
# the true path's own state at the window's start and are weighted by the
# transcript part of the incremental weight alone, the metabolite term left
# out. The weighted proposal then targets p(window path | start, transcripts),
# which the true window path is a draw from. So the reference's weighted rank of
# each panel pool at 60m is uniform when the proposal and weight are right.
#
# Ranks are pooled over windows, and cells for M0, per pool; a Kolmogorov–Smirnov
# distance against the uniform is reported with its asymptotic p-value. Pools are
# correlated within a window, so the 17 tests are not independent.
#
# Task 16a.7b (§12 2026-10-02) runs it on five independent full-scale cells,
# 30001 to 30005, one per array task, and passes each pool whose pooled ranks are
# uniform by KS at p > 0.01 after Bonferroni over the 17 pools.
#
# Usage: julia --project dev/scripts/csmc_rank_16a7.jl run <scale> <cell> <K>
#        julia --project dev/scripts/csmc_rank_16a7.jl merge <scale> <K>
# Regenerate: sbatch --array=1-5 dev/scripts/csmc_rank_16a7.slurm run full 50, then
# sbatch dev/scripts/csmc_rank_16a7.slurm merge full 50.

include(joinpath(@__DIR__, "csmc_variants_16a7.jl"))

latent(d, rows) = [d.ode.u[i] * d.factor + d.rounding.remainders[i] for i in rows]

# Asymptotic Kolmogorov p-value for distance D on n points.
function ks_p(D, n)
    λ = (sqrt(n) + 0.12 + 0.11 / sqrt(n)) * D
    return clamp(2 * sum((-1)^(k - 1) * exp(-2 * k^2 * λ^2) for k in 1:100), 0.0, 1.0)
end

function ks(u)
    s = sort(u)
    n = length(s)
    return maximum(max(i / n - s[i], s[i] - (i - 1) / n) for i in 1:n)
end

# The weighted ranks of one cell's windows: pools × windows.
function cell_ranks(scale, c, K)
    out = Vector{Float64}[]                  # per window, per pool rank
    x = cell(scale, c)
    rows = x.spec.panel_rows
    met(d, m) = InferCell._metabolite_loglik(d, x.spec, view(x.data.metabolites, :, m))
    d = restore(x.base)
    y0 = [d.jump.u[s] for s in x.spec.tmap.states]
    @printf("%s cell %d: %d windows, %d proposals each\n", scale, x.seed, x.T, K)
    flush(stdout)
    for m in 1:x.T
        yprev = m == 1 ? y0 : x.data.transcripts[:, m - 1]
        start = snapshot(d)
        advance_window!(Xoshiro(1), d, x.spec, yprev, x.data.transcripts[:, m],
                        view(x.data.metabolites, :, m); fixed = x.ref[m])
        ref = latent(d, rows)
        lw = Float64[]
        vals = Vector{Float64}[]
        for k in 1:K
            e = restore(start)
            l, _ = advance_window!(Xoshiro(16_077_000 + 100_000c + 1000m + k), e, x.spec, yprev,
                                   x.data.transcripts[:, m], view(x.data.metabolites, :, m))
            isfinite(l) || continue
            push!(lw, l - met(e, m))          # the transcript part alone
            push!(vals, latent(e, rows))
        end
        w = exp.(lw .- maximum(lw))
        w ./= sum(w)
        # Weighted rank, ties split evenly, so a pool that never moves gives 0.5.
        r = [sum(w[k] * ((v[j] < ref[j]) + 0.5 * (v[j] == ref[j])) for (k, v) in enumerate(vals))
             for j in eachindex(rows)]
        push!(out, r)
        m % 10 == 0 && (@printf("  window %d done\n", m); flush(stdout))
    end
    R = reduce(hcat, out)
    serialize(joinpath(DIR, "ranks_$(scale)_c$(c)_K$(K).jls"), R)
    return R
end

function report(R, K; bonferroni = false)
    n = size(R, 1)
    @printf("\n%d windows; weighted rank of the true path among %d proposals, per pool:\n", size(R, 2), K)
    ok = true
    for (j, s) in enumerate(M0_PANEL)
        u = R[j, :]
        D = ks(u)
        p = ks_p(D, length(u))
        pass = !bonferroni || p * n > 0.01
        ok &= pass
        @printf("  %-14s mean %.3f; below 0.1: %d, above 0.9: %d; KS D = %.3f, p = %.3g%s\n",
                s, mean(u), count(<(0.1), u), count(>(0.9), u), D, p,
                bonferroni ? @sprintf(", Bonferroni p = %.3g%s", min(1, p * n), pass ? "" : "  FAIL") : "")
    end
    bonferroni && println(ok ? "16a.7b passes" : "16a.7b FAILS")
end

if ARGS[1] == "run"
    report(cell_ranks(ARGS[2], parse(Int, ARGS[3]), parse(Int, ARGS[4])), parse(Int, ARGS[4]))
elseif ARGS[1] == "merge"
    scale, K = ARGS[2], parse(Int, ARGS[3])
    files = sort(filter(f -> startswith(f, "ranks_$(scale)_c") && endswith(f, "_K$(K).jls"), readdir(DIR)))
    println("merging ", join(files, ", "))
    report(reduce(hcat, [deserialize(joinpath(DIR, f)) for f in files]), K; bonferroni = true)
else
    error("stage must be run or merge")
end
