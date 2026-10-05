# spec/phases/16-recovery.md task 16a.7a: inside one annealed window.
#
# For the first W windows of an M0 cell, from the true path's own state at each
# window's start, run C forward AIS chains at each K, and the reference's
# backward chain once. Reports, per K:
# - the MH acceptance rate by fifth of the schedule;
# - each chain's final log weight w against the reference's, which is where the
#   chain should land if it reaches the posterior;
# - the AIS log weights against the reference's. Under particle Gibbs with N
#   particles, a fresh particle replaces the reference with probability about
#   (N − 1)·mean(exp(logW))/(exp(logW_ref) + (N − 1)·mean(exp(logW))).
#
# Usage: julia --project dev/scripts/csmc_anneal_diag_16a7.jl <cell> <W> <C> <K...>
# Regenerate: sbatch dev/scripts/csmc_anneal_diag_16a7.slurm <cell> <W> <C> <K...>

include(joinpath(@__DIR__, "csmc_variants_16a7.jl"))

# The forward chain, as annealed_window runs it, with each stage's outcome kept.
function forward(rng, d, spec, yprev, yend, ymet, sched)
    e = restore(d)
    lwx, x = advance_window!(rng, e, spec, yprev, yend, ymet)
    isfinite(lwx) || return -Inf, -Inf, Bool[]
    logw, βprev = 0.0, 0.0
    acc = Bool[]
    for β in sched
        logw += (β - βprev) * lwx
        x2, lw2, e2 = InferCell._stage_move(rng, d, spec, yprev, yend, ymet, x, lwx, β)
        push!(acc, e2 !== nothing)
        x, lwx = x2, lw2
        βprev = β
    end
    return logw, lwx, acc
end

function diag(c, W, C, Ks)
    x = cell("m0", c)
    d = restore(x.base)
    y0 = [d.jump.u[s] for s in x.spec.tmap.states]
    @printf("m0 cell %d, σ = %.3g; %d chains per K; schedules geometric from β1 = %g\n",
            x.seed, x.σ, C, BETA1)
    for m in 1:W
        yprev = m == 1 ? y0 : x.data.transcripts[:, m - 1]
        yend, ymet = x.data.transcripts[:, m], view(x.data.metabolites, :, m)
        start = snapshot(d)
        lw_ref, _ = advance_window!(Xoshiro(1), d, x.spec, yprev, yend, ymet; fixed = x.ref[m])
        @printf("\nwindow %d: reference log w = %.1f\n", m, lw_ref)
        for K in Ks
            sched = geometric_schedule(K; β1 = BETA1)
            rs = [forward(Xoshiro(16_078_000 + 10_000m + 100K + i), start, x.spec, yprev, yend, ymet, sched)
                  for i in 1:C]
            acc = reduce(hcat, [r[3] for r in rs if !isempty(r[3])])
            fifths = [mean(acc[(1 + (j - 1) * K ÷ 5):(j * K ÷ 5), :]) for j in 1:5]
            lwref_b, _, _ = annealed_window(Xoshiro(16_079_000 + 10_000m + K), start, x.spec, yprev,
                                            yend, ymet, sched; fixed = x.ref[m])
            fin = [r[2] for r in rs]
            lws = [r[1] for r in rs]
            mx = max(lwref_b, maximum(lws))
            share(N) = (N - 1) * mean(exp.(lws .- mx)) /
                       (exp(lwref_b - mx) + (N - 1) * mean(exp.(lws .- mx)))
            @printf("  K = %4d: acceptance by fifth %s\n", K, join([@sprintf("%.2f", f) for f in fifths], " "))
            @printf("            final log w − reference: median %+.1f, best %+.1f\n",
                    median(fin .- lw_ref), maximum(fin .- lw_ref))
            @printf("            AIS log weight − reference's: median %+.1f, best %+.1f; replacement chance at N = 2, 5: %.3f, %.3f\n",
                    median(lws .- lwref_b), maximum(lws .- lwref_b), share(2), share(5))
            flush(stdout)
        end
    end
end

diag(parse(Int, ARGS[1]), parse(Int, ARGS[2]), parse(Int, ARGS[3]), parse.(Int, ARGS[4:end]))
