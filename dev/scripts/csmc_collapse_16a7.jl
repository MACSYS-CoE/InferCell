# spec/phases/16-recovery.md task 16a.7: why the path update's weights collapse.
#
# The variant smoke runs changed no window. A first diagnosis found the
# reference particle takes nearly all the weight, and that three panel pools,
# M_pi_c, M_3pg_c and M_2pg_c, carry the metabolite likelihood's gap. This asks
# which genes' translation explains those pools' deviations.
#
# For each of the first W windows, K proposals start from the reference's own
# state at the window's start, so they differ from the reference only within the
# window. For each proposal, at 60m:
# - each gene's protein count minus the reference's, the window's translation
#   difference;
# - each watched pool's signed deviation from the reference, in σ units of log;
# - its metabolite log-likelihood minus the reference's.
# Each response is regressed on the protein differences: the single best gene's
# correlation and R², and all genes' R² by least squares with an intercept.
#
# Usage: julia --project dev/scripts/csmc_collapse_16a7.jl <scale> <W> <K>
# Regenerate: sbatch dev/scripts/csmc_collapse_16a7.slurm <scale> <W> <K>

include(joinpath(@__DIR__, "csmc_variants_16a7.jl"))
using LinearAlgebra

const POOLS = [:M_pi_c, :M_3pg_c, :M_2pg_c]

latent(d, rows) = [d.ode.u[i] * d.factor + d.rounding.remainders[i] for i in rows]

function r2(X, y)
    A = hcat(ones(length(y)), X)
    β = A \ y
    r = y .- A * β
    return 1 - sum(abs2, r) / sum(abs2, y .- mean(y))
end

function collapse(scale, W, K)
    x = cell(scale, 1)
    tl = only(m for m in x.ms if m isa CoreATranslation)
    jnames = InferCell._block_names(x.ms, :jump)
    prot = [findfirst(==(protein_state(g.locus)), jnames) for g in tl.genes]
    loci = [g.locus for g in tl.genes]
    any(isnothing, prot) && error("a translated gene's protein state is missing")
    prow = [x.spec.panel_rows[findfirst(==(s), M0_PANEL)] for s in POOLS]
    met(d, m) = InferCell._metabolite_loglik(d, x.spec, view(x.data.metabolites, :, m))
    @printf("%s cell %d, σ = %.3g; %d translated genes; %d proposals per window from the reference's state\n",
            scale, x.seed, x.σ, length(loci), K)
    d = restore(x.base)
    y0 = [d.jump.u[s] for s in x.spec.tmap.states]
    for m in 1:W
        yprev = m == 1 ? y0 : x.data.transcripts[:, m - 1]
        start = snapshot(d)
        advance_window!(Xoshiro(1), d, x.spec, yprev, x.data.transcripts[:, m],
                        view(x.data.metabolites, :, m); fixed = x.ref[m])
        ref_p = [d.jump.u[i] for i in prot]
        ref_x = log.(max.(latent(d, prow), 1.0))
        ref_met = met(d, m)
        dp = zeros(K, length(prot))
        dx = zeros(K, length(POOLS))
        dmet = zeros(K)
        ok = trues(K)
        for k in 1:K
            e = restore(start)
            lw, _ = advance_window!(Xoshiro(16_076_000 + 1000m + k), e, x.spec, yprev,
                                    x.data.transcripts[:, m], view(x.data.metabolites, :, m))
            ok[k] = isfinite(lw)
            ok[k] || continue
            dp[k, :] = [e.jump.u[i] for i in prot] .- ref_p
            dx[k, :] = (log.(max.(latent(e, prow), 1.0)) .- ref_x) ./ x.σ
            dmet[k] = met(e, m) - ref_met
        end
        dp, dx, dmet = dp[ok, :], dx[ok, :], dmet[ok]
        moved = [j for j in axes(dp, 2) if std(dp[:, j]) > 0]
        @printf("\nwindow %d: %d of %d proposals finite; genes whose protein count varies: %d\n",
                m, count(ok), K, length(moved))
        for j in moved
            @printf("  %-16s protein − reference: mean %+.2f, sd %.2f\n", loci[j], mean(dp[:, j]), std(dp[:, j]))
        end
        responses = vcat([(string(s, " Δlog/σ"), dx[:, i]) for (i, s) in enumerate(POOLS)],
                         [("metabolite loglik − reference", dmet)])
        for (name, y) in responses
            cs = [cor(dp[:, j], y) for j in moved]
            b = argmax(abs.(cs))
            @printf("  %-32s mean %+8.2f, sd %7.2f; best gene %s, r = %+.3f (R² %.3f); all genes R² %.3f\n",
                    name, mean(y), std(y), loci[moved[b]], cs[b], cs[b]^2, r2(dp[:, moved], y))
        end
    end
end

collapse(ARGS[1], parse(Int, ARGS[2]), parse(Int, ARGS[3]))
