# spec/phases/16-recovery.md task 16a.5: block 1's pieces, on short M0 runs.
# V4 itself, the chain against a direct-replay grid on three full M0 cells, is
# dev/scripts/block1_16a5.jl.

using Random
using Statistics
using Test
using InferCell
using Distributions: LogNormal, Normal, logpdf, params

# Record `n` M0 cells over `H` s at `truth`, with the panel's carry-inclusive
# counts at every 60 s save, and a snapshot at t = 0 with only the stochastic
# parameters written.
function _m0_cells(ms, truth, seeds, H)
    names = InferCell._block_names(ms, :ode)
    rows = [findfirst(==(s), names) for s in M0_PANEL]
    cells = map(seeds) do s
        Random.seed!(s)
        d = build_m0(tspan = (0.0, H))
        set_parameters!(d, ms, truth.values)
        Random.seed!(s)
        record_path!(d)
        latent = zeros(length(rows), round(Int, H) ÷ 60)
        for k in axes(latent, 2)
            for _ in 1:60
                handshake_step!(d)
            end
            latent[:, k] = [d.ode.u[i] * d.factor + d.rounding.remainders[i] for i in rows]
        end
        (path = recorded_path(d), latent = latent)
    end
    base = build_m0(tspan = (0.0, H))
    cme = [p for p in truth.values if !(first(p) in (:kcatF_R_ENO, :kcatR_R_ENO,
                                                      :kcatF_R_FBA, :kcatR_R_FBA))]
    set_parameters!(base, ms, cme)
    return base, cells
end

_batch_se(x; nb = 50) = (L = length(x) ÷ nb;
                         std([mean(@view x[(b - 1) * L + 1:b * L]) for b in 1:nb]) / sqrt(nb))

@testset "16a.5: block 1's pieces" begin
    ms = m0_models()
    truth = draw_truth(Xoshiro(1605), ms; names = M0_TARGETS, purpose = :calibration)
    tv = Dict(truth.values)
    fw = [:kcatF_R_ENO, :kcatF_R_FBA]

    @testset "reverse constants are derived as the truth's were" begin
        w = Dict(derived_ode_values(ms, [f => tv[f] for f in fw]))
        @test w[:kcatR_R_ENO] ≈ tv[:kcatR_R_ENO] rtol = 1e-12
        @test w[:kcatR_R_FBA] ≈ tv[:kcatR_R_FBA] rtol = 1e-12
        @test all(r -> r.relative <= r.bound, assert_haldane(ms, collect(w)))
    end

    H = 120.0
    base, rec = _m0_cells(ms, truth, [16051, 16052], H)
    cells = [Block1Cell(snapshot(base), c.path, c.latent) for c in rec]
    b = Block1(ms, cells; forwards = fw, panel = M0_PANEL)
    u_true = [log(tv[f]) for f in fw]

    @testset "the prior is on the forward constants only, with its Jacobian" begin
        @test_throws ArgumentError Block1(ms, cells; forwards = [:kcatF_R_ENO, :kcatR_R_ENO],
                                          panel = M0_PANEL)
        priors = [InferCell._parameter_by_name(ms, f).prior for f in fw]
        for u in (u_true, u_true .+ 0.3, u_true .- [0.1, 0.5])
            @test block1_log_prior(b, u) ≈
                  sum(logpdf(p, exp(x)) + x for (p, x) in zip(priors, u)) rtol = 1e-12
        end
    end

    @testset "the likelihood is assembled from the replay" begin
        # Observed equal to the replayed latent leaves no residual at the truth.
        lp, ss, n = block1_replay(b, u_true)
        @test ss == 0.0 && n == 2 * length(M0_PANEL) * 2
        # The path term is the replay's own density at the full truth.
        ref = sum(rec) do c
            d = build_m0(tspan = (0.0, H))
            set_parameters!(d, ms, truth.values)
            path_logdensity(replay!(d, c.path, round(Int, H); density = true))
        end
        @test lp ≈ ref rtol = 1e-12
        # Away from the truth the residuals appear.
        _, ss2, _ = block1_replay(b, u_true .+ [0.3, 0.0])
        @test ss2 > 0
        σ = 0.1
        @test block1_logtarget(b, u_true, σ) ≈
              block1_log_prior(b, u_true) + lp - n * log(σ) - n * log(2π) / 2 rtol = 1e-12
    end

    @testset "the σ step samples its exact conditional" begin
        f = sigma_log_conditional(b, 2.0, 200)
        grid = range(log(0.02), log(0.5); length = 4001)
        lw = f.(grid)
        w = exp.(lw .- maximum(lw))
        cdf = cumsum(w) ./ sum(w)
        rng = Xoshiro(16053)
        v = log(0.1)
        draws = [v = slice_step(rng, f, v) for _ in 1:40_000]
        for p in (0.05, 0.25, 0.5, 0.75, 0.95)
            q = grid[findfirst(>=(p), cdf)]
            ind = Float64.(draws .<= q)
            @test abs(mean(ind) - p) <= 3 * _batch_se(ind) + 1e-3
        end
    end

    @testset "one sweep runs, keeps Haldane, and can hold σ" begin
        st = Block1State(b, u_true, 0.1)
        block1_update!(Xoshiro(16054), b, st; update_sigma = false)
        @test st.σ == 0.1 && all(isfinite, st.u)
        block1_update!(Xoshiro(16055), b, st)
        @test st.σ != 0.1
    end
end
