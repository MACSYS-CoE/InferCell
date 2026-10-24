# spec/phases/approximate-method.md task B5: the mixed likelihood, its verification
# clause (§12 2026-10-23, "B5 as planned"), and the Turing model over it.

using Test
using InferCell
using InferCell: _mixed_solve, _mixed_writes, AbstractSubModel
using ForwardDiff
using Distributions
using Random
using Statistics: mean, cov
using Turing: logjoint, sample, NUTS

@testset "B5: the mixed likelihood" begin
    H = 180.0
    ds = m0_dataset(16085; ncells = 3, horizon = H)
    s = mixed_setup(ds)
    truth = Dict(ds.truth.values)
    x = vcat([log(truth[n]) for n in s.names], log(ds.sigma))
    θ = exp.(x[1:length(s.names)])

    @testset "the setup" begin
        @test mixed_parameter_names(s) ==
              [:log_S_JCVISYN3A_0607, :log_S_JCVISYN3A_0779, :log_krnadeg,
               :log_kcatF_R_ENO, :log_kcatF_R_FBA, :log_sigma_b]
        @test s.switches == [0.0, 59.0, 119.0, 179.0]
        @test s.saves == [60.0, 120.0, 180.0]
        @test length(s.factors) == 3
        @test Set(s.forwards) == Set([:kcatF_R_ENO, :kcatF_R_FBA])
        # The record drops t = 0; every cell starts from the same fixed count.
        @test all(f -> f.y0 == s.factors[1].y0, s.factors)
        @test all(f -> f.kmax ≈ InferCell.TURNOVER_CEILING ./ (s.lengths .- 1), s.factors)
    end

    @testset "Haldane: the dual-safe reverse constants" begin
        ms = s.models
        fw = [:kcatF_R_ENO => 1.3 * truth[:kcatF_R_ENO]]
        plain = derived_ode_values(ms, fw)
        @test eltype(plain) == Pair{Symbol, Float64}
        g = ForwardDiff.gradient(z -> sum(last, derived_ode_values(ms, [:kcatF_R_ENO => z[1]])),
                                 [last(fw[1])])
        dual = derived_ode_values(ms, [:kcatF_R_ENO => ForwardDiff.Dual(last(fw[1]), 1.0)])
        @test first.(dual) == first.(plain)
        @test [ForwardDiff.value(last(p)) for p in dual] == last.(plain)
        @test g[1] > 0
        # Each reverse is linear in its forward.
        @test all(ForwardDiff.partials(last(p))[1] ≈ last(q) / last(fw[1])
                  for (p, q) in zip(dual, plain))
    end

    @testset "clock identity: the factor reads the constants the route ran" begin
        r = _mixed_solve(s, θ)
        k = route_transcription_rates(s, r.sol, r.p)
        @test size(k) == (2, 4)
        @test all(>(0), k)
        # Each column against the constants a solve ending at that switch's save,
        # 60m, left in its parameter vector, read from the route's own rebuild. A
        # route rebuilds at 60m − 1 only when its horizon reaches 60m. The last
        # column comes from the same solve, so it agrees to roundoff.
        slots(rt) = [rt.rebuilds[i].fill_idxs[j] for (i, j) in s.tx_slots]
        @test k[:, end] ≈ r.sol.prob.p[slots(s.route)] rtol = 1e-12
        for (c, t) in enumerate(s.switches)
            rt = build_mean_field_problem(s.models; tspan = (0.0, t + 1.0), complete = true)
            p = mean_field_parameters(rt, _mixed_writes(s, θ))
            sol = solve_mean_field(rt; p)
            @test k[:, c] ≈ sol.prob.p[slots(rt)] rtol = 1e-6
        end
        # And the constants switch: the windows do not all run at R₀.
        @test maximum(abs, k[:, end] ./ k[:, 1] .- 1) > 1e-6
    end

    @testset "bulk reduction: δ = 0 and S = 0 is block 1's lognormal" begin
        r = _mixed_solve(s, θ)
        m = mean_field_log_bulk(s, r.sol)
        @test length(m) == length(M0_PANEL) * 3
        σ = ds.sigma
        expected = sum(logpdf(Normal(m[i], σ), s.logz[i]) for i in eachindex(m))
        @test bulk_loglik(s, m, σ) ≈ expected atol = 1e-10
        # Against csmc's form, from the route's particles directly.
        rt = s.route
        row = 0
        alt = 0.0
        for t in s.saves
            u = r.sol.u[findfirst(==(t), r.sol.t)]
            f = InferCell._route_geometry(rt, u[rt.n_idx]).factor
            for (j, i) in enumerate(s.rows)
                alt += logpdf(Normal(log(max(u[i] * f, 1.0)), σ),
                              log(ds.bulk[j, findfirst(==(t), s.saves)]))
            end
        end
        @test bulk_loglik(s, m, σ) ≈ alt atol = 1e-10
    end

    @testset "delta-method S" begin
        rng = Xoshiro(3)
        cells = 10 .+ 5 .* rand(rng, 40, 2, 3)
        m_ref = log.(fill(11.0, 6))
        out = delta_and_S(cells, m_ref; C = 7)
        X = reshape(cells, 40, 6)
        x̄ = vec(mean(X; dims = 1))
        @test out.delta ≈ log.(x̄) .- m_ref
        for a in 1:6, b in 1:6
            @test out.S[a, b] ≈ cov(X[:, a], X[:, b]) / (x̄[a] * x̄[b] * 7)
        end
        # Pools fastest, as the bulk record's vec.
        @test X[:, 2] == cells[:, 2, 1] && X[:, 3] == cells[:, 1, 2]
        @test_throws DimensionMismatch delta_and_S(cells, m_ref[1:5]; C = 7)
    end

    @testset "AV3 on the whole log density" begin
        # Tight tolerances, so central differences resolve 1e-6 relative.
        n = length(M0_PANEL) * 3
        rng = Xoshiro(11)
        A = randn(rng, n, 4)
        tight = with_reference(mixed_setup(ds; abstol = 1e-14, reltol = 1e-12),
                               0.02 .* randn(rng, n), 1e-3 .* (A * A'))
        ℓ(z) = mixed_logdensity(tight, z)
        @test isfinite(ℓ(x))
        g = ForwardDiff.gradient(ℓ, x)
        h = 1e-4
        fd = map(eachindex(x)) do j
            e = zeros(length(x))
            e[j] = h
            (ℓ(x + e) - ℓ(x - e)) / 2h
        end
        @test all(isfinite, g)
        @test all(!iszero, g)
        for j in eachindex(x)
            @test g[j] ≈ fd[j] rtol = 1e-6
        end
    end

    @testset "−∞ outside the support, not an error" begin
        j = findfirst(==(:krnadeg), s.names)
        lo = copy(x)
        lo[j] = s.prior_mu[j] - 6 * s.prior_sd[j]
        @test mixed_logdensity(s, lo) == -Inf
        @test mixed_loglik(s, lo) == -Inf
        @test ForwardDiff.gradient(z -> mixed_logdensity(s, z), lo) isa Vector
        # The support box: outside it, −∞ without a solve.
        lb, ub = mixed_box(s)
        @test lb ≈ mixed_prior_median(s) .- 6 .* vcat(s.prior_sd, s.sigma_prior.σ)
        far = copy(x)
        far[4] = ub[4] + 1.0
        t_far = @elapsed v_far = mixed_logdensity(s, far)
        @test t_far < 1.0
        @test v_far == -Inf
        # A solve that reaches its step cap is −∞ too.
        capped = mixed_setup(ds; maxiters = 10)
        @test mixed_loglik(capped, x) == -Inf
    end

    @testset "the Turing model is the same density" begin
        model = mixed_model(s)
        for z in (x, mixed_prior_median(s))
            @test logjoint(model, (x = z,)) ≈ mixed_logdensity(s, z) rtol = 1e-12
        end
        v = mixed_theta_values(s, x)
        @test Dict(v)[:krnadeg] ≈ truth[:krnadeg]
        @test Dict(v)[:kcatR_R_ENO] ≈ truth[:kcatR_R_ENO] rtol = 1e-12
    end

    # B7.3 (§12 2026-10-24, "B7 as planned"): one stream at a time.
    @testset "streams: each factor alone" begin
        @test s.streams == [:tx, :bulk]
        @test_throws ArgumentError mixed_setup(ds; streams = ())
        @test_throws ArgumentError mixed_setup(ds; streams = (:tx, :tx))
        @test_throws ArgumentError mixed_setup(ds; streams = (:cells,))
        rng = Xoshiro(12)
        n = length(M0_PANEL) * 3
        A = randn(rng, n, 4)
        δ, Σ = 0.02 .* randn(rng, n), 1e-3 .* (A * A')
        both = with_reference(mixed_setup(ds; streams = (:bulk, :tx)), δ, Σ)
        txo = with_reference(mixed_setup(ds; streams = (:tx,)), δ, Σ)
        bulko = with_reference(mixed_setup(ds; streams = (:bulk,)), δ, Σ)
        @test txo.streams == [:tx] && bulko.streams == [:bulk]
        r = _mixed_solve(s, θ)
        k = route_transcription_rates(s, r.sol, r.p)
        μ = θ[findfirst(==(:krnadeg), s.names)] ./ s.lengths
        ltx = transcript_loglik(s.factors, k, μ)
        lbulk = bulk_loglik(both, mean_field_log_bulk(s, r.sol), ds.sigma)
        @test mixed_logdensity(both, x) === mixed_logdensity(with_reference(s, δ, Σ), x)
        @test mixed_loglik(both, x) === ltx + lbulk
        @test mixed_loglik(txo, x) === ltx
        @test mixed_loglik(bulko, x) === lbulk
        @test mixed_logdensity(txo, x) == mixed_log_prior(s, x) + ltx
        @test mixed_logdensity(bulko, x) == mixed_log_prior(s, x) + lbulk
        # AV3 on each stream alone, at tight tolerances as above.
        for st in ((:tx,), (:bulk,))
            tight = with_reference(mixed_setup(ds; abstol = 1e-14, reltol = 1e-12, streams = st),
                                   δ, Σ)
            ℓ(z) = mixed_logdensity(tight, z)
            g = ForwardDiff.gradient(ℓ, x)
            h = 1e-4
            fd = map(eachindex(x)) do j
                e = zeros(length(x))
                e[j] = h
                (ℓ(x + e) - ℓ(x - e)) / 2h
            end
            @test all(isfinite, g)
            for j in eachindex(x)
                @test g[j] ≈ fd[j] rtol = 1e-6
            end
        end
    end

    if get(ENV, "INFERCELL_INTEGRATION_TESTS", "false") == "true"
        @testset "integration: the reference θ and a few NUTS steps" begin
            ref = reference_theta(s; maxiters = 20)
            @test ref.method in (:map, :prior_median)
            @test isfinite(ref.logdensity)
            @test length(ref.x) == length(x)
            chain = sample(Xoshiro(1), mixed_model(s), NUTS(5, 0.65; max_depth = 5), 5; progress = false)
            @test size(chain, 1) == 5
            @test all(isfinite, chain[:lp])
        end
    end
end
