# Spec §11 task 15.4 and §3 check 9: per-modality observation noise on a log
# scale, and recovery of a known noise scale.

using Test
using InferCell
using OrdinaryDiffEq: Tsit5, solve
using Random
using Distributions: Normal, MvNormal, LogNormal, logpdf, quantile
using Turing: logjoint

@testset "Per-modality noise (spec §11 task 15.4)" begin
    @testset "the smallest pools now contribute" begin
        # One large pool and one 1,000 times smaller, as Core A′'s span from
        # 0.0098 to 17.8 mM. Double the small pool's prediction and ask how much
        # the log-likelihood notices.
        data = ObservedData([1.0], reshape([10.0, 0.01], 2, 1), [:big, :small])
        pred = [10.0, 0.01]
        moved = [10.0, 0.02]
        # One additive σ sized for the large pool: the small pool is invisible.
        σ_add = 0.1 * 10.0
        Δ_add = logpdf(MvNormal(pred, σ_add), data.observations[:, 1]) -
                logpdf(MvNormal(moved, σ_add), data.observations[:, 1])
        @test abs(Δ_add) < 1e-3
        # On a log scale it moves the log-likelihood by (ln 2)²/2σ².
        nm = NoiseModel(Modality(:metabolite, :lognormal, [:big, :small]))
        Δ_log = observation_loglik(nm, data, i -> pred, [0.2]) -
                observation_loglik(nm, data, i -> moved, [0.2])
        @test Δ_log ≈ log(2)^2 / (2 * 0.2^2)
        @test Δ_log > 5
    end

    @testset "counts, partitions and refusals" begin
        data = ObservedData([1.0, 2.0], [3.0 5.0; 1.2 0.8], [:mRNA, :atp])
        nm = NoiseModel(Modality(:transcript, :poisson, [:mRNA]),
                        Modality(:metabolite, :lognormal, [:atp]))
        @test noise_scales(nm) == [(:sigma_metabolite, LogNormal(log(0.2), 1.0))]
        ll = observation_loglik(nm, data, i -> [4.0, 1.0], [0.1])
        @test ll ≈ sum(logpdf.(Ref(InferCell.Poisson(4.0)), [3, 5])) +
                   sum(logpdf.(Ref(Normal(0.0, 0.1)), log.([1.2, 0.8])))
        @test_throws ArgumentError observation_loglik(
            NoiseModel(Modality(:m, :lognormal, [:atp])), data, i -> [4.0, 1.0], [0.1])
        @test_throws ArgumentError NoiseModel(Modality(:a, :lognormal, [:atp]),
                                              Modality(:b, :poisson, [:atp]))
        @test_throws ArgumentError Modality(:a, :additive, [:atp])
        # A count modality needs counts, and a log modality positive values.
        bad = ObservedData([1.0], reshape([2.5, -1.0], 2, 1), [:mRNA, :atp])
        @test_throws ArgumentError observation_loglik(nm, bad, i -> [4.0, 1.0], [0.1])
    end

    txl = TranscriptionTranslation()
    prob = build_problem(txl; tspan = (0.0, 50.0))
    sol = solve(prob, Tsit5(); abstol = 1e-10, reltol = 1e-10)
    times = 5.0:5.0:50.0
    nm = NoiseModel(Modality(:all, :lognormal, [:mRNA, :protein]))
    σ_true = 0.15

    @testset "the Turing model scores the noise model" begin
        data = observe(sol, times, [txl], nm; scales = Dict(:sigma_all => σ_true),
                       rng = Xoshiro(4))
        tm, names = build_turing_model([txl], data, prob; noise = nm)
        @test names == [:k_tx, :k_tl, :gamma_mRNA, :gamma_protein, :sigma_all]
        θ = [p.value for p in parameters(txl) if !p.fixed && p.role === :rate]
        priors = [p.prior for p in parameters(txl) if !p.fixed && p.role === :rate]
        lj = logjoint(tm, (theta = vcat(θ, σ_true),))
        pred(i) = solve(prob, Tsit5(); saveat = times).u[i]
        expect = sum(logpdf.(priors, θ)) + logpdf(LogNormal(log(0.2), 1.0), σ_true) +
                 observation_loglik(nm, data, pred, [σ_true])
        @test lj ≈ expect rtol = 1e-6
    end

    # Check 9: generate at a known σ, infer it, and require the 90% interval to
    # cover the truth at its nominal rate over repeated datasets. The rates are
    # held at truth, so the posterior in σ is one-dimensional and is computed on
    # a grid of the same likelihood the sampler uses.
    @testset "check 9: a known noise scale is recovered (spec §3)" begin
        grid = range(0.02, 0.6; length = 2000)
        prior = LogNormal(log(0.2), 1.0)
        pred(i) = sol(times[i])
        function interval(data; misread = false)
            lp = [logpdf(prior, s) +
                  observation_loglik(nm, data, pred, [misread ? s^2 : s]) for s in grid]
            w = exp.(lp .- maximum(lp))
            c = cumsum(w) ./ sum(w)
            return grid[findfirst(>=(0.05), c)], grid[findfirst(>=(0.95), c)],
                   grid[findfirst(>=(0.5), c)]
        end
        rng = Xoshiro(909)
        n = 300
        covered = 0
        covered_misread = 0
        medians = Float64[]
        for _ in 1:n
            data = observe(sol, times, [txl], nm; scales = Dict(:sigma_all => σ_true), rng)
            lo, hi, med = interval(data)
            covered += lo <= σ_true <= hi
            push!(medians, med)
            lo, hi, _ = interval(data; misread = true)
            covered_misread += lo <= σ_true <= hi
        end
        # Binomial band for 90% at n = 300: 3 SD is ±5.2 points.
        @test 0.848 <= covered / n <= 0.952
        @test abs(sum(medians) / n / σ_true - 1) < 0.05
        # The mutation: reading the scale as a variance, the misreading R13 names,
        # breaks coverage.
        @test covered_misread / n < 0.5
    end
end
