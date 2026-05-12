using Test
using InferCell
using Distributions
using Random
using MCMCChains: Chains

@testset "Boundary Conditioning" begin

    @testset "KDEPrior basics" begin
        Random.seed!(42)
        samples = randn(500) .+ 3.0  # Normal(3, 1) samples
        kp = KDEPrior(samples)

        @test kp isa KDEPrior
        @test kp.bandwidth > 0

        # rand produces values in reasonable range
        draws = [rand(Random.default_rng(), kp) for _ in 1:1000]
        @test 1.0 < mean(draws) < 5.0
        @test 0.5 < std(draws) < 3.0

        # pdf is positive near the mean and near-zero far away
        @test pdf(kp, 3.0) > 0.1
        @test pdf(kp, 100.0) == 0.0
        @test pdf(kp, -100.0) == 0.0

        # logpdf
        @test logpdf(kp, 3.0) < 0
        @test logpdf(kp, 100.0) == -Inf

        # Support bounds for Bijectors / Turing (NUTS needs these to derive a transform)
        @test isfinite(minimum(kp))
        @test isfinite(maximum(kp))
        @test minimum(kp) < mean(samples) < maximum(kp)
        @test Distributions.insupport(kp, 3.0)
        @test !Distributions.insupport(kp, 100.0)
    end

    @testset "boundary_condition with :kde method" begin
        Random.seed!(42)

        # Build a synthetic chain with 4 parameters matching StochasticGeneExpression
        n_samples = 200
        chain_data = hcat(
            rand(LogNormal(0, 0.3), n_samples),   # k_tx ~ 1.0
            rand(LogNormal(0.7, 0.3), n_samples),  # k_tl ~ 2.0
            rand(LogNormal(-0.7, 0.3), n_samples),  # gamma_mRNA ~ 0.5
            rand(LogNormal(-2.3, 0.3), n_samples),  # gamma_protein ~ 0.1
            abs.(rand(Normal(0.3, 0.05), n_samples)),  # sigma_obs
        )
        chain_data = reshape(chain_data, n_samples, 5, 1)
        chain = Chains(chain_data, ["k_tx", "k_tl", "gamma_mRNA", "gamma_protein", "sigma_obs"])

        sge = StochasticGeneExpression()
        bc = boundary_condition(chain, sge; method=:kde)

        @test length(bc.priors) == 4
        @test bc.param_names == [:k_tx, :k_tl, :gamma_mRNA, :gamma_protein]
        @test all(bc.conditioned)
        @test all(p -> p isa KDEPrior, bc.priors)
    end

    @testset "boundary_condition with :normal method" begin
        Random.seed!(42)

        n_samples = 200
        chain_data = hcat(
            rand(LogNormal(0, 0.3), n_samples),
            rand(LogNormal(0.7, 0.3), n_samples),
            rand(LogNormal(-0.7, 0.3), n_samples),
            rand(LogNormal(-2.3, 0.3), n_samples),
            abs.(rand(Normal(0.3, 0.05), n_samples)),
        )
        chain_data = reshape(chain_data, n_samples, 5, 1)
        chain = Chains(chain_data, ["k_tx", "k_tl", "gamma_mRNA", "gamma_protein", "sigma_obs"])

        sge = StochasticGeneExpression()
        bc = boundary_condition(chain, sge; method=:normal)

        @test length(bc.priors) == 4
        @test all(bc.conditioned)
        @test all(p -> p isa Normal, bc.priors)

        # Normal fit should have reasonable mean for k_tx (centered around 1)
        @test 0.5 < mean(bc.priors[1]) < 2.0
    end

    @testset "partial conditioning (unmatched params keep original prior)" begin
        Random.seed!(42)

        # Chain with only 2 of the 4 SGE parameters
        n_samples = 200
        chain_data = hcat(
            rand(LogNormal(0, 0.3), n_samples),
            rand(LogNormal(0.7, 0.3), n_samples),
        )
        chain_data = reshape(chain_data, n_samples, 2, 1)
        chain = Chains(chain_data, ["k_tx", "k_tl"])

        sge = StochasticGeneExpression()
        bc = boundary_condition(chain, sge; method=:kde)

        @test length(bc.priors) == 4
        @test bc.conditioned == [true, true, false, false]
        @test bc.priors[1] isa KDEPrior    # k_tx: conditioned
        @test bc.priors[2] isa KDEPrior    # k_tl: conditioned
        @test bc.priors[3] isa LogNormal   # gamma_mRNA: original
        @test bc.priors[4] isa LogNormal   # gamma_protein: original
    end

    @testset "KDEPrior pdf support check compatibility" begin
        # abc_smc.jl line 51: all(pdf(priors[j], candidate[j]) > 0 for j in 1:n_params)
        samples = rand(LogNormal(0, 0.5), 500)
        kp = KDEPrior(samples)

        # Positive values near the samples should have pdf > 0
        @test pdf(kp, median(samples)) > 0
        # Extreme values should have pdf == 0
        @test pdf(kp, 1e6) == 0.0
    end

    @testset "boundary_condition on ABCPosterior (SSA -> ODE backward direction)" begin
        Random.seed!(42)
        # Hand-build an ABCPosterior over k_tx, k_tl with clearly-located particles
        n_particles = 200
        particles = zeros(2, n_particles)
        particles[1, :] = rand(LogNormal(0.0, 0.2), n_particles)   # k_tx
        particles[2, :] = rand(LogNormal(0.7, 0.2), n_particles)   # k_tl
        weights = fill(1 / n_particles, n_particles)
        post = InferCell.ABCPosterior(particles, weights, [:k_tx, :k_tl], 0.1, 5)

        # ODE block that shares both names
        txl = TranscriptionTranslation()
        bc = boundary_condition(post, txl; method=:kde, n_samples=300)

        @test length(bc.priors) == 4
        @test bc.param_names == [:k_tx, :k_tl, :gamma_mRNA, :gamma_protein]
        @test bc.conditioned == [true, true, false, false]
        @test bc.priors[1] isa KDEPrior
        @test bc.priors[2] isa KDEPrior
        @test bc.priors[3] isa LogNormal  # unconditioned, original prior preserved
        @test bc.priors[4] isa LogNormal

        # Sanity: the KDE prior for k_tx is centred near the particle mean
        kp_draws = [rand(Random.default_rng(), bc.priors[1]) for _ in 1:1000]
        @test 0.6 < mean(kp_draws) < 1.6  # particles centred near exp(0)=1
    end

    @testset "kl_divergence basics" begin
        Random.seed!(42)
        s1 = randn(500)
        s2 = randn(500)
        @test kl_divergence(s1, s1) ≈ 0.0 atol=0.05  # KL(p||p) = 0
        @test kl_divergence(s1, s2) >= 0.0
        # Shifted distributions: KL grows
        s_shift = randn(500) .+ 3.0
        @test kl_divergence(s1, s_shift) > kl_divergence(s1, s1)
    end

    @testset "chain_kl over named parameters" begin
        Random.seed!(42)
        n = 200
        # Two near-identical "chains" → low KL on shared params
        data_a = reshape(hcat(randn(n), randn(n) .+ 1.0), n, 2, 1)
        data_b = reshape(hcat(randn(n), randn(n) .+ 1.0), n, 2, 1)
        ch_a = Chains(data_a, ["k_tx", "k_tl"])
        ch_b = Chains(data_b, ["k_tx", "k_tl"])
        kl_close = chain_kl(ch_a, ch_b, [:k_tx, :k_tl])
        @test kl_close < 1.0  # similar distributions
        # Now shift one chain by a lot
        data_c = reshape(hcat(randn(n) .+ 5.0, randn(n) .+ 6.0), n, 2, 1)
        ch_c = Chains(data_c, ["k_tx", "k_tl"])
        kl_far = chain_kl(ch_a, ch_c, [:k_tx, :k_tl])
        @test kl_far > kl_close
    end
end
