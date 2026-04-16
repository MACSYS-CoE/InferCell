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
end
