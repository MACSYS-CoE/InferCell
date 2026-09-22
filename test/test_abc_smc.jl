using Random
using Statistics

@testset "ABC-SMC" begin
    @testset "Toy problem: infer Gaussian mean" begin
        # True parameter: mu = 3.0
        # Observed data summary: mean of 100 samples from Normal(3.0, 1.0)
        true_mu = 3.0
        observed_stats = [true_mu]  # The ideal summary stat

        priors = [Normal(0.0, 10.0)]
        param_names = [:mu]

        function simulate_gaussian(theta)
            samples = theta[1] .+ randn(100)
            return [mean(samples)]
        end

        result = abc_smc(simulate_gaussian, observed_stats, priors, param_names;
                         n_particles=200, n_populations=5, alpha=0.5)

        @test result isa ABCPosterior
        @test size(result.particles, 1) == 1
        @test size(result.particles, 2) == 200
        @test length(result.weights) == 200
        @test result.n_populations == 5

        # Posterior mean should be close to true value
        posterior_mean = sum(result.weights .* result.particles[1, :])
        @test abs(posterior_mean - true_mu) < 1.5
    end

    @testset "Returned particles lie within the reported tolerance" begin
        # A deterministic simulator, so each returned particle's distance can
        # be recomputed exactly. Before the fix the returned tolerance was the
        # next population's threshold, and about half the particles exceeded it.
        simulate(theta) = [theta[1]]
        observed = [3.0]
        for n_pop in (2, 5)
            post = abc_smc(simulate, observed, [Normal(0.0, 10.0)], [:mu];
                           n_particles=200, n_populations=n_pop, alpha=0.5,
                           rng=Random.Xoshiro(20260923))
            d = [summary_distance(simulate(post.particles[:, i]), observed)
                 for i in axes(post.particles, 2)]
            @test all(d .<= post.tolerance)
            @test maximum(d) > post.tolerance / 2   # the label is tight, not inflated
        end
        # Population 1 is the unfiltered prior: refused, not returned.
        @test_throws ArgumentError abc_smc(simulate, observed, [Normal(0.0, 10.0)],
                                           [:mu]; n_particles=50, n_populations=1)
    end

    @testset "ABCPosterior fields" begin
        particles = randn(2, 50)
        weights = fill(1.0/50, 50)
        result = ABCPosterior(particles, weights, [:a, :b], 0.1, 5)
        @test result.param_names == [:a, :b]
        @test result.tolerance == 0.1
        @test result.n_populations == 5
    end
end
