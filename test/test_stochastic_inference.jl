using Statistics
using Random

@testset "Stochastic Inference (Integration)" begin
    # Seed explicitly. Without this the file inherits whatever RNG state the
    # preceding test file happened to leave behind, so editing an unrelated
    # test upstream silently changes the trajectories drawn here and, with
    # them, whether the recovery bounds below hold.
    Random.seed!(2718)

    @testset "ABC-SMC parameter recovery on StochasticGeneExpression" begin
        # True parameters
        true_params = (k_tx=1.0, k_tl=2.0, gamma_mRNA=0.5, gamma_protein=0.1)
        model = StochasticGeneExpression(; true_params...)

        # Generate synthetic data: many replicate trajectories
        tspan = (0.0, 100.0)
        times = collect(0.0:10.0:100.0)
        prob = build_problem(model; tspan=tspan)
        n_data_replicates = 200
        trajectories = [solve(prob, SSAStepper(); saveat=times) for _ in 1:n_data_replicates]
        data = observe(trajectories, times, model)

        # Run ABC-SMC with modest settings for test speed
        result = infer(model, data;
                       n_samples=200, n_populations=6, alpha=0.5,
                       n_replicates=50, tspan=tspan)

        @test result isa ABCPosterior
        @test size(result.particles, 1) == 4  # 4 free rate params

        # Check posterior means are in the right ballpark
        # Wide tolerance (factor of 10) is appropriate for 200 particles / 6 populations.
        # gamma_protein in particular is weakly identifiable from mean trajectory statistics.
        for (j, (name, true_val)) in enumerate(pairs(true_params))
            posterior_mean = sum(result.weights .* result.particles[j, :])
            @test 0.1 * true_val < posterior_mean < 10.0 * true_val
        end
    end
end
