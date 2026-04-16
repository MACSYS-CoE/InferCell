@testset "Summary Statistics" begin
    @testset "Computation" begin
        model = StochasticGeneExpression()
        prob = build_problem(model; tspan=(0.0, 50.0))
        trajectories = [solve(prob, SSAStepper(); saveat=[50.0]) for _ in 1:50]

        stats = compute_summary_stats(trajectories, states(model))
        # 2 species x 3 stats (mean, var, fano) = 6
        @test length(stats) == 6
        @test all(isfinite, stats)
    end

    @testset "Computation at timepoints" begin
        model = StochasticGeneExpression()
        prob = build_problem(model; tspan=(0.0, 50.0))
        times = [10.0, 25.0, 50.0]
        trajectories = [solve(prob, SSAStepper(); saveat=times) for _ in 1:50]

        stats = compute_summary_stats(trajectories, states(model); times=times)
        # 2 species x 3 timepoints = 6
        @test length(stats) == 6
        @test all(isfinite, stats)
    end

    @testset "Distance properties" begin
        s1 = [1.0, 2.0, 3.0]
        s2 = [1.0, 2.0, 3.0]
        @test summary_distance(s1, s2) == 0.0
        @test summary_distance(s1, [2.0, 3.0, 4.0]) > 0.0
        # Symmetry
        @test summary_distance(s1, [4.0, 5.0, 6.0]) == summary_distance([4.0, 5.0, 6.0], s1)
    end
end
