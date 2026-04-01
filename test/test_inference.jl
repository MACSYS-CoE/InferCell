using Test
using InferCell
using OrdinaryDiffEq: Tsit5, solve
using Random

@testset "Inference" begin
    @testset "observe" begin
        Random.seed!(123)
        txl = TranscriptionTranslation()
        prob = build_problem(txl; tspan=(0.0, 50.0))
        sol = solve(prob, Tsit5())
        times = 0.0:5.0:50.0

        data = observe(sol, times, txl; sigma=0.1)

        @test data isa ObservedData
        @test length(data.times) == 11
        @test size(data.observations) == (2, 11)
        @test data.species == [:mRNA, :protein]

        # Observations should differ from true values (noise added)
        pred = Array(sol(times))
        @test !isapprox(data.observations, pred; atol=1e-10)
    end

    @testset "build_turing_model" begin
        Random.seed!(42)
        txl = TranscriptionTranslation()
        prob = build_problem(txl; tspan=(0.0, 50.0))
        sol = solve(prob, Tsit5())
        data = observe(sol, 0.0:5.0:50.0, txl; sigma=0.3)

        turing_model, param_names = build_turing_model([txl], data, prob)

        @test length(param_names) == 5  # 4 rate + 1 obs noise
        @test :k_tx in param_names
        @test :k_tl in param_names
        @test :gamma_mRNA in param_names
        @test :gamma_protein in param_names
        @test :sigma_obs in param_names
    end

    @testset "check_identifiability" begin
        txl = TranscriptionTranslation()
        prob = build_problem(txl; tspan=(0.0, 50.0))
        times = collect(0.0:2.5:50.0)

        result = check_identifiability(txl, prob, times)

        @test result.n_params == 4
        @test result.n_obs == 2 * length(times)  # 2 species x 21 times
        @test result.full_rank == true
        @test result.rank >= 4
    end
end
