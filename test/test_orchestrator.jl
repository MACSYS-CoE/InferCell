using Test
using InferCell
using OrdinaryDiffEq: Tsit5, solve
using SciMLBase: ReturnCode

@testset "Orchestrator" begin
    @testset "build_problem with TranscriptionTranslation" begin
        txl = TranscriptionTranslation()
        prob = build_problem([txl]; tspan=(0.0, 50.0))

        # Check problem structure
        @test prob.tspan == (0.0, 50.0)
        @test length(prob.u0) == 2       # mRNA, protein
        @test prob.u0[1] == 0.0          # mRNA0
        @test prob.u0[2] == 0.0          # protein0
        @test length(prob.p) == 4        # k_tx, k_tl, gamma_mRNA, gamma_protein
    end

    @testset "forward solve" begin
        txl = TranscriptionTranslation()
        prob = build_problem(txl; tspan=(0.0, 50.0))
        sol = solve(prob, Tsit5(); saveat=0.0:2.5:50.0)

        @test sol.retcode == ReturnCode.Success
        @test size(Array(sol), 1) == 2   # 2 species
        @test size(Array(sol), 2) == 21  # 21 timepoints

        # mRNA should reach steady state k_tx / gamma_mRNA = 1.0 / 0.5 = 2.0
        @test Array(sol)[1, end] ≈ 2.0 atol=0.1

        # protein steady state: k_tl * mRNA_ss / gamma_protein = 2.0 * 2.0 / 0.1 = 40.0
        @test Array(sol)[2, end] ≈ 40.0 atol=1.0
    end

    @testset "single model convenience" begin
        txl = TranscriptionTranslation()
        prob1 = build_problem([txl]; tspan=(0.0, 10.0))
        prob2 = build_problem(txl; tspan=(0.0, 10.0))
        @test prob1.u0 == prob2.u0
        @test prob1.p == prob2.p
        @test prob1.tspan == prob2.tspan
    end

    @testset "formalism validation" begin
        @test_throws ArgumentError InferCell.validate_formalism(:ODE)
        @test_throws ArgumentError InferCell.validate_formalism(:invalid)
        @test InferCell.validate_formalism(:ode) === true
    end
end
