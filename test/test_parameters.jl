using Test
using InferCell
using Distributions

@testset "InferParameter" begin
    @testset "construction" begin
        p = InferParameter(1.0, LogNormal(0, 1), false, :k_tx, :txl, :rate)
        @test p.value == 1.0
        @test p.fixed == false
        @test p.name == :k_tx
        @test p.module_id == :txl
        @test p.role == :rate
    end

    @testset "filtering" begin
        params = [
            InferParameter(1.0, LogNormal(0, 1), false, :k_tx,    :txl, :rate),
            InferParameter(2.0, LogNormal(0, 1), false, :k_tl,    :txl, :rate),
            InferParameter(0.0, Normal(0, 1),    true,  :mRNA0,   :txl, :initial_condition),
            InferParameter(0.0, Normal(0, 1),    true,  :protein0,:txl, :initial_condition),
            InferParameter(0.3, truncated(Normal(0,1); lower=0), false, :sigma, :txl, :observation),
        ]

        @test length(free_params(params)) == 3
        @test length(rate_params(params)) == 2
        @test length(ic_params(params)) == 2
        @test length(obs_params(params)) == 1
        @test length(model_free_params(params)) == 2
        @test length(obs_free_params(params)) == 1

        # All fixed
        fixed_params = [InferParameter(1.0, Normal(0,1), true, :x, :m, :rate)]
        @test length(free_params(fixed_params)) == 0
        @test length(model_free_params(fixed_params)) == 0

        # Empty
        @test length(free_params(InferParameter[])) == 0
    end
end
