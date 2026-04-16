using InferCell
using Test
using StaticArrays

@testset "LightMetabolism" begin
    @testset "Construction" begin
        m = LightMetabolism()
        @test states(m) == [:ATP, :NTP, :AA]
        @test length(parameters(m)) == 9  # 5 free rate + 3 fixed IC + 1 obs
        @test length(model_free_params(parameters(m))) == 5
        @test length(obs_free_params(parameters(m))) == 1
        @test formalism(m) == :ode
        @test inference_mode(m) == :differentiable
        @test inputs(m) == [:mRNA]
    end

    @testset "Custom parameters" begin
        m = LightMetabolism(k_atp=1.0, ATP_max=20.0)
        mfp = model_free_params(parameters(m))
        @test mfp[1].value == 1.0  # k_atp
        @test m.ATP_max == 20.0
    end

    @testset "Dynamics direct call" begin
        m = LightMetabolism()
        u_local = SA[5.0, 2.0, 2.0]    # ATP, NTP, AA
        p_local = [1.0, 0.2, 0.4, 1.0, 2.0]  # k_atp, k_ntp, k_aa, k_tx, k_tl
        u_inputs = [2.0]  # mRNA = 2.0

        du = dynamics(u_local, p_local, 0.0, m, u_inputs)
        @test length(du) == 3

        # ATP production = 1.0*(10-5) = 5.0
        # ATP costs: e_tx*k_tx = 0.1*1.0 = 0.1, e_tl*k_tl*mRNA = 0.05*2.0*2.0 = 0.2
        # NTP synth = 0.2*5 = 1.0, AA synth = 0.4*5 = 2.0
        # dATP = 5.0 - 0.1 - 0.2 - 1.0 - 2.0 = 1.7
        @test du[1] ≈ 1.7

        # dNTP = 1.0 - 0.5*1.0 - 0.1*2.0 = 0.3
        @test du[2] ≈ 0.3

        # dAA = 2.0 - 0.5*2.0*2.0 - 0.1*2.0 = 2.0 - 2.0 - 0.2 = -0.2
        @test du[3] ≈ -0.2
    end
end
