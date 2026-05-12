using InferCell
using Test
using StaticArrays
using OrdinaryDiffEq: Tsit5, solve, ODEProblem
using SciMLBase: ReturnCode

@testset "LightMetabolism" begin
    @testset "Construction" begin
        m = LightMetabolism()
        @test states(m) == [:ATP, :NTP, :AA]
        @test length(parameters(m)) == 9  # 5 free rate + 3 fixed IC + 1 obs
        @test length(model_free_params(parameters(m))) == 5
        @test length(obs_free_params(parameters(m))) == 1
        @test formalism(m) == :ode
        @test inference_mode(m) == :differentiable
        @test inputs(m) == [:mRNA, :protein]
        @test m.K_M_enzyme == 5.0
    end

    @testset "Custom parameters" begin
        m = LightMetabolism(k_atp=1.0, ATP_max=20.0, K_M_enzyme=2.0)
        mfp = model_free_params(parameters(m))
        @test mfp[1].value == 1.0  # k_atp
        @test m.ATP_max == 20.0
        @test m.K_M_enzyme == 2.0
    end

    @testset "Dynamics direct call (enzyme=0 preserves Step-3 baseline)" begin
        m = LightMetabolism()
        u_local = SA[5.0, 2.0, 2.0]    # ATP, NTP, AA
        p_local = [1.0, 0.2, 0.4, 1.0, 2.0]  # k_atp, k_ntp, k_aa, k_tx, k_tl
        u_inputs = [2.0, 0.0]  # mRNA = 2.0, enzyme = 0.0

        du = dynamics(u_local, p_local, 0.0, m, u_inputs)
        @test length(du) == 3

        # enzyme=0 ⇒ enzyme_factor=0 ⇒ v_atp_prod = k_atp * 1 * (ATP_max - ATP)
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

    @testset "Enzyme feedback: bounded saturation" begin
        m = LightMetabolism()  # K_M_enzyme = 5.0
        u_local = SA[5.0, 2.0, 2.0]
        p_local = [1.0, 0.2, 0.4, 1.0, 2.0]
        mRNA = 2.0

        # Baseline (no enzyme)
        du_zero = dynamics(u_local, p_local, 0.0, m, [mRNA, 0.0])

        # Mid enzyme: factor = 5/(5+5) = 0.5 ⇒ v_atp_prod multiplied by 1.5
        du_mid = dynamics(u_local, p_local, 0.0, m, [mRNA, 5.0])

        # Large enzyme: factor → 1 ⇒ v_atp_prod multiplied by ~2
        du_high = dynamics(u_local, p_local, 0.0, m, [mRNA, 1e6])

        # Only the ATP production term changes; the difference is bounded
        @test du_mid[1] > du_zero[1]
        @test du_high[1] > du_mid[1]
        # In the enzyme→∞ limit, ATP production at most doubles (5 → 10), so dATP gains exactly k_atp*(ATP_max - ATP) = 5
        @test du_high[1] - du_zero[1] ≈ 5.0 atol=1e-3
        # NTP and AA dynamics depend only on ATP, k_*, mRNA — unchanged across enzyme levels
        @test du_zero[2] ≈ du_mid[2] ≈ du_high[2]
        @test du_zero[3] ≈ du_mid[3] ≈ du_high[3]
    end

    @testset "Closed-loop forward sim does not blow up" begin
        m = LightMetabolism()
        u0 = SA[5.0, 2.0, 2.0]
        # Hold mRNA, enzyme constant at high values; verify ATP stays bounded
        rhs(u, p, t) = dynamics(u, p, t, m, [5.0, 100.0])
        prob = ODEProblem{false}(rhs, u0, (0.0, 500.0), [1.0, 0.2, 0.4, 1.0, 2.0])
        sol = solve(prob, Tsit5(); saveat=10.0)
        @test sol.retcode == ReturnCode.Success
        # ATP must stay below ATP_max (production saturates the relaxation toward ATP_max,
        # and consumption keeps the steady state below it)
        @test all(s[1] < m.ATP_max for s in sol.u)
        @test all(s[1] > 0 for s in sol.u)
    end
end
