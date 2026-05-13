using OrdinaryDiffEq: Tsit5, solve
using SciMLBase: ReturnCode

@testset "TierBMetabolism (generator only)" begin
    @testset "Construction" begin
        m = TierBMetabolism()
        @test states(m) == [:ATP, :NTP, :AA]
        @test inputs(m) == [:mRNA, :protein]
        @test formalism(m) == :ode
        @test m.K_M_enzyme == 5.0
        @test m.k_atp_leak == 0.2
        @test m.n_hill == 2.0
    end

    @testset "Differs from LightMetabolism at intermediate enzyme levels" begin
        # Both models reduce to similar forms at enzyme=0 (modulo the ATP leak).
        # At intermediate enzyme, the Hill (n=2) sharpening + ATP leak produces
        # measurably different ATP dynamics — the source of mismatch Tier B exists
        # to surface.
        light = LightMetabolism()
        tierb = TierBMetabolism()
        u = SA[5.0, 2.0, 2.0]
        p = [1.0, 0.2, 0.4, 1.0, 2.0]
        inputs_vec = [2.0, 3.0]  # mRNA, intermediate enzyme

        du_light = dynamics(u, p, 0.0, light, inputs_vec)
        du_tierb = dynamics(u, p, 0.0, tierb, inputs_vec)

        # ATP dynamics differ; the leak alone gives -0.2, plus Hill-vs-MM difference
        @test du_light[1] != du_tierb[1]
        @test du_tierb[1] < du_light[1]  # leak dominates the difference

        # NTP / AA dynamics share the same form; should match
        @test du_light[2] ≈ du_tierb[2]
        @test du_light[3] ≈ du_tierb[3]
    end

    @testset "Closed-loop forward sim does not blow up (high enzyme)" begin
        m = TierBMetabolism()
        u0 = SA[5.0, 2.0, 2.0]
        rhs(u, p, t) = dynamics(u, p, t, m, [5.0, 100.0])
        prob = ODEProblem{false}(rhs, u0, (0.0, 500.0),
                                  [1.0, 0.2, 0.4, 1.0, 2.0])
        sol = solve(prob, Tsit5(); saveat=10.0)
        @test sol.retcode == ReturnCode.Success
        @test all(s[1] > 0 for s in sol.u)
        @test all(s[1] < m.ATP_max for s in sol.u)
    end

    @testset "Composable with Block 2" begin
        txl = TranscriptionTranslation()
        gen = TierBMetabolism()
        prob = build_problem([txl, gen]; tspan=(0.0, 100.0))
        @test length(prob.u0) == 5
        sol = solve(prob, Tsit5(); saveat=0:10.0:100.0)
        @test sol.retcode == ReturnCode.Success
        # The two metabolism flavours produce visibly different ATP trajectories.
        # (Direction isn't fixed: Hill boost and ATP leak push in opposite directions
        # at steady state. The point of Tier B is the mismatch, not its sign.)
        txl2 = TranscriptionTranslation()
        light = LightMetabolism()
        prob_light = build_problem([txl2, light]; tspan=(0.0, 100.0))
        sol_light = solve(prob_light, Tsit5(); saveat=0:10.0:100.0)
        diff = abs(sol.u[end][3] - sol_light.u[end][3])
        @test diff > 0.01  # measurable divergence in ATP at steady state
    end
end
