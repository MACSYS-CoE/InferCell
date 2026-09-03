using Test
using InferCell
using JumpProcesses: SSAStepper

# Jump composition (spec §11 phase 2). Doubles in jump_test_models.jl.

@testset "Jump composition" begin
    # Task 2.1: pin the bug before fixing it. Module B's only jump is written in
    # the pre-phase-2 style — `integrator.u[1] += 1`, `p[1]` — and B sits second
    # in the composition, so today its firings land on A's state and its rate
    # reads A's parameter. Over 50 s at rate 1 the chance of no firing is e^-50.
    @testset "2.1 witness: a second jump module's affect writes the first module's state" begin
        a = LegacyOwner(x0 = 3)
        b = LegacyWriter(b = 1.0)
        prob = build_problem([a, b]; tspan = (0.0, 50.0))
        sol = solve(prob, SSAStepper(); saveat = [50.0])
        @test sol[1, end] > 3      # A's :legacy_x moved, though A fires nothing
        @test sol[2, end] == 0     # B's own :legacy_y never moved
    end

    @testset "2.1 witness: a second jump module's rate reads the first module's parameter" begin
        a = LegacyOwner(x0 = 3, k = 0.0)   # p[1] is A's k_dead = 0
        b = LegacyWriter(b = 1.0)          # B's rate reads p[1], intending b_legacy = 1
        prob = build_problem([a, b]; tspan = (0.0, 50.0))
        @test prob.prob.p == [0.0, 1.0]
        sol = solve(prob, SSAStepper(); saveat = [50.0])
        @test sol[1, end] == 3     # nothing fired: B's rate was A's zero
        @test sol[2, end] == 0
    end
end
