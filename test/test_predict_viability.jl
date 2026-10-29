# KIM-30: viability-conditioned completion, on a toy rule V and made-up cells.
# Failed, starved, below-floor, clipped and clean stay distinct, and a rejected
# parameter set keeps its reason.

using Test
using InferCell

# A made-up run: `clips` as (counter, drains clipped) of `drains`.
_toy_cell(; failure = :none, drains = 100, clips = [], gtp = 1e4) =
    CellRun(1, failure === :none ? 1.0 : NaN, failure, "", drains,
            NamedTuple[(counter = c, drains = n, first = 1.0, max_deficit = 1.0)
                       for (c, n) in clips],
            Dict(:M_gtp_c => gtp, :M_atp_c => 1e5), nothing)

@testset "KIM-30: viability-conditioned completion" begin
    V = ViabilityRule(; label = "toy", allow_clipped = false, min_viable_fraction = 0.5,
                      floors = [:M_gtp_c => 500.0])
    clean = _toy_cell()
    clipped = _toy_cell(; clips = [(:GTP_translat, 3)])
    starved = _toy_cell(; clips = [(:GTP_translat, 3), (:tRNA_translat, 51)])
    low = _toy_cell(; clips = [(:GTP_translat, 3)], gtp = 120.0)
    failed = _toy_cell(; failure = :not_advanced)
    overdraw = _toy_cell(; failure = :overdraw)

    @testset "one category per cell, in order of precedence" begin
        @test cell_status(clean, V) === :clean
        @test cell_status(clipped, V) === :clipped
        @test cell_status(starved, V) === :starved
        @test cell_status(low, V) === :below_floor      # below the floor outranks clipping
        @test cell_status(failed, V) === :failed
        @test cell_status(overdraw, V) === :starved     # the mean field's starvation
        # 50 of 100 drains is not more than half; 14c's rule is strictly more.
        @test cell_status(_toy_cell(; clips = [(:x, 50)]), V) === :clipped
        # Without a rule: no floors, so the low cell reads as clipped.
        @test cell_status(low) === :clipped && cell_status(starved) === :starved
        @test cell_status(starved; starved_fraction = 0.9) === :clipped
        bad = ViabilityRule(; label = "bad", allow_clipped = true, min_viable_fraction = 1,
                            floors = [:M_nope_c => 1.0])
        @test_throws ArgumentError cell_status(clean, bad)
    end

    @testset "V is a property of θ: the fraction of its cells that are viable" begin
        a = assess_viability([clean, clean, clipped, starved], V)
        @test a.viable && a.fraction == 0.5
        @test a.counts[:clean] == 2 && a.counts[:clipped] == 1 && a.counts[:starved] == 1
        r = assess_viability([clean, low, failed, starved], V)
        @test !r.viable && r.fraction == 0.25
        # The rejected set keeps its reason: what each cell was.
        @test r.counts == Dict(:clean => 1, :below_floor => 1, :failed => 1, :starved => 1,
                               :clipped => 0)
        lenient = ViabilityRule(; label = "lenient", allow_clipped = true,
                                min_viable_fraction = 0.75)
        @test assess_viability([clean, clipped, clipped, starved], lenient).viable
        @test !assess_viability([clean, clipped, clipped, starved], V).viable
        @test_throws ArgumentError assess_viability(CellRun[], V)
        @test_throws ArgumentError ViabilityRule(; label = "x", allow_clipped = true,
                                                 min_viable_fraction = 0)
    end

    @testset "both variants as weights on one crossed matrix" begin
        v = Bool[1 1 0 0;     # row 1: 2 of 4 viable
                 0 0 0 0;     # row 2: none; (b) fails here
                 1 1 1 1;     # row 3: all
                 0 0 0 1]     # row 4: 1 of 4
        w = viability_weights(v)
        @test w.acceptance == 7 / 16
        @test w.row_acceptance == [0.5, 0.0, 1.0, 0.25]
        @test w.failed_rows == [2]
        # (b): each surviving row weighs 1/3, shared over its viable cells.
        @test vec(sum(w.per_theta_I; dims = 2)) ≈ [1 / 3, 0, 1 / 3, 1 / 3]
        @test w.per_theta_I[4, 4] ≈ 1 / 3 && w.per_theta_I[3, 1] ≈ 1 / 12
        # (a): every viable cell weighs the same, so a row weighs its acceptance.
        @test sum(w.joint) ≈ 1 && all(w.joint[v] .≈ 1 / 7) && all(w.joint[.!v] .== 0)
        @test w.ess_joint ≈ 7
        @test w.ess_per_theta_I ≈ 1 / sum(abs2, w.per_theta_I)
        @test w.ess_per_theta_I < w.ess_joint
        none = viability_weights(falses(2, 3))
        @test none.acceptance == 0 && none.ess_joint == 0 && none.ess_per_theta_I == 0
        @test none.failed_rows == [1, 2]
    end
end
