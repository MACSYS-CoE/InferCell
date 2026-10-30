# KIM-68: `predict` against the frozen reference. The M4 freeze (spec §12
# 2026-10-30) requires any change to `predict` after it to reproduce the frozen
# code's results on the same seeds. test/fixtures/predict_reference.jl holds them:
# the calls in predict_reference_calls.jl, run at the frozen commit by
# dev/scripts/make_predict_reference.jl. This reruns the calls at this tree, serial
# and threaded. Integers, symbols and seeds must match exactly. Floats must match to
# rtol 1e-6 and atol 1e-9. Roundoff differs between machines, and the ODE solves at
# reltol 1e-8, so the step sizes it picks can move a result by about that much: on
# CI a minimum of 0.0013 moved by 1.2e-8 relative (PR #129's first run). The bitwise
# check, on one node, is dev/scripts/make_predict_reference.slurm, run on the frozen
# tree and on this one.

using Test

include(joinpath(@__DIR__, "fixtures", "predict_reference_calls.jl"))
include(joinpath(@__DIR__, "fixtures", "predict_reference.jl"))

@testset "KIM-68: predict against the frozen reference ($PREDICT_REFERENCE_COMMIT)" begin
    @test isequal(predict_reference_posterior(), PREDICT_REFERENCE_POSTERIOR)
    @test Set(first.(predict_reference_calls())) == Set(keys(PREDICT_REFERENCE))
    for (name, run) in predict_reference_calls(), threads in (false, true)
        d = _pref_diff(run(threads), PREDICT_REFERENCE[name])
        d === nothing || @info "KIM-68 reference: call $name, threads = $threads, first difference" d
        @test d === nothing
    end
end
