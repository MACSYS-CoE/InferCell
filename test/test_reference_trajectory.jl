using Test
using InferCell
using OrdinaryDiffEq: Tsit5, solve

# Bit-for-bit reference for the composed ODE path (spec §11 task 1.2).
#
# The fixture was generated on the tree before phase 1 touched the orchestrator
# (commit recorded in the file). Any change to `_build_rhs`, `_build_u0` or
# `_build_p0` that alters even the sign of a zero in the trajectory of the
# simplest model fails here. Tolerances are read from the fixture so the test
# and the generator cannot drift apart.
include(joinpath(@__DIR__, "fixtures", "txl_reference.jl"))

_bits(v) = reinterpret(UInt64, collect(Float64, v))

@testset "Reference trajectory is byte-identical" begin
    prob = build_problem(TranscriptionTranslation(); tspan = (0.0, 50.0))
    sol = solve(prob, Tsit5(); saveat = TXL_REFERENCE_SAVEAT,
                abstol = TXL_REFERENCE_ABSTOL, reltol = TXL_REFERENCE_RELTOL)

    @test _bits(prob.u0) == TXL_REFERENCE_U0_BITS
    @test _bits(prob.p) == TXL_REFERENCE_P_BITS
    @test _bits(sol.t) == TXL_REFERENCE_T_BITS
    @test length(sol.u) == length(TXL_REFERENCE_U_BITS)
    # One assertion per save point so a failure names the time it first differs.
    for (i, u) in enumerate(sol.u)
        @test _bits(u) == TXL_REFERENCE_U_BITS[i]
    end
end
