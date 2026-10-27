# Issue #89 B6: the dispatcher's identity checks, which test/test_dispatch.jl gates
# behind INFERCELL_INTEGRATION_TESTS (spec/phases/approximate-method.md §12
# 2026-10-24, "B6 as planned"). Loads the test doubles as test/runtests.jl does.
#
# Usage: INFERCELL_INTEGRATION_TESTS=true julia --project --threads=4 dev/scripts/dispatch_89_b6.jl
# Regenerate: sbatch dev/scripts/dispatch_89_b6.slurm, from a worktree at the commit.

using Test, InferCell

const T = joinpath(@__DIR__, "..", "..", "test")
for f in ("corea_test_models.jl", "contribution_test_models.jl", "jump_test_models.jl",
          "hybrid_test_models.jl")
    include(joinpath(T, f))
end
t = @elapsed include(joinpath(T, "test_dispatch.jl"))
println("test_dispatch.jl took $(round(t; digits = 1)) s")
