using Test
using InferCell
using JumpProcesses
using Aqua

@testset "InferCell" begin
    @testset "Aqua quality" begin
        # `ambiguities` and `piracies` are disabled initially: extending
        # SciML / Turing types produces many false positives. Re-enable
        # after an audit (tracked as a follow-up to issue #14).
        Aqua.test_all(
            InferCell;
            ambiguities = false,
            piracies = false,
        )
    end

    include("test_parameters.jl")
    include("test_orchestrator.jl")
    include("test_inference.jl")
    include("test_stochastic_ge.jl")
    include("test_bursty_gene_expression.jl")
    include("test_summary_stats.jl")
    include("test_abc_smc.jl")
    include("test_light_metabolism.jl")
    include("test_composition.jl")
    include("test_multi_gene_txl.jl")
    include("test_tier_b.jl")
    include("test_boundary.jl")

    if get(ENV, "INFERCELL_INTEGRATION_TESTS", "false") == "true"
        include("test_txl.jl")
        include("test_stochastic_inference.jl")
        include("test_joint_inference.jl")
        include("test_sequential_inference.jl")
    end
end
