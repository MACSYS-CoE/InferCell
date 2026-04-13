using Test
using InferCell
using JumpProcesses

@testset "InferCell" begin
    include("test_parameters.jl")
    include("test_orchestrator.jl")
    include("test_inference.jl")
    include("test_stochastic_ge.jl")
    include("test_summary_stats.jl")
    include("test_abc_smc.jl")

    if get(ENV, "INFERCELL_INTEGRATION_TESTS", "false") == "true"
        include("test_txl.jl")
        include("test_stochastic_inference.jl")
    end
end
