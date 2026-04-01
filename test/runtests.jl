using Test
using InferCell

@testset "InferCell" begin
    include("test_parameters.jl")
    include("test_orchestrator.jl")
    include("test_inference.jl")

    if get(ENV, "INFERCELL_INTEGRATION_TESTS", "false") == "true"
        include("test_txl.jl")
    end
end
