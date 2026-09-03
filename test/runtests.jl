using Test
using InferCell
using JumpProcesses
using Aqua

@testset "InferCell" begin
    @testset "Aqua quality" begin
        # `ambiguities` and `piracies` are disabled initially: extending
        # SciML / Turing types produces many false positives. Re-enable
        # after an audit (tracked as a follow-up to issue #14).
        #
        # `persistent_tasks` builds a fresh package depending on InferCell and
        # resolves it from Project.toml alone, which needs the registry. Cluster
        # compute nodes have no network, so it cannot run there — set
        # INFERCELL_OFFLINE_TESTS=true to skip it. It stays on by default, so CI
        # still runs it.
        offline = get(ENV, "INFERCELL_OFFLINE_TESTS", "false") == "true"
        Aqua.test_all(
            InferCell;
            ambiguities = false,
            piracies = false,
            persistent_tasks = !offline,
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
    include("test_reference_trajectory.jl")

    # The composition framework (edges, resolver, loader, labels) and the one
    # organism that uses it. `corea_test_models.jl` defines the sub-model
    # doubles the resolver tests compose, so it must come first.
    include("corea_test_models.jl")
    include("test_corea_registry.jl")
    include("test_edges.jl")
    include("test_resolver.jl")
    include("test_loader.jl")
    include("test_labels.jl")
    # Executable doubles for the contribution channel, then its tests.
    include("contribution_test_models.jl")
    include("test_contributions.jl")

    if get(ENV, "INFERCELL_INTEGRATION_TESTS", "false") == "true"
        include("test_txl.jl")
        include("test_stochastic_inference.jl")
        include("test_joint_inference.jl")
        include("test_sequential_inference.jl")
    end
end
