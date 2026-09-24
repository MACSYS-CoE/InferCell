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
        # `persistent_tasks` is off, and not because it fails — it *throws*
        # before it can run. Aqua walks the manifest to locate each dependency's
        # directory and errors with `Unable to locate MathOptInterface, a
        # dependency of Optim`, because MathOptInterface is a weak dependency
        # behind an Optim extension and is not installed unless its trigger
        # loads. That is a limitation of Aqua's manifest walker (seen on 0.8.17),
        # not a property of InferCell: the check never gets far enough to look
        # for a task. It went unnoticed while the exhausted Actions quota kept CI
        # from running at all, and surfaced on 2026-09-18 when the repo went
        # public. Worth retrying on a newer Aqua — a green suite should mean the
        # suite passed, not that one check errors every time.
        Aqua.test_all(
            InferCell;
            ambiguities = false,
            piracies = false,
            persistent_tasks = false,
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
    # Executable jump doubles, then the jump-composition tests (spec phase 2).
    include("jump_test_models.jl")
    include("test_jump_composition.jl")
    # Executable hybrid doubles, then the handshake tests (spec phase 3).
    include("hybrid_test_models.jl")
    include("test_hybrid_handshake.jl")
    # The 60 s rebuild (spec phase 4), on the same doubles file.
    include("test_rebuild.jl")
    # Growth and volume (spec phase 5), on the same doubles file.
    include("test_growth.jl")
    # The volume channel's inbound half (spec phase 5b), same doubles file.
    include("test_volume_inbound.jl")
    # Central glycolysis (spec phase 6), the first Core A′ module. Its
    # energy double is local to the file, not in the shared doubles.
    include("test_corea_central_glycolysis.jl")
    # Phosphotransferase transport and lactate export (spec phase 7). Needs
    # corea_test_models.jl for `caught` and contribution_test_models.jl for
    # `_gidx`, both included above.
    include("test_corea_pts_transport.jl")
    # Nucleotide recycling (spec phase 8), which owns the pools the other
    # three ODE modules route energy through. Its doubles — the charging
    # drain and the held glycolytic pools — are shared with the full-cycle
    # driver, so they live in their own file.
    include("nucleotide_test_models.jl")
    include("test_corea_nucleotide_recycling.jl")
    # tRNA charging (spec phase 9). Its doubles compose with phase 8's, so
    # they are included after them.
    include("trna_test_models.jl")
    include("test_corea_trna_charging.jl")
    # Transcription (spec phase 10): the decay double first, then the tests.
    include("corea_transcription_doubles.jl")
    include("test_corea_transcription.jl")
    # Transcript decay (spec phase 12): its transcript-source double, then the tests.
    include("corea_decay_doubles.jl")
    include("test_corea_transcript_decay.jl")
    # Catalytic edges on jump-owned protein counts (spec phase 11a). Uses the
    # hybrid doubles and `CoreAStub`, and phase 7's file is not needed.
    include("test_catalytic_slots.jl")
    # Phase 11 reuses phase 8's `HeldGlycolytic`, so it comes after it.
    include("test_corea_translation.jl")
    # The framework fixes assembly needs (spec phase 13a), on the hybrid doubles.
    include("test_chemostat_and_products.jl")
    # The assembled model (spec phase 13b) composes all seven modules.
    include("test_corea_assembly.jl")

    if get(ENV, "INFERCELL_INTEGRATION_TESTS", "false") == "true"
        include("test_txl.jl")
        include("test_stochastic_inference.jl")
        include("test_joint_inference.jl")
        include("test_sequential_inference.jl")
    end
end
