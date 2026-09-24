using Test
using InferCell

# Spec §11 phase 13b — the assembled Core A′: all seven modules, every counter
# wired, built by `corea_models()` / `build_corea()`.

_without(id) = AbstractSubModel[m for m in corea_models() if module_id(m) !== id]

@testset "Phase 13b: the assembled Core A′" begin

    @testset "13.1 completeness mode passes on the assembly" begin
        graph = resolve_coupling(corea_models(); complete = true)
        @test isempty(graph.unowned_states)
        @test isempty(graph.dead_ends)
        # The four PTS carriers are translated and never degraded: their only
        # sink is growth dilution, so they accumulate rather than dead-end.
        @test sort([d.species for d in graph.accumulating]) ==
              [:M_crr_c, :M_ptsg_c, :M_ptsh_c, :M_ptsi_c]
        @test all(d -> d.missing_role === :consumer, graph.accumulating)
        @test all(d -> d.modules == [:CoreATranslation], graph.accumulating)
    end

    @testset "13.1 completeness mode fails, naming what is missing" begin
        # The default stays a report: a partial composition resolves.
        partial = _without(:NucleotideRecycling)
        @test !isempty(resolve_coupling(partial).unowned_states)

        err = caught(() -> resolve_coupling(partial; complete = true))
        @test err isa IncompleteComposition
        for s in (:M_atp_c, :M_gtp_c, :M_gdp_c, :M_amp_c)
            @test s in err.unowned
        end
        msg = sprint(showerror, err)
        @test occursin(":M_gtp_c", msg)
        @test occursin("adenylate moiety", msg)
        @test occursin("guanylate moiety", msg)

        err = caught(() -> resolve_coupling(_without(:TrnaCharging); complete = true))
        @test err isa IncompleteComposition
        @test Set(err.unowned) == Set([:M_trna_c, :M_trna_chg_c])
        stranded = Dict(d.species => d.missing_role for d in err.dead_ends)
        @test stranded[:M_trna_chg_c] === :producer   # translation draws it
        @test stranded[:M_trna_c] === :consumer       # translation returns it
        @test occursin("trna moiety", sprint(showerror, err))

        # The same through the hybrid build, which is where the assembly meets it.
        @test caught(() -> build_problem(_without(:TrnaCharging);
                                         tspan = (0.0, 60.0), complete = true)) isa
              IncompleteComposition
    end

    @testset "13.1 a protein consumed with no producer is still a dead end" begin
        # The accumulating class is the produced-with-no-consumer half only.
        drawer = CoreAStub(:Drawer; ins = [:M_ptsg_c], contribs = [:M_ptsg_c],
                           edges = [MassEdge(species = :M_ptsg_c, direction = :in)])
        graph = resolve_coupling([drawer])
        @test only(graph.dead_ends).species === :M_ptsg_c
        @test only(graph.dead_ends).missing_role === :producer
        @test isempty(graph.accumulating)
    end
end
