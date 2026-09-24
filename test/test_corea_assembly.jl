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

    @testset "13.2 every declared edge of the assembly is executed" begin
        ms = corea_models()
        d = build_corea(tspan = (0.0, 60.0))   # builds, so the assertion held
        @test isempty(unexecuted_edges(ms, d))

        graph = resolve_coupling(ms)
        kinds = Dict(k => count(r -> r.kind === k, graph.edges)
                     for k in unique(r.kind for r in graph.edges))
        @test kinds == Dict(:mass => 17, :currency => 16, :deferred_counter => 30,
                            :catalytic => 15, :rate_constant => 5, :volume => 3,
                            :clamped => 3)
        @test length(d.debits) == 30 && length(d.catalytic) == 15
        @test length(d.rebuilds) == 2 && length(d.growth) == 2
        @test length(d.geometry) == 1

        # Task 13.5: the solver and tolerances are recorded on the model, and
        # the assembly runs. The full 6,300 s cycle is the Slurm driver's
        # (dev/scripts/corea_full_cycle.jl); this is the smoke test.
        @test solver_settings(d) == (solver = :Rodas5P, abstol = 1e-10, reltol = 1e-8)
        run_handshake!(d, 60)
        @test d.ode.t == 60.0
        @test d.n_handshakes == 60

        # Task 13.10's last clause: every product edge the four product-bearing
        # counters declare is credited by the driver.
        for c in (:ATP_trsc, :ATP_mRNAdeg, :GTP_translat, :ATP_transloc)
            declared = Set(e.species for m in ms for e in coupling(m)
                           if e isa DeferredCounterEdge && e.counter === c &&
                              e.direction === :out)
            credited = Set(b.species for b in d.debits if b.counter === c && b.sign > 0)
            @test !isempty(declared)
            @test declared == credited
        end
    end

    @testset "13.2 an inert edge fails the build, naming it" begin
        # A jump module declaring an outbound currency edge on a peer's pool
        # with no written_states entry: it resolves, and before task 13.2 it
        # built and never ran.
        inert = ToyExpression(edges = CouplingEdge[
            DeferredCounterEdge(species = :M_atp_c, direction = :in, counter = :atp_cost),
            CurrencyEdge(species = :M_adp_c, direction = :out)])
        @test resolve_coupling(AbstractSubModel[ToyPool(), inert]) isa CouplingGraph
        err = caught(() -> build_problem(AbstractSubModel[ToyPool(), inert];
                                         tspan = (0.0, 10.0)))
        @test err isa ArgumentError
        @test occursin("nothing in this composition executes", err.msg)
        @test occursin("currency on :M_adp_c (out) declared by ToyExpression", err.msg)

        # The same pair without the extra edge builds, and every edge runs.
        ok = AbstractSubModel[ToyPool(), ToyExpression()]
        @test isempty(unexecuted_edges(ok, build_problem(ok; tspan = (0.0, 10.0))))
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
