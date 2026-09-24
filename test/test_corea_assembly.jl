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

    @testset "13.5 a full cycle at published parameters" begin
        # The whole 6,300 s, never shortened (spec §3). Warm, it costs about
        # 32 s of stepping (job 17298242), so it stays in the default suite
        # rather than behind a gate (task 13.8). The hook throws on any
        # non-success solver retcode, so completing is the no-failure check.
        d = build_corea()
        @test solver_settings(d) == (solver = :Rodas5P, abstol = 1e-10, reltol = 1e-8)
        out = run_handshake!(d, round(Int, COREA_CYCLE_S))
        @test d.ode.t == COREA_CYCLE_S
        @test last(out.t) == COREA_CYCLE_S
        @test d.n_handshakes == 6300
        @test all(isfinite, d.ode.u)
        @test growth_report(d).fractional > 1
    end

    @testset "13.6 what is ours, across the whole composition" begin
        ms = corea_models()
        d = build_corea(tspan = (0.0, 60.0))
        labels = reduction_declarations(ms, d)
        of(c) = [l for l in labels if l.category === c]
        has(c, subject, text = "") = any(l -> l.subject === subject &&
                                              occursin(text, l.description), of(c))

        # The lumping once, and its deterministic formalism as a second,
        # separate declaration. Translation's per-amino-acid note mentions the
        # lumping and is a model note, not a second lumping.
        @test length(of(:lumping)) == 1 && has(:lumping, :TrnaCharging)
        @test length(of(:formalism)) == 1 && has(:formalism, :TrnaCharging, "deterministic")

        # The fourteen asserted priors an inference might free (spec §4 D7).
        fourteen = [Symbol(k, :_glcpts, i) for i in 0:4 for k in (:kf, :kr)]
        append!(fourteen, [:p_lact2r, :k_chg, :trna_pool_mM, :trna_charged_fraction])
        @test length(fourteen) == 14
        @test issubset(fourteen, [l.subject for l in of(:asserted_prior)])

        # Asserted priors on slots the driver overwrites, reported apart: the
        # fifteen translated-enzyme slots and the cell radius the inbound
        # volume channel fills.
        written = Set(w.param_slot for w in driver_written_params(d))
        @test all(l -> l.subject in written, of(:discarded_prior))
        @test !any(l -> l.subject in written, of(:asserted_prior))
        @test count(l -> startswith(string(l.subject), "enz_"), of(:discarded_prior)) == 15
        @test has(:discarded_prior, :r_cell_nm)
        @test length(of(:discarded_prior)) == 16

        # The chemostatted pools: CTP and UTP clamped by transcription, and the
        # amino-acid pool, which nothing reads.
        @test has(:clamp, :M_ctp_c) && has(:clamp, :M_utp_c)
        @test has(:model_note, :CoreATranscription, "chemostat")
        @test has(:model_note, :TrnaCharging, "M_aa_pool_c")

        # The module notes the list names, each from the module that owns it.
        @test has(:model_note, :CoreATranscription, "mapping")
        @test has(:model_note, :CentralGlycolysis, "balanced")      # Km column
        @test has(:model_note, :PtsTransport, "medium-to-cell")     # lactate ratio

        # The driver's: rounding and drain granularity on every report, the
        # capped geometry, and exogenous membrane growth.
        @test has(:driver_policy, :fractional_carry)
        @test has(:driver_policy, :drain_interval)
        @test has(:capped_rate_law_geometry, :inbound_volume_channel)
        @test has(:exogenous_growth, :membrane_area_baseline)

        # The report counts departures and leaves the policy records out.
        @test occursin("$(length(labels) - 2) declaration(s)", reduction_report(ms, d))
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
