using Test
using InferCell
using Distributions

@testset "Core A′ coupling resolver" begin

    @testset "A consistent composition resolves" begin
        producer = CoreAStub(:Central;
            st = [:M_atp_c, :M_g6p_c],
            edges = [MassEdge(species=:M_atp_c, direction=:out, peer=:Expression)])
        consumer = CoreAStub(:Expression;
            edges = [MassEdge(species=:M_atp_c, direction=:in, peer=:Central)],
            ins = [:M_atp_c])

        graph = resolve_coupling([producer, consumer])

        @test length(graph.edges) == 2
        @test all(r -> r.kind == :mass, graph.edges)
        @test Set(r.declared_by for r in graph.edges) == Set([:Central, :Expression])
        @test all(r -> r.species == :M_atp_c, graph.edges)

        # Registry positions come back with the edge.
        @test all(r -> r.state_index == species_index(:M_atp_c), graph.edges)

        # ATP has both a producer and a consumer, so it is not a dead end.
        @test isempty(graph.dead_ends)
    end

    @testset "Resolution needs no ODEProblem and works on one module" begin
        lone = CoreAStub(:Central;
            st = [:M_atp_c],
            edges = [MassEdge(species=:M_atp_c, direction=:out)])
        graph = resolve_coupling(lone)          # single-model method
        @test length(graph.edges) == 1
        @test graph.edges[1].peer === nothing
    end

    @testset "An unregistered species is rejected, naming it" begin
        bad = CoreAStub(:Central;
            edges = [MassEdge(species=:M_not_a_species_c, direction=:out)])
        err = caught(() -> resolve_coupling([bad]))
        @test err isa ArgumentError
        @test occursin("M_not_a_species_c", err.msg)
        @test occursin("Central", err.msg)
    end

    @testset "A missing peer is rejected, and distinguished from a lone module" begin
        orphan = CoreAStub(:Central;
            st = [:M_atp_c],
            edges = [MassEdge(species=:M_atp_c, direction=:out, peer=:NotHere)])
        err = caught(() -> resolve_coupling([orphan]))
        @test err isa ArgumentError
        @test occursin("NotHere", err.msg)

        # The same module with an unnamed peer is a deliberate partial
        # composition and resolves cleanly.
        lone = CoreAStub(:Central;
            st = [:M_atp_c],
            edges = [MassEdge(species=:M_atp_c, direction=:out)])
        @test length(resolve_coupling([lone]).edges) == 1
    end

    @testset "Two modules disagreeing on the transport description is rejected" begin
        # Mass and currency are two spellings of one continuous crossing —
        # direct shared state versus the same state routed via a pool — so one
        # (species, direction) admits only one of them.
        as_mass = CoreAStub(:Central;
            st = [:M_atp_c],
            edges = [MassEdge(species=:M_atp_c, direction=:in)])
        as_currency = CoreAStub(:Expression;
            edges = [CurrencyEdge(species=:M_atp_c, direction=:in)],
            ins = [:M_atp_c])

        err = caught(() -> resolve_coupling([as_mass, as_currency]))
        @test err isa ArgumentError
        # Names both modules, the species, and both kinds.
        @test occursin("Central", err.msg)
        @test occursin("Expression", err.msg)
        @test occursin("M_atp_c", err.msg)
        @test occursin("mass", err.msg)
        @test occursin("currency", err.msg)
        # And distinguishes the semantics rather than reporting a generic clash.
        @test occursin("continuous", err.msg)
        @test occursin("pool", err.msg)
    end

    @testset "Distinct mechanisms coexist on one species and direction" begin
        # The published model routes ATP through a currency pool, a
        # deferred-counter debit and a rate-constant rebuild simultaneously —
        # three kinds, two directions. The kind-agreement check is scoped to
        # the transport descriptions (mass, currency) precisely so this
        # flagship boundary is declarable.
        owner = CoreAStub(:Nucleotide;
            st = [:M_atp_c],
            edges = [CurrencyEdge(species=:M_atp_c, direction=:out)])
        expression = CoreAStub(:Expression;
            edges = [CurrencyEdge(species=:M_atp_c, direction=:in),
                     DeferredCounterEdge(species=:M_atp_c, direction=:in,
                                         counter=:ATP_trsc),
                     RateConstantEdge(species=:M_atp_c, direction=:in)],
            ins = [:M_atp_c])

        graph = resolve_coupling([owner, expression])
        @test length(graph.edges) == 4
        @test Set(r.kind for r in graph.edges) ==
              Set([:currency, :deferred_counter, :rate_constant])
        # The currency return path covers the debit, so nothing is stranded.
        @test isempty(graph.dead_ends)

        # A mass consumer beside a deferred-counter debit is likewise two
        # mechanisms, not a disagreement: the ODE block consumes continuously
        # while the hook debits the same pool a step later.
        continuous = CoreAStub(:Central;
            st = [:M_atp_c],
            edges = [MassEdge(species=:M_atp_c, direction=:in)])
        debit = CoreAStub(:Hook;
            edges = [DeferredCounterEdge(species=:M_atp_c, direction=:in,
                                         counter=:ATP_translat)])
        @test length(resolve_coupling([continuous, debit]).edges) == 2
    end

    @testset "One module describing one crossing two ways is also rejected" begin
        both = CoreAStub(:Central;
            st = [:M_atp_c],
            edges = [MassEdge(species=:M_atp_c, direction=:in),
                     CurrencyEdge(species=:M_atp_c, direction=:in)],
            ins = [:M_atp_c])
        err = caught(() -> resolve_coupling([both]))
        @test err isa ArgumentError
        # Named as one module describing one crossing, not as two modules
        # disagreeing with themselves.
        @test occursin("Module Central describes", err.msg)
        @test occursin("mass", err.msg)
        @test occursin("currency", err.msg)
    end

    @testset "A module may not integrate a chemostat" begin
        bad = CoreAStub(:Central; st = [:M_ctp_c])
        err = caught(() -> resolve_coupling([bad]))
        @test err isa ArgumentError
        @test occursin("M_ctp_c", err.msg)
        @test occursin("Central", err.msg)
    end

    @testset "A dynamic state may not be owned twice" begin
        one = CoreAStub(:Central; st = [:M_atp_c])
        two = CoreAStub(:Nucleotide; st = [:M_atp_c])
        err = caught(() -> resolve_coupling([one, two]))
        @test err isa ArgumentError
        @test occursin("M_atp_c", err.msg)
        @test occursin("Central", err.msg)
        @test occursin("Nucleotide", err.msg)
    end

    @testset "Unowned states are reported, not thrown" begin
        # A single module under test legitimately leaves most states unowned.
        lone = CoreAStub(:Central; st = [:M_atp_c])
        graph = resolve_coupling([lone])
        @test :M_atp_c ∉ graph.unowned_states
        @test :M_g6p_c ∈ graph.unowned_states
        @test length(graph.unowned_states) == n_dynamic_states() - 1
        # Chemostats are never "unowned" — nothing is supposed to integrate them.
        @test :M_ctp_c ∉ graph.unowned_states
    end

    @testset "A cost with no paying state is reported, not thrown" begin
        # Nothing integrates AMP and the registry does not chemostat it. A
        # module validated alone legitimately imports species whose owners are
        # absent — the standalone-validation constraint — so this resolves and
        # reports rather than failing.
        spender = CoreAStub(:Charging;
            edges = [DeferredCounterEdge(species=:M_amp_c, direction=:in,
                                         counter=:AMP_charging)])
        graph = resolve_coupling([spender])
        @test :M_amp_c in graph.unowned_states
        @test length(graph.dead_ends) == 1
        @test graph.dead_ends[1].species == :M_amp_c
        @test graph.dead_ends[1].missing_role == :producer
        @test graph.dead_ends[1].modules == [:Charging]
    end

    @testset "A consumed species with no producer is reported with its moiety" begin
        # The adenylate dead end: charging converts ATP to AMP and, without
        # ADK1, nothing returns it. Same class of error as stranded GMP.
        adenylate = CoreAStub(:Charging;
            st = [:M_amp_c],
            edges = [MassEdge(species=:M_amp_c, direction=:in)])
        graph = resolve_coupling([adenylate])
        @test length(graph.dead_ends) == 1
        @test graph.dead_ends[1].species == :M_amp_c
        @test graph.dead_ends[1].moiety == :adenylate
        @test graph.dead_ends[1].missing_role == :producer
        @test graph.dead_ends[1].modules == [:Charging]
        @test occursin("adenylate", dead_end_report(graph))

        guanylate = CoreAStub(:Transcription;
            st = [:M_gmp_c],
            edges = [MassEdge(species=:M_gmp_c, direction=:in)])
        g2 = resolve_coupling([guanylate])
        @test g2.dead_ends[1].species == :M_gmp_c
        @test g2.dead_ends[1].moiety == :guanylate
        @test occursin("guanylate", dead_end_report(g2))
    end

    @testset "A produced species with no consumer is the mirror dead end" begin
        # The spec requires detecting both halves: no producer, and no consumer.
        producer = CoreAStub(:Nucleotide;
            edges = [MassEdge(species=:M_amp_c, direction=:out)])
        owner = CoreAStub(:Charging; st = [:M_amp_c])
        graph = resolve_coupling([producer, owner])
        @test length(graph.dead_ends) == 1
        @test graph.dead_ends[1].species == :M_amp_c
        @test graph.dead_ends[1].missing_role == :consumer
        @test graph.dead_ends[1].modules == [:Nucleotide]
        @test occursin("nothing drawing it down", dead_end_report(graph))
    end

    @testset "Adding the producer closes the dead end" begin
        # ADK1 returns AMP to ADP, so the moiety is no longer stranded.
        charging = CoreAStub(:Charging;
            st = [:M_amp_c],
            edges = [MassEdge(species=:M_amp_c, direction=:in)])
        adk1 = CoreAStub(:Nucleotide;
            edges = [MassEdge(species=:M_amp_c, direction=:out)])
        graph = resolve_coupling([charging, adk1])
        @test isempty(graph.dead_ends)
        @test occursin("No dead ends", dead_end_report(graph))
    end

    @testset "Chemostatted species are exempt, and the exemption is enumerable" begin
        # CTP is chemostatted, so the chemostat is the return path and no
        # dead end is raised. Recorded so a later change making it live
        # reinstates the check.
        user = CoreAStub(:Transcription;
            edges = [MassEdge(species=:M_ctp_c, direction=:in)])
        graph = resolve_coupling([user])
        @test isempty(graph.dead_ends)
        @test :M_ctp_c in graph.chemostat_exemptions

        # The producer side is exempted and recorded the same way: the
        # chemostat absorbs what is produced into it.
        feeder = CoreAStub(:Nucleotide;
            edges = [MassEdge(species=:M_utp_c, direction=:out)])
        g2 = resolve_coupling([feeder])
        @test isempty(g2.dead_ends)
        @test :M_utp_c in g2.chemostat_exemptions
    end

    @testset "Catalytic edges are excluded from the dead-end accounting" begin
        # Counts entering a rate law strand no moiety, so a catalytic edge with
        # no counterpart producer is not a dead end.
        reader = CoreAStub(:Central;
            st = [:M_ptsg_c],
            edges = [CatalyticEdge(species=:M_ptsg_c, direction=:in,
                                   param_slot=:enzyme_conc)])
        graph = resolve_coupling([reader])
        @test isempty(graph.dead_ends)
        @test length(graph.edges) == 1
        @test !carries_mass(graph.edges[1].edge)
    end

    @testset "inputs() must not drift from the typed declaration" begin
        drifted = CoreAStub(:Metabolism;
            st = [:M_atp_c],
            edges = [MassEdge(species=:M_atp_c, direction=:in)],
            ins = [:M_gtp_c])          # declared as an input, but no inbound edge
        err = caught(() -> resolve_coupling([drifted]))
        @test err isa ArgumentError
        @test occursin("M_gtp_c", err.msg)
        @test occursin("Metabolism", err.msg)

        consistent = CoreAStub(:Metabolism;
            st = [:M_atp_c],
            edges = [MassEdge(species=:M_atp_c, direction=:in)],
            ins = [:M_atp_c])
        @test length(resolve_coupling([consistent]).edges) == 1
    end

    @testset "A legacy non-registry input coexists with typed coupling" begin
        # A hybrid module bridging Core A′ and the pre-contract blocks reads
        # legacy state (mRNA) through inputs(), which no edge kind can name.
        # Only registry species are held to the typed contract, so the two
        # channels are not mutually exclusive.
        hybrid = CoreAStub(:Expression;
            st = [:M_atp_c],
            edges = [RateConstantEdge(species=:M_gtp_c, direction=:in)],
            ins = [:mRNA])
        graph = resolve_coupling([hybrid])
        @test length(graph.edges) == 1
    end

    @testset "An inbound mass or currency edge must be listed in inputs()" begin
        # The converse drift: only inputs() wires a state into dynamics, so an
        # inbound mass or currency edge on a state the module does not
        # integrate would silently never arrive.
        silent = CoreAStub(:Expression;
            edges = [MassEdge(species=:M_atp_c, direction=:in)])
        err = caught(() -> resolve_coupling([silent]))
        @test err isa ArgumentError
        @test occursin("M_atp_c", err.msg)
        @test occursin("inputs()", err.msg)
        @test occursin("Expression", err.msg)

        # A currency consumer is continuous shared state exactly like mass, so
        # the same omission would silently starve it of its declared coupling.
        silent_currency = CoreAStub(:Expression;
            edges = [CurrencyEdge(species=:M_atp_c, direction=:in)])
        err2 = caught(() -> resolve_coupling([silent_currency]))
        @test err2 isa ArgumentError
        @test occursin("currency", err2.msg)
        @test occursin("inputs()", err2.msg)

        # A module integrating the species itself needs no input for it, and a
        # deferred counter is not wired through inputs() at all — the hook
        # debits it, the RHS never reads it.
        own = CoreAStub(:Central;
            st = [:M_atp_c],
            edges = [MassEdge(species=:M_atp_c, direction=:in)])
        @test length(resolve_coupling([own]).edges) == 1

        counter = CoreAStub(:Expression;
            st = [:M_atp_c],
            edges = [DeferredCounterEdge(species=:M_atp_c, direction=:in,
                                         counter=:ATP_trsc)])
        @test length(resolve_coupling([counter]).edges) == 1
    end

    @testset "A currency pool must be a registry species" begin
        typo = CoreAStub(:Central;
            st = [:M_atp_c],
            edges = [CurrencyEdge(species=:M_atp_c, direction=:out,
                                  pool=:M_atp_typo_c)])
        err = caught(() -> resolve_coupling([typo]))
        @test err isa ArgumentError
        @test occursin("M_atp_typo_c", err.msg)
        @test occursin("pool", err.msg)
    end

    @testset "A clamp must hold a chemostat at the registry's value" begin
        # Glucose is chemostatted at 40 mM. Two modules clamping the medium
        # differently would each pass alone and silently disagree at the join.
        wrong = CoreAStub(:Transport;
            edges = [ClampedEdge(species=:M_glc__D_e, direction=:in,
                                 origin=:published, held_value=20.0)])
        err = caught(() -> resolve_coupling([wrong]))
        @test err isa ArgumentError
        @test occursin("20.0", err.msg)
        @test occursin("40.0", err.msg)
        @test occursin("M_glc__D_e", err.msg)

        right = CoreAStub(:Transport;
            edges = [ClampedEdge(species=:M_glc__D_e, direction=:in,
                                 origin=:published, held_value=40.0)])
        @test length(resolve_coupling([right]).edges) == 1

        # Agreement is up to the registry's 4-decimal transcription — the same
        # tolerance as the loader — so a full-precision table value passes.
        full_precision = CoreAStub(:Transport;
            edges = [ClampedEdge(species=:M_glc__D_e, direction=:in,
                                 origin=:published, held_value=40.00002)])
        @test length(resolve_coupling([full_precision]).edges) == 1

        # A chemostat with no imported value constrains nothing yet.
        open_value = CoreAStub(:Transcription;
            edges = [ClampedEdge(species=:M_ctp_c, direction=:in,
                                 origin=:ours, held_value=1.0)])
        @test length(resolve_coupling([open_value]).edges) == 1
    end

    @testset "A clamp on a species another module integrates is rejected" begin
        # A dependence cannot be both replaced by a constant and live in one
        # composition: the clamped module would run against a frozen value
        # while the owner evolves it.
        owner = CoreAStub(:Nucleotide; st = [:M_gtp_c])
        clamper = CoreAStub(:Transcription;
            edges = [ClampedEdge(species=:M_gtp_c, direction=:in,
                                 origin=:ours, held_value=1.6627)])
        err = caught(() -> resolve_coupling([owner, clamper]))
        @test err isa ArgumentError
        @test occursin("M_gtp_c", err.msg)
        @test occursin("Transcription", err.msg)
        @test occursin("Nucleotide", err.msg)

        # The clamp alone is the standalone case — the owner is absent by
        # construction, so a module author validating alone still passes.
        @test length(resolve_coupling([clamper]).edges) == 1
    end

    @testset "Two coupled modules sharing one module_id are rejected" begin
        one = CoreAStub(:Central; st = [:M_atp_c])
        other = CoreAStub(:Central; st = [:M_g6p_c])
        err = caught(() -> resolve_coupling([one, other]))
        @test err isa ArgumentError
        @test occursin("Central", err.msg)
        @test occursin("module_id", err.msg)

        # Two instances of a legacy sub-model — no edges, no registry states —
        # compose exactly as they did before the contract.
        graph = resolve_coupling([TranscriptionTranslation(),
                                  TranscriptionTranslation()])
        @test isempty(graph.edges)
    end

    @testset "Two clamps at different held values are rejected" begin
        # CTP has no imported registry value, so the edge-vs-registry check
        # cannot fire; the two edges must still agree with each other, or
        # composition order silently decides the held concentration.
        one = CoreAStub(:Transcription;
            edges = [ClampedEdge(species=:M_ctp_c, direction=:in,
                                 origin=:ours, held_value=1.0)])
        other = CoreAStub(:Nucleotide;
            edges = [ClampedEdge(species=:M_ctp_c, direction=:in,
                                 origin=:ours, held_value=2.0)])
        err = caught(() -> resolve_coupling([one, other]))
        @test err isa ArgumentError
        @test occursin("M_ctp_c", err.msg)
        @test occursin("1.0", err.msg)
        @test occursin("2.0", err.msg)
        @test occursin("Transcription", err.msg)
        @test occursin("Nucleotide", err.msg)

        # Agreeing clamps resolve cleanly.
        agreeing = CoreAStub(:Nucleotide;
            edges = [ClampedEdge(species=:M_ctp_c, direction=:in,
                                 origin=:ours, held_value=1.0)])
        @test length(resolve_coupling([one, agreeing]).edges) == 2
    end

    @testset "Gradient obstructions are collected and reported" begin
        clamped = CoreAStub(:Expression;
            st = [:M_atp_c],
            edges = [DeferredCounterEdge(species=:M_atp_c, direction=:in,
                                         counter=:ATP_trsc)],
            mode = :simulation)
        graph = resolve_coupling([clamped])
        @test length(graph.gradient_obstructions) == 1
        @test occursin("max(0, ·)", gradient_report(graph))

        smoothed = CoreAStub(:Expression;
            st = [:M_atp_c],
            edges = [DeferredCounterEdge(species=:M_atp_c, direction=:in,
                                         counter=:ATP_trsc, clip=:smoothed,
                                         smoothing=0.05)])
        @test isempty(resolve_coupling([smoothed]).gradient_obstructions)
        @test occursin("No gradient obstructions",
                       gradient_report(resolve_coupling([smoothed])))
    end

    @testset "A differentiable composition warns about a clamped boundary" begin
        clamped = CoreAStub(:Expression;
            st = [:M_atp_c],
            edges = [DeferredCounterEdge(species=:M_atp_c, direction=:in,
                                         counter=:ATP_trsc)],
            mode = :differentiable)
        obstructions = @test_logs (:warn,) match_mode=:any check_gradient_safety([clamped])
        @test length(obstructions) == 1
        # A warning, not a throw: a clamped interface is a legitimate model,
        # just not a differentiable one.
        @test obstructions[1].species == :M_atp_c
    end

    @testset "Existing compositions are unaffected" begin
        # Sub-models outside Core A′ name no registry species and declare no
        # typed edges, so resolution is a no-op for them.
        txl = TranscriptionTranslation()
        metab = LightMetabolism()
        graph = resolve_coupling([txl, metab])
        @test isempty(graph.edges)
        @test isempty(graph.dead_ends)
        @test isempty(graph.gradient_obstructions)
        # Every registry state is unowned, because these models own none of them.
        @test length(graph.unowned_states) == n_dynamic_states()

        # And build_problem, which now runs the resolver, still works.
        prob = build_problem([txl, metab])
        @test length(prob.u0) == length(states(txl)) + length(states(metab))
    end

    @testset "A chemostatted input fails the contract with the actual remedy" begin
        # inputs() resolves against integrated states only, and nothing may
        # integrate a chemostat, so the wiring can never deliver one. The
        # resolver rejects the declaration outright — build_problem would
        # otherwise fail later whatever edges the module declared — and the
        # error names the fix: drop the input, keep the ClampedEdge, carry the
        # value as a fixed parameter until wave 2 executes coupling.
        reader = CoreAStub(:Transcription;
            st = [:M_atp_c],
            edges = [ClampedEdge(species=:M_ctp_c, direction=:in,
                                 origin=:ours, held_value=1.0)],
            ins = [:M_ctp_c])
        err = caught(() -> resolve_coupling([reader]))
        @test err isa ArgumentError
        @test occursin("M_ctp_c", err.msg)
        @test occursin("chemostat", err.msg)
        @test occursin("Remove :M_ctp_c from inputs()", err.msg)
        @test occursin("ClampedEdge", err.msg)

        # A module outside the typed contract — no coupling declared — skips
        # the resolver's inputs checks, so the orchestrator keeps the same
        # diagnostic as a backstop.
        legacy = CoreAStub(:Transcription;
            st = [:M_atp_c],
            ins = [:M_ctp_c])
        err2 = caught(() -> build_problem([legacy]))
        @test err2 isa ErrorException
        @test occursin("M_ctp_c", err2.msg)
        @test occursin("chemostatted", err2.msg)
        @test occursin("Remove :M_ctp_c from inputs()", err2.msg)
        @test occursin("ClampedEdge", err2.msg)
    end

    @testset "The jump path validates shared parameters too" begin
        # Two :jump modules sharing a parameter name with different values must
        # fail the same way an ODE composition does; the check runs before any
        # reactions are assembled.
        p_a = InferParameter(1.0, Normal(0, 1), false, :k_shared, :A, :rate)
        p_b = InferParameter(2.0, Normal(0, 1), false, :k_shared, :B, :rate)
        a = CoreAStub(:A; params = [p_a], form = :jump)
        b = CoreAStub(:B; params = [p_b], form = :jump)
        err = caught(() -> build_problem([a, b]))
        @test err isa ErrorException
        @test occursin("k_shared", err.msg)

        # And the cross-file provenance report fires for agreeing values from
        # different source files.
        tagged_a = InferParameter(1.0, Normal(0, 1), false, :k_shared, :A, :rate,
                                  ParameterSource("central_balanced"))
        tagged_b = InferParameter(1.0, Normal(0, 1), false, :k_shared, :B, :rate,
                                  ParameterSource("nucleotide_balanced"))
        @test_logs (:warn, r"more than one source file") InferCell._validate_shared_params(
            [CoreAStub(:A; params = [tagged_a]), CoreAStub(:B; params = [tagged_b])])
    end
end
