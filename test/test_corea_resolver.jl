using Test
using InferCell

@testset "Core A′ coupling resolver" begin

    @testset "A consistent composition resolves" begin
        producer = CoreAStub(:Central;
            st = [:M_atp_c, :M_g6p_c],
            edges = [MassEdge(species=:M_atp_c, direction=:out, peer=:Expression)])
        consumer = CoreAStub(:Expression;
            edges = [MassEdge(species=:M_atp_c, direction=:in, peer=:Central)])

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
        err = try
            resolve_coupling([bad])
        catch e
            e
        end
        @test err isa ArgumentError
        @test occursin("M_not_a_species_c", err.msg)
        @test occursin("Central", err.msg)
    end

    @testset "A missing peer is rejected, and distinguished from a lone module" begin
        orphan = CoreAStub(:Central;
            st = [:M_atp_c],
            edges = [MassEdge(species=:M_atp_c, direction=:out, peer=:NotHere)])
        err = try
            resolve_coupling([orphan])
        catch e
            e
        end
        @test err isa ArgumentError
        @test occursin("NotHere", err.msg)

        # The same module with an unnamed peer is a deliberate partial
        # composition and resolves cleanly.
        lone = CoreAStub(:Central;
            st = [:M_atp_c],
            edges = [MassEdge(species=:M_atp_c, direction=:out)])
        @test length(resolve_coupling([lone]).edges) == 1
    end

    @testset "Two modules disagreeing on an edge kind is rejected" begin
        as_mass = CoreAStub(:Central;
            st = [:M_atp_c],
            edges = [MassEdge(species=:M_atp_c, direction=:in)])
        as_counter = CoreAStub(:Expression;
            edges = [DeferredCounterEdge(species=:M_atp_c, direction=:in,
                                         counter=:ATP_trsc)])

        err = try
            resolve_coupling([as_mass, as_counter])
        catch e
            e
        end
        @test err isa ArgumentError
        # Names both modules, the species, and both kinds.
        @test occursin("Central", err.msg)
        @test occursin("Expression", err.msg)
        @test occursin("M_atp_c", err.msg)
        @test occursin("mass", err.msg)
        @test occursin("deferred_counter", err.msg)
        # And distinguishes the semantics rather than reporting a generic clash.
        @test occursin("continuous", err.msg)
        @test occursin("debited a step later", err.msg)
    end

    @testset "A module may not integrate a chemostat" begin
        bad = CoreAStub(:Central; st = [:M_ctp_c])
        err = try
            resolve_coupling([bad])
        catch e
            e
        end
        @test err isa ArgumentError
        @test occursin("M_ctp_c", err.msg)
        @test occursin("Central", err.msg)
    end

    @testset "A dynamic state may not be owned twice" begin
        one = CoreAStub(:Central; st = [:M_atp_c])
        two = CoreAStub(:Nucleotide; st = [:M_atp_c])
        err = try
            resolve_coupling([one, two])
        catch e
            e
        end
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

    @testset "A cost with no paying state is an error" begin
        # Nothing integrates AMP and the registry does not chemostat it.
        spender = CoreAStub(:Charging;
            edges = [DeferredCounterEdge(species=:M_amp_c, direction=:in,
                                         counter=:AMP_charging)])
        err = try
            resolve_coupling([spender])
        catch e
            e
        end
        @test err isa ArgumentError
        @test occursin("M_amp_c", err.msg)
        @test occursin("Charging", err.msg)
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
        @test graph.dead_ends[1].consumers == [:Charging]
        @test occursin("adenylate", dead_end_report(graph))

        guanylate = CoreAStub(:Transcription;
            st = [:M_gmp_c],
            edges = [MassEdge(species=:M_gmp_c, direction=:in)])
        g2 = resolve_coupling([guanylate])
        @test g2.dead_ends[1].species == :M_gmp_c
        @test g2.dead_ends[1].moiety == :guanylate
        @test occursin("guanylate", dead_end_report(g2))
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
        err = try
            resolve_coupling([drifted])
        catch e
            e
        end
        @test err isa ArgumentError
        @test occursin("M_gtp_c", err.msg)
        @test occursin("Metabolism", err.msg)

        consistent = CoreAStub(:Metabolism;
            st = [:M_atp_c],
            edges = [MassEdge(species=:M_atp_c, direction=:in)],
            ins = [:M_atp_c])
        @test length(resolve_coupling([consistent]).edges) == 1
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
                                         counter=:ATP_trsc, clip=:smoothed)])
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
end
