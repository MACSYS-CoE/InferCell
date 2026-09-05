using Test
using InferCell
using Distributions

@testset "Labelling what is ours" begin

    @testset "All four categories of deviation are enumerated" begin
        pts_constant = load_parameter(
            [read_source_table(joinpath(@__DIR__, "fixtures", "central_balanced.tsv");
                               file = "central_balanced")],
            "k_pts_GLCpts0_fwd";
            name = :k_GLCpts0_fwd, module_id = :Transport,
            prior = LogNormal(0.0, 2.0))

        composition = [
            # A clamp this reduction introduced: CTP is chemostatted because the
            # nucleotide module is out of scope, not because the model does it.
            CoreAStub(:Transcription;
                edges = [ClampedEdge(species=:M_ctp_c, direction=:in,
                                     origin=:ours, held_value=1.0)]),
            # A smoothed expression-cost drain, replacing the published max(0, ·).
            CoreAStub(:Expression;
                st = [:M_atp_c],
                edges = [DeferredCounterEdge(species=:M_atp_c, direction=:in,
                                             counter=:ATP_trsc, clip=:smoothed,
                                             smoothing=0.05),
                         MassEdge(species=:M_atp_c, direction=:out)]),
            # A continuous rebuild, where the published model holds rate
            # constants piecewise-constant between 60 s refreshes. The module
            # whose rate constants are rebuilt declares the edge inbound.
            CoreAStub(:Rebuild;
                edges = [RateConstantEdge(species=:M_gtp_c, direction=:in,
                                          cadence=:continuous)]),
            # A parameter whose prior we asserted rather than inherited.
            CoreAStub(:Transport; params = [pts_constant]),
        ]

        labels = reduction_declarations(composition)
        categories = Set(l.category for l in labels)
        @test :clamp in categories
        @test :smoothed_counter in categories
        @test :continuous_rebuild in categories
        @test :asserted_prior in categories
        @test length(categories) == 4

        @test any(l -> l.subject == :M_ctp_c, labels)
        @test any(l -> l.subject == :k_GLCpts0_fwd, labels)

        report = reduction_report(composition)
        @test occursin("clamp", report)
        @test occursin("asserted_prior", report)
        @test occursin("M_ctp_c", report)
    end

    @testset "The lumped tRNA charging step registers as a lumping" begin
        # Wave 1's add-lumped-trna-charging needs somewhere to declare that its
        # one effective charging step is ours, not the published model's.
        note = "one lumped tRNA charging step replaces the 20 per-amino-acid " *
               "synthetase chains; this lumping is ours, not the published model's"
        charging = CoreAStub(:Charging;
            st = [:M_trna_chg_c],
            notes = [note])

        labels = reduction_declarations([charging])
        lumpings = filter(l -> l.category === :lumping, labels)
        @test length(lumpings) == 1
        @test lumpings[1].subject == :Charging
        @test occursin("lumped tRNA charging", lumpings[1].description)
        @test occursin("lumping", reduction_report([charging]))
    end

    @testset "A composition that follows the published model is labelled clean" begin
        published = CoreAStub(:Central;
            st = [:M_atp_c],
            edges = [MassEdge(species=:M_atp_c, direction=:out),
                     # The published clipped drain and the published 60 s
                     # rebuild are not deviations.
                     DeferredCounterEdge(species=:M_atp_c, direction=:in,
                                         counter=:ATP_trsc),
                     RateConstantEdge(species=:M_gtp_c, direction=:in)])
        @test isempty(reduction_declarations([published]))
        @test occursin("Nothing in this composition departs",
                       reduction_report([published]))
    end

    @testset "An unclamped counter is a labelled deviation" begin
        # The pool may go negative, which the published model never allows —
        # a result depending on that must not ship unlabelled.
        unclamped = CoreAStub(:Expression;
            st = [:M_atp_c],
            edges = [DeferredCounterEdge(species=:M_atp_c, direction=:in,
                                         counter=:ATP_trsc, clip=:unclamped)])
        labels = reduction_declarations([unclamped])
        @test length(labels) == 1
        @test labels[1].category === :unclamped_counter
        @test labels[1].subject === :M_atp_c
        @test occursin("unclamped_counter", reduction_report([unclamped]))
    end

    @testset "A label category outside the vocabulary is rejected" begin
        @test_throws ArgumentError ReductionLabel(:clmap, :M_atp_c, "typo")
        # Pinned, because a category outside this tuple is counted by the
        # report's header and dropped from its body. The last two are the
        # *driver's* policy rather than any module's declaration, and are
        # enumerated by `driver_declarations` (spec §11 tasks 4.6 and 3.3);
        # they share the vocabulary because they share the report.
        @test REDUCTION_CATEGORIES == (:clamp, :smoothed_counter,
                                       :unclamped_counter, :continuous_rebuild,
                                       :coarse_drain, :rounding_policy,
                                       :asserted_prior, :lumping)
    end

    @testset "Asserted-prior labelling does not depend on composition order" begin
        # unique_params prefers the copy that carries provenance, so a module
        # constructing the same-named parameter directly cannot shadow the
        # imported copy just by being composed first.
        prior = LogNormal(0.0, 2.0)
        imported = InferParameter(1.0, prior, false, :k_pts, :Transport, :rate,
                                  ParameterSource("central_balanced";
                                                  informedness = :asserted))
        direct = InferParameter(1.0, prior, false, :k_pts, :Other, :rate)
        a = CoreAStub(:Transport; params = [imported])
        b = CoreAStub(:Other; params = [direct])

        for composition in ([a, b], [b, a])
            labels = reduction_declarations(composition)
            @test count(l -> l.category === :asserted_prior &&
                             l.subject === :k_pts, labels) == 1
        end
    end

    @testset "A parameter shared by two modules is one asserted prior, not two" begin
        shared = InferParameter(1.0, LogNormal(0.0, 2.0), false, :k_shared, :A, :rate,
                                ParameterSource("central_balanced";
                                                informedness = :asserted))
        composition = [CoreAStub(:A; params = [shared]),
                       CoreAStub(:B; params = [shared])]
        labels = reduction_declarations(composition)
        @test count(l -> l.subject === :k_shared, labels) == 1
    end

    @testset "Existing sub-models carry no reduction labels" begin
        @test isempty(reduction_declarations([TranscriptionTranslation()]))
        @test isempty(reduction_declarations(LightMetabolism()))
    end
end
