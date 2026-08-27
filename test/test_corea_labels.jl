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
            # constants piecewise-constant between 60 s refreshes.
            CoreAStub(:Rebuild;
                edges = [RateConstantEdge(species=:M_gtp_c, direction=:out,
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
                     RateConstantEdge(species=:M_gtp_c, direction=:out)])
        @test isempty(reduction_declarations([published]))
        @test occursin("Nothing in this composition departs",
                       reduction_report([published]))
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
