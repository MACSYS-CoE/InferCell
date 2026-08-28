using Test
using InferCell

@testset "Core A′ coupling edge kinds" begin

    @testset "All seven kinds construct and name themselves" begin
        edges = [
            MassEdge(species=:M_atp_c, direction=:out),
            CurrencyEdge(species=:M_atp_c, direction=:in),
            DeferredCounterEdge(species=:M_atp_c, direction=:in, counter=:ATP_trsc),
            CatalyticEdge(species=:M_ptsg_c, direction=:in, param_slot=:enzyme_conc),
            RateConstantEdge(species=:M_gtp_c, direction=:out),
            VolumeEdge(species=:M_ptsg_c, direction=:out),
            ClampedEdge(species=:M_ctp_c, direction=:in, origin=:ours),
        ]
        @test [edge_kind(e) for e in edges] == collect(EDGE_KINDS)
        @test length(EDGE_KINDS) == 7
        @test all(e -> e isa CouplingEdge, edges)
    end

    @testset "A missing required field fails at construction, naming field and kind" begin
        # Not at composition, and not at run time.
        err = caught(() -> DeferredCounterEdge(species=:M_atp_c, direction=:in))
        @test err isa ArgumentError
        @test occursin("counter", err.msg)
        @test occursin("DeferredCounterEdge", err.msg)

        err2 = caught(() -> CatalyticEdge(species=:M_ptsg_c, direction=:in))
        @test err2 isa ArgumentError
        @test occursin("param_slot", err2.msg)
        @test occursin("CatalyticEdge", err2.msg)

        @test_throws ArgumentError MassEdge(direction=:in)               # no species
        @test_throws ArgumentError MassEdge(species=:M_atp_c)            # no direction
        @test_throws ArgumentError ClampedEdge(species=:M_ctp_c, direction=:in)  # no origin
    end

    @testset "Positional construction cannot bypass validation" begin
        # Validation lives in the inner constructors, so the positional form is
        # checked exactly like the keyword one.
        @test_throws ArgumentError MassEdge(:M_atp_c, :sideways, nothing)
        @test_throws ArgumentError DeferredCounterEdge(:M_atp_c, :in, nothing,
                                                       :c, :hope, nothing)
        @test_throws ArgumentError RateConstantEdge(:M_gtp_c, :out, nothing,
                                                    :continuous, 60.0)
        @test_throws ArgumentError ClampedEdge(:M_ctp_c, :in, nothing, 1.0, :borrowed)
    end

    @testset "Vocabulary is checked" begin
        @test_throws ArgumentError MassEdge(species=:M_atp_c, direction=:sideways)
        @test_throws ArgumentError DeferredCounterEdge(
            species=:M_atp_c, direction=:in, counter=:c, clip=:hope_for_the_best)
        @test_throws ArgumentError RateConstantEdge(
            species=:M_gtp_c, direction=:out, cadence=:whenever)
        @test_throws ArgumentError ClampedEdge(
            species=:M_ctp_c, direction=:in, origin=:borrowed)
    end

    @testset "Deferred counters default to the published clipped drain" begin
        e = DeferredCounterEdge(species=:M_atp_c, direction=:in, counter=:ATP_trsc)
        @test e.clip == :clamped_deficit_carried
        @test obstructs_gradients(e)              # the max(0, ·) on the interface
        @test !deviates_from_published(e)
        @test deviation_reason(e) === nothing
    end

    @testset "Unclamped and smoothed policies are declarable, and both deviate" begin
        # Either departure from the published clamped drain is a labelled
        # deviation: an unclamped pool can go negative, which the published
        # model never allows, so a result depending on it must carry a label.
        unclamped = DeferredCounterEdge(species=:M_atp_c, direction=:in,
                                        counter=:ATP_trsc, clip=:unclamped)
        @test !obstructs_gradients(unclamped)
        @test deviates_from_published(unclamped)
        @test occursin("negative", deviation_reason(unclamped))
        @test deviation_category(unclamped) === :unclamped_counter

        smoothed = DeferredCounterEdge(species=:M_atp_c, direction=:in,
                                       counter=:ATP_trsc, clip=:smoothed,
                                       smoothing=0.05)
        @test !obstructs_gradients(smoothed)
        @test smoothed.smoothing == 0.05
        @test deviates_from_published(smoothed)   # differs from the published model
        @test occursin("smoothed", deviation_reason(smoothed))
        @test occursin("0.05", deviation_reason(smoothed))
        @test deviation_category(smoothed) === :smoothed_counter
    end

    @testset "The smoothing parameter is exposed, not hidden" begin
        # A smoothed clip without its width is underdeclared, and the other
        # policies take no width at all.
        err = caught() do
            DeferredCounterEdge(species=:M_atp_c, direction=:in,
                                counter=:ATP_trsc, clip=:smoothed)
        end
        @test err isa ArgumentError
        @test occursin("smoothing", err.msg)

        @test_throws ArgumentError DeferredCounterEdge(
            species=:M_atp_c, direction=:in, counter=:ATP_trsc,
            clip=:smoothed, smoothing=-1.0)
        @test_throws ArgumentError DeferredCounterEdge(
            species=:M_atp_c, direction=:in, counter=:ATP_trsc,
            smoothing=0.05)   # clamped policy takes no smoothing

        published = DeferredCounterEdge(species=:M_atp_c, direction=:in,
                                        counter=:ATP_trsc)
        @test published.smoothing === nothing
    end

    @testset "Rate-constant edges default to the published 60 s rebuild" begin
        e = RateConstantEdge(species=:M_gtp_c, direction=:out)
        @test e.cadence == :piecewise_constant
        @test e.interval == 60.0
        @test !deviates_from_published(e)

        custom = RateConstantEdge(species=:M_gtp_c, direction=:out, interval=30.0)
        @test custom.interval == 30.0
        @test !deviates_from_published(custom)

        @test_throws ArgumentError RateConstantEdge(species=:M_gtp_c, direction=:out,
                                                    interval=nothing)
        @test_throws ArgumentError RateConstantEdge(species=:M_gtp_c, direction=:out,
                                                    interval=-1.0)
    end

    @testset "A continuous cadence is marked as a deviation and carries no interval" begin
        e = RateConstantEdge(species=:M_gtp_c, direction=:out, cadence=:continuous)
        @test e.cadence == :continuous
        @test e.interval === nothing      # carrying one would imply a cadence it lacks
        @test deviates_from_published(e)
        @test occursin("continuous", deviation_reason(e))

        # An interval passed explicitly is validated, not silently discarded.
        @test_throws ArgumentError RateConstantEdge(species=:M_gtp_c, direction=:out,
                                                    cadence=:continuous, interval=30.0)
    end

    @testset "Clamped edges record whether the clamp is ours" begin
        ours = ClampedEdge(species=:M_ctp_c, direction=:in, origin=:ours, held_value=1.0)
        @test ours.origin == :ours
        @test ours.held_value == 1.0
        @test deviates_from_published(ours)
        @test occursin("this reduction", deviation_reason(ours))

        published = ClampedEdge(species=:M_glc__D_e, direction=:in,
                                origin=:published, held_value=40.0)
        @test !deviates_from_published(published)
        @test deviation_reason(published) === nothing
    end

    @testset "Direction gives producer and consumer" begin
        producer = MassEdge(species=:M_atp_c, direction=:out)
        consumer = MassEdge(species=:M_atp_c, direction=:in)
        @test is_producer(producer) && !is_consumer(producer)
        @test is_consumer(consumer) && !is_producer(consumer)
        @test mass_contribution(producer) == 1
        @test mass_contribution(consumer) == -1
    end

    @testset "Only three kinds carry mass" begin
        @test carries_mass(MassEdge(species=:M_atp_c, direction=:in))
        @test carries_mass(CurrencyEdge(species=:M_atp_c, direction=:in))
        @test carries_mass(DeferredCounterEdge(species=:M_atp_c, direction=:in,
                                               counter=:c))
        @test !carries_mass(CatalyticEdge(species=:M_ptsg_c, direction=:in,
                                          param_slot=:e))
        @test !carries_mass(RateConstantEdge(species=:M_gtp_c, direction=:out))
        @test !carries_mass(VolumeEdge(species=:M_ptsg_c, direction=:out))
        @test !carries_mass(ClampedEdge(species=:M_ctp_c, direction=:in, origin=:ours))
    end

    @testset "A catalytic edge is rejected by a conservation check, not counted as zero" begin
        # Counts entering a rate law move no matter. Returning a silent zero
        # would let the category error pass unnoticed.
        cat = CatalyticEdge(species=:M_ptsg_c, direction=:in, param_slot=:enzyme_conc)
        err = caught(() -> mass_contribution(cat))
        @test err isa ArgumentError
        @test occursin("carries no mass", err.msg)
        @test occursin("M_ptsg_c", err.msg)

        @test_throws ArgumentError mass_contribution(
            RateConstantEdge(species=:M_gtp_c, direction=:out))
        @test_throws ArgumentError mass_contribution(
            VolumeEdge(species=:M_ptsg_c, direction=:out))
    end

    @testset "The currency pool" begin
        e = CurrencyEdge(species=:M_atp_c, direction=:in, peer=:Central)
        @test e.species == :M_atp_c
        @test e.direction == :in
        @test e.peer == :Central
        @test e.pool == :M_atp_c            # defaults to the species itself
        @test CurrencyEdge(species=:M_atp_c, direction=:in, pool=:M_adp_c).pool == :M_adp_c

        @test MassEdge(species=:M_atp_c, direction=:in).peer === nothing
    end

    @testset "Existing sub-models declare no coupling" begin
        # The contract is additive: every sub-model written before it keeps
        # working and reports an empty declaration.
        for m in (TranscriptionTranslation(), StochasticGeneExpression(),
                  BurstyGeneExpression(), LightMetabolism(), TierBMetabolism())
            @test coupling(m) == CouplingEdge[]
            @test isempty(coupling(m))
            @test reduction_notes(m) == String[]
        end

        # And their existing declarations are untouched.
        @test inputs(LightMetabolism()) == [:mRNA, :protein]
        @test inputs(TranscriptionTranslation()) == Symbol[]
        @test module_id(LightMetabolism()) == :LightMetabolism
    end
end
