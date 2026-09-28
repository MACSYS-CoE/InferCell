# Spec §11 phase 14c: the smoothed drain.
#
# Task 14c.4: every consumer counter declares `clip = :smoothed` with its width,
# the width's choice is recorded with its reason, `reduction_report` carries the
# smoothing as ours (T2), and `obstructs_gradients` reports no edge.

using Test
using InferCell

# The consumer counters that can clip: inbound debits on a pool a module
# integrates. A product credit or a chemostat debit never clips.
_consumers(ms) = [e for m in ms for e in coupling(m)
                  if e isa DeferredCounterEdge && e.direction === :in &&
                     !(is_registered(e.species) && is_chemostatted(e.species))]

@testset "The smoothed drain (spec §11 task 14c.4)" begin
    clamped = corea_models()
    smoothed = corea_models(smoothing = COREA_SMOOTHING_WIDTH)

    @testset "every consumer counter is smoothed, at the recorded width" begin
        @test COREA_SMOOTHING_WIDTH == 1.0
        # The width's reason is recorded where it is defined.
        doc = string(@doc COREA_SMOOTHING_WIDTH)
        @test occursin("round trip", doc) && occursin("ours", doc)

        cs = _consumers(smoothed)
        # The seven counters on live pools of 14b's census, plus the two
        # CTP/UTP chemostat debits, which are excluded above.
        @test Set(e.counter for e in cs) ==
              Set([:ATP_mRNA, :ATP_mRNAdeg, :ATP_transloc, :ATP_trsc, :GTP_mRNA,
                   :GTP_translat, :tRNA_translat])
        @test all(e -> e.clip === :smoothed && e.smoothing == 1.0, cs)
        @test all(e -> e.clip === :clamped_deficit_carried, _consumers(clamped))
        # The default is the published clamp, and nothing else moved.
        @test length(_consumers(clamped)) == length(cs)
    end

    @testset "no edge obstructs gradients, and the published model's do" begin
        @test !any(obstructs_gradients, _consumers(smoothed))
        @test isempty(resolve_coupling(smoothed).gradient_obstructions)
        @test Set(r.species for r in resolve_coupling(clamped).gradient_obstructions) ==
              Set(e.species for e in _consumers(clamped))
        @test isempty(check_gradient_safety(smoothed))
    end

    @testset "the smoothing is labelled ours (T2)" begin
        labels = reduction_declarations(smoothed)
        sm = [l for l in labels if l.category === :smoothed_counter]
        @test Set(l.subject for l in sm) == Set(e.species for e in _consumers(smoothed))
        @test all(l -> occursin("smoothed clip of width 1.0", l.description), sm)
        @test !any(l -> l.category === :smoothed_counter, reduction_declarations(clamped))
        d = build_corea(; smoothing = COREA_SMOOTHING_WIDTH, tspan = (0.0, 10.0))
        @test occursin("smoothed clip of width 1.0", reduction_report(smoothed, d))
    end
end
