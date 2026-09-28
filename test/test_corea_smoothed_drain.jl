# Spec §11 phase 14c: the smoothed drain.
#
# Task 14c.4: every consumer counter declares `clip = :smoothed` with its width,
# the width's choice is recorded with its reason, `reduction_report` carries the
# smoothing as ours (T2), and `obstructs_gradients` reports no edge.

using Test
using InferCell
using Random

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

# Task 14c.5, the two halves that are tests. The agreement of the candidate
# observables over 14b's seeds, and the same continuity on the assembled model,
# are `dev/scripts/corea_smoothing_14c.jl`'s, whose result file is the record.

# One handshake of the toy with transcription and translation off, so the debit
# is a fixed `accrued` and the ATP pool before it is a smooth function of
# `kcat` through the ODE step. The drain is then the only kink. Returns the pool
# after the debit, in particles, exactly (the carried remainder included).
function _toy_pool_after(kcat; clip, accrued = 400.0, n0 = 1000, protein = 500)
    kw = clip === :smoothed ? (; smoothing = COREA_SMOOTHING_WIDTH) : (;)
    f = corea_particles_per_mM()
    d = build_problem([ToyPool(kcat = kcat, atp0 = n0 / f),
                       ToyExpression(k_tx = 0.0, k_tl = 0.0, protein0 = protein,
                           edges = [DeferredCounterEdge(; species = :M_atp_c,
                               direction = :in, counter = :atp_cost, clip, kw...)])];
                      tspan = (0.0, 10.0), abstol = 1e-14, reltol = 1e-12)
    d.jump.u[d.counters[1].counter_idx] = accrued
    handshake_step!(d)
    return d.ode.u[1] * d.factor + d.rounding.remainders[1], d.debits[1].clipped
end

@testset "The smoothed drain against the clamped one (spec §11 task 14c.5)" begin
    @testset "the derivative is continuous across a clip, and the clamp's jumps" begin
        # The clip: the kcat at which the pool before the debit equals the
        # accrual, bisected on the clamped debit's own flag.
        lo, hi = 10.0, 100.0
        @test !last(_toy_pool_after(lo; clip = :clamped_deficit_carried))
        @test last(_toy_pool_after(hi; clip = :clamped_deficit_carried))
        for _ in 1:60
            mid = (lo + hi) / 2
            last(_toy_pool_after(mid; clip = :clamped_deficit_carried)) ? (hi = mid) : (lo = mid)
        end
        θc = (lo + hi) / 2

        δ = 1e-4
        slope(θ, clip) = (first(_toy_pool_after(θ + δ; clip)) -
                          first(_toy_pool_after(θ - δ; clip))) / 2δ
        # The pool's own sensitivity, well away from the clip, sets the scale.
        dP = slope(θc - 1.0, :clamped_deficit_carried)
        @test dP < 0
        jump(k, clip) = abs(slope(θc + k * δ, clip) - slope(θc - k * δ, clip)) / abs(dP)
        for k in (30, 3)
            # The clamp: one side pays in full and moves with the pool, the other
            # floors it at zero, so the slope jumps by the whole sensitivity.
            @test jump(k, :clamped_deficit_carried) > 0.9
        end
        # The smoothed drain: the slope difference across the clip shrinks with
        # the spacing, which is what a continuous derivative does.
        far, near = jump(30, :smoothed), jump(3, :smoothed)
        @test far < 0.05
        @test near < far / 5
        @info "14c.5 toy: slope jump across the clip, as a fraction of the pool's sensitivity" θc far near clamped = jump(3, :clamped_deficit_carried)
    end

    @testset "14a's closures hold under smoothing, at their own gates" begin
        # validation_models is the metered assembly; smoothing all three
        # stochastic modules is corea_models(smoothing = w) with the meters.
        drain = (clip = :smoothed, smoothing = COREA_SMOOTHING_WIDTH)
        ms = validation_models(transcription = CoreATranscription(; drain...),
                               decay = CoreATranscriptDecay(; drain...),
                               translation = CoreATranslation(; drain...))
        @test typeof.(ms) == typeof.(corea_models(metered = true,
                                                  smoothing = COREA_SMOOTHING_WIDTH))
        n = round(Int, COREA_CYCLE_S)
        function vrun(abstol, reltol)
            Random.seed!(1410)
            d = build_problem(ms; tspan = (0.0, Float64(n)), complete = true, abstol, reltol)
            Random.seed!(1410)
            return validation_run!(ms, d, n)
        end
        run, tight = vrun(1e-10, 1e-8), vrun(1e-11, 1e-9)
        maxres(r, name) = maximum(abs, closure_residual(r, name))
        # The same gate as 14a's tolerance-principle testset: every closure at
        # least three orders below its tol_C at both rungs, and flat.
        for m in run.moieties
            for r in (run, tight)
                @test maxres(r, m.name) < 1e-3 * moiety_bound(r, m.name)
            end
            a, b = maxres(run, m.name), maxres(tight, m.name)
            @test max(a, b) <= 100 * max(min(a, b), 1e-12)
        end
        @info "14c.5: smoothed closures, max residual / tol_C at the pinned pair" [
            m.name => maxres(run, m.name) / moiety_bound(run, m.name) for m in run.moieties]
    end
end
