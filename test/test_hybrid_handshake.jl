using Test
using InferCell
using StaticArrays
using Statistics: mean, std
using Random

# Spec §11 phase 3 — the 1 s handshake on a two-module toy. The kill phase.
#
# Doubles are in `hybrid_test_models.jl`, included first by `runtests.jl`.
# `caught` comes from `corea_test_models.jl`.

const _F = corea_particles_per_mM()          # particles per mM at the published radius

# A pool value that is an exact whole number of particles, so a round trip
# through the conversion has nothing to quantise on the first hook and task
# 3.4's assertion can be `==` rather than "unchanged after the first".
_exact_mM(n) = n / _F

@testset "Phase 3 — the 1 s handshake" begin

    @testset "3.2 a mixed composition returns a driver rather than throwing" begin
        d = build_problem([ToyPool(), ToyExpression()]; tspan = (0.0, 600.0))
        @test d isa HandshakeDriver
        @test d.interval == 1.0
        @test length(d.catalytic) == 1
        @test length(d.debits) == 1
        # Both blocks were built, and each block's own layout is its own.
        @test length(d.ode.u) == 2          # M_atp_c, M_adp_c
        @test length(d.jump.u) == 3         # toy_mrna, M_ptsg_c, atp_cost
        @test length(d.counters) == 1                 # one counter feeds one pool

        # The refusal the composition used to give is gone for a genuine
        # hybrid, but the two shipped models still cannot compose — they both
        # own :mRNA and :protein, one on each side of the boundary. The reason
        # changed; the refusal did not. Note that no existing check caught
        # this: `_check_state_ownership` does span both blocks — it
        # iterates every model regardless of formalism — but skips any name the
        # registry does not know, and :mRNA is one.
        err = caught(() -> build_problem([TranscriptionTranslation(),
                                          StochasticGeneExpression()]))
        @test err !== nothing
        msg = sprint(showerror, err)
        @test occursin("owned in both blocks", msg)
        @test occursin("mRNA", msg)
    end

    @testset "3.2 a missing exchange declaration throws, named" begin
        # A catalytic edge whose parameter slot is not one of the module's own.
        bad = ToyPool(edges = [CatalyticEdge(species = :M_ptsg_c, direction = :in,
                                             param_slot = :no_such_slot)])
        err = caught(() -> build_problem([bad, ToyExpression()]))
        @test err isa ArgumentError
        msg = sprint(showerror, err)
        @test occursin("no_such_slot", msg) && occursin("ToyPool", msg)

        # A catalytic edge on a species no jump module owns.
        bad2 = ToyPool(edges = [CatalyticEdge(species = :M_crr_c, direction = :in,
                                              param_slot = :enzyme_conc)])
        err2 = caught(() -> build_problem([bad2, ToyExpression()]))
        @test err2 isa ArgumentError
        @test occursin("M_crr_c", sprint(showerror, err2))

        # A deferred counter debiting a pool no ODE module integrates.
        bad3 = ToyExpression(edges = [DeferredCounterEdge(species = :M_gtp_c,
                                                          direction = :in,
                                                          counter = :atp_cost)])
        err3 = caught(() -> build_problem([ToyPool(), bad3]))
        @test err3 isa ArgumentError
        @test occursin("M_gtp_c", sprint(showerror, err3))

        # A deferred counter accruing into a state no jump module owns.
        bad4 = ToyExpression(edges = [DeferredCounterEdge(species = :M_atp_c,
                                                          direction = :in,
                                                          counter = :not_a_state)])
        err4 = caught(() -> build_problem([ToyPool(), bad4]))
        @test err4 isa ArgumentError
        @test occursin("not_a_state", sprint(showerror, err4))

        # A counter that is a registry species. The hook clears its counters
        # every handshake, so accruing into a modelled pool would destroy that
        # pool once a second — silently, since clearing looks like a debit.
        bad5 = ToyExpression(edges = [DeferredCounterEdge(species = :M_atp_c,
                                                          direction = :in,
                                                          counter = :M_ptsg_c)])
        err5 = caught(() -> build_problem([ToyPool(), bad5]))
        @test err5 isa ArgumentError
        @test occursin("M_ptsg_c", sprint(showerror, err5))
        @test occursin("registry", sprint(showerror, err5))

        # A catalytic slot whose name two ODE modules both declare. Parameters
        # dedup by name into one slot, so the count would drive both rate laws.
        err6 = caught(() -> build_problem([ToyPool(), ToyPool2(), ToyExpression()]))
        @test err6 isa ArgumentError
        @test occursin("enzyme_conc", sprint(showerror, err6))

        # A policy keyword on a homogeneous composition is refused rather than
        # silently ignored.
        err5 = caught(() -> build_problem([ToyPool()]; rounding = :deterministic))
        @test err5 !== nothing
        @test occursin("rounding", sprint(showerror, err5))

        # An :sde sub-model has no build path in any composition. It used to
        # fall through to the mixed-formalism refusal and be reported as a
        # mixed composition, which it is not (spec §2, G1).
        sde = CoreAStub(:SdeBlock; form = :sde, st = [:M_gtp_c])
        err6 = caught(() -> build_problem([sde]))
        @test err6 !== nothing
        @test occursin("sde", sprint(showerror, err6))
        err7 = caught(() -> build_problem([ToyPool(), sde]))
        @test err7 !== nothing
        @test occursin("sde", sprint(showerror, err7))
    end

    @testset "3.2 a named peer across the boundary is accepted, not refused" begin
        # The edge kinds invite naming the counterpart module. Each block is
        # built from its own modules, so a per-block re-resolution would see the
        # peer as absent and refuse a composition that is in fact correct — the
        # failure phase 4's real modules would hit on their first declaration.
        pool = ToyPool(edges = [CatalyticEdge(species = :M_ptsg_c, direction = :in,
                                              param_slot = :enzyme_conc,
                                              peer = :ToyExpression)])
        expr = ToyExpression(edges = [DeferredCounterEdge(species = :M_atp_c,
                                          direction = :in, counter = :atp_cost,
                                          peer = :ToyPool)])
        d = build_problem([pool, expr]; tspan = (0.0, 10.0))
        @test d isa HandshakeDriver
        @test length(d.catalytic) == 1 && length(d.debits) == 1
    end

    @testset "3.2 a cross-block inputs() is refused, and says why" begin
        # `inputs()` is wired within a block, from that block's own state
        # vector, so it cannot reach the other side. Left to the block builders
        # this fails as "not owned by any sub-model" — true of that block, and
        # no help at all in diagnosing a boundary crossing.
        reader = CoreAStub(:JumpReader; form = :jump, st = [:jr_counter],
                           ins = [:M_atp_c],
                           edges = [CurrencyEdge(species = :M_atp_c, direction = :in)])
        err = caught(() -> build_problem([ToyPool(), reader]; tspan = (0.0, 10.0)))
        @test err !== nothing
        msg = sprint(showerror, err)
        @test occursin("M_atp_c", msg) && occursin("JumpReader", msg)
        @test occursin("inputs()", msg) && occursin("DeferredCounterEdge", msg)
    end

    @testset "3.3 the conversion is derived, not transcribed" begin
        # The scoping note records 1 mM = 20,180 particles at a 200 nm radius,
        # and 3.35e-17 L. Both must fall out of Avogadro and a sphere.
        @test round(Int, corea_particles_per_mM()) == 20180
        @test cell_volume_litres(200.0) ≈ 3.351e-17 rtol = 1e-3
        # Phase 5 grows the radius, so the factor must depend on it.
        @test corea_particles_per_mM(400.0) ≈ 8 * corea_particles_per_mM(200.0)
    end

    @testset "3.3 check 0: the round trip is exactly zero under fractional carry" begin
        st = RoundingState(:fractional_carry; nspecies = 1)
        n0 = 73717
        x = _exact_mM(n0)
        stable = true
        for _ in 1:6300
            n = round_to_counts!(st, 1, x * _F)
            stable &= (n == n0)
            x = counts_to_mM(n, _F)
        end
        # One assertion, not 6,300: a per-iteration @test would add six
        # thousand cases to the suite count and say nothing more.
        @test stable
        @test x == _exact_mM(n0)          # exact, not approximate
    end

    @testset "3.3 check 0: the three policies have three signatures" begin
        # A pool losing a fractional number of particles per handshake is what
        # separates the policies: the fraction is either carried, gambled or
        # discarded. Deterministic rounding discards it the same way every
        # time, which is why its error is a bias and not a wobble.
        function drift(policy, n_hooks; start = 73717.0, per_hook = 10.4, seed = 7)
            Random.seed!(seed)
            st = RoundingState(policy; nspecies = 1)
            pool = start
            for _ in 1:n_hooks
                pool = float(round_to_counts!(st, 1, pool - per_hook))
            end
            return pool - (start - n_hooks * per_hook)
        end

        # Fractional carry: bounded by half a particle and *not growing*.
        @test abs(drift(:fractional_carry, 630)) <= 0.5
        @test abs(drift(:fractional_carry, 6300)) <= 0.5

        # Deterministic: linear. Ten times the handshakes, ten times the error.
        d630, d6300 = drift(:deterministic, 630), drift(:deterministic, 6300)
        @test d6300 / d630 ≈ 10.0 rtol = 1e-9
        @test abs(d6300) > 2000          # ~0.4 particle per hook, biased one way

        # Stochastic: square root. Averaged over seeds, since one realisation of
        # a random walk says little. The ratio of two 60-seed RMS estimates has
        # a sampling standard deviation of about 0.39 (measured over 2,000
        # independent blocks), i.e. ~12.4% of sqrt(10), so four standard errors
        # is rtol ≈ 0.50. An earlier 0.35 was called "four standard errors" and
        # was in fact 2.8, which would flake on roughly one RNG stream in 160.
        rms(policy, n) = sqrt(mean(drift(policy, n; seed = s)^2 for s in 1:60))
        r630, r6300 = rms(:stochastic, 630), rms(:stochastic, 6300)
        @test r6300 / r630 ≈ sqrt(10) rtol = 0.50
        @test r6300 < 200                # far below deterministic's ~2500
        @info "3.3 check 0 signatures" carry_630 = drift(:fractional_carry, 630) carry_6300 = drift(:fractional_carry, 6300) det_630 = d630 det_6300 = d6300 stoch_rms_630 = r630 stoch_rms_6300 = r6300
    end

    @testset "3.3 the rounding policy is a field on the driver" begin
        d = build_problem([ToyPool(), ToyExpression()]; tspan = (0.0, 10.0))
        @test d.rounding.policy === :fractional_carry        # the default
        d2 = build_problem([ToyPool(), ToyExpression()]; tspan = (0.0, 10.0),
                           rounding = :stochastic)
        @test d2.rounding.policy === :stochastic
        err = caught(() -> build_problem([ToyPool(), ToyExpression()];
                                         tspan = (0.0, 10.0), rounding = :nearest))
        @test err isa ArgumentError
        @test occursin("nearest", sprint(showerror, err))
    end

    @testset "3.3 the policy keywords reach the driver" begin
        # Five of the six documented keywords were never passed through
        # `build_problem`, so nothing pinned that they arrive where the docs say.
        # Phase 5 changes `radius_nm` at a call site, and phase 4 the interval.
        d = build_problem([ToyPool(), ToyExpression()]; tspan = (0.0, 10.0),
                          interval = 0.5, radius_nm = 400.0,
                          ode_solver = Tsit5(), abstol = 1e-8, reltol = 1e-6)
        @test d.interval == 0.5
        @test d.factor == corea_particles_per_mM(400.0)     # eight times the default
        handshake_step!(d)
        @test d.ode.t == 0.5                                 # the interval is honoured

        err = caught(() -> build_problem([ToyPool(), ToyExpression()];
                                         tspan = (0.0, 10.0), interval = 0.0))
        @test err isa ArgumentError
        @test occursin("interval", sprint(showerror, err))
    end

    @testset "3.4 a run past the declared tspan is refused rather than silently taken" begin
        d = build_problem([ToyPool(kcat = 0.0), ToyExpression()]; tspan = (0.0, 10.0))
        err = caught(() -> run_handshake!(d, 20))
        @test err !== nothing
        @test occursin("tspan", sprint(showerror, err))
    end

    @testset "3.4 exchange error is isolated from integration error" begin
        # Zero derivative in the ODE block, no reactions in the jump block: the
        # only thing that can move a number is the exchange itself.
        d = build_problem([ToyPool(kcat = 0.0, atp0 = _exact_mM(73717),
                                   adp0 = _exact_mM(4400)),
                           ToyExpression(k_tx = 0.0, k_tl = 0.0)];
                          tspan = (0.0, 600.0))
        u0_ode = collect(Float64, d.ode.u)
        u0_jump = collect(Int, d.jump.u)
        rec = run_handshake!(d, 600)
        @test rec.t[end] == 600.0
        # Not merely at the end: no handshake moved anything.
        @test all(u -> u == u0_ode, rec.ode)
        @test all(u -> u == u0_jump, rec.jump)
        # The carry never grows: it holds the sub-particle difference between a
        # pool and its own round trip, which is a float rounding error here and
        # not an accumulating one.
        @test all(r -> abs(r) < 1e-6, d.rounding.remainders)
    end

    @testset "3.5 the deferred debit clamps, carries the deficit, and repays it" begin
        pool0 = 100
        d = build_problem([ToyPool(kcat = 0.0, atp0 = _exact_mM(pool0)),
                           ToyExpression(k_tx = 0.0, k_tl = 0.0)];
                          tspan = (0.0, 10.0))
        counter = d.debits[1].counter_idx
        pool = d.debits[1].pool_idx
        @test d.debits[1].clip === :clamped_deficit_carried    # the published policy

        # A cost larger than the pool: the pool floors at zero and the shortfall
        # is carried exactly.
        d.jump.u[counter] = 250
        handshake_step!(d)
        @test d.ode.u[pool] == 0.0
        @test d.debits[1].deficit == 150.0
        @test d.n_clipped == 1

        # Refill the pool; the carried deficit is debited at the next hook.
        InferCell._set_ode_state!(d.ode, pool, _exact_mM(200))
        d.jump.u[counter] = 0
        handshake_step!(d)
        @test d.ode.u[pool] == _exact_mM(50)
        @test d.debits[1].deficit == 0.0
        @test d.n_clipped == 1                                  # no new clip
    end

    @testset "3.5 one counter feeding two pools debits both from the same accrual" begin
        # Spec §3's charged-tRNA transfer is one accrual debiting a pool and
        # crediting another. Reading the counter inside the per-pool loop makes
        # the second pool see zero, which conserves nothing and reports nothing.
        cost = 10
        d = build_problem([ToyPool(kcat = 0.0, atp0 = _exact_mM(5000),
                                   adp0 = _exact_mM(1000)),
                           ToyExpression(k_tx = 0.0, k_tl = 0.0, cost = cost,
                               edges = [DeferredCounterEdge(species = :M_atp_c,
                                            direction = :in, counter = :atp_cost),
                                        DeferredCounterEdge(species = :M_adp_c,
                                            direction = :out, counter = :atp_cost)])];
                          tspan = (0.0, 10.0))
        @test length(d.debits) == 2
        @test length(d.counters) == 1            # both debits share one accrual

        counter = d.counters[1].counter_idx
        d.jump.u[counter] = 250
        handshake_step!(d)

        # Every particle taken from ATP arrives in ADP: the moiety is conserved
        # across the pair, which is the whole point of a shared counter.
        @test d.ode.u[1] == _exact_mM(5000 - 250)
        @test d.ode.u[2] == _exact_mM(1000 + 250)
        @test d.jump.u[counter] == 0
    end

    @testset "3.5 a credit on a shared counter matches the debit, not the demand" begin
        # If the drawn pool clips, only what actually left it may arrive
        # anywhere else. Crediting the raw accrual would mint the shortfall —
        # and on the charged-tRNA transfer this shape exists for, it would do so
        # every time the charged pool ran dry.
        d = build_problem([ToyPool(kcat = 0.0, atp0 = _exact_mM(100),
                                   adp0 = _exact_mM(1000)),
                           ToyExpression(k_tx = 0.0, k_tl = 0.0,
                               edges = [DeferredCounterEdge(species = :M_atp_c,
                                            direction = :in, counter = :atp_cost),
                                        DeferredCounterEdge(species = :M_adp_c,
                                            direction = :out, counter = :atp_cost)])];
                          tspan = (0.0, 10.0))
        d.jump.u[d.counters[1].counter_idx] = 250
        handshake_step!(d)

        # 100 particles were all ATP had, so 100 is all ADP may receive.
        @test d.ode.u[1] == 0.0
        @test d.ode.u[2] == _exact_mM(1000 + 100)
        # The unpaid 150 is carried on the consumer, not conjured on the producer.
        consumer = only(b for b in d.debits if b.sign < 0)
        producer = only(b for b in d.debits if b.sign > 0)
        @test consumer.deficit == 150.0
        @test producer.deficit == 0.0
        @test d.n_clipped == 1
    end

    @testset "3.5 a channel the counter's owner never declares is refused" begin
        # Lowering from the accruing side means a declaration by the pool's
        # owner alone would be skipped. Silently: the counter would grow without
        # bound, the pool would never be debited, and the census would report
        # nothing wrong.
        pool_only = ToyPool(kcat = 0.0,
                            edges = [DeferredCounterEdge(species = :M_atp_c,
                                         direction = :out, counter = :atp_cost)])
        err = caught(() -> build_problem([pool_only,
                                          ToyExpression(edges = CouplingEdge[])];
                                         tspan = (0.0, 10.0)))
        @test err isa ArgumentError
        msg = sprint(showerror, err)
        @test occursin("atp_cost", msg) && occursin("ToyExpression", msg)

        # And a catalytic edge declared only by the side that owns the count,
        # which names no rate law and so would write nothing.
        jump_only = ToyExpression(edges = [CatalyticEdge(species = :M_ptsg_c,
                                               direction = :out,
                                               param_slot = :enzyme_conc)])
        err2 = caught(() -> build_problem([ToyPool(edges = CouplingEdge[]), jump_only];
                                          tspan = (0.0, 10.0)))
        @test err2 isa ArgumentError
        @test occursin("M_ptsg_c", sprint(showerror, err2))
    end

    @testset "3.5 two pools credited from one counter is refused" begin
        # Consumers share what a counter holds; producers would each be credited
        # the whole of it. One accrual, two credits, matter from nothing — the
        # defect this phase already shipped once and had to fix.
        two_credits = ToyExpression(k_tx = 0.0, k_tl = 0.0,
            edges = [DeferredCounterEdge(species = :M_atp_c, direction = :out,
                                         counter = :atp_cost),
                     DeferredCounterEdge(species = :M_adp_c, direction = :out,
                                         counter = :atp_cost)])
        err = caught(() -> build_problem([ToyPool(kcat = 0.0), two_credits];
                                         tspan = (0.0, 10.0)))
        @test err isa ArgumentError
        @test occursin("atp_cost", sprint(showerror, err))
    end

    @testset "3.5 the sign of a debit is the counter owner's, not the vector order's" begin
        # Both blocks may name one channel from their own side, with opposite
        # directions. If the sign came from whichever edge was seen first, the
        # same two models in the other order would credit instead of debit.
        function final_atp(order)
            expr = ToyExpression(k_tx = 0.0, k_tl = 0.0,
                       edges = [DeferredCounterEdge(species = :M_atp_c,
                                    direction = :in, counter = :atp_cost)])
            # The pool's owner names the same channel from its own side.
            pool = ToyPool(kcat = 0.0, atp0 = _exact_mM(5000),
                       edges = [CatalyticEdge(species = :M_ptsg_c, direction = :in,
                                              param_slot = :enzyme_conc),
                                DeferredCounterEdge(species = :M_atp_c,
                                    direction = :out, counter = :atp_cost)])
            d = build_problem(order ? [pool, expr] : [expr, pool]; tspan = (0.0, 10.0))
            d.jump.u[d.counters[1].counter_idx] = 250
            handshake_step!(d)
            return d.ode.u[d.debits[1].pool_idx]
        end
        # One channel, so one debit however many times it is declared.
        @test final_atp(true) == final_atp(false) == _exact_mM(5000 - 250)
    end

    @testset "3.5 the two labelled clip policies run, and the census judges them" begin
        function clip_run(clip; smoothing = nothing, pool = 100, accrued = 250)
            kw = smoothing === nothing ? (;) : (; smoothing = smoothing)
            d = build_problem([ToyPool(kcat = 0.0, atp0 = _exact_mM(pool)),
                               ToyExpression(k_tx = 0.0, k_tl = 0.0,
                                   edges = [DeferredCounterEdge(; species = :M_atp_c,
                                       direction = :in, counter = :atp_cost,
                                       clip = clip, kw...)])]; tspan = (0.0, 10.0))
            d.jump.u[d.counters[1].counter_idx] = accrued
            handshake_step!(d)
            return d
        end

        # :unclamped — a labelled departure: the pool is allowed to go negative.
        u = clip_run(:unclamped)
        @test u.ode.u[1] == _exact_mM(100 - 250)
        @test u.ode.u[1] < 0
        @test u.debits[1].deficit == 0.0
        @test u.n_clipped == 0                    # nothing was withheld

        # :smoothed — the other labelled departure, and the one whose whole
        # purpose is to be a *differentiable* approximation of the clamp. Two
        # things are worth pinning, and neither is "the pool stays positive":
        # that is true only while the smoothing width is comparable to the
        # shortfall.
        #
        # First, a cost far inside a deep pool must NOT be scored as clipping.
        # The softplus leaves a strictly positive residue at every step however
        # deep the pool, so a bare `deficit > 0` test would report 100% clipping
        # and make the census that scores K5 useless exactly where :smoothed is
        # adopted to escape K5.
        deep = clip_run(:smoothed; smoothing = 1.0, pool = 73717, accrued = 100)
        @test deep.ode.u[1] > 0
        @test deep.n_clipped == 0

        # Second, the approximation converges to the policy it approximates. At
        # a width far below the shortfall the softplus is the hard max to well
        # inside a particle, so :smoothed reproduces the clamped answer exactly
        # — pool on the floor, deficit the full shortfall.
        narrow = clip_run(:smoothed; smoothing = 1.0)
        @test narrow.ode.u[1] == 0.0
        @test narrow.debits[1].deficit ≈ 150.0
        @test narrow.n_clipped == 1

        # At a width comparable to the shortfall it visibly departs: the pool is
        # held off the floor, at the price of withholding more than the hard
        # shortfall. That gap is the deviation the label is for.
        wide = clip_run(:smoothed; smoothing = 200.0)
        @test wide.ode.u[1] > 0
        @test wide.debits[1].deficit > 150.0
        @test wide.n_clipped == 1
    end

    @testset "3.5 the clamped counter is reported as a gradient obstruction" begin
        obstructions = check_gradient_safety([ToyPool(), ToyExpression()])
        @test length(obstructions) == 1
        @test obstructions[1].species === :M_atp_c
        @test obstructions[1].kind === :deferred_counter
        graph = resolve_coupling([ToyPool(), ToyExpression()])
        @test occursin("M_atp_c", gradient_report(graph))
    end

    @testset "3.6 the enzyme-concentration channel scales the flux by the count ratio" begin
        function after_one_handshake(protein)
            d = build_problem([ToyPool(), ToyExpression(k_tx = 0.0, k_tl = 0.0,
                                                        protein0 = protein)];
                              tspan = (0.0, 10.0))
            handshake_step!(d)
            return collect(Float64, d.ode.p)
        end
        n = 400
        p1, p2 = after_one_handshake(n), after_one_handshake(2n)
        slot = 3                                   # :enzyme_conc, third free parameter

        # The count arrives as a concentration, exactly.
        @test p1[slot] == counts_to_mM(n, _F)
        @test p2[slot] == counts_to_mM(2n, _F)

        # The flux scales by the ratio of counts, at a common state.
        u = SA[3.6529, 0.2178]
        @test toy_flux(p2, u) / toy_flux(p1, u) == 2.0
        # And nothing else moved: every other parameter is untouched.
        @test p1[[1, 2]] == p2[[1, 2]]
        @test length(p1) == length(p2) == 3
    end

    @testset "3.4 the metabolic block really integrates its own rate law" begin
        # Every other assertion in this file would still pass if `step!` moved
        # the clock without solving anything: the debit arithmetic is the
        # driver's, the catalytic write is a parameter, and task 3.4's toy has
        # a zero derivative by construction. So check the one thing they do not
        # — that a handshake's ODE step reproduces an independent solve of the
        # same problem over the same second.
        n = 400
        d = build_problem([ToyPool(kcat = 30.0, atp0 = _exact_mM(73717)),
                           ToyExpression(k_tx = 0.0, k_tl = 0.0, protein0 = n)];
                          tspan = (0.0, 10.0))
        u0 = collect(Float64, d.ode.u)
        handshake_step!(d)
        @test d.ode.t == 1.0                       # the clock advanced by the interval

        # The same problem, solved on its own for one second at the enzyme
        # concentration the handshake wrote.
        ref = build_problem([ToyPool(kcat = 30.0, enzyme = counts_to_mM(n, _F),
                                     atp0 = _exact_mM(73717))]; tspan = (0.0, 1.0))
        sol = solve(ref, Rodas5P(); abstol = 1e-10, reltol = 1e-8)

        # ADP is not named by any deferred counter, so nothing quantises it and
        # it must match the reference to solver tolerance.
        @test d.ode.u[2] ≈ sol.u[end][2] rtol = 1e-6

        # ATP *is* the debited pool, so the hook writes it back through whole
        # particles and it is quantised by construction — the two cannot agree
        # to solver tolerance, and an assertion that they do would be asserting
        # the rounding policy is inert. The right bound is the conversion's own
        # resolution: one particle is 1/20180.4 mM, and the step agrees with an
        # independent solve to better than that.
        one_particle = 1.0 / _F
        @test abs(d.ode.u[1] - sol.u[end][1]) < one_particle
        @test d.ode.u[1] * _F ≈ round(d.ode.u[1] * _F) atol = 1e-6   # it is whole
        @info "3.4 integration against an independent solve" atp_handshake = d.ode.u[1] atp_reference = sol.u[end][1] gap_in_particles = abs(d.ode.u[1] - sol.u[end][1]) * _F

        # And it actually moved: a no-op step would leave the pool untouched.
        @test d.ode.u[1] < u0[1]
        @test d.ode.u[2] > u0[2]
    end

    @testset "3.7 check 7: the clipping census" begin
        # At nominal parameters the pool covers the cost many times over, so no
        # handshake may carry a deficit. K5 itself is scored on the assembled
        # model in phase 14, not here — see §12, 2026-09-05.
        Random.seed!(3001)
        d = build_problem([ToyPool(kcat = 0.0), ToyExpression()]; tspan = (0.0, 600.0))
        rec = run_handshake!(d, 600)
        @test rec.census.handshakes == 600
        @test rec.census.clipped == 0
        @test all(iszero, rec.census.deficits)

        # And across 200 draws from the modules' own priors, the fraction that
        # clips at any handshake. Above 5% the non-smooth drain is in the
        # operating regime and K5 fires (spec §3 check 7, §8 K5).
        n_draws, horizon = 200, 60
        clipped = 0
        Random.seed!(3002)
        for i in 1:n_draws
            k_tx = rand(LogNormal(log(2.0), 0.5))
            k_tl = rand(LogNormal(log(1.0), 0.5))
            gamma_m = rand(LogNormal(log(0.5), 0.5))
            di = build_problem([ToyPool(kcat = 0.0),
                                ToyExpression(k_tx = k_tx, k_tl = k_tl,
                                              gamma_m = gamma_m)];
                               tspan = (0.0, Float64(horizon)))
            run_handshake!(di, horizon)
            di.n_clipped > 0 && (clipped += 1)
        end
        fraction = clipped / n_draws
        @info "3.7 clipping census" nominal_clipped = rec.census.clipped prior_draws = n_draws prior_clipped = clipped prior_fraction = fraction
        @test fraction <= 0.05
    end

    @testset "Done when: 600 s of handshake, both channels live" begin
        Random.seed!(4242)
        atp0_particles = 73717
        cost = 10
        d = build_problem([ToyPool(kcat = 0.0, atp0 = _exact_mM(atp0_particles)),
                           ToyExpression(cost = cost)]; tspan = (0.0, 600.0))
        rec = run_handshake!(d, 600)

        protein = d.jump.u[2]
        counter = d.jump.u[3]
        @test protein > 0                                  # translation ran

        # The protein count sets the ODE rate law — as of the last exchange,
        # not as of now. The count the rate law carries is the one read at the
        # start of the final handshake; the translations that fired during the
        # stochastic step after it are not visible to the ODE block until the
        # next exchange. That lag *is* the published piecewise-constant
        # coupling, so it is asserted rather than tolerated.
        @test d.ode.p[3] == counts_to_mM(rec.jump[end - 1][2], _F)
        @test protein >= rec.jump[end - 1][2]

        # The pool paid the deferred cost, and paid exactly it. Every particle
        # debited was accrued by a translation event, and every accrued particle
        # was either debited or is still sitting in the counter.
        debited = atp0_particles - d.ode.u[1] * _F
        accrued = cost * protein - counter
        @test abs(debited - accrued) < 1.0                 # within one rounding
        @test rec.census.clipped == 0

        @info "phase 3 done-when" handshakes = rec.census.handshakes protein = protein debited_particles = debited accrued_particles = accrued final_atp_mM = d.ode.u[1]
    end
end
