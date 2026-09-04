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
        @test d.ode_ids == [:ToyPool] && d.jump_ids == [:ToyExpression]

        # The refusal the composition used to give is gone for a genuine
        # hybrid, but the two shipped models still cannot compose — they both
        # own :mRNA and :protein, one on each side of the boundary. The reason
        # changed; the refusal did not. Note that no existing check caught
        # this: `_check_state_ownership` only sees duplicates within a block,
        # and only on registry species, and :mRNA is neither.
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
        # a random walk says little. The tolerance is four standard errors of
        # the mean of |drift|, whose scale is sqrt(n_hooks * p(1-p)).
        rms(policy, n) = sqrt(mean(drift(policy, n; seed = s)^2 for s in 1:60))
        r630, r6300 = rms(:stochastic, 630), rms(:stochastic, 6300)
        @test r6300 / r630 ≈ sqrt(10) rtol = 0.35
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
        @test rec.ode[end] == u0_ode       # exactly, after 600 handshakes
        @test rec.jump[end] == u0_jump
        # And not merely at the end: no handshake moved anything.
        @test all(u -> u == u0_ode, rec.ode)
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

    @testset "3.7 check 7: the clipping census" begin
        # At nominal parameters the pool covers the cost many times over, so no
        # handshake may carry a deficit. K5 fires here or nowhere.
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
