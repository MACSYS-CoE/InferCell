using Test
using InferCell
using StaticArrays
using Random

# Spec §11 phase 4 — the 60 s rebuild: live ODE pools become the stochastic
# block's rate constants, held piecewise-constant between refreshes.
#
# Doubles are in `hybrid_test_models.jl`, included first by `runtests.jl`.
# `caught` comes from `corea_test_models.jl`.

# How many handshakes changed a given slot of the jump parameter vector. This is
# what makes a refresh count a *measurement* rather than a restatement of the
# schedule that produced it: it is read off the recorded trajectory, and it
# would still be right if the schedule were wrong.
function _slot_changes(jump_p, slot)
    n = 0
    for i in 2:length(jump_p)
        jump_p[i][slot] == jump_p[i - 1][slot] || (n += 1)
    end
    return n
end

# A metabolic half that *moves* without running dry. Phase 3's long-horizon
# tests freeze the pool (`kcat = 0.0`) so that any drift is the exchange's;
# phase 4 needs the opposite — an upstream pool that changes at every handshake,
# so that "the rate constant is unchanged between refreshes" is a statement
# about the coupling rather than about a constant pool. A slow enzyme and a
# large pool give both: ATP declines steadily and ends around two thirds full,
# so nothing clips and check 7's census stays clean.
_slow_pool(; kwargs...) = ToyPool(; kcat = 0.02, atp0 = 20.0, kwargs...)

# The global slot a named free parameter occupies in a block's parameter vector.
function _slot_of(models, name)
    seen = Symbol[]
    for m in models, q in InferCell.model_free_params(parameters(m))
        q.name in seen || push!(seen, q.name)
    end
    return findfirst(==(name), seen)
end

@testset "Phase 4 — the 60 s rebuild" begin

    @testset "4.1 a rebuilt constant lives in the parameter vector, not on the struct" begin
        m = ToyRebuiltExpression()
        d = build_problem([_slow_pool(), m]; tspan = (0.0, 600.0))
        @test d isa HandshakeDriver
        @test length(d.rebuilds) == 1
        r = only(d.rebuilds)
        @test r.names == [:k_tx_rb]
        @test r.species == [:M_atp_c]
        @test r.steps == 60

        # The point of the mechanism: the rebuilt slot is still a *named free
        # parameter*, so `remake(prob; p = θ)` and `model_free_params` reach it.
        # A constant held as mutable state on the sub-model would be invisible to
        # both, and so could never appear in a posterior, a provenance table or
        # an identifiability Jacobian.
        slot = _slot_of([m], :k_tx_rb)
        @test slot !== nothing
        @test r.fill_idxs == [slot]
        @test d.jump.p[slot] == m.params[1].value
        @test :k_tx_rb in [q.name for q in InferCell.model_free_params(parameters(m))]
    end

    @testset "4.1 a composition declaring no rebuild is untouched" begin
        d = build_problem([ToyPool(), ToyExpression()]; tspan = (0.0, 100.0))
        @test isempty(d.rebuilds)
        @test d.drain_interval == d.interval
        @test d.steps_per_drain == 1
        @test isempty(driver_declarations(d))

        # In trajectory form: with nothing rebuilding, the jump block's
        # parameter vector is bitwise constant for the whole run.
        Random.seed!(4001)
        out = run_handshake!(d, 100)
        @test all(p -> p == out.jump_p[1], out.jump_p)
        @test isempty(out.rebuilds)

        # And the protocol's default really is empty, so a module that rebuilds
        # nothing needs no method.
        @test isempty(rebuilt_params(ToyPool()))
        @test length(rate_constants(Float64[], 0.0, ToyPool(), SA[1.0])) == 0
    end

    @testset "4.2 the refresh count is the edge's interval, read from the trajectory" begin
        Random.seed!(4002)
        m60 = ToyRebuiltExpression(interval = 60.0)
        d60 = build_problem([_slow_pool(), m60]; tspan = (0.0, 600.0))
        slot = only(d60.rebuilds).fill_idxs[1]
        out60 = run_handshake!(d60, 600)
        # 600 s at 60 s is ten refreshes. Counted from the recorded parameter
        # vector, not from the schedule: the first refresh at t = 60 changes the
        # slot away from its declared nominal value, so ten changes.
        @test _slot_changes(out60.jump_p, slot) == 10
        @test only(out60.rebuilds).refreshes == 10
        @test only(out60.rebuilds).interval == 60.0

        Random.seed!(4002)
        m30 = ToyRebuiltExpression(interval = 30.0)
        d30 = build_problem([_slow_pool(), m30]; tspan = (0.0, 600.0))
        slot30 = only(d30.rebuilds).fill_idxs[1]
        out30 = run_handshake!(d30, 600)
        @test _slot_changes(out30.jump_p, slot30) == 20
        @test only(out30.rebuilds).refreshes == 20
    end

    @testset "4.2 the value written is the module's own law at the live pool" begin
        Random.seed!(4003)
        m = ToyRebuiltExpression(interval = 60.0)
        d = build_problem([_slow_pool(), m]; tspan = (0.0, 120.0))
        slot = only(d.rebuilds).fill_idxs[1]
        out = run_handshake!(d, 60)
        # The rebuild reads the pools the hook has just left behind, so the
        # value at the refresh handshake is the law evaluated at that same
        # handshake's recorded ATP concentration — exactly, not approximately.
        atp = out.ode[60][1]
        @test out.jump_p[60][slot] == toy_rebuilt_k(m, [atp])
        @test out.jump_p[60][slot] != m.params[1].value
    end

    @testset "4.3 the coupling is piecewise-constant, not accidentally continuous" begin
        Random.seed!(4004)
        m = ToyRebuiltExpression(interval = 60.0)
        d = build_problem([_slow_pool(), m]; tspan = (0.0, 180.0))
        slot = only(d.rebuilds).fill_idxs[1]
        out = run_handshake!(d, 180)

        # The upstream pool moves at (very nearly) every handshake …
        @test length(unique(o[1] for o in out.ode)) > 150
        # … and the rate constant does not: bitwise unchanged between refreshes.
        for window in (61:119, 121:179)
            @test all(p -> p[slot] === out.jump_p[first(window)][slot],
                      out.jump_p[window])
        end
        # It does change at the boundaries, so the assertion above is not
        # vacuously true of a channel that never fired.
        @test out.jump_p[120][slot] !== out.jump_p[119][slot]

        # Propensities themselves move with the jump state, so the statement
        # "the propensities change only at interval boundaries" is made at a
        # *fixed reference state*: the module's own rate function there is
        # bitwise unchanged across the same handshakes.
        rate = first(reactions(m)).rate
        ref_u, ref_pidx = SA[3, 7, 0], only(d.rebuilds).p_idxs
        props = [rate(ref_u, view(p, ref_pidx), 0.0, nothing) for p in out.jump_p]
        @test all(x -> x === props[61], props[61:119])
        @test props[120] !== props[119]
    end

    @testset "4.4 a continuous cadence is refused, by name" begin
        bad = ToyRebuiltExpression(cadence = :continuous)
        err = caught(() -> build_problem([ToyPool(), bad]; tspan = (0.0, 60.0)))
        @test err isa ArgumentError
        msg = sprint(showerror, err)
        @test occursin("ToyRebuiltExpression", msg)
        @test occursin("M_atp_c", msg)
        @test occursin(":continuous", msg)
        @test occursin("piecewise_constant", msg)

        # The label is unchanged where the edge is declared but not executed:
        # refusing to run it does not stop it being a declared deviation.
        labels = reduction_declarations([bad])
        @test :continuous_rebuild in [l.category for l in labels]
    end

    @testset "4.2 the refusals, each naming what is wrong" begin
        # A rebuilt name that is not one of the module's own free parameters.
        err = caught(() -> build_problem([ToyPool(),
                                          ToyRebuiltExpression(rebuilt = [:not_mine])];
                                         tspan = (0.0, 60.0)))
        @test err isa ArgumentError
        @test occursin("not_mine", sprint(showerror, err))

        # A name a second jump module also declares: parameters deduplicate by
        # name, so the rebuild would drive that module's propensities too.
        err = caught(() -> build_problem([ToyPool(), ToyRebuiltExpression(),
                                          ToyRebuildClash()]; tspan = (0.0, 60.0)))
        @test err isa ArgumentError
        msg = sprint(showerror, err)
        @test occursin("k_tx_rb", msg) && occursin("ToyRebuildClash", msg)

        # An inbound edge on a pool no ODE module integrates.
        err = caught(() -> build_problem([ToyPool(),
                                          ToyRebuiltExpression(pools = [:M_gtp_c])];
                                         tspan = (0.0, 60.0)))
        @test err isa ArgumentError
        @test occursin("M_gtp_c", sprint(showerror, err))

        # rebuilt_params with no edge to rebuild from.
        noedge = ToyRebuiltExpression(
            edges = [DeferredCounterEdge(species = :M_atp_c, direction = :in,
                                         counter = :atp_cost)])
        err = caught(() -> build_problem([ToyPool(), noedge]; tspan = (0.0, 60.0)))
        @test err isa ArgumentError
        @test occursin("no inbound RateConstantEdge", sprint(showerror, err))

        # An edge with nothing to rebuild: declared and never executed.
        err = caught(() -> build_problem([ToyPool(),
                                          ToyRebuiltExpression(rebuilt = Symbol[])];
                                         tspan = (0.0, 60.0)))
        @test err isa ArgumentError
        @test occursin("never", sprint(showerror, err))

        # An ODE module declaring rebuilt_params: the channel rebuilds the
        # stochastic block's constants, and an ODE module reads a pool directly.
        err = caught(() -> build_problem([ToyPool(), ToyPoolRebuilder(),
                                          ToyExpression()]; tspan = (0.0, 60.0)))
        @test err isa ArgumentError
        @test occursin("ToyPoolRebuilder", sprint(showerror, err))

        # An interval that is not a whole number of handshakes is unreachable,
        # and rounding it silently would run a cadence nobody declared.
        err = caught(() -> build_problem([ToyPool(),
                                          ToyRebuiltExpression(interval = 1.5)];
                                         tspan = (0.0, 60.0)))
        @test err isa ArgumentError
        @test occursin("whole number", sprint(showerror, err))

        # Two intervals on one module: one module rebuilds on one schedule.
        two = ToyRebuiltExpression(
            edges = [DeferredCounterEdge(species = :M_atp_c, direction = :in,
                                         counter = :atp_cost),
                     RateConstantEdge(species = :M_atp_c, direction = :in,
                                      interval = 60.0),
                     RateConstantEdge(species = :M_adp_c, direction = :in,
                                      interval = 30.0)])
        err = caught(() -> build_problem([ToyPool(), two]; tspan = (0.0, 60.0)))
        @test err isa ArgumentError
        @test occursin("different intervals", sprint(showerror, err))

        # The mirror trap: the pool's owner names a channel nothing consumes.
        err = caught(() -> build_problem([ToyPool(), ToyPoolRebuildPeer(),
                                          ToyRebuiltExpression()];
                                         tspan = (0.0, 60.0)))
        @test err isa ArgumentError
        msg = sprint(showerror, err)
        @test occursin("ToyPoolRebuildPeer", msg) && occursin("M_gtp_c", msg)
    end

    @testset "4.2 several pools feed one rebuild, and each is read" begin
        Random.seed!(4005)
        m = ToyRebuiltExpression(pools = [:M_atp_c, :M_adp_c], interval = 60.0)
        d = build_problem([_slow_pool(), m]; tspan = (0.0, 60.0))
        r = only(d.rebuilds)
        @test r.species == [:M_atp_c, :M_adp_c]
        @test length(r.pool_idxs) == 2
        out = run_handshake!(d, 60)
        @test out.jump_p[60][r.fill_idxs[1]] ==
              toy_rebuilt_k(m, [out.ode[60][1], out.ode[60][2]])
    end

    @testset "4.5 the channel's gain, against its closed form" begin
        # The rebuild law is Michaelis-shaped in each pool, so the elasticity
        # d ln k / d ln pool has the closed form km / (km + pool) and the
        # finite-difference diagnostic is checked rather than merely reported.
        for km in (1.0, 0.19)
            m = ToyRebuiltExpression(km_tx = km)
            d = build_problem([ToyPool(), m]; tspan = (0.0, 60.0))
            rows = rate_constant_elasticity(d)
            row = only(rows)
            @test row.module_id === :ToyRebuiltExpression
            @test row.param === :k_tx_rb
            @test row.pool === :M_atp_c
            @test isapprox(row.elasticity, toy_rebuilt_elasticity(m, row.concentration);
                           rtol = 1e-6)
        end

        # A gain in the band the real transcription channel is measured at
        # (0.044–0.051, `dev/notes/reduced-syn3a-scoping.md`), so the diagnostic
        # is exercised where phase 10 will read it.
        m = ToyRebuiltExpression(km_tx = 0.175)
        d = build_problem([ToyPool(), m]; tspan = (0.0, 60.0))
        eps = only(rate_constant_elasticity(d)).elasticity
        @test 0.044 <= eps <= 0.051
    end

    @testset "4.6 the drain interval is a declared, labelled policy" begin
        d = build_problem([ToyPool(), ToyExpression()];
                          tspan = (0.0, 600.0), drain_interval = 60.0)
        @test d.drain_interval == 60.0
        @test d.steps_per_drain == 60

        labels = driver_declarations(d)
        @test :coarse_drain in [l.category for l in labels]

        # A non-carry rounding policy is a departure too, and phase 3's own
        # docstring said so while nothing labelled it.
        dd = build_problem([ToyPool(), ToyExpression()];
                           tspan = (0.0, 60.0), rounding = :deterministic)
        @test :rounding_policy in [l.category for l in driver_declarations(dd)]

        # A drain between handshakes is unreachable.
        err = caught(() -> build_problem([ToyPool(), ToyExpression()];
                                         tspan = (0.0, 60.0), drain_interval = 1.5))
        @test err isa ArgumentError
        @test occursin("whole number", sprint(showerror, err))
    end

    @testset "4.6 a coarser drain accrues and pays in one lump" begin
        Random.seed!(4006)
        d = build_problem([ToyPool(kcat = 0.0), ToyExpression()];
                          tspan = (0.0, 60.0), drain_interval = 10.0)
        out = run_handshake!(d, 20)
        # Between drains the counter carries the accrual rather than the pool
        # paying it, so the pool is bitwise unchanged and the counter is not.
        atp = [o[1] for o in out.ode]
        @test all(x -> x === atp[11], atp[11:19])
        @test atp[20] !== atp[19]
        # Twenty handshakes, two drains: the pool takes exactly two values after
        # its initial one, wherever the SSA happened to fire in between.
        @test length(unique(atp)) == 3
        # The counter carries the accrual between drains. What it holds at a
        # drain handshake is one second's worth, not none: the debit clears it
        # and the SSA step that follows immediately starts refilling it — the
        # same one-exchange lag phase 3's done-when test asserts on the protein
        # count.
        @test out.jump[19][3] > out.jump[20][3] > 0
    end

    @testset "Done when: all three executed channels run together" begin
        Random.seed!(4007)
        m = ToyRebuiltExpression(interval = 60.0)
        d = build_problem([_slow_pool(), m]; tspan = (0.0, 600.0))
        out = run_handshake!(d, 600)
        @test d.n_handshakes == 600
        @test only(out.rebuilds).refreshes == 10          # the reverse channel ran
        @test out.jump[600][2] > 0                        # protein, the catalytic count
        @test out.ode[600][2] > out.ode[1][2]             # ADP rose: the ODE ran
        @test out.census.clipped == 0                     # check 7 on this toy
        @test only(rate_constant_elasticity(d)).elasticity > 0
    end
end
