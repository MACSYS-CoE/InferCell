using Test
using InferCell
using StaticArrays
using Random

# Spec §11 phase 13a — the framework fixes assembly needs. Tasks 13a.1 (task
# 13.9, the ownerless path for chemostatted pools) and 13a.2 (task 13.10,
# per-product stoichiometry on deferred counters), on the phase 3 and 4
# doubles in `hybrid_test_models.jl`. The Core A′ halves of both tasks are in
# the transcription, decay and translation test files, where the pinned throws
# they replace used to be.

_F13 = corea_particles_per_mM()
_mM13(n) = n / _F13

# A phosphate pool that does nothing between hooks, so every particle of change
# in it is the hook's.
struct FrozenPi <: AbstractSubModel
    params::Vector{InferParameter}
end
FrozenPi(; pi0 = 1000) = FrozenPi([_jic(_mM13(pi0), :M_pi_c0, :FrozenPi)])
InferCell.states(::FrozenPi) = [:M_pi_c]
InferCell.parameters(m::FrozenPi) = m.params
InferCell.formalism(::FrozenPi) = :ode
InferCell.coupling(::FrozenPi) = CouplingEdge[]
InferCell.dynamics(u, p, t, ::FrozenPi) = SA[0.0]

_dc(species, direction; stoichiometry = 1.0) =
    DeferredCounterEdge(species = species, direction = direction,
                        counter = :atp_cost, stoichiometry = stoichiometry)

# A frozen ATP/ADP pool, and a jump block whose only cross-boundary costs are
# the edges given.
_frozen(; atp = 5000, adp = 1000) =
    ToyPool(kcat = 0.0, atp0 = _mM13(atp), adp0 = _mM13(adp))
_costs(edges) = ToyExpression(k_tx = 0.0, k_tl = 0.0, edges = edges)

@testset "Phase 13a — chemostats and products" begin

    @testset "13a.1 a chemostat pays a debit with no owner behind it" begin
        d = build_problem([_frozen(), _costs([_dc(:M_ctp_c, :in)])];
                          tspan = (0.0, 10.0))
        b = only(d.debits)
        @test b.pool_idx == 0
        d.jump.u[b.counter_idx] = 250
        u0 = copy(d.ode.u)
        handshake_step!(d)

        # Paid in full, nothing carried, nothing clipped, no pool written.
        @test b.deficit == 0.0
        @test d.n_clipped == 0
        @test d.jump.u[b.counter_idx] == 0
        @test d.ode.u == u0
        row = only(chemostat_census(d))
        @test (row.species, row.sign, row.particles) == (:M_ctp_c, -1, 250.0)
    end

    @testset "13a.1 a chemostat's payment is what a product is credited" begin
        # CTP → PPi's shape, with ADP standing in for the product.
        d = build_problem([_frozen(), _costs([_dc(:M_ctp_c, :in), _dc(:M_adp_c, :out)])];
                          tspan = (0.0, 10.0))
        d.jump.u[d.counters[1].counter_idx] = 250
        handshake_step!(d)
        @test d.ode.u[2] == _mM13(1000 + 250)
    end

    @testset "13a.1 a chemostat absorbs a credit, matched to what was paid" begin
        # ATP holds 100 against an accrual of 250, so it clips; the chemostat
        # credited from it may take only the 100 that left.
        d = build_problem([_frozen(atp = 100), _costs([_dc(:M_atp_c, :in),
                                                       _dc(:M_utp_c, :out)])];
                          tspan = (0.0, 10.0))
        d.jump.u[d.counters[1].counter_idx] = 250
        handshake_step!(d)
        @test d.ode.u[1] == 0.0
        @test d.n_clipped == 1
        row = only(chemostat_census(d))
        @test (row.species, row.sign, row.particles) == (:M_utp_c, 1, 100.0)
    end

    @testset "13a.1 a dynamic pool with no owner is still refused" begin
        # The ownerless path is for chemostats only. GTP is dynamic, so a debit
        # on it with nobody integrating it is a missing module, as before.
        err = caught(() -> build_problem([_frozen(), _costs([_dc(:M_gtp_c, :in)])];
                                         tspan = (0.0, 10.0)))
        @test err isa ArgumentError
        @test occursin("no ODE module in this composition integrates",
                       sprint(showerror, err))
    end

    @testset "13a.1 a rebuild reads a chemostat's clamped value" begin
        pools = [:M_atp_c, :M_ctp_c]
        base = ToyRebuiltExpression(pools = pools)
        clamp = ClampedEdge(species = :M_ctp_c, direction = :in, origin = :ours,
                            held_value = 0.6874)
        m = ToyRebuiltExpression(pools = pools, edges = [base.edges; clamp])
        d = build_problem([toy_slow_pool(), m]; tspan = (0.0, 120.0))
        r = only(d.rebuilds)
        @test r.pool_idxs[2] == 0
        @test r.pool_held[2] == 0.6874
        Random.seed!(7)
        rec = run_handshake!(d, 60)
        # The refresh at 60 s read live ATP and the held CTP.
        k = rec.jump_p[end][r.fill_idxs[1]]
        @test k ≈ toy_rebuilt_k(m, [rec.ode[end][1], 0.6874]) rtol = 1e-14

        # With no clamp the registry imports no CTP value, so nothing says what
        # the pool is held at, and the rebuild is refused rather than guessed.
        err = caught(() -> build_problem([toy_slow_pool(), base]; tspan = (0.0, 120.0)))
        @test err isa ArgumentError
        @test occursin("nothing gives it a value", sprint(showerror, err))
    end

    @testset "13a.2 products are credited per unit paid, at their stoichiometry" begin
        # ATP → ADP + Pi, with ATP clipping: 100 paid of 250 accrued.
        d = build_problem([_frozen(atp = 100), FrozenPi(pi0 = 1000),
                           _costs([_dc(:M_atp_c, :in), _dc(:M_adp_c, :out),
                                   _dc(:M_pi_c, :out)])];
                          tspan = (0.0, 10.0))
        @test length(d.counters) == 1
        aden0 = (d.ode.u[1] + d.ode.u[2]) * _F13
        phos0 = (3d.ode.u[1] + 2d.ode.u[2] + d.ode.u[3]) * _F13
        d.jump.u[d.counters[1].counter_idx] = 250
        handshake_step!(d)
        @test d.n_clipped == 1
        @test d.ode.u[1] == 0.0
        @test d.ode.u[2] == _mM13(1000 + 100)
        @test d.ode.u[3] == _mM13(1000 + 100)
        # Both moieties balance across the clip: nothing minted from the 150
        # that ATP could not pay, which stays carried on the consumer.
        @test (d.ode.u[1] + d.ode.u[2]) * _F13 ≈ aden0 atol = 1e-6
        @test (3d.ode.u[1] + 2d.ode.u[2] + d.ode.u[3]) * _F13 ≈ phos0 atol = 1e-6
        @test only(b for b in d.debits if b.sign < 0).deficit == 150.0

        # A non-unit stoichiometry scales the payment, not the demand.
        d2 = build_problem([_frozen(atp = 100),
                            _costs([_dc(:M_atp_c, :in),
                                    _dc(:M_adp_c, :out; stoichiometry = 2.0)])];
                           tspan = (0.0, 10.0))
        d2.jump.u[d2.counters[1].counter_idx] = 250
        handshake_step!(d2)
        @test d2.ode.u[2] == _mM13(1000 + 200)
    end

    @testset "13a.2 two products with no consumer are refused" begin
        # Nothing was paid, so each would be credited the raw accrual: one event
        # counted twice rather than a reaction with two products.
        err = caught(() -> build_problem([_frozen(), FrozenPi(),
                                          _costs([_dc(:M_adp_c, :out),
                                                  _dc(:M_pi_c, :out)])];
                                         tspan = (0.0, 10.0)))
        @test err isa ArgumentError
        @test occursin("debits none", sprint(showerror, err))
    end

    @testset "13a.2 a debit is one unit per accrued unit" begin
        @test_throws ArgumentError DeferredCounterEdge(species = :M_atp_c,
            direction = :in, counter = :c, stoichiometry = 2.0)
        @test_throws ArgumentError DeferredCounterEdge(species = :M_adp_c,
            direction = :out, counter = :c, stoichiometry = 0.0)
        @test DeferredCounterEdge(species = :M_adp_c, direction = :out,
                                  counter = :c).stoichiometry == 1.0
    end
end
