# spec/phases/16-recovery.md task 16a.2, check V1: the path log-density a replay
# accumulates. Two hand checks on the phase 3 and phase 4 toys, the three
# mutations the second must catch, and the importance identity at suite size.
# The identity at 10⁴ paths is dev/scripts/path_density_16a2.jl.

using Random
using Statistics
using Test
using InferCell

_label_index(d, l) = findfirst(==(l), d.events.labels)

@testset "16a.2: the path log-density (V1)" begin

    @testset "hand check: one gene, frozen rates" begin
        m = ToyExpression()                       # k_tx 2, gamma 0.5, k_tl 1
        d = build_problem([ToyPool(kcat = 0.0), m]; tspan = (0.0, 5.0))
        tx, dec, tl = (_label_index(d, Symbol(:ToyExpression_, j)) for j in 1:3)
        path = JumpPath([0.3, 1.7, 2.2, 3.4, 4.1], [tx, tx, tl, dec, tl], d.events.labels)
        r = replay!(d, path, 5; density = true)
        k, γ, kl = 2.0, 0.5, 1.0
        hand = log(k) + log(k) + log(kl * 2) + log(γ * 2) + log(kl * 1) -
               (0.3 * k + 1.4 * (k + 1.5) + 1.7 * (k + 3.0) + 1.6 * (k + 1.5))
        @test abs(path_logdensity(r) - hand) <= 1e-12
        @test r.counts[tx] == 2 && r.counts[dec] == 1 && r.counts[tl] == 2
        @test r.exposure[tx] ≈ 5k atol = 1e-12
    end

    @testset "hand check: rates 2× apart across a rebuild at 60 − 1 s" begin
        m = ToyRebuiltExpression()                # k_tx rebuilt every 60 s from ATP
        d = build_problem([toy_slow_pool(), m]; tspan = (0.0, 120.0))
        tx, dec = _label_index(d, :ToyRebuiltExpression_1), _label_index(d, :ToyRebuiltExpression_2)
        path = JumpPath([30.0, 58.5, 59.5, 70.2], [tx, tx, tx, dec], d.events.labels)
        slot = only(only(d.rebuilds).fill_idxs)
        iatp = findfirst(==(:M_atp_c), InferCell._block_names([toy_slow_pool(), m], :ode))
        r = PathReplay(path, d; density = true)
        for _ in 1:59
            replay_step!(d, r)
        end
        k_old = d.jump.p[slot]
        # Drop the pool so the rebuild at handshake 60 moves the constant by far
        # more than 2×.
        InferCell._set_ode_state!(d.ode, iatp, 0.1)
        replay_step!(d, r)
        k_new = d.jump.p[slot]
        @test max(k_old / k_new, k_new / k_old) >= 2
        for _ in 61:120
            replay_step!(d, r)
        end
        k_end = d.jump.p[slot]                    # rebuilt at handshake 120, on [119, 120)
        γ, kl = 0.5, 1.0
        occupancy = 1 * 28.5 + 2 * 1.0 + 3 * 10.7 + 2 * 49.8   # ∫ m dt over [0, 120)
        logs = log(k_old) + log(k_old) + log(k_new) + log(γ * 3)
        hand = logs - (k_old * 59 + k_new * 60 + k_end * 1 + (γ + kl) * occupancy)
        @test abs(path_logdensity(r) - hand) <= 1e-12 * max(1, abs(hand))

        # The three mutations, each a different number by far more than roundoff.
        no_exposure = logs
        previous = log(k_old) + log(k_old) + log(k_old) + log(γ * 3) -
                   (k_old * 60 + k_new * 59 + k_end * 1 + (γ + kl) * occupancy)
        at_60m = log(k_old) + log(k_old) + log(k_old) + log(γ * 3) -
                 (k_old * 60 + k_new * 60 + (γ + kl) * occupancy)
        for mutant in (no_exposure, previous, at_60m)
            @test abs(mutant - hand) > 1e-3
        end
    end

    @testset "importance identity, suite size" begin
        # E_{X ~ p(·|θ0)}[p(X|θ1) / p(X|θ0)] = 1, for a CME-only shift and an
        # ODE-only shift. 10⁴ paths are the Slurm run; 400 here.
        # The declared k_tx must equal the rebuild law at the initial pool.
        # set_parameters! recomputes a module's rebuilt constants from the
        # initial pools when any of its parameters changes, so with the default
        # declared 2.0 against the law's 5.71 a decay shift would also move
        # transcription over [0, 59), and the two densities would describe two
        # different models. Core A′'s constructors compute their initial
        # constants from the same law.
        k0 = toy_rebuilt_k(ToyRebuiltExpression(), [toy_slow_pool().params[4].value])
        models() = [toy_slow_pool(), ToyRebuiltExpression(k_tx = k0)]
        H = 180.0
        function logp(path, values)
            ms = models()
            d = build_problem(ms; tspan = (0.0, H))
            isempty(values) || set_parameters!(d, ms, values)
            return path_logdensity(replay!(d, path, round(Int, H); density = true))
        end
        rng = Xoshiro(1602)
        for shift in ([:gamma_m_rb => 0.5 * 1.1], [:kcat_toy => 0.02 * 2.0])
            w = map(1:400) do i
                Random.seed!(rand(rng, UInt32))
                ms = models()
                d = build_problem(ms; tspan = (0.0, H))
                record_path!(d)
                for _ in 1:round(Int, H)
                    handshake_step!(d)
                end
                path = recorded_path(d)
                exp(logp(path, shift) - logp(path, Pair{Symbol, Float64}[]))
            end
            se = std(w) / sqrt(length(w))
            @test abs(mean(w) - 1) <= 3se
            @test std(w) > 0                       # the shift did move the density
        end
    end
end
