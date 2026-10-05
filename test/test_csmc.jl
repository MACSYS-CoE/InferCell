# spec/phases/16-recovery.md task 16a.7: the conditional SMC path update's
# pieces, its streams, snapshot independence (V9), and a smoke sweep on M0.
# V5, against brute force, is dev/scripts/csmc_v5_16a7.jl.

using Random
using Test
using InferCell
using Statistics
using Distributions: Normal, logpdf

# The toy of phase 4: one gene whose transcription constant is rebuilt from ATP.
function _toy_csmc(; horizon = 180.0, k_tx = nothing)
    m = ToyRebuiltExpression()
    k0 = k_tx === nothing ? toy_rebuilt_k(m, [toy_slow_pool().params[4].value]) : k_tx
    ms = [toy_slow_pool(), ToyRebuiltExpression(k_tx = k0)]
    d = build_problem(ms; tspan = (0.0, horizon))
    names = InferCell._block_names(ms, :jump)
    lab(s) = findfirst(==(s), d.events.labels)
    tm = TranscriptMap([findfirst(==(:toy_mrna), names)], [lab(:ToyRebuiltExpression_1)],
                       [lab(:ToyRebuiltExpression_2)])
    iatp = findfirst(==(:M_atp_c), InferCell._block_names(ms, :ode))
    return ms, d, tm, iatp
end

# Record a toy cell and observe it: transcripts exactly, ATP lognormally.
function _toy_data(d, tm, iatp, T, σ, seed)
    Random.seed!(seed)
    record_path!(d)
    tx = zeros(Int, 1, T)
    met = zeros(1, T)
    rng = Xoshiro(seed + 1)
    for m in 1:T
        for _ in 1:60
            handshake_step!(d)
        end
        tx[1, m] = d.jump.u[tm.states[1]]
        x = d.ode.u[iatp] * d.factor + d.rounding.remainders[iatp]
        met[1, m] = max(x, 1.0) * exp(σ * randn(rng))
    end
    return CellData(tx, met), recorded_path(d)
end

@testset "16a.7: the conditional SMC path update" begin

    @testset "the weight's pieces match matrix exponentials, with rates 2× apart" begin
        ms, d, tm, iatp = _toy_csmc()
        spec = CSMCSpec(tm, [iatp]; sigma = 0.1, cap_floor = bridge_cap(4))
        slot = only(only(d.rebuilds).fill_idxs)
        γ = 0.5
        u = copy(d.jump.u)
        u[tm.states[1]] = 2
        for k in (0.8, 1.6), (a, b, τ) in ((2, 4, 60.0), (3, 2, 1.0), (2, 3, 1.0))
            d.jump.p[slot] = k
            kk, μ = InferCell._birth_death(d, tm, 1, u)
            @test kk == k
            @test μ ≈ γ rtol = 1e-14
            cap = adequate_cap(k, γ, τ, a; floor = spec.cap_floor)
            exact = log(exp(InferCell.generator(TranscriptChain(k, γ, cap)) * τ)[a + 1, b + 1])
            @test abs(InferCell._log_transitions(d, spec, u, [a], [b], τ) - exact) <= 1e-12
        end
        @test window_log_weight(1.0, 2.0, 0.5, -3.0) == 1.0 + 2.0 - 0.5 - 3.0
    end

    T, σ = 3, 0.1
    ms, d, tm, iatp = _toy_csmc()
    base = snapshot(d)
    data, path = _toy_data(d, tm, iatp, T, σ, 16071)
    ref = window_events(path, T)
    spec = CSMCSpec(tm, [iatp]; sigma = σ, cap_floor = bridge_cap(maximum(data.transcripts)))

    @testset "a reference particle replays its window bitwise" begin
        e = restore(base)
        lw, ev = advance_window!(Xoshiro(1), e, spec, [base.jump.u[tm.states[1]]],
                                 data.transcripts[:, 1], data.metabolites[:, 1]; fixed = ref[1])
        f = restore(base)
        replay!(f, path, 60)
        @test e.ode.u == f.ode.u && e.jump.u == f.jump.u && e.rounding.remainders == f.rounding.remainders
        @test ev.times == ref[1].times && ev.reactions == ref[1].reactions
        @test isfinite(lw)
    end

    @testset "per-particle streams: alike seeds agree, different seeds differ" begin
        y0 = [base.jump.u[tm.states[1]]]
        run(seed) = advance_window!(Xoshiro(seed), restore(base), spec, y0,
                                    data.transcripts[:, 1], data.metabolites[:, 1])
        a, b, c = run(7), run(7), run(8)
        @test a[1] == b[1] && a[2].times == b[2].times && a[2].reactions == b[2].reactions
        @test a[2].times != c[2].times
        # A proposal hits the observed count exactly.
        g = restore(base)
        advance_window!(Xoshiro(9), g, spec, y0, data.transcripts[:, 1], data.metabolites[:, 1])
        @test g.jump.u[tm.states[1]] == data.transcripts[1, 1]
    end

    @testset "V9: particles share no state" begin
        x, y = restore(base), restore(base)
        x.debits[1].deficit = 5.0
        x.jump.u .+= 1
        x.rounding.remainders .+= 0.25
        InferCell._set_ode_state!(x.ode, iatp, 2 * x.ode.u[iatp])
        @test y.debits[1].deficit == base.debits[1].deficit
        @test y.jump.u == base.jump.u && y.rounding.remainders == base.rounding.remainders
        @test y.ode.u == base.ode.u
        @test x.ode !== y.ode && x.jump !== y.jump
    end

    @testset "a sweep returns a path that reproduces the data" begin
        for lag in (0, 2)
            new, changed = csmc_sweep(Xoshiro(16072 + lag), base, spec, data, ref; N = 4, lag)
            @test length(new) == T && length(changed) == T
            h = restore(base)
            r = PathReplay(join_path(new, h.events.labels), h)
            for m in 1:T
                for _ in 1:60
                    replay_step!(h, r)
                end
                @test h.jump.u[tm.states[1]] == data.transcripts[1, m]
            end
        end
    end

    @testset "annealing (§12 2026-10-02): the schedule" begin
        β = geometric_schedule(5; β1 = 1e-3)
        @test length(β) == 5 && β[1] ≈ 1e-3 && β[end] == 1.0 && issorted(β)
        @test geometric_schedule(1) == [1.0]
        y0 = [base.jump.u[tm.states[1]]]
        args = (Xoshiro(1), restore(base), spec, y0, data.transcripts[:, 1], data.metabolites[:, 1])
        @test_throws ArgumentError annealed_window(args..., [0.5, 0.4, 1.0])
        @test_throws ArgumentError annealed_window(args..., [0.5, 0.9])
        @test_throws ArgumentError csmc_sweep(Xoshiro(1), base, spec, data, ref; N = 2, lag = 1,
                                              schedule = [1.0])
    end

    @testset "annealing: a head keeps the path before its cut and redraws the rest" begin
        y0 = [base.jump.u[tm.states[1]]]
        s = 25.0
        for seed in 1:5
            e = restore(base)
            lw, ev = advance_window!(Xoshiro(seed), e, spec, y0, data.transcripts[:, 1],
                                     data.metabolites[:, 1]; head = (ref[1], s))
            k = searchsortedfirst(ref[1].times, s)
            @test ev.times[1:(k - 1)] == ref[1].times[1:(k - 1)]
            @test ev.reactions[1:(k - 1)] == ref[1].reactions[1:(k - 1)]
            @test all(>=(s), ev.times[k:end])
            @test isfinite(lw) && e.jump.u[tm.states[1]] == data.transcripts[1, 1]
            # Its weight is the weight of its own path, replayed.
            lw2, _ = advance_window!(Xoshiro(0), restore(base), spec, y0, data.transcripts[:, 1],
                                     data.metabolites[:, 1]; fixed = ev)
            @test lw == lw2
        end
        @test_throws ArgumentError advance_window!(Xoshiro(1), restore(base), spec, y0,
                                                   data.transcripts[:, 1], data.metabolites[:, 1];
                                                   head = (ref[1], 59.5))
    end

    @testset "annealing: the reference keeps its path and its replayed state" begin
        y0 = [base.jump.u[tm.states[1]]]
        lw, ev, e = annealed_window(Xoshiro(4), base, spec, y0, data.transcripts[:, 1],
                                    data.metabolites[:, 1], geometric_schedule(4); fixed = ref[1])
        f = restore(base)
        replay!(f, path, 60)
        @test ev === ref[1] && isfinite(lw)
        @test e.ode.u == f.ode.u && e.jump.u == f.jump.u
        # The start is not changed.
        @test base.jump.t == 0.0
    end

    @testset "annealing: the AIS weight is unbiased for ∫ q·w (3 SE)" begin
        # Window 1 of the toy at σ = 0.3, where plain importance sampling is
        # precise enough to compare against.
        sp = CSMCSpec(tm, [iatp]; sigma = 0.3, cap_floor = spec.cap_floor)
        y0 = [base.jump.u[tm.states[1]]]
        n = 300
        plain = [advance_window!(Xoshiro(20_000 + i), restore(base), sp, y0, data.transcripts[:, 1],
                                 data.metabolites[:, 1])[1] for i in 1:n]
        ais = [annealed_window(Xoshiro(30_000 + i), base, sp, y0, data.transcripts[:, 1],
                               data.metabolites[:, 1], geometric_schedule(3; β1 = 0.1))[1] for i in 1:n]
        c = max(maximum(plain), maximum(ais))
        a, b = exp.(plain .- c), exp.(ais .- c)
        se = sqrt(var(a) / n + var(b) / n)
        @test abs(mean(a) - mean(b)) <= 3se
    end

    @testset "annealing: a sweep returns a path that reproduces the data" begin
        new, changed = csmc_sweep(Xoshiro(16074), base, spec, data, ref; N = 3,
                                  schedule = geometric_schedule(3))
        @test length(new) == T && length(changed) == T
        h = restore(base)
        r = PathReplay(join_path(new, h.events.labels), h)
        for m in 1:T
            for _ in 1:60
                replay_step!(h, r)
            end
            @test h.jump.u[tm.states[1]] == data.transcripts[1, m]
        end
    end

    @testset "bulk (§12 2026-10-05): the weight scores the population mean with this particle in it" begin
        y0 = [base.jump.u[tm.states[1]]]
        others = [1234.5, 2000.0]            # two windows' other-cell sums of ATP
        z = [700.0, 650.0]
        bulk = bulk_window_observations(reshape(z, 1, 2), reshape(others, 1, 2), 3)
        e = restore(base)
        lw_b, _ = advance_window!(Xoshiro(1), e, spec, y0, data.transcripts[:, 1], bulk[1];
                                  fixed = ref[1])
        f = restore(base)
        lw_p, _ = advance_window!(Xoshiro(1), f, spec, y0, data.transcripts[:, 1],
                                  data.metabolites[:, 1]; fixed = ref[1])
        x = f.ode.u[iatp] * f.factor + f.rounding.remainders[iatp]
        met_p = logpdf(Normal(log(max(x, 1.0)), σ), log(data.metabolites[1, 1]))
        met_b = logpdf(Normal(log(max((others[1] + x) / 3, 1.0)), σ), log(z[1]))
        @test lw_b - lw_p ≈ met_b - met_p rtol = 1e-12
        # A sweep returns the chosen path's latent, which a replay reproduces.
        zz = reshape([700.0, 650.0, 600.0], 1, 3)
        oo = reshape([1234.5, 2000.0, 1500.0], 1, 3)
        new, changed, lat = csmc_sweep(Xoshiro(16075), base, spec, data, ref; N = 3,
                                       bulk = bulk_window_observations(zz, oo, 3))
        h = restore(base)
        r = PathReplay(join_path(new, h.events.labels), h)
        for m in 1:T
            for _ in 1:60
                replay_step!(h, r)
            end
            @test lat[1, m] == h.ode.u[iatp] * h.factor + h.rounding.remainders[iatp]
        end
        @test_throws ArgumentError csmc_sweep(Xoshiro(1), base, spec, data, ref; N = 2, lag = 1,
                                              bulk = bulk_window_observations(zz, oo, 3))
    end

    @testset "M0: transcript map and one window" begin
        mm = m0_models()
        dm = build_m0(tspan = (0.0, 120.0))
        tmm = transcript_map(dm, mm)
        @test length(tmm.states) == 2
        @test [dm.events.labels[k] for k in tmm.births] == [:CoreATranscription_1, :CoreATranscription_2]
        rows = [findfirst(==(s), InferCell._block_names(mm, :ode)) for s in M0_PANEL]
        Random.seed!(16073)
        dd = restore(dm)
        record_path!(dd)
        for _ in 1:60
            handshake_step!(dd)
        end
        y1 = [dd.jump.u[s] for s in tmm.states]
        ymet = [dd.ode.u[i] * dd.factor + dd.rounding.remainders[i] for i in rows]
        sp = CSMCSpec(tmm, rows; sigma = 0.1, cap_floor = bridge_cap(maximum(y1)))
        e = restore(dm)
        lw, ev = advance_window!(Xoshiro(3), e, sp, [dm.jump.u[s] for s in tmm.states], y1, ymet)
        @test isfinite(lw) && [e.jump.u[s] for s in tmm.states] == y1
    end
end
