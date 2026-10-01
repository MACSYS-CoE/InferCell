# spec/phases/16-recovery.md task 16a.7: the conditional SMC path update's
# pieces, its streams, snapshot independence (V9), and a smoke sweep on M0.
# V5, against brute force, is dev/scripts/csmc_v5_16a7.jl.

using Random
using Test
using InferCell
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
        x.ode.u[iatp] *= 2
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
