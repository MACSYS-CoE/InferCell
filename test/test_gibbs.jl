# spec/phases/16-recovery.md task 16a.9a: the three blocks composed into one
# Gibbs chain, on a three-cell, 180 s M0 dataset. The pilot that measures its
# mixing is dev/scripts/gibbs_pilot_16a9a.jl.

using Random
using Serialization
using Test
using InferCell
using Distributions: Normal, logpdf, params

# Replay `path` from `base`: the transcripts and panel latent at each window's
# end, and the path's log density.
function _gibbs_replay(g, base, path)
    d = restore(base)
    r = PathReplay(join_path(path, d.events.labels), d; density = true)
    T = size(g.bulk, 2)
    tx = zeros(Int, length(g.tmap.states), T)
    lat = zeros(length(g.rows), T)
    for m in 1:T
        for _ in 1:g.window
            replay_step!(d, r)
        end
        tx[:, m] = [d.jump.u[s] for s in g.tmap.states]
        lat[:, m] = [d.ode.u[i] * d.factor + d.rounding.remainders[i] for i in g.rows]
    end
    return tx, lat, path_logdensity(r)
end

# Everything a chain's draws depend on, for bitwise comparison.
_gibbs_key(s) = (s.cme, s.u, s.σ, [[(w.times, w.reactions) for w in p] for p in s.paths], s.latents,
                 [(t.cme, t.u, t.σ, t.changed) for t in s.trace])

@testset "16a.9a: the Gibbs chain" begin
    H = 180.0
    ds = m0_dataset(16085; ncells = 3, horizon = H)
    g = gibbs_setup(ds; N = 4, horizon = H)
    tv = Dict(ds.truth.values)
    cme_true = [tv[n] for n in g.cme]
    u_true = [log(tv[f]) for f in g.forwards]
    @test g.cme == [:S_JCVISYN3A_0607, :S_JCVISYN3A_0779, :krnadeg] ||
          g.cme == [:S_JCVISYN3A_0779, :S_JCVISYN3A_0607, :krnadeg]

    # A path recorded at the truth, from each cell's base.
    bases = InferCell._gibbs_bases(g, cme_true, u_true; threads = false)
    T = size(g.bulk, 2)
    true_paths = map(eachindex(bases)) do c
        d = restore(bases[c])
        Random.seed!(ds.seeds[c])
        record_path!(d)
        for _ in 1:round(Int, H)
            handshake_step!(d)
        end
        window_events(recorded_path(d), T)
    end

    @testset "the bases reproduce the dataset" begin
        for c in eachindex(bases)
            tx, lat, _ = _gibbs_replay(g, bases[c], true_paths[c])
            @test tx == g.transcripts[c]
            @test isapprox(lat, ds.latent[c, [findfirst(==(s), ds.species) for s in g.panel], :];
                           rtol = 1e-10)
        end
    end

    @testset "a base written at once replays as one written in steps" begin
        cme_a, u_a = cme_true .* [1.4, 0.7, 1.2], u_true .+ [0.3, -0.2]
        writes_b = gibbs_writes(g, cme_true, u_true)
        stepped = restore(g.pristine[1])
        set_parameters!(stepped, g.models, gibbs_writes(g, cme_a, u_a))
        set_parameters!(stepped, g.models, writes_b)
        once = restore(g.pristine[1])
        set_parameters!(once, g.models, writes_b)
        tx1, lat1, lp1 = _gibbs_replay(g, stepped, true_paths[1])
        tx2, lat2, lp2 = _gibbs_replay(g, once, true_paths[1])
        @test tx1 == tx2
        @test lat1 == lat2
        @test lp1 == lp2
    end

    # Block 1's replay over every cell, at a θ.
    function replay_all(cme, u)
        bs = InferCell._gibbs_bases(g, cme, u; threads = false)
        labels = first(bs).events.labels
        cells = [Block1Cell(bs[c], join_path(true_paths[c], labels), zeros(0, 0))
                 for c in eachindex(bs)]
        b = Block1(g.models, cells; forwards = g.forwards, panel = g.panel,
                   save_every = g.window, bulk = g.bulk)
        return block1_replay(b, u; threads = false)
    end

    @testset "block 2's conditional is exact on the hybrid" begin
        lp0, _, _, counts, exposure = replay_all(cme_true, u_true)
        for k in eachindex(g.cme)
            js = k < length(g.cme) ? [g.tmap.births[k]] : g.tmap.deaths
            n, A = sum(counts[js]), sum(exposure[js]) / cme_true[k]
            f = rate_log_conditional(n, A, g.cme_priors[k])
            μ, sd = params(g.cme_priors[k])
            for factor in (0.6, 1.7)
                cme1 = copy(cme_true)
                cme1[k] *= factor
                lp1 = first(replay_all(cme1, u_true))
                u0, u1 = log(cme_true[k]), log(cme1[k])
                expected = (lp1 - lp0) + logpdf(Normal(μ, sd), u1) - logpdf(Normal(μ, sd), u0)
                @test abs((f(u1) - f(u0)) - expected) <= 1e-10 * max(1.0, abs(expected))
            end
        end
    end

    @testset "threads change nothing in block 1's replay" begin
        bs = InferCell._gibbs_bases(g, cme_true, u_true; threads = false)
        labels = first(bs).events.labels
        cells = [Block1Cell(bs[c], join_path(true_paths[c], labels), zeros(0, 0))
                 for c in eachindex(bs)]
        b = Block1(g.models, cells; forwards = g.forwards, panel = g.panel,
                   save_every = g.window, bulk = g.bulk)
        @test block1_replay(b, u_true; threads = false) == block1_replay(b, u_true; threads = true)
    end

    start = (cme = cme_true .* [1.5, 0.6, 1.3], u = u_true .+ [0.4, -0.5], σ = 0.3)

    @testset "a chain starts on a path that meets the data" begin
        s = gibbs_init(Xoshiro(1609), g; start..., threads = false)
        @test s.sweep == 0
        for c in eachindex(s.paths)
            tx, lat, lp = _gibbs_replay(g, s.bases[c], s.paths[c])
            @test tx == g.transcripts[c]
            @test lat == s.latents[c]
            @test isfinite(lp)
        end
    end

    @testset "a sweep keeps the chain's invariants" begin
        s = gibbs_init(Xoshiro(1609), g; start..., threads = false)
        for _ in 1:2
            gibbs_sweep!(s, g; threads = false)
        end
        @test s.sweep == 2 && length(s.trace) == 2
        @test size(s.trace[end].changed) == (T, length(s.paths))
        for c in eachindex(s.paths)
            tx, lat, lp = _gibbs_replay(g, s.bases[c], s.paths[c])
            @test tx == g.transcripts[c]
            @test lat == s.latents[c]
            @test isfinite(lp)
        end
        ode = derived_ode_values(g.models, [f => exp(x) for (f, x) in zip(g.forwards, s.u)])
        @test all(r -> r.relative <= r.bound, assert_haldane(g.models, ode))
        @test all(>(0), s.cme) && s.σ > 0
        # Something moved: two sweeps from an overdispersed start.
        @test s.u != start.u && s.cme != start.cme
    end

    @testset "threads change nothing in the chain" begin
        a = gibbs_init(Xoshiro(1610), g; start..., threads = false)
        b = gibbs_init(Xoshiro(1610), g; start..., threads = true)
        for _ in 1:2
            gibbs_sweep!(a, g; threads = false)
            gibbs_sweep!(b, g; threads = true)
        end
        @info "16a.9a threads test ran on $(Threads.nthreads()) threads"
        @test _gibbs_key(a) == _gibbs_key(b)
    end

    @testset "a resumed chain continues as the uninterrupted one" begin
        a = gibbs_init(Xoshiro(1611), g; start..., threads = false)
        b = gibbs_init(Xoshiro(1611), g; start..., threads = false)
        for _ in 1:3
            gibbs_sweep!(a, g; threads = false)
        end
        for _ in 1:2
            gibbs_sweep!(b, g; threads = false)
        end
        file = tempname()
        serialize(file, gibbs_checkpoint(b))
        b = gibbs_resume(deserialize(file), g; threads = false)
        gibbs_sweep!(b, g; threads = false)
        @test _gibbs_key(a) == _gibbs_key(b)
    end
end
