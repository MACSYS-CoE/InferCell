# KIM-28: predict's posterior ingest. One file of each kind is written in the
# layout its producer writes and read back (dev/notes/2026-10-29-archived-posteriors.md).

using Test
using InferCell
using InferCell: _log_names
using Random
using Serialization
using MCMCChains: Chains

const M0_NAMES = [:S_JCVISYN3A_0607, :S_JCVISYN3A_0779, :krnadeg, :kcatF_R_ENO, :kcatF_R_FBA]

@testset "KIM-28: posterior ingest" begin
    rng = Xoshiro(28)
    dir = mktempdir()

    @testset "a V7 replicate file" begin
        # sbc_v7_16a9b.jl's layout: the CME parameters natural, the forward
        # constants as u = ln θ, σ natural, in meta.names order.
        nsweep, burn = 12, 4
        trace = [(sweep = k, cme = rand(rng, 3) .+ 0.5, u = randn(rng, 2), σ = 0.2 + rand(rng),
                  changed = falses(2, 3), seconds = (0.1, 0.1, 0.1, 0.0)) for k in 1:nsweep]
        truth = randn(rng, 6)
        meta = (; replicate = 7, seed = 16_200_007, names = vcat(M0_NAMES, :sigma_b), truth,
                start = nothing, init_seconds = 1.0, ncells = 50, N = 20, threads = 1)
        path = joinpath(dir, "replicate7.jls")
        serialize(path, (; state = (; trace, sweep = nsweep), meta))

        pd = posterior_draws(path; burn)
        @test pd isa PosteriorDraws
        @test pd.names == M0_NAMES
        @test length(pd) == nsweep - burn
        @test pd.values[1, :] ≈ vcat(trace[burn + 1].cme, exp.(trace[burn + 1].u))
        @test pd.values[end, 4:5] ≈ exp.(trace[end].u)
        @test pd.sigma_b == [t.σ for t in trace[(burn + 1):end]]
        @test pd.iteration == collect((burn + 1):nsweep)
        @test pd.truth[:krnadeg] ≈ exp(truth[3]) && pd.truth[:sigma_b] ≈ exp(truth[6])
        @test pd.source.kind === :v7_replicate && pd.source.replicate == 7
        @test !(:sigma_b in pd.names)
        @test draw_values(pd, 2) == [n => pd.values[2, j] for (j, n) in enumerate(M0_NAMES)]
        # V7's default burn-in is its own 100 sweeps, more than this file holds.
        @test_throws ArgumentError posterior_draws(path)
        thinned = posterior_draws(path; burn, thin = 3)
        @test thinned.iteration == collect((burn + 1):3:nsweep)
    end

    @testset "a B7 chain file, with its fit's prep" begin
        # av6_89_b7.jl's layout: x[1]…x[6] on the log scale in the mixed route's
        # order, NUTS internals beside them, and the dataset in prep.jls.
        nit = 5
        X = randn(rng, nit, 6)
        lp = randn(rng, nit)
        ch = Chains(reshape(hcat(X, lp), nit, 7, 1),
                    vcat([Symbol("x[$i]") for i in 1:6], :lp), (internals = [:lp],))
        fit = joinpath(dir, "r3")
        mkpath(fit)
        serialize(joinpath(fit, "chain_2.jls"),
                  (; fit = "r3", c = 2, seed = 9, start = nothing, start_tries = 1, chain = ch,
                   seconds = 1.0, nadapt = 0, ndraw = nit, commit = "abc1234", job = "1"))
        tv = vcat([n => 0.1 * i for (i, n) in enumerate(M0_NAMES)], :kcatR_R_ENO => 9.0)
        serialize(joinpath(fit, "prep.jls"),
                  (; fit = "r3", ds = (; truth = (; values = tv), sigma = 0.3)))

        pd = posterior_draws(joinpath(fit, "chain_2.jls"))
        @test pd.names == M0_NAMES
        @test pd.values ≈ exp.(X[:, 1:5])
        @test pd.sigma_b ≈ exp.(X[:, 6])
        @test all(==(2), pd.chain)
        @test pd.truth == Dict(vcat(tv[1:5], :sigma_b => 0.3))
        @test !haskey(pd.truth, :kcatR_R_ENO)   # reverse constants are re-derived, never read
        @test pd.source.kind === :b7_chain && pd.source.commit == "abc1234"

        # Two chain files of one fit merge; a different truth is refused.
        serialize(joinpath(fit, "chain_1.jls"),
                  (; fit = "r3", c = 1, seed = 8, start = nothing, start_tries = 1, chain = ch,
                   seconds = 1.0, nadapt = 0, ndraw = nit, commit = "abc1234", job = "1"))
        both = posterior_draws(joinpath.(fit, ["chain_1.jls", "chain_2.jls"]))
        @test length(both) == 2nit && both.chain == vcat(fill(1, nit), fill(2, nit))
        @test both.source.kind === :merged && both.truth == pd.truth
        other = PosteriorDraws(pd.names, pd.values, pd.sigma_b, pd.chain, pd.iteration,
                               nothing, Dict(:krnadeg => 1.0), pd.source)
        @test_throws ArgumentError vcat(pd, other)
        # Another fit is refused even where neither records a truth, as most of
        # the B7 archive's fits do not.
        ks = joinpath(dir, "ks_low1")
        mkpath(ks)
        serialize(joinpath(ks, "chain_1.jls"),
                  (; fit = "ks_low1", c = 1, seed = 8, start = nothing, start_tries = 1,
                   chain = ch, seconds = 1.0, nadapt = 0, ndraw = nit, commit = "abc1234",
                   job = "1"))
        @test posterior_draws(joinpath(ks, "chain_1.jls")).truth === nothing
        @test_throws ArgumentError posterior_draws([joinpath(fit, "chain_1.jls"),
                                                    joinpath(ks, "chain_1.jls")])
        # burn and thin are V7's, and refused for a B7 file.
        @test_throws ArgumentError posterior_draws(joinpath(fit, "chain_1.jls"); thin = 2)
    end

    @testset "a reduced-model pilot's chain and extension files (KIM-94)" begin
        # corea_pilot_kim92.jl's layout: each stage's record under `result`;
        # x[1]…x[21] on the log scale, COREA_TARGETS then σ_b; a chain file holds
        # its adaptation then its kept draws, an extension its own draws numbered
        # from 1; the dataset, with the targets and truth, beside them.
        P = length(COREA_TARGETS) + 1
        nadapt, ndraw = 3, 4
        wrap(result, stage) = (; result, commit = "def5678", commit_full = "def5678" ^ 5,
                               job = "42", stage, smoke = false)
        mkch(X) = Chains(reshape(hcat(X, randn(rng, size(X, 1))), size(X, 1), P + 1, 1),
                         vcat([Symbol("x[$i]") for i in 1:P], :lp), (internals = [:lp],))
        fit = joinpath(dir, "pilot")
        mkpath(fit)
        tv = vcat([n => 0.5 + i for (i, n) in enumerate(COREA_TARGETS)],
                  :kcatR_R_ENO => 9.0, :kcatR_R_FBA => 8.0)
        serialize(joinpath(fit, "dataset.jls"),
                  wrap((; truth = (; values = tv), sigma = 0.35, targets = COREA_TARGETS),
                       "dataset"))
        Xc = Dict(c => randn(rng, nadapt + ndraw, P) for c in 1:2)
        Xe = Dict(c => randn(rng, ndraw, P) for c in 1:2)
        for c in 1:2
            serialize(joinpath(fit, "chain_$c.jls"),
                      wrap((; c, seed = 10 + c, start = nothing, start_tries = 1, chain = mkch(Xc[c]),
                            seconds = 1.0, eps0 = 0.1, metric = "warm", nadapt, ndraw), "chain"))
            serialize(joinpath(fit, "extend_$c.jls"),
                      wrap((; c, seed = 10 + c, chain = mkch(Xe[c]), seconds = 1.0, ndraw), "extend"))
        end

        ch1 = posterior_draws(joinpath(fit, "chain_1.jls"))
        @test ch1.names == COREA_TARGETS
        # The adaptation is dropped; the kept draws are exp(x) to the bit.
        @test length(ch1) == ndraw
        @test ch1.values == exp.(Xc[1][(nadapt + 1):end, 1:(P - 1)])
        @test ch1.sigma_b == exp.(Xc[1][(nadapt + 1):end, P])
        @test ch1.iteration == collect((nadapt + 1):(nadapt + ndraw)) && all(==(1), ch1.chain)
        # The truth is the dataset's targets and σ, never the derived reverse constants.
        @test ch1.truth == Dict(vcat(tv[1:(P - 1)], :sigma_b => 0.35))
        @test ch1.source.kind === :pilot_chain && ch1.source.stage == "chain" &&
              ch1.source.commit == "def5678"

        # An extension is numbered after its chain's run.
        ex1 = posterior_draws(joinpath(fit, "extend_1.jls"))
        @test ex1.values == exp.(Xe[1][:, 1:(P - 1)])
        @test ex1.iteration == collect((nadapt + ndraw + 1):(nadapt + 2ndraw))

        # The pooled posterior, in the pilot's own order: each chain's kept draws,
        # then its extension's, chain by chain (pilot_draws.tsv's rows).
        paths = [joinpath(fit, f) for c in 1:2 for f in ("chain_$c.jls", "extend_$c.jls")]
        pd = posterior_draws(paths)
        flat = reduce(vcat, (vcat(Xc[c][(nadapt + 1):end, :], Xe[c]) for c in 1:2))
        @test length(pd) == 2 * 2ndraw
        @test pd.values == exp.(flat[:, 1:(P - 1)]) && pd.sigma_b == exp.(flat[:, P])
        @test pd.chain == repeat(1:2; inner = 2ndraw)
        @test allunique(zip(pd.chain, pd.iteration))
        @test InferCell._route.(pd.source.parts) == fill(:mixed, 4)
        # A posterior over every target leaves nothing for completion to draw.
        hp = held_parameters(d11_models(), pd)
        @test isempty(hp.held)

        # The files beside it are needed: without the dataset there are no names,
        # and without its chain file an extension has no place in the run.
        lone = joinpath(dir, "lone")
        mkpath(lone)
        cp(joinpath(fit, "extend_1.jls"), joinpath(lone, "extend_1.jls"))
        @test_throws ArgumentError posterior_draws(joinpath(lone, "extend_1.jls"))
        cp(joinpath(fit, "dataset.jls"), joinpath(lone, "dataset.jls"))
        @test_throws ArgumentError posterior_draws(joinpath(lone, "extend_1.jls"))
        # Two pilot fits are two posteriors, even over the same targets.
        cp(joinpath(fit, "chain_1.jls"), joinpath(lone, "chain_1.jls"))
        @test_throws ArgumentError vcat(ch1, posterior_draws(joinpath(lone, "chain_1.jls")))
        # An extension beside another chain's file is refused, not misnumbered; so
        # is a chain file that does not hold its adaptation and kept draws.
        odd = joinpath(dir, "odd")
        mkpath(odd)
        for f in ("dataset.jls", "extend_1.jls")
            cp(joinpath(fit, f), joinpath(odd, f))
        end
        cp(joinpath(fit, "chain_2.jls"), joinpath(odd, "chain_1.jls"))
        @test_throws ArgumentError posterior_draws(joinpath(odd, "extend_1.jls"))
        serialize(joinpath(odd, "chain_3.jls"),
                  wrap((; c = 3, seed = 13, start = nothing, start_tries = 1,
                        chain = mkch(randn(rng, nadapt + ndraw - 1, P)), seconds = 1.0,
                        eps0 = 0.1, metric = "warm", nadapt, ndraw), "chain"))
        @test_throws ArgumentError posterior_draws(joinpath(odd, "chain_3.jls"))
    end

    @testset "infer's Chains, mapped by name" begin
        # `infer(models, ds; route)` writes log_<name> columns; the mixed route
        # keeps Turing's internals. Columns in another order map the same way.
        nit, nch = 4, 2
        X = randn(rng, nit, 6, nch)
        names = _log_names(M0_NAMES)
        ch = Chains(X, names; info = (route = :exact,))
        perm = [3, 6, 1, 5, 2, 4]
        shuffled = Chains(X[:, perm, :], names[perm])
        a, b = posterior_draws(ch), posterior_draws(shuffled)
        @test a.names == M0_NAMES
        @test b.values[:, indexin(a.names, b.names)] ≈ a.values
        @test b.sigma_b ≈ a.sigma_b ≈ exp.(vcat(X[:, 6, 1], X[:, 6, 2]))
        # Every chain is read, stacked draw-wise.
        @test length(a) == nit * nch && a.chain == vcat(fill(1, nit), fill(2, nit))
        @test a.values[nit + 1, 1] ≈ exp(X[1, 1, 2])
        @test a.source.route === :exact

        # A saved Chains file reads the same.
        path = joinpath(dir, "chains.jls")
        serialize(path, ch)
        @test posterior_draws(path).values == a.values
        @test_throws ArgumentError posterior_draws(path; burn = 1)
        # Two routes are two likelihoods, and one file read twice is not two chains.
        mixed = Chains(X, names; info = (route = :mixed,))
        @test_throws ArgumentError vcat(posterior_draws(ch), posterior_draws(mixed))
        @test_throws ArgumentError posterior_draws([path, path])

        # Separate runs each number their chains from 1; merged, each keeps a place
        # of its own, and the sources stay one flat list.
        single = Chains(X[:, :, 1:1], names)
        three = vcat(posterior_draws(single), posterior_draws(single), posterior_draws(single))
        @test three.chain == vcat(fill(1, nit), fill(2, nit), fill(3, nit))
        @test allunique(zip(three.chain, three.iteration))
        @test length(three.source.parts) == 3 && all(s -> s.kind === :chains, three.source.parts)

        # The toy models' NUTS path writes natural-scale names.
        nat = Chains(rand(rng, nit, 2, 1), [:k_tx, :sigma_obs])
        p = posterior_draws(nat)
        @test p.names == [:k_tx, :sigma_obs] && p.sigma_b === nothing
        @test p.values ≈ nat.value.data[:, :, 1]

        # x[i] columns need their names.
        xs = Chains(randn(rng, nit, 2, 1), [Symbol("x[1]"), Symbol("x[2]")])
        @test_throws ArgumentError posterior_draws(xs; names = [:log_a, :log_b, :log_c])
        @test posterior_draws(xs; names = [:log_a, :log_sigma_b]).names == [:a]
    end

    @testset "an ABC posterior keeps its weights" begin
        r = ABCPosterior(rand(rng, 2, 10), collect(1.0:10.0), [:a, :b], 0.1, 3)
        pd = posterior_draws(r)
        @test pd.values == permutedims(r.particles)
        @test pd.weights ≈ r.weights ./ sum(r.weights) && sum(pd.weights) ≈ 1
        @test_throws ArgumentError vcat(pd, pd)
    end

    @testset "refusals" begin
        bad = joinpath(dir, "bad.jls")
        serialize(bad, (; x = 1))
        @test_throws ArgumentError posterior_draws(bad)
        @test_throws ArgumentError PosteriorDraws([:a, :sigma_b], ones(2, 2), nothing,
                                                  [1, 1], [1, 2], nothing, nothing, (kind = :x,))
        @test_throws ArgumentError PosteriorDraws([:a], fill(NaN, 2, 1), nothing,
                                                  [1, 1], [1, 2], nothing, nothing, (kind = :x,))
    end

    # The archived posteriors themselves (KIM-18), where this machine has them.
    v7 = "/fred/oz022/tkimpson/infercell-16a9b-archive-2026-10-17/v7/sbc_v7_16a9b/replicate1.jls"
    b7 = "/fred/oz022/tkimpson/infercell-89-b7-archive-2026-10-27/run/r1/chain_1.jls"
    if get(ENV, "INFERCELL_INTEGRATION_TESTS", "false") == "true" && isfile(v7) && isfile(b7)
        @testset "the archived V7 and B7 posteriors" begin
            e, m = posterior_draws(v7), posterior_draws(b7)
            @test e.names == m.names == M0_NAMES
            @test length(e) == 556 && length(m) == 200
            @test keys(e.truth) == keys(m.truth)
            @test all(e.truth[k] ≈ m.truth[k] for k in keys(e.truth))
        end
    end

    # KIM-92's posterior, as M9 predicts from it (KIM-93's Go, 2026-11-17): the
    # pooled draws equal the pilot's own pilot_draws.tsv, written at %.17g, to the bit.
    pilot = "/fred/oz022/tkimpson/infercell-kim92-archive-2026-11-17"
    if get(ENV, "INFERCELL_INTEGRATION_TESTS", "false") == "true" &&
       isfile(joinpath(pilot, "pilot_draws.tsv"))
        @testset "the archived reduced-model pilot posterior" begin
            paths = [joinpath(pilot, f) for c in 1:4 for f in ("chain_$c.jls", "extend_$c.jls")]
            pd = posterior_draws(paths)
            lines = readlines(joinpath(pilot, "pilot_draws.tsv"))
            @test split(lines[1], '\t') == String.(_log_names(COREA_TARGETS))
            X = permutedims(reduce(hcat, [parse.(Float64, split(l, '\t')) for l in lines[2:end]]))
            @test pd.names == COREA_TARGETS && length(pd) == 1600 == size(X, 1)
            @test pd.values == exp.(X[:, 1:20]) && pd.sigma_b == exp.(X[:, 21])
            @test pd.chain == repeat(1:4; inner = 400)
            @test pd.iteration == repeat(201:600, 4)
            @test isempty(held_parameters(d11_models(), pd).held)
        end
    end
end
