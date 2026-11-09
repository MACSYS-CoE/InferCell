# corea_dataset and the mixed likelihood on the reduced syn3A model. The dataset
# tests were built for B8.2 and the setup, AV3 and stream tests for B8.3
# (spec/phases/approximate-method.md task B8, kept at b92af35 when B8 was parked;
# spec/spec.md §12 2026-10-29, "v0.3.0 launch decisions"). KIM-91 revives them
# for the reduced-model pilot (Linear KIM-89, KIM-92).

using Test
using InferCell
using InferCell: _mixed_solve
using ForwardDiff
using Random
using Distributions: LogNormal

@testset "corea_dataset and the mixed likelihood on the reduced model" begin
    # B8.2: an injection, drawn from one recorded seed.
    @testset "corea_dataset: reproducible, drawn from the prior, its seeds disjoint" begin
        a = corea_dataset(17_200_001; ncells = 2, horizon = 120.0)
        @test first.(a.truth.values)[1:length(D11_TARGETS)] == D11_TARGETS
        @test a.truth.purpose === :calibration
        check_truth(d11_models(), a.truth)
        @test a.sigma > 0
        @test a.times == [60.0, 120.0]
        @test length(a.seeds) == 2 && allunique(a.seeds)
        genes = [g.locus for g in read_transcription_genes()]
        @test a.genes == genes
        @test size(a.transcripts) == (2, length(genes), 2)
        @test size(a.bulk) == (length(M0_PANEL), 2)
        @test a.bulk == observe_bulk(a.latent, a.species, M0_PANEL; sigma = a.sigma,
                                     noise_seed = a.bulk_seed)
        b = corea_dataset(17_200_001; ncells = 2, horizon = 120.0)
        @test b.latent == a.latent && b.bulk == a.bulk && b.transcripts == a.transcripts
        @test b.truth.values == a.truth.values && b.sigma == a.sigma && b.seeds == a.seeds
        # The dataset's draws are corea_draws', in its order.
        d = corea_draws(17_200_001; ncells = 2)
        @test d.truth.values == a.truth.values && d.σ == a.sigma && d.cells == a.seeds
        @test d.noise_seed == a.noise_seed && d.bulk_seed == a.bulk_seed
        # σ_b is SIGMA_MET_PRIOR's draw, made after the truth's.
        rng = Xoshiro(17_200_001)
        draw_truth(rng, d11_models(); names = D11_TARGETS, purpose = :calibration)
        @test a.sigma == rand(rng, SIGMA_MET_PRIOR)
        c = corea_draws(17_200_002; ncells = 2)
        @test c.σ != a.sigma && c.truth.values != a.truth.values
    end

    # KIM-89's free set (Spec D1): twenty names, the dataset's truth over all of
    # them, and a different stream from D11's six at the same seed.
    @testset "COREA_TARGETS and a truth over it" begin
        genes = [g.locus for g in read_transcription_genes()]
        @test COREA_TARGETS == vcat(promoter_param.(genes),
                                    [:krnadeg, :kcatF_R_ENO, :kcatF_R_FBA])
        @test length(COREA_TARGETS) == 20 && allunique(COREA_TARGETS)
        @test D11_TARGETS ⊆ COREA_TARGETS
        ms = d11_models()
        @test all(n -> InferCell._parameter_by_name(ms, n).prior isa LogNormal, COREA_TARGETS)
        a = corea_dataset(17_200_001; ncells = 2, horizon = 120.0, targets = COREA_TARGETS)
        @test a.targets == COREA_TARGETS
        @test first.(a.truth.values)[1:20] == COREA_TARGETS
        check_truth(ms, a.truth; names = COREA_TARGETS)
        nominal = nominal_parameter_values(ms)
        @test all(nv -> nv[2] != nominal[nv[1]], a.truth.values[1:20])
        d6 = corea_draws(17_200_001; ncells = 2)
        @test a.sigma != d6.σ && a.seeds != d6.cells
        @test corea_dataset(17_200_001; ncells = 2, horizon = 120.0).targets == D11_TARGETS
    end

    # KIM-89's fit window (Spec D2): a cut of the full run, so the window's saves
    # are the full record's first ones at the same seed.
    @testset "corea_dataset's fit window: a prefix of the full run" begin
        kw = (; ncells = 2, horizon = 300.0, targets = COREA_TARGETS)
        full = corea_dataset(17_200_001; kw...)
        win = corea_dataset(17_200_001; kw..., fit_to = 150.0)   # off the 60 s grid
        @test full.times == 60.0 .* (1:5) && full.fit_to == 300.0
        @test win.times == [60.0, 120.0] && win.fit_to == 150.0
        k = 1:2
        @test win.latent == full.latent[:, :, k]
        @test win.transcripts == full.transcripts[:, :, k]
        @test win.observed == full.observed[:, :, k]
        @test win.bulk == full.bulk[:, k]
        @test win.truth.values == full.truth.values && win.sigma == full.sigma
        @test win.seeds == full.seeds && win.bulk_seed == full.bulk_seed
        @test_throws ArgumentError corea_dataset(17_200_001; kw..., fit_to = 301.0)
        @test_throws ArgumentError corea_dataset(17_200_001; kw..., fit_to = 59.0)
        # mixed_setup takes the windowed record: its saves are every 60 s.
        sw = mixed_setup(win; models = d11_models(), targets = COREA_TARGETS)
        @test sw.saves == [60.0, 120.0] && length(sw.logz) == 2 * length(M0_PANEL)
        # At production size: the window keeps 52 of the cycle's 105 saves, the
        # last at 3120 s.
        @test COREA_FIT_S == 3150.0 == COREA_CYCLE_S / 2
        grid = 60.0 .* (1:round(Int, COREA_CYCLE_S / 60))
        @test length(grid) == 105 && count(<=(COREA_FIT_S), grid) == 52 && grid[52] == 3120.0
    end

    @testset "the 200 injections' seeds, at production size" begin
        reps = 17_200_000 .+ (1:200)
        # Replicate seeds disjoint from V7's and A2's.
        @test isdisjoint(reps, 16_200_000 .+ (1:50))
        @test isdisjoint(reps, 20_000 .+ (1:1470))
        # Within an injection, cells, reference cells and chains never meet
        # (reference_seeds and mixed_chain_seeds refuse it); chains never meet
        # across injections either. Reference and cell seeds are UInt32 draws, so
        # across injections they may coincide, at different truths. The draws
        # are the pilot's, over COREA_TARGETS (KIM-92).
        chains = Int[]
        for r in reps
            d = corea_draws(r; targets = COREA_TARGETS)
            ds = (; replicate_seed = r, seeds = d.cells)
            refs = reference_seeds(ds, 800)
            ch = reduce(vcat, [mixed_chain_seeds(ds, 16; arm) for arm in 0:3])
            @test isdisjoint(refs, d.cells) && isdisjoint(ch, d.cells) && isdisjoint(ch, refs)
            append!(chains, ch)
        end
        @test allunique(chains)
        # And from V7's 50 M0 replicates' chain seeds.
        v7 = reduce(vcat, [mixed_chain_seeds((; replicate_seed = 16_200_000 + r, seeds = Int[]),
                                             16; arm) for r in 1:50 for arm in 0:3])
        @test isdisjoint(chains, v7)
    end

    # B8.3: the reduced model in mixed_setup, AV3 and B7.3's stream identities,
    # on KIM-89's free set: the pilot's own setup. One tight route (ii) solve
    # takes minutes here (Slurm 20358816: 211 s at 180 s; the gradient 598 s), so
    # AV3 and the per-stream sweeps, about 5 h each, run only with
    # INFERCELL_SLOW_TESTS=true (KIM-91). The unit suite keeps the setup and the
    # stream identities, at the default tolerances and a 120 s horizon. To run the
    # slow tier alone (about 10 h; KIM-91's run was Slurm 20378723), from the tree
    # to test:
    #   sbatch --time=12:00:00 --export=ALL,INFERCELL_SLOW_TESTS=true \
    #       dev/scripts/testfiles.slurm test/test_corea_mixed.jl
    function corea_case(H)
        ds = corea_dataset(17_200_001; ncells = 2, horizon = H, targets = COREA_TARGETS)
        corea(; kw...) = mixed_setup(ds; models = d11_models(), targets = COREA_TARGETS, kw...)
        s = corea()
        truth = Dict(ds.truth.values)
        xt = vcat([log(truth[n]) for n in s.names], log(ds.sigma))
        n = length(M0_PANEL) * length(s.saves)
        rng = Xoshiro(13)
        A = randn(rng, n, 4)
        return (; ds, corea, s, xt, xs = [xt, mixed_prior_median(s)], n,
                δ = 0.02 .* randn(rng, n), Σ = 1e-3 .* (A * A'))
    end
    # Gradient against central differences, coordinate by coordinate.
    function av3_check(ℓ, x; nonzero = true)
        g = ForwardDiff.gradient(ℓ, x)
        @test all(isfinite, g)
        nonzero && @test all(!iszero, g)
        for j in eachindex(x)
            e = zeros(length(x))
            e[j] = 1e-4
            @test g[j] ≈ (ℓ(x + e) - ℓ(x - e)) / 2e-4 rtol = 1e-6
        end
    end

    c = corea_case(120.0)
    @testset "the setup on the reduced model" begin
        s = c.s
        @test mixed_parameter_names(s) == vcat([Symbol(:log_, t) for t in COREA_TARGETS],
                                               :log_sigma_b)
        @test length(mixed_parameter_names(s)) == 21
        @test s.loci == [g.locus for g in read_transcription_genes()]
        @test length(s.factors) == 2
        @test s.saves == [60.0, 120.0] && s.switches == [0.0, 59.0, 119.0]
        @test Set(s.forwards) == Set([:kcatF_R_ENO, :kcatF_R_FBA])
        @test length(s.logz) == c.n
        @test all(x -> isfinite(mixed_logdensity(s, x)), c.xs)
    end

    @testset "streams on the reduced model: each factor alone" begin
        (; s, corea, δ, Σ) = c
        both = with_reference(corea(; streams = (:bulk, :tx)), δ, Σ)
        txo = with_reference(corea(; streams = (:tx,)), δ, Σ)
        bulko = with_reference(corea(; streams = (:bulk,)), δ, Σ)
        x = c.xt
        θ = exp.(x[1:length(s.names)])
        r = _mixed_solve(s, θ)
        k = route_transcription_rates(s, r.sol, r.p)
        μ = θ[findfirst(==(:krnadeg), s.names)] ./ s.lengths
        ltx = transcript_loglik(s.factors, k, μ)
        lbulk = bulk_loglik(both, mean_field_log_bulk(s, r.sol), c.ds.sigma)
        @test mixed_logdensity(both, x) === mixed_logdensity(with_reference(s, δ, Σ), x)
        @test mixed_loglik(both, x) === ltx + lbulk
        @test mixed_loglik(txo, x) === ltx
        @test mixed_loglik(bulko, x) === lbulk
        @test mixed_logdensity(txo, x) == mixed_log_prior(s, x) + ltx
        @test mixed_logdensity(bulko, x) == mixed_log_prior(s, x) + lbulk
    end

    if get(ENV, "INFERCELL_SLOW_TESTS", "false") == "true"
        cs = corea_case(180.0)
        tight(; kw...) = with_reference(cs.corea(; abstol = 1e-14, reltol = 1e-12, kw...),
                                        cs.δ, cs.Σ)
        @testset "AV3 on the reduced model at two x, 21 coordinates" begin
            @test cs.s.saves == [60.0, 120.0, 180.0]
            t = tight()
            for x in cs.xs
                av3_check(z -> mixed_logdensity(t, z), x)
            end
        end
        @testset "AV3 on the reduced model, each stream alone" begin
            for st in ((:tx,), (:bulk,))
                t = tight(; streams = st)
                av3_check(z -> mixed_logdensity(t, z), cs.xt; nonzero = false)
            end
        end
    else
        @info "test_corea_mixed.jl: AV3 on the reduced model skipped; set " *
              "INFERCELL_SLOW_TESTS=true (about 10 h: sbatch --time=12:00:00 " *
              "--export=ALL,INFERCELL_SLOW_TESTS=true dev/scripts/testfiles.slurm " *
              "test/test_corea_mixed.jl)"
    end
end
