# spec/phases/approximate-method.md task B8: the mixed likelihood on Core A′
# (§12 2026-10-29, "B8 as planned").

using Test
using InferCell
using InferCell: _mixed_solve
using ForwardDiff
using Random

@testset "B8: the mixed likelihood on Core A′" begin
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

    @testset "the 200 injections' seeds, at production size" begin
        reps = 17_200_000 .+ (1:200)
        # Replicate seeds disjoint from V7's and A2's.
        @test isdisjoint(reps, 16_200_000 .+ (1:50))
        @test isdisjoint(reps, 20_000 .+ (1:1470))
        # Within an injection, cells, reference cells and chains never meet
        # (reference_seeds and mixed_chain_seeds refuse it); chains never meet
        # across injections either. Reference and cell seeds are UInt32 draws, so
        # across injections they may coincide, at different truths.
        chains = Int[]
        for r in reps
            d = corea_draws(r)
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

    # B8.3: Core A′ in mixed_setup, AV3 and B7.3's stream identities.
    H = 180.0
    ds = corea_dataset(17_200_001; ncells = 2, horizon = H)
    corea(; kw...) = mixed_setup(ds; models = d11_models(), targets = D11_TARGETS, kw...)
    s = corea()
    truth = Dict(ds.truth.values)
    xt = vcat([log(truth[n]) for n in s.names], log(ds.sigma))
    xs = [xt, mixed_prior_median(s)]
    n = length(M0_PANEL) * 3
    rng = Xoshiro(13)
    A = randn(rng, n, 4)
    δ, Σ = 0.02 .* randn(rng, n), 1e-3 .* (A * A')

    @testset "the setup on Core A′" begin
        @test mixed_parameter_names(s) == vcat([Symbol(:log_, t) for t in D11_TARGETS],
                                               :log_sigma_b)
        @test s.loci == [g.locus for g in read_transcription_genes()]
        @test length(s.factors) == 2
        @test s.saves == [60.0, 120.0, 180.0] && s.switches == [0.0, 59.0, 119.0, 179.0]
        @test Set(s.forwards) == Set([:kcatF_R_ENO, :kcatF_R_FBA])
        @test length(s.logz) == n
        @test all(x -> isfinite(mixed_logdensity(s, x)), xs)
    end

    @testset "AV3 on Core A′ at two x" begin
        tight = with_reference(corea(; abstol = 1e-14, reltol = 1e-12), δ, Σ)
        ℓ(z) = mixed_logdensity(tight, z)
        for x in xs
            g = ForwardDiff.gradient(ℓ, x)
            @test all(isfinite, g) && all(!iszero, g)
            for j in eachindex(x)
                e = zeros(length(x))
                e[j] = 1e-4
                @test g[j] ≈ (ℓ(x + e) - ℓ(x - e)) / 2e-4 rtol = 1e-6
            end
        end
    end

    @testset "streams on Core A′: each factor alone" begin
        both = with_reference(corea(; streams = (:bulk, :tx)), δ, Σ)
        txo = with_reference(corea(; streams = (:tx,)), δ, Σ)
        bulko = with_reference(corea(; streams = (:bulk,)), δ, Σ)
        x = xt
        θ = exp.(x[1:length(s.names)])
        r = _mixed_solve(s, θ)
        k = route_transcription_rates(s, r.sol, r.p)
        μ = θ[findfirst(==(:krnadeg), s.names)] ./ s.lengths
        ltx = transcript_loglik(s.factors, k, μ)
        lbulk = bulk_loglik(both, mean_field_log_bulk(s, r.sol), ds.sigma)
        @test mixed_logdensity(both, x) === mixed_logdensity(with_reference(s, δ, Σ), x)
        @test mixed_loglik(both, x) === ltx + lbulk
        @test mixed_loglik(txo, x) === ltx
        @test mixed_loglik(bulko, x) === lbulk
        @test mixed_logdensity(txo, x) == mixed_log_prior(s, x) + ltx
        @test mixed_logdensity(bulko, x) == mixed_log_prior(s, x) + lbulk
        for st in ((:tx,), (:bulk,))
            tight = with_reference(corea(; abstol = 1e-14, reltol = 1e-12, streams = st), δ, Σ)
            ℓ(z) = mixed_logdensity(tight, z)
            g = ForwardDiff.gradient(ℓ, x)
            for j in eachindex(x)
                e = zeros(length(x))
                e[j] = 1e-4
                @test g[j] ≈ (ℓ(x + e) - ℓ(x - e)) / 2e-4 rtol = 1e-6
            end
        end
    end
end
