# corea_dataset, the Core A′ injections built for B8.2. B8 is parked; the generator
# is kept for a later scored recovery test on the reduced model (spec/spec.md §12
# 2026-10-29, "v0.3.0 launch decisions").

using Test
using InferCell
using Random

@testset "corea_dataset: Core A′ injections" begin
    # An injection, drawn from one recorded seed.
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

end
