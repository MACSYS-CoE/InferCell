# spec/phases/16-recovery.md task 16a.8: M0, the two-gene reference system, and
# its datasets.

using Random
using Test
using InferCell

@testset "16a.8: M0 and its datasets" begin
    ms = m0_models()

    @testset "M0 builds in completeness mode, with two genes and 15 held proteins" begin
        d = build_m0()                                   # complete = true
        @test d.t_end == M0_HORIZON_S
        tx = only(m for m in ms if m isa CoreATranscription)
        @test [g.locus for g in tx.genes] == M0_GENES
        held = only(m for m in ms if m isa HeldEnzymeCounts)
        @test length(held.genes) == 12
        @test isempty(intersect([g.locus for g in held.genes], M0_GENES))
        rep = reduction_report(ms, d)
        @test all(g -> occursin(String(g.locus), rep), held.genes)
        @test occursin("ptsI, ptsH and crr", rep)
        # The held counts stay at their proteomics copy numbers, and GAPD's is live.
        dm = build_m0(tspan = (0.0, 120.0))
        names = InferCell._block_names(ms, :jump)
        idx(s) = findfirst(==(s), names)
        before = [dm.jump.u[idx(protein_state(g.locus))] for g in held.genes]
        @test before == [g.ptn_count for g in held.genes]
        Random.seed!(1608)
        for _ in 1:120
            handshake_step!(dm)
        end
        @test [dm.jump.u[idx(protein_state(g.locus))] for g in held.genes] == before
    end

    @testset "an M0 replicate carries its truth, σ, seeds and report, and reruns identically" begin
        a = m0_dataset(16081; ncells = 2, horizon = 120.0)
        @test first.(a.truth.values)[1:length(M0_TARGETS)] == M0_TARGETS
        @test a.truth.purpose === :calibration
        @test a.sigma > 0 && a.replicate_seed == 16081
        @test length(a.seeds) == 2 && allunique(a.seeds)
        @test a.times == [60.0, 120.0]
        @test size(a.transcripts) == (2, 2, 2)
        @test a.observed_species == vcat(M0_PANEL, [transcript_state(g) for g in M0_GENES])
        @test occursin("M0 holds 12 enzyme counts", a.report)
        # Transcripts are exact.
        for (j, g) in enumerate(M0_GENES)
            @test a.observed[:, length(M0_PANEL) + j, :] == a.transcripts[:, j, :]
        end
        b = m0_dataset(16081; ncells = 2, horizon = 120.0)
        @test b.latent == a.latent && b.observed == a.observed
        @test b.truth.values == a.truth.values && b.sigma == a.sigma && b.seeds == a.seeds
        # σ is drawn per replicate.
        c = m0_dataset(16082; ncells = 1, horizon = 60.0)
        @test c.sigma != a.sigma && c.truth.values != a.truth.values
    end
end
