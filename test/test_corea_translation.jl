using Test
using InferCell
using Random
using Statistics: mean, median

import InferCell: states, parameters, reactions, inputs, written_states,
                  coupling, rebuilt_params, reduction_notes

# Spec §11 phase 11 — translation: one reaction per transcript plus ptsG
# translocation, with the rate constant read from the lumped charged-tRNA pool.
#
# `TranscriptSource` is in `corea_decay_doubles.jl`; `HeldGlycolytic` in
# `nucleotide_test_models.jl`; `caught` in `corea_test_models.jl`. All are
# included first by `runtests.jl`.

@testset "Phase 11 — translation" begin
    genes = read_transcription_genes()

    @testset "11.1 residues per gene, stop codon excluded" begin
        res = read_translation_residues()
        @test length(res) == 17
        # The scoping note's full-proteome figure, and the one k_chg was
        # calibrated on: Σ copies × residues.
        @test sum(g.ptn_count * res[g.locus] for g in genes) == 3_484_518
        # ptsG, ptsI, Crr, ptsH. The recorded 746/574/155/90 counted the stop.
        @test [res[l] for l in (:JCVISYN3A_0779, :JCVISYN3A_0233,
                                :JCVISYN3A_0234, :JCVISYN3A_0694)] == [745, 573, 154, 89]
        @test all(res[g.locus] == g.length ÷ 3 - 1 for g in genes)
    end

    @testset "11.3 the restart law, with the lumped pool's per-aa share" begin
        res = read_translation_residues()
        ppm = corea_particles_per_mM()
        ribo = RIBOSOME_COPIES / ppm
        @test RIBO_KCAT == 12 && RIBO_KD == 1e-3 && RIBO_K0 == 1e-4

        # The polysome: min(15, max(1, round(L/125 - 1))). ptsH (270 nt) has
        # one ribosome; ptsG (2238 nt) sixteen, capped at fifteen.
        @test ribosomes_per_transcript(270) == 1
        @test ribosomes_per_transcript(1017) == 7
        @test ribosomes_per_transcript(2238) == 15
        @test translation_kcat(270) == 0.45 * 12
        @test translation_kcat(1017) == (0.25 * 7 + 0.2) * 12

        # Hand-checked at GAPD and the nominal 0.2 mM charged pool: each of
        # the twenty-one tRNA concentrations reads 0.2/20 = 0.01 mM.
        g = genes[5]
        @test g.locus === :JCVISYN3A_0607
        c = 0.01
        hand = 23.4 / ((1 + 1e-4 / ribo) * 1e-6 / c^2 + 338 * 1e-3 / c + 338)
        @test translation_rate_constant(g.length, res[g.locus], 0.2) ≈ hand rtol = 1e-14

        # At the per-aa share upstream's own pools hold (150 copies each), the
        # lumped law is the restart law: nothing but the concentration moved.
        up = 150 / ppm
        kcat = translation_kcat(g.length)
        upstream = kcat / ((1 + RIBO_K0 / ribo) * RIBO_KD^2 / up^2 +
                           res[g.locus] * RIBO_KD / up + res[g.locus])
        @test translation_rate_constant(g.length, res[g.locus], 20 * up) ≈ upstream rtol = 1e-14

        # The elasticity is the analytic derivative: check it by a central
        # difference in log space.
        for gg in genes
            n, r = gg.length, res[gg.locus]
            e = translation_elasticity(n, r, 0.2)
            h = 1e-6
            fd = (log(translation_rate_constant(n, r, 0.2 * exp(h))) -
                  log(translation_rate_constant(n, r, 0.2 * exp(-h)))) / 2h
            @test e ≈ fd rtol = 1e-6
        end

        # An exhausted pool is floored at one particle per share, as upstream
        # floors each pool with max(1, count), rather than dividing by zero.
        k0 = translation_rate_constant(g.length, res[g.locus], 0.0)
        @test k0 ≈ translation_rate_constant(g.length, res[g.locus], 20 / ppm) rtol = 1e-14
        @test 0 < k0 < translation_rate_constant(g.length, res[g.locus], 0.2)
    end

    tl = CoreATranslation()
    res = last.(residue_counts(tl))
    n_state = length(states(tl))
    cyto = 18
    ix(s) = findfirst(==(s), states(tl))

    @testset "11.2 seventeen jumps, catalytic in the transcript" begin
        @test length(reactions(tl)) == 18                # 17 + translocation
        @test states(tl)[1:17] == [protein_state(g.locus) for g in genes]
        @test states(tl)[cyto] === TL_PTSG_CYTO
        @test !any(is_registered, states(tl))
        # Transcripts are read, as a peer's state, and never written.
        @test inputs(tl) == [transcript_state(g.locus) for g in genes]
        @test isempty(written_states(tl))

        free = filter(p -> !p.fixed, parameters(tl))
        @test [q.name for q in free] == [translation_rate_param(g.locus) for g in genes]
        @test rebuilt_params(tl) == [q.name for q in free]
        ks = [q.value for q in free]

        gtp, trna = ix(:GTP_translat), ix(:tRNA_translat)
        iptsg = findfirst(g -> g.locus === TL_PTSG, genes)
        for i in 1:17
            u = zeros(Int, n_state); w = fill(4, 17)
            r = reactions(tl)[i]
            # First order in the transcript: the propensity is k·mRNA, and
            # zero at zero copies.
            @test r.rate(u, ks, 0.0, w) == ks[i] * 4
            @test r.rate(u, ks, 0.0, [j == i ? 0 : 4 for j in 1:17]) == 0.0
            r.affect!(u, w)
            @test w == fill(4, 17)                        # the transcript is unchanged
            target = i == iptsg ? cyto : i
            @test u[target] == 1
            @test all(u[j] == 0 for j in 1:cyto if j != target)
            @test u[gtp] == 2 * res[i]                    # two GTP per residue
            @test u[trna] == res[i]                       # one charged tRNA per residue
        end

        # It composes with the transcripts' real owner, reading them as a
        # phase-2 peer state, and makes protein only where there is mRNA.
        Random.seed!(1102)
        tx = CoreATranscription()
        ms = AbstractSubModel[tx, tl]
        sol = solve(build_problem(ms; tspan = (0.0, 600.0)), SSAStepper(); saveat = 600.0)
        off = length(states(tx))
        p0, p1 = sol.u[1][off .+ (1:17)], sol.u[end][off .+ (1:17)]
        @test all(p1 .>= p0)
        @test sum(p1 .- p0) + sol.u[end][off + cyto] > 0
    end
end
