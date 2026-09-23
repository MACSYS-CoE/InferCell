using Test
using InferCell
using Random
using Statistics: mean, median

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
end
