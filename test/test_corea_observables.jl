# Spec §11 phase 15: observables and synthetic data.

using Test
using InferCell
using Random

@testset "Observables and synthetic data (spec §11 phase 15)" begin
    ms = corea_models()
    promoters = [promoter_param(g.locus) for g in read_transcription_genes()]
    metabolites = [:M_atp_c, :M_gtp_c, :M_fdp_c, :M_13dpg_c]
    transcripts = [transcript_state(g.locus) for g in read_transcription_genes()]

    @testset "15.1: the circularity guard" begin
        @test length(promoters) == 17
        proteins = protein_count_states(ms)
        # Seventeen translated proteins, cytosolic ptsG, and the eight carriers.
        @test length(proteins) == 17 + 1 + 8
        @test :M_ptsg_P_c in proteins && protein_state(:JCVISYN3A_0607) in proteins

        # Metabolites, transcripts and the volume pass.
        @test assert_no_circularity(ms, vcat(metabolites, transcripts, [:volume_litres])) ===
              nothing

        # A protein count is rejected, and the error names all seventeen promoters.
        err = try
            assert_no_circularity(ms, vcat(metabolites, [protein_state(:JCVISYN3A_0607)]))
            nothing
        catch e
            e
        end
        @test err isa ArgumentError
        @test all(p -> occursin(string(p), err.msg), promoters)
        @test occursin("observes protein counts", err.msg)
        # So is a carrier.
        @test_throws ArgumentError assert_no_circularity(ms, [:M_atp_c, :M_ptsi_c])

        # Freeing a protein initial condition is rejected, and names them too.
        ic = Symbol(protein_state(:JCVISYN3A_0445), "0")
        err = try
            assert_no_circularity(ms, metabolites; free = [ic])
            nothing
        catch e
            e
        end
        @test err isa ArgumentError && occursin("frees protein initial conditions", err.msg)
        @test all(p -> occursin(string(p), err.msg), promoters)
        # A freed non-protein initial condition is not the guard's business.
        @test assert_no_circularity(ms, metabolites; free = [:M_atp_c0]) === nothing
    end

    @testset "15.2: truths are drawn, and nominal is a smoke test" begin
        dm = d11_models()
        nominal = nominal_parameter_values(dm)
        t = draw_truth(Xoshiro(15), dm; purpose = :coverage)
        @test first.(t.values)[1:6] == D11_TARGETS
        @test all(n -> Dict(t.values)[n] != nominal[n], D11_TARGETS)
        @test t.label.category === :model_note && occursin("Mode", t.label.description)
        # Haldane-consistent: the draw keeps every equilibrium constant.
        @test all(r -> r.relative <= r.bound, assert_haldane(dm, t.values))
        # Two draws at one seed are identical.
        @test draw_truth(Xoshiro(15), dm).values == t.values

        # A coverage or calibration run refuses a truth at the nominal values.
        smoke = nominal_truth(dm)
        @test smoke.purpose === :smoke && check_truth(dm, smoke) === nothing
        for purpose in (:coverage, :calibration)
            @test_throws ArgumentError nominal_truth(dm; purpose)
            @test_throws ArgumentError check_truth(dm, (; smoke..., purpose))
        end
        # One target set to nominal among drawn ones is refused too.
        v = copy(t.values)
        v[1] = D11_TARGETS[1] => nominal[D11_TARGETS[1]]
        @test_throws ArgumentError check_truth(dm, (; t..., values = v))
        # And a draw cannot be labelled a smoke test.
        @test_throws ArgumentError draw_truth(Xoshiro(1), dm; purpose = :smoke)
        # A reverse constant drawn independently breaks Haldane and is refused.
        bad = draw_parameters(Xoshiro(3), dm, [:kcatF_R_ENO, :kcatR_R_ENO];
                              reverse = :independent)
        @test_throws ArgumentError check_truth(dm, (values = bad, purpose = :coverage,
                                                    label = truth_label());
                                               names = [:kcatF_R_ENO])
    end
end
