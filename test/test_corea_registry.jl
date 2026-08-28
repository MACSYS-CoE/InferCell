using Test
using InferCell

@testset "Core A′ species registry" begin

    @testset "Cardinalities match the scoping note" begin
        # dev/notes/reduced-syn3a-scoping.md, "State list". If that table changes
        # and the registry does not, this is the test that says so.
        @test n_dynamic_states() == 32
        @test n_chemostats() == 5
        @test length(COREA_SPECIES) == 37

        @test length(species_in_group(:glycolytic)) == 11
        @test length(species_in_group(:adenylate)) == 4
        @test length(species_in_group(:guanylate)) == 3
        @test length(species_in_group(:redox)) == 2
        @test length(species_in_group(:other)) == 2
        @test length(species_in_group(:trna)) == 2
        @test length(species_in_group(:pts)) == 8
        @test length(species_in_group(:chemostat)) == 5

        # Every group in the vocabulary is accounted for above.
        @test sum(length(species_in_group(g)) for g in InferCell.SPECIES_GROUPS) ==
              length(COREA_SPECIES)
    end

    @testset "H2O and H+ are absent" begin
        # Core A′ derives from the NoH2O model with hydrogen-ion accounting
        # removed, which is why PPA is carried as PPi -> 2 Pi.
        @test !is_registered(:M_h2o_c)
        @test !is_registered(:M_h_c)
    end

    @testset "Lookup round-trips in both directions" begin
        for (i, e) in enumerate(COREA_SPECIES)
            @test species_index(e.name) == i
            @test species_entry(i).name == e.name
            @test species_entry(e.name).name == e.name
        end
    end

    @testset "Ordering is stable" begin
        first_read = [e.name for e in COREA_SPECIES]
        second_read = [e.name for e in COREA_SPECIES]
        @test first_read == second_read
        # Canonical order is registry order, and chemostats come last.
        @test first_read[1:n_dynamic_states()] == dynamic_species()
        @test first_read[(n_dynamic_states() + 1):end] == chemostat_species()
    end

    @testset "An unregistered species is rejected by name" begin
        err = caught(() -> species_index(:M_not_a_species_c))
        @test err isa ArgumentError
        @test occursin("M_not_a_species_c", err.msg)

        @test !is_registered(:M_not_a_species_c)
        @test_throws ArgumentError species_entry(:M_not_a_species_c)
        @test_throws ArgumentError species_group(:M_not_a_species_c)
    end

    @testset "Nucleotide species resolve to the nucleotide file" begin
        # The cross-file trap, in the concentration table. The central file shows
        # these at the 0.1 mM prior default because they are out of its module's
        # scope; the nucleotide file governs and has real balanced values.
        gtp = species_entry(:M_gtp_c)
        @test gtp.initial_value == 1.6627
        @test gtp.gstd == 1.57
        @test gtp.source_file == InferCell.NUCLEOTIDE_FILE
        @test gtp.source_file != InferCell.CENTRAL_FILE
        @test gtp.initial_value != 0.1

        for (name, value) in ((:M_gdp_c, 0.2981), (:M_amp_c, 0.0832), (:M_gmp_c, 0.0117))
            e = species_entry(name)
            @test e.initial_value == value
            @test e.source_file == InferCell.NUCLEOTIDE_FILE
        end
    end

    @testset "Central-file species resolve to the central file" begin
        atp = species_entry(:M_atp_c)
        @test atp.initial_value == 3.6529   # ~73,700 particles, not the 1.04 mM
        @test atp.source_file == InferCell.CENTRAL_FILE
        @test species_entry(:M_pi_c).initial_value == 17.8185
        @test species_entry(:M_g6p_c).initial_value == 3.7076
    end

    @testset "Uninformed initial conditions are flagged" begin
        # G3P, PPi and cytosolic lactate sit at the prior median with a geometric
        # standard deviation of 10 in both files. Nothing informed them.
        for name in (:M_g3p_c, :M_ppi_c, :M_lac__L_c)
            e = species_entry(name)
            @test e.informedness == :prior_default
            @test e.initial_value == 0.1
            @test e.gstd == 10.0
            @test !is_informed(name)
        end

        for name in (:M_atp_c, :M_gtp_c, :M_g6p_c, :M_nad_c)
            @test is_informed(name)
            @test species_entry(name).informedness == :balanced
        end

        u = uninformed_species()
        @test :M_g3p_c in u
        @test :M_ppi_c in u
        @test :M_lac__L_c in u
        @test !(:M_atp_c in u)
    end

    @testset "Chemostats are marked and carry their held value" begin
        @test is_chemostatted(:M_glc__D_e)
        @test is_chemostatted(:M_ctp_c)
        @test is_chemostatted(:M_utp_c)
        @test is_chemostatted(:M_aa_pool_c)
        @test is_chemostatted(:M_o2_c)
        @test !is_dynamic(:M_glc__D_e)

        @test held_value(:M_glc__D_e) == 40.0   # medium
        @test held_value(:M_aa_pool_c) == 0.1   # medium

        # CTP, UTP and O2 are chemostatted but nothing has been imported for
        # them yet. `nothing` is the honest answer; a number would not be.
        @test held_value(:M_ctp_c) === nothing
        @test held_value(:M_utp_c) === nothing
        @test held_value(:M_o2_c) === nothing

        @test chemostat_species() == [:M_glc__D_e, :M_ctp_c, :M_utp_c,
                                      :M_aa_pool_c, :M_o2_c]
    end

    @testset "Asking a dynamic state for a held value is an error" begin
        @test is_dynamic(:M_atp_c)
        @test !is_chemostatted(:M_atp_c)
        err = caught(() -> held_value(:M_atp_c))
        @test err isa ArgumentError
        @test occursin("M_atp_c", err.msg)
    end

    @testset "Groups carry the conserved moiety" begin
        # This is what lets the dead-end report say "adenylate" for AMP with no
        # route back and "guanylate" for stranded GMP — one class of error.
        @test species_group(:M_amp_c) == :adenylate
        @test species_group(:M_gmp_c) == :guanylate
        @test species_group(:M_nad_c) == :redox
        @test species_group(:M_ptsg_P_c) == :pts
        @test species_group(:M_ctp_c) == :chemostat
    end

    @testset "Copy-number regimes are recorded" begin
        # Three regimes spanning nearly five orders of magnitude is why Core A′
        # needs several formalisms at once.
        @test species_entry(:M_atp_c).regime == :metabolite
        @test species_entry(:M_ptsg_c).regime == :protein
        @test all(e -> e.regime in InferCell.SPECIES_REGIMES, COREA_SPECIES)
    end

    @testset "Entry validation rejects bad vocabulary" begin
        @test_throws ArgumentError SpeciesEntry(:M_test_c, :not_a_group, :dynamic,
                                                :metabolite, nothing, nothing, nothing, :balanced)
        @test_throws ArgumentError SpeciesEntry(:M_test_c, :redox, :held,
                                                :metabolite, nothing, nothing, nothing, :balanced)
        @test_throws ArgumentError SpeciesEntry(:M_test_c, :redox, :dynamic,
                                                :enormous, nothing, nothing, nothing, :balanced)
        @test_throws ArgumentError SpeciesEntry(:M_test_c, :redox, :dynamic,
                                                :metabolite, nothing, nothing, nothing, :vibes)
    end

    @testset "The prior-default rule is enforced per entry" begin
        # A gstd at or above the prior width is :prior_default and nothing else
        # is — the loader's rule, so a hand-transcribed row cannot drift from it.
        @test_throws ArgumentError SpeciesEntry(:M_test_c, :redox, :dynamic,
                                                :metabolite, 0.1, 10.0, nothing, :balanced)
        @test_throws ArgumentError SpeciesEntry(:M_test_c, :redox, :dynamic,
                                                :metabolite, 0.5, 2.0, nothing, :prior_default)
        # And the consistent forms construct.
        @test SpeciesEntry(:M_test_c, :redox, :dynamic, :metabolite,
                           0.1, 10.0, nothing, :prior_default).informedness == :prior_default
        @test SpeciesEntry(:M_test_c, :redox, :dynamic, :metabolite,
                           0.5, 2.0, nothing, :balanced).informedness == :balanced
    end
end
