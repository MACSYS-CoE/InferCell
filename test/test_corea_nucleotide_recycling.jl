using Test
using InferCell

# Spec §11 phase 8 — nucleotide recycling: the five reactions that make GTP and
# close the adenylate, guanylate and phosphate moieties.

const RECYCLING_DATA = joinpath(@__DIR__, "..", "src", "organisms", "coreA", "data")

recycling_tables() = [
    read_source_table(joinpath(RECYCLING_DATA, "nucleotide_recycling.tsv");
                      file = "nucleotide_balanced"),
    read_source_table(joinpath(RECYCLING_DATA, "nucleotide_recycling_central.tsv");
                      file = "central_balanced"),
]

@testset "Core A′ nucleotide recycling" begin

    @testset "8.1 the two vendored extracts" begin
        nucleo, central = recycling_tables()

        # Row classes and counts. Ten catalytic constants, seventeen Michaelis
        # constants and eight concentrations — seventeen and not nineteen
        # because ADK1 and PPA each name one product species twice under one
        # Michaelis constant (§3, and design.md D1).
        for t in (nucleo, central)
            ids = collect(keys(t.values))
            @test length(ids) == 35
            @test count(startswith("kcatF_R_"), ids) == 5
            @test count(startswith("kcatR_R_"), ids) == 5
            @test count(startswith("km_R_"), ids) == 17
            @test count(startswith("conc_M_"), ids) == 8
        end

        # Every one of the 35 is held by *both*. This is the assertion the
        # phase's provenance machinery rests on: with one holder a `governing`
        # declaration passes without being exercised, so an extract that
        # quietly lost a row would disarm the ambiguity check rather than fail.
        @test keys(nucleo.values) == keys(central.values)

        report = ambiguity_report([nucleo, central])
        @test length(report) == 35
        @test all(a -> length(a.values) == 2, report)

        # Spot values, upstream against the extract. The two that matter most
        # are the ones D2 quotes: pyrophosphatase differs across the files by
        # 902x and guanylate kinase by 54x.
        @test nucleo.values["kcatF_R_PPA"] == 646.7271
        @test central.values["kcatF_R_PPA"] == 583611.6071
        @test central.values["kcatF_R_PPA"] / nucleo.values["kcatF_R_PPA"] ≈ 902.4 atol = 0.1
        @test nucleo.values["kcatF_R_GK1"] == 410.2268
        @test central.values["kcatF_R_GK1"] == 7.5991
        @test nucleo.values["kcatF_R_GK1"] / central.values["kcatF_R_GK1"] ≈ 53.98 atol = 0.01
        @test nucleo.values["conc_M_gtp_c"] == 1.6627
        @test nucleo.values["conc_M_ppi_c"] == 0.1

        # The central file holds no information about these reactions at all:
        # every one of its seventeen Michaelis constants sits at the prior
        # median, at prior width. That is what makes it a balancing artefact
        # rather than a rival estimate.
        km_ids = filter(startswith("km_R_"), collect(keys(central.values)))
        @test all(id -> central.values[id] == 0.1, km_ids)
        @test all(id -> central.informedness[id] === :prior_default, km_ids)

        # Informedness is derived from the geometric standard deviation, never
        # declared: the nucleotide file's forward constants are balanced at
        # gstd 1.0513, and the central file's run from 1e33 to 1e63.
        @test all(r -> nucleo.informedness["kcatF_$r"] === :balanced,
                  ["R_PGK3", "R_PYK3", "R_ADK1", "R_GK1", "R_PPA"])
        @test all(r -> central.informedness["kcatF_$r"] === :prior_default,
                  ["R_PGK3", "R_PYK3", "R_ADK1", "R_GK1", "R_PPA"])

        # Three reverse constants are prior-default even in the governing file
        # (design.md D4). Not asserted — a distribution exists — but not
        # evidence either, and PPi's steady state depends on one of them.
        @test nucleo.informedness["kcatR_R_PYK3"] === :prior_default
        @test nucleo.informedness["kcatR_R_GK1"] === :prior_default
        @test nucleo.informedness["kcatR_R_PPA"] === :prior_default
        @test nucleo.informedness["kcatR_R_PGK3"] === :balanced
        @test nucleo.informedness["kcatR_R_ADK1"] === :balanced

        # The concentration trap runs the other way and is the reason the
        # nucleotide file governs the guanylate species: reading them from the
        # central file would put the pool at 0.1 mM, an order of magnitude low.
        agreeing = ["conc_M_atp_c", "conc_M_adp_c", "conc_M_pi_c", "conc_M_ppi_c"]
        disagreeing = ["conc_M_amp_c", "conc_M_gtp_c", "conc_M_gdp_c", "conc_M_gmp_c"]
        by_id = Dict(a.identifier => a for a in report)
        @test all(id -> by_id[id].agrees, agreeing)
        @test all(id -> !by_id[id].agrees, disagreeing)
        @test all(id -> !by_id["kcatF_$id"].agrees,
                  ["R_PGK3", "R_PYK3", "R_ADK1", "R_GK1", "R_PPA"])
        @test all(id -> central.values[id] == 0.1, disagreeing)

        @test length(disagreements([nucleo, central])) == 31

        # The extracts agree with the registry on all eight initial conditions,
        # each read from the file the registry names as its source. The loader
        # checks this on every import; asserting it here says the extract is
        # what makes that check meaningful rather than vacuous.
        for s in [:M_atp_c, :M_adp_c, :M_pi_c, :M_ppi_c]
            @test central.values["conc_$s"] ≈ species_entry(s).initial_value atol = 5e-5
            @test species_entry(s).source_file == "central_balanced"
        end
        for s in [:M_amp_c, :M_gtp_c, :M_gdp_c, :M_gmp_c]
            @test nucleo.values["conc_$s"] ≈ species_entry(s).initial_value atol = 5e-5
            @test species_entry(s).source_file == "nucleotide_balanced"
        end

        @info "8.1 the cross-file ratios" PPA = central.values["kcatF_R_PPA"] / nucleo.values["kcatF_R_PPA"] GK1 = nucleo.values["kcatF_R_GK1"] / central.values["kcatF_R_GK1"] ADK1 = nucleo.values["kcatF_R_ADK1"] / central.values["kcatF_R_ADK1"] PGK3 = central.values["kcatF_R_PGK3"] / nucleo.values["kcatF_R_PGK3"] PYK3 = central.values["kcatF_R_PYK3"] / nucleo.values["kcatF_R_PYK3"]
    end

end
