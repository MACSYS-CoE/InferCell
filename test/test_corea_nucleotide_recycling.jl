using Test
using InferCell
using Distributions: LogNormal

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


    @testset "8.2 a governing file for every value" begin
        m = NucleotideRecycling()
        ps = parameters(m)
        @test length(ps) == 40                       # 27 kinetic + 5 enzyme + 8 IC

        # Every imported value names the file it chose and the file it rejected.
        choices = governing_choices(ps)
        @test length(choices) == 35
        @test all(c -> length(c[3]) == 1, choices)

        by_name = Dict(p.name => p for p in ps)
        ppa = by_name[:kcatF_R_PPA]
        @test ppa.value == 646.7271
        @test source_file(ppa) == "nucleotide_balanced"
        @test provenance_of(ppa).alternatives == ["central_balanced" => 583611.6071]
        @test informedness(ppa) === :balanced

        # The five reactions are governed by the nucleotide file without
        # exception; the eight initial conditions are not, and assuming they
        # were would be wrong for four of them.
        for r in RECYCLING_REACTIONS, pre in ("kcatF_", "kcatR_")
            @test source_file(by_name[Symbol(pre, r)]) == "nucleotide_balanced"
        end
        for s in RECYCLING_STATES
            @test source_file(by_name[Symbol(s, "0")]) == species_entry(s).source_file
            @test recycling_governing_file("conc_$s") == species_entry(s).source_file
        end
        @test recycling_governing_file("conc_M_atp_c") == "central_balanced"
        @test recycling_governing_file("conc_M_gtp_c") == "nucleotide_balanced"

        # Dropping a governing declaration is a load-time failure naming the
        # identifier and both values — not a silent choice.
        tables = recycling_tables()
        err = caught(() -> load_parameter(tables, "kcatF_R_PPA";
                                          name = :kcatF_R_PPA,
                                          module_id = :NucleotideRecycling,
                                          prior = LogNormal(log(646.7271), log(1.0513))))
        @test err isa ArgumentError
        @test occursin("kcatF_R_PPA", err.msg)
        @test occursin("646.7271", err.msg)
        @test occursin("583611.6071", err.msg)

        # Declaring the wrong file does not import from it quietly: the
        # registry-agreement check catches a guanylate concentration read from
        # the file that holds it only as a prior default.
        err2 = caught(() -> load_parameter(tables, "conc_M_gtp_c";
                                           name = :M_gtp_c0,
                                           module_id = :NucleotideRecycling,
                                           prior = LogNormal(log(0.1), log(10.0)),
                                           role = :initial_condition,
                                           governing = "central_balanced"))
        @test err2 isa ArgumentError
        @test occursin("M_gtp_c", err2.msg)

        # Not one of the 35 imported values is an asserted prior: every one
        # carries a real balancing distribution. Phase 6 labels a
        # Michaelis-constant column choice and this phase does not, and that
        # asymmetry is a property of the source files rather than an oversight
        # (§4 D3) — so this module declares no column deviation and no lumping.
        imported = [by_name[Symbol(id)] for id in RECYCLING_KINETIC_IDS
                    if !startswith(String(id), "enz_")]
        @test isempty(asserted_prior_params(imported))
        @test isempty(reduction_notes(m))

        # The five enzyme concentrations are the module's only asserted priors,
        # and calling them anything else would be a friendlier label than the
        # source supports: the proteomics table states a count and no width.
        labels = reduction_declarations(m)
        @test length(labels) == 5
        @test all(l -> l.category === :asserted_prior, labels)
        @test Set(l.subject for l in labels) ==
              Set(Symbol("enz_", r) for r in RECYCLING_REACTIONS)

        # Three reverse constants are prior-default even in the governing file:
        # a distribution exists, but nothing informed it.
        @test informedness(by_name[:kcatR_R_PPA]) === :prior_default
        @test informedness(by_name[:kcatR_R_PYK3]) === :prior_default
        @test informedness(by_name[:kcatR_R_GK1]) === :prior_default
        @test informedness(by_name[:M_ppi_c0]) === :prior_default
    end

    @testset "8.3 the five reactions, with repeated stoichiometry expanded" begin
        m = NucleotideRecycling()
        ids = collect(RECYCLING_KINETIC_IDS)
        @test count(id -> startswith(String(id), "km_"), ids) == 17
        @test count(id -> startswith(String(id), "kcat"), ids) == 10
        @test count(id -> startswith(String(id), "enz_"), ids) == 5

        # Seventeen and not nineteen: ADK1 names ADP twice and PPA names
        # phosphate twice, each under one Michaelis constant.
        @test :km_R_ADK1_M_adp_c in ids
        @test count(id -> startswith(String(id), "km_R_ADK1_"), ids) == 3
        @test count(id -> startswith(String(id), "km_R_PPA_"), ids) == 2

        # No rate law names water or a hydrogen ion: this is the NoH2O model,
        # and PYK3's published H+ is not a state.
        @test !any(id -> occursin("M_h_c", String(id)) || occursin("M_h2o_c", String(id)), ids)
        @test !(:M_h_c in states(m)) && !(:M_h_c in inputs(m))

        p = Float64[]
        zero_in = [0.0, 0.0, 0.0, 0.0]

        # The squaring, in closed form and with nothing else left to vary. With
        # both substrates at zero the forward term vanishes and the denominator
        # collapses to the product bracket alone, so
        #     -v / (E·kcatR) = (r/(1+r))^2   for a doubled product,
        # while a reaction with two distinct products gives the plain product of
        # two such ratios. That distinguishes "squared" from "twice" exactly.
        E_adk1 = 213 / 20180
        kcatR_adk1, km_adp = 783.2866, 0.2669
        for adp in (0.05, 0.2669, 1.7)
            u = [0.0, adp, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0]
            v = recycling_fluxes(u, p, 0.0, m, zero_in)[3]
            r = adp / km_adp
            @test -v / (E_adk1 * kcatR_adk1) ≈ (r / (1 + r))^2 rtol = 1e-12
        end

        E_ppa = 190 / 20180
        kcatR_ppa, km_pi = 0.1874, 0.0976
        for pin in (0.01, 0.0976, 17.8185)
            u = [0.0, 0.0, 0.0, pin, 0.0, 0.0, 0.0, 0.0]
            v = recycling_fluxes(u, p, 0.0, m, zero_in)[5]
            q = pin / km_pi
            @test -v / (E_ppa * kcatR_ppa) ≈ (q / (1 + q))^2 rtol = 1e-12
        end

        # GK1 has two *distinct* products, so its reverse term is a product of
        # two ratios and not a square — the control that says the test above is
        # measuring the stoichiometry and not the algebra.
        E_gk1, kcatR_gk1 = 186 / 20180, 127.921
        km_gk1_adp, km_gk1_gdp = 0.0676, 0.0236
        u = [0.0, 0.3, 0.0, 0.0, 0.0, 0.11, 0.0, 0.0]
        v = recycling_fluxes(u, p, 0.0, m, zero_in)[4]
        ra, rg = 0.3 / km_gk1_adp, 0.11 / km_gk1_gdp
        @test -v / (E_gk1 * kcatR_gk1) ≈ (ra * rg) / ((1 + ra) * (1 + rg)) rtol = 1e-12

        # Zero at balance, and it changes sign: the flux is a difference of two
        # terms rather than a one-way rate.
        kcatF_adk1, km_amp, km_atp = 319.1556, 0.096, 0.1426
        amp, atp = 0.0832, 3.6529
        # kcatF·(amp/Kamp)(atp/Katp) = kcatR·(adp/Kadp)^2
        adp_bal = km_adp * sqrt(kcatF_adk1 * (amp / km_amp) * (atp / km_atp) / kcatR_adk1)
        balanced = [atp, adp_bal, amp, 0.0, 0.0, 0.0, 0.0, 0.0]
        @test abs(recycling_fluxes(balanced, p, 0.0, m, zero_in)[3]) < 1e-12
        excess = copy(balanced); excess[2] = 2 * adp_bal
        @test recycling_fluxes(excess, p, 0.0, m, zero_in)[3] < 0
        starved = copy(balanced); starved[2] = 0.5 * adp_bal
        @test recycling_fluxes(starved, p, 0.0, m, zero_in)[3] > 0

        # The pyrophosphatase is carried reversibly. At the registry's initial
        # conditions the forward term still wins, which is what bounds PPi, but
        # the reverse term is not negligible: phosphate sits at 17.8 mM against
        # a Michaelis constant of 0.0976.
        u0 = Float64[species_entry(s).initial_value for s in RECYCLING_STATES]
        v_ppa = recycling_fluxes(u0, p, 0.0, m, zero_in)[5]
        @test v_ppa > 0
        rev_only = copy(u0); rev_only[8] = 0.0
        @test recycling_fluxes(rev_only, p, 0.0, m, zero_in)[5] < 0

        # A deactivated reaction is exactly zero and every other rate is
        # untouched, so a removal configuration is the same module at a
        # different setting rather than a second model.
        no_kinase = NucleotideRecycling(reactions = (:R_PGK3, :R_PYK3, :R_GK1, :R_PPA))
        full = recycling_fluxes(u0, p, 0.0, m, [0.0098, 0.0409, 1.1015, 3.3660])
        cut = recycling_fluxes(u0, p, 0.0, no_kinase, [0.0098, 0.0409, 1.1015, 3.3660])
        @test cut[3] == 0.0
        @test cut[[1, 2, 4, 5]] == full[[1, 2, 4, 5]]
        @test caught(() -> NucleotideRecycling(reactions = (:R_NOPE,))) isa ArgumentError
    end

    @testset "8.4 enzyme concentrations, and the two shared with glycolysis" begin
        enz = recycling_enzymes()
        @test length(enz) == 5
        expected = Dict(:R_ADK1 => 0.010555, :R_PPA => 0.009415, :R_GK1 => 0.009217,
                        :R_PGK3 => 0.020367, :R_PYK3 => 0.027304)
        for e in enz
            @test round(e.concentration, digits = 6) == expected[e.reaction]
            @test e.concentration == e.copies / 20180
        end

        loci = Dict(e.reaction => e.locus for e in enz)
        @test loci[:R_PGK3] == "JCVISYN3A_0606"
        @test loci[:R_PYK3] == "JCVISYN3A_0221"
        @test loci[:R_ADK1] == "JCVISYN3A_0651"
        @test loci[:R_GK1] == "JCVISYN3A_0203"
        @test loci[:R_PPA] == "JCVISYN3A_0344"

        # PGK3 and PYK3 add reactions but no genes, so two modules carry a
        # nominal for each of those two enzymes.
        shared = Set(e.reaction for e in enz if e.shared_with_glycolysis)
        @test shared == Set([:R_PGK3, :R_PYK3])

        m = NucleotideRecycling()
        by_name = Dict(p.name => p for p in parameters(m))
        for e in enz
            @test by_name[Symbol("enz_", e.reaction)].value == e.concentration
            @test by_name[Symbol("enz_", e.reaction)].fixed
        end

        # Not independent of phase 6. Composing central glycolysis and asserting
        # that neither PGK/PGK3 nor PYK/PYK3 is run at two different
        # concentrations belongs to whichever pull request lands second;
        # stubbing glycolysis to make it pass here would assert nothing.
        @test_skip "one enzyme, one concentration across modules — needs spec §11 phase 6"
    end

    @testset "8.5 the boundary" begin
        m = NucleotideRecycling()
        edges = coupling(m)
        @test length(edges) == 11

        declared = Set((edge_kind(e), e.species, e.direction) for e in edges)
        @test declared == Set([
            (:mass, :M_13dpg_c, :in), (:mass, :M_pep_c, :in),
            (:mass, :M_3pg_c, :out), (:mass, :M_3pg_c, :in),
            (:mass, :M_pyr_c, :out), (:mass, :M_pyr_c, :in),
            (:currency, :M_amp_c, :in), (:currency, :M_gmp_c, :in),
            (:currency, :M_ppi_c, :in),
            (:currency, :M_pi_c, :out), (:currency, :M_gtp_c, :out),
        ])

        # ATP and ADP deliberately carry no edge. This module is neither their
        # principal producer nor their principal consumer, so claiming the
        # crossing here would put two principals on one pool. The rule said
        # *sole* until D13 gave the charging step AMP and pyrophosphate too.
        @test !any(e -> e.species in (:M_atp_c, :M_adp_c), edges)

        # Only mass and currency: nothing here is debited on a clock, fills a
        # rate-law slot, is rebuilt at 60 s, reads the geometry, or is clamped.
        @test all(e -> e isa MassEdge || e isa CurrencyEdge, edges)

        # The two products are read as well as written, because the published
        # reverse terms and denominators name them. contributed_states still
        # names each species once and carries the net signed rate.
        @test inputs(m) == [:M_13dpg_c, :M_pep_c, :M_3pg_c, :M_pyr_c]
        @test Set(contributed_states(m)) == Set(inputs(m))
        @test length(contributed_states(m)) == 4

        # The eight owned states, in strictly increasing registry order.
        @test length(states(m)) == 8
        @test Set(states(m)) ==
              union(Set(species_in_group(:adenylate)), Set(species_in_group(:guanylate)),
                    Set([:M_ppi_c]))
        idx = species_index.(states(m))
        @test all(idx[i] < idx[i + 1] for i in 1:(length(idx) - 1))
        @test formalism(m) === :ode

        # Standalone resolution reports rather than fails: the four glycolytic
        # species are unowned, which is a deliberate partial composition and not
        # an incomplete declaration.
        graph = resolve_coupling(m)
        @test :M_13dpg_c in graph.unowned_states
        @test :M_pep_c in graph.unowned_states
        @test !(:M_atp_c in graph.unowned_states)
        @test length(graph.edges) == 11

        # The contribution channel and the edges are held to each other in both
        # directions: dropping either half is drift, and drift is what this
        # check exists to catch.
        drifted = NucleotideRecycling()
        empty!(drifted.ins)
        append!(drifted.ins, [:M_pep_c, :M_3pg_c, :M_pyr_c])
        err = caught(() -> resolve_coupling(drifted))
        @test err isa ArgumentError
        @test occursin("M_13dpg_c", err.msg)

        dropped = NucleotideRecycling()
        empty!(dropped.contribs)
        append!(dropped.contribs, [:M_13dpg_c, :M_pep_c, :M_3pg_c])
        err2 = caught(() -> resolve_coupling(dropped))
        @test err2 isa ArgumentError
        @test occursin("M_pyr_c", err2.msg)

        # Not independent of phases 6, 7 and 9. The four-module assertion — that
        # no species and direction is described as both mass and currency across
        # every ODE module — belongs to the last of them to land, or to phase 13.
        @test_skip "no species is both mass and currency across all four ODE modules — needs spec §11 phases 6, 7, 9"
    end

end
