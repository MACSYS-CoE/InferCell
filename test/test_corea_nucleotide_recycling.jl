using Test
using InferCell
using Distributions: LogNormal
using SciMLBase: ReturnCode

# Spec §11 phase 8 — nucleotide recycling: the five reactions that make GTP and
# close the adenylate, guanylate and phosphate moieties.

# The composition, its state indices, the three moiety sums and the drift
# helpers all live in `nucleotide_test_models.jl`, which `test/runtests.jl`
# includes first and which `dev/scripts/full_cycle_recycling.jl` also includes.
# One text for one measurement: the artefact records what this file asserts.

# The all-five 6,300 s trajectory, solved once. Testsets 8.6, 8.6b and 8.7 all
# read it, and the ladder's fifth rung is this pair, so binding it here removes
# four bitwise-identical solves.
const FULL_SOL = recycling_solve(recycling_models())
const TIGHT_SOL = recycling_solve(recycling_models();
                                  abstol = ABSTOL_R / 10, reltol = RELTOL_R / 10)

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
        E_adk1 = 213 / corea_particles_per_mM()
        kcatR_adk1, km_adp = 783.2866, 0.2669
        for adp in (0.05, 0.2669, 1.7)
            u = [0.0, adp, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0]
            v = recycling_fluxes(u, p, 0.0, m, zero_in)[3]
            r = adp / km_adp
            @test -v / (E_adk1 * kcatR_adk1) ≈ (r / (1 + r))^2 rtol = 1e-12
        end

        E_ppa = 190 / corea_particles_per_mM()
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
        E_gk1, kcatR_gk1 = 186 / corea_particles_per_mM(), 127.921
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
        # Six decimals, against `corea_particles_per_mM()` = 20,180.39. PGK3 is
        # 0.020366 and not the 0.020367 a rounded 20,180 gives — the sixth
        # decimal is exactly where the two factors part, and it is the decimal
        # the cross-module equality below is decided on.
        expected = Dict(:R_ADK1 => 0.010555, :R_PPA => 0.009415, :R_GK1 => 0.009217,
                        :R_PGK3 => 0.020366, :R_PYK3 => 0.027304)
        for e in enz
            @test round(e.concentration, digits = 6) == expected[e.reaction]
            @test e.concentration == e.copies / corea_particles_per_mM()
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

        # One enzyme, one concentration, across modules. Spec task 8.4 assigns
        # this to whichever pull request lands second; phase 6 landed first as
        # PR #50, so it is this one's, and both halves are in the tree. It is
        # an equality rather than a tolerance on purpose: a transcribed 20,180
        # against the handshake's derived 20,180.39 fails it at 1.9e-5, which
        # is the whole reason the factor is derived in both modules.
        glycolytic = Dict(String(r.locus) => r for r in GLYCOLYTIC_REACTIONS)
        shared_checked = 0
        for e in enz
            e.shared_with_glycolysis || continue
            g = glycolytic[e.locus]
            @test g.copies == e.copies
            @test g.copies / corea_particles_per_mM() == e.concentration
            shared_checked += 1
        end
        @test shared_checked == 2
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
        # every ODE module — belongs to the last of them to land, which is phase
        # 9. The resolver throws on exactly that conflict (`_check_kind_agreement`),
        # so resolving the four is the assertion.
        @test resolve_coupling(AbstractSubModel[CentralGlycolysis(), PtsTransport(),
                                                NucleotideRecycling(), TrnaCharging()]) isa CouplingGraph
    end


    @testset "8.6 a full cycle against the charging drain" begin
        sol = FULL_SOL
        @test sol.retcode == ReturnCode.Success
        @test sol.t[end] == CYCLE_S
        u0 = sol.u[1]

        # Both moieties conserved over the full cycle. These assert the SECOND
        # gate of §3's exact-conservation exception, not the fivefold fall —
        # see spec §12, 2026-09-11, and testset "the second gate" below for the
        # per-evaluation assertion and the ladder that license it.
        bound = tolerance_bound(sol, ADENYLATE_COEFFS)
        @test max_drift(sol, adenylate_of) < 1e-10
        @test max_drift(sol, guanylate_of) < 1e-10
        @test max_drift(sol, adenylate_of) < bound / 100
        @test adenylate_of(u0) ≈ 3.9539 atol = 1e-4
        @test guanylate_of(u0) ≈ 1.9725 atol = 1e-4

        # The glycolytic double holds its four pools. Asserted, because it did
        # not: with a zero derivative the pools took `NucleotideRecycling`'s
        # contributions straight into their own `du` and 13DPG was 99.8% gone
        # by t = 0.5 s, which killed the GTP branch inside the first save
        # interval while every docstring said the pools were clamped.
        for (i, sp) in zip((DPG_I, PG3_I, PEP_I, PYR_I), HELD_GLYCOLYTIC_VALUES[1:4])
            @test maximum(abs(u[i] - sp) for u in sol.u) / sp < 0.01
        end

        # And with the pools held, what stops the GTP branch is guanylate,
        # which is the structural limit: nothing in Core A′ consumes GTP until
        # phase 9's translation, so PGK3 and PYK3 run until GDP is spent.
        @test sol.u[end][GTP_I] > 0.99 * guanylate_of(u0)
        @test sol.u[end][GDP_I] < 0.01 * u0[GDP_I]

        # Tightening the solver tenfold does not shrink it, and that is the
        # claim rather than a shortfall: both settings sit on the same floor.
        tight = TIGHT_SOL
        @test max_drift(tight, adenylate_of) < 1e-10
        @test max_drift(tight, guanylate_of) < 1e-10

        # ATP stays positive, and pyrophosphate settles rather than climbing.
        # The reverse term is what puts it there: a pyrophosphatase written
        # irreversibly would not settle at a finite pool.
        @test minimum(u[ATP_I] for u in sol.u) > 3.4
        ppi_late = [u[PPI_I] for u in sol.u[(end - 10):end]]
        @test maximum(ppi_late) - minimum(ppi_late) < 1e-9
        @test 0.2 < sol.u[end][PPI_I] < 0.6

        # The drain delivers the published demand. Read off the trajectory, not
        # assumed: the double integrates what it charged.
        delivered = sol.u[end][DRAIN_CUM_I] / CYCLE_S
        @test delivered ≈ RECYCLING_DRAIN_MM_PER_S rtol = 5e-3

        # Adenylate kinase removed: the threshold crossing, on a grid fine
        # enough to tell 120 s from 180 s.
        no_kinase = recycling_solve(recycling_models(reactions = without(:R_ADK1));
                                    horizon = 600.0, saveat = 1.0)
        t_cross = crossing_time(no_kinase, 0.01)
        @test t_cross !== nothing
        # Two constant-drain figures, and the difference matters because the
        # spec quotes one of them. The scoping note computes the FULL pool at
        # the published demand, 3.9539 / 0.027408 = 144.3 s. The archived
        # design subtracts the AMP already present, (3.9539 − 0.0832) /
        # 0.027408 = 141.2 s, and labels that "the module's own arithmetic".
        # Spec §3 check 4 and task 8.6 both ask for agreement with the note's
        # 144 s within an order of magnitude, which is what is asserted; the
        # closer agreement with 141.2 s is recorded, not required, because it
        # is a property of this composition's saturating drain and phase 9's
        # mass-action module will not reproduce it.
        @test 14.4 < t_cross < 1443.0
        @test t_cross ≈ 144.3 rtol = 0.05      # the scoping note's figure
        @test t_cross ≈ 141.2 rtol = 0.05      # the archived design's

        # A *crossing*, not exhaustion: the drain saturates as ATP falls, so ATP
        # decays toward zero rather than through it, and phase 9's mass-action
        # module will do the same for a different reason. Asserted against the
        # state's own integrator bound, which is spec §3 check 1's rule — a
        # decaying state may undershoot by the local error the solver is allowed
        # and a tighter floor would be asserting something about round-off.
        atp_bound = max(ABSTOL_R, RELTOL_R * no_kinase.u[1][ATP_I])
        atp_floor = minimum(u[ATP_I] for u in no_kinase.u)
        @test atp_floor > -atp_bound
        @test no_kinase.u[end][ATP_I] < 1e-6
        # The pool is stranded as AMP rather than lost, which is what says the
        # kinase is the missing return path and not a leak.
        @test max_drift(no_kinase, adenylate_of) < 1e-10
        @test no_kinase.u[end][AMP_I] > 0.9 * adenylate_of(no_kinase.u[1])

        # Pyrophosphatase removed. It strands the phosphate moiety rather than
        # letting the concentration diverge, and it cannot do otherwise: this
        # composition's phosphate is closed. The scoping note's 173 mM is
        # open-pool arithmetic that assumes charging runs the whole cycle
        # (spec §11 task 8.6, annotated 2026-09-10).
        no_ppa = recycling_solve(
            recycling_models(reactions = (:R_PGK3, :R_PYK3, :R_ADK1, :R_GK1)))
        budget = phosphate_of(u0)
        @test no_ppa.u[end][PPI_I] > 100 * no_ppa.u[1][PPI_I]
        @test 2 * no_ppa.u[end][PPI_I] > 0.5 * budget
        @test no_ppa.u[end][ATP_I] < 0.01 * u0[ATP_I]
        # ... and it flattens, because charging stops once ATP is gone.
        tail = [u[PPI_I] for u in no_ppa.u[(end - 5):end]]
        @test maximum(tail) - minimum(tail) < 1e-6
        # Against the enzyme being present, which is the comparison that says
        # the reaction is required rather than merely present.
        @test no_ppa.u[end][PPI_I] > 20 * sol.u[end][PPI_I]

        @info "8.6 full cycle" atp_min = minimum(u[ATP_I] for u in sol.u) atp_floor_no_kinase = atp_floor ppi_settled = sol.u[end][PPI_I] t_cross = t_cross adenylate_drift = max_drift(sol, adenylate_of) guanylate_drift = max_drift(sol, guanylate_of) ppi_no_ppa = no_ppa.u[end][PPI_I] phosphate_budget = budget
    end

    # ------------------------------------------------------------------
    # Spec §3's exact-conservation exception, second gate (§12, 2026-09-11).
    # This is what licenses tasks 8.6 and 8.7 to assert a floor rather than
    # the fivefold fall, so it is asserted rather than argued.
    # ------------------------------------------------------------------
    @testset "8.6b the second gate: one ulp per evaluation, and the ladder" begin
        prob = build_problem(recycling_models(); tspan = (0.0, CYCLE_S))
        rhs(u) = prob.f(u, prob.p, 0.0)

        # The three moiety vectors, as weights on the composed state.
        aden_w(du) = du[ATP_I] + du[ADP_I] + du[AMP_I]
        guan_w(du) = du[GTP_I] + du[GDP_I] + du[GMP_I]
        # Flux-corrected phosphate. Standalone the inbound flux is exactly the
        # GTP gain, so the correction is a state difference and the corrected
        # closure stays a linear functional with GTP weighted 3 − 1 = 2.
        phos_w(du) = du[PI_I] + 3du[ATP_I] + 2du[ADP_I] + du[AMP_I] +
                     2du[GTP_I] + 2du[GDP_I] + du[GMP_I] + 2du[PPI_I]

        sol = FULL_SOL
        u0 = sol.u[1]
        sums = (adenylate_of(u0), guanylate_of(u0), phosphate_of(u0))

        # One ulp of the conserved sum, at every state on the trajectory. Not
        # an `≈`: one ulp is the tightest non-zero bound there is.
        for (w, total, name) in zip((aden_w, guan_w, phos_w), sums,
                                    ("adenylate", "guanylate", "phosphate"))
            worst = maximum(abs(w(rhs(u))) for u in sol.u)
            @test worst <= eps(total)
            @info "second gate, $name" worst ulps = worst / eps(total)
        end

        # It is a real gate: a mutated stoichiometry misses it by orders of
        # magnitude, not by a factor.
        mprob = build_problem(recycling_models(mutation = :adk1_adp_coefficient);
                              tspan = (0.0, CYCLE_S))
        mworst = maximum(abs(aden_w(mprob.f(u, mprob.p, 0.0))) for u in sol.u)
        @test mworst > 1e6 * eps(sums[1])

        # And the ladder, over six decades. The bound is each rung's own
        # `tol_C`, not a flat ulp count: the residual neither falls nor
        # accumulates — over nine decades it wanders between 128 and 1,119
        # ulps of the sum while the evaluation count grows 36-fold — so a flat
        # ulp bound would be a number fitted to this one composition, whereas
        # `tol_C` tightens with the solver and cannot be passed by loosening.
        rungs = [(1e-6, 1e-4), (1e-7, 1e-5), (1e-8, 1e-6), (1e-9, 1e-7)]
        ladder = [recycling_solve(recycling_models(); abstol = a, reltol = r)
                  for (a, r) in rungs]
        # The last two rungs are the pair 8.6 already solved.
        push!(rungs, (ABSTOL_R, RELTOL_R));         push!(ladder, FULL_SOL)
        push!(rungs, (ABSTOL_R / 10, RELTOL_R / 10)); push!(ladder, TIGHT_SOL)

        residuals = Float64[]
        for ((a, r), l) in zip(rungs, ladder)
            push!(residuals, max_drift(l, adenylate_of))
            @test residuals[end] <
                  tolerance_bound(l, ADENYLATE_COEFFS; abstol = a, reltol = r) / 1e3
        end
        @test length(residuals) == 6
        @test maximum(residuals) < 100 * minimum(residuals)
        # The point of the ladder: it does not fall. If it ever starts to, the
        # exception no longer applies and the fivefold fall is the right rule.
        @test maximum(residuals) / minimum(residuals) > 2
        @info "second gate ladder, adenylate" residuals
    end

    @testset "8.7 phosphate closure, both forms, and the mutations" begin
        # Exact: with the GTP branch inactive and no phosphorylation, nothing
        # carries phosphate across the boundary. The charging drain still runs
        # and closes internally — three phosphates in ATP become one in AMP and
        # two in pyrophosphate — so this tests the module rather than the
        # absence of traffic.
        exact = recycling_solve(recycling_models(reactions = (:R_ADK1, :R_GK1, :R_PPA),
                                                 k_slp = 0.0))
        @test max_drift(exact, phosphate_of) < 1e-10

        # Flux-corrected: all five active, and the inbound flux *subtracted*
        # rather than the bound relaxed. That the uncorrected drift is twelve
        # orders of magnitude larger is what says the correction does work.
        # The inbound flux is now bounded by *guanylate* — GDP is driven to GTP
        # and the branch stops — so it is ~0.309 mM, the whole GDP + GMP pool,
        # rather than the 0.0505 mM a draining 13DPG pool used to cut it off at.
        full = FULL_SOL
        uncorrected = max_drift(full, phosphate_of)
        @test corrected_phosphate_drift(full) < 1e-10
        @test uncorrected > 1e6 * corrected_phosphate_drift(full)
        @test uncorrected ≈ full.u[end][GTP_I] - full.u[1][GTP_I] rtol = 1e-6
        @test uncorrected ≈ full.u[1][GDP_I] + full.u[1][GMP_I] rtol = 0.02

        # The substrate-level phosphorylation crosses nothing: it takes every
        # phosphate from the free pool. A first version of the double drew them
        # from the held 13DPG pool instead, which let 345 mM of phosphate into a
        # closed moiety and carried pyrophosphate to 51 mM.
        @test full.u[end][SLP_CUM_I] > 100.0

        # A conservation test that cannot fail is not evidence. Each mutation
        # breaks exactly one moiety and leaves the other passing, so a failure
        # localises; the phosphate sum spans both groups and a phosphorylated
        # species going wrong shows up there too, which is why it is asserted
        # separately rather than as a third independent moiety.
        # `Tsit5` here, not `Rodas5P`: these three assert gross violations at
        # 600 s, and specialising the stiff solver a second time on the
        # mutant's type costs about fifteen seconds of suite compilation for
        # values that agree to ten significant figures.
        adk1 = recycling_solve(recycling_models(mutation = :adk1_adp_coefficient);
                               horizon = 600.0, alg = Tsit5())
        @test max_drift(adk1, adenylate_of) > 1e-3
        @test max_drift(adk1, guanylate_of) < 1e-10

        gk1 = recycling_solve(recycling_models(mutation = :gk1_gdp_created);
                              horizon = 600.0, alg = Tsit5())
        @test max_drift(gk1, guanylate_of) > 1e-3
        @test max_drift(gk1, adenylate_of) < 1e-10

        ppa = recycling_solve(recycling_models(mutation = :ppa_phosphate_coefficient);
                              horizon = 600.0, alg = Tsit5())
        @test corrected_phosphate_drift(ppa) > 1e-3
        @test max_drift(ppa, adenylate_of) < 1e-10
        @test max_drift(ppa, guanylate_of) < 1e-10

        @test caught(() -> MutatedRecycling(:not_a_mutation)) isa ArgumentError

        @info "8.7 phosphate closure" exact = max_drift(exact, phosphate_of) corrected = corrected_phosphate_drift(full) uncorrected = uncorrected mutated_adenylate = max_drift(adk1, adenylate_of) mutated_guanylate = max_drift(gk1, guanylate_of) mutated_phosphate = corrected_phosphate_drift(ppa)
    end

    @testset "8.8 the drain's stoichiometry and rate, recorded for phase 9" begin
        # 3,484,518 residues over a 6,300 s cycle. Derived, and derived here
        # rather than typed: phase 9 must *meet* this demand, and a value it may
        # calibrate against and then re-check would verify arithmetic (§4 D14).
        @test 3_484_518 / 6300 ≈ RECYCLING_DRAIN_PER_S rtol = 1e-4
        @test RECYCLING_DRAIN_MM_PER_S == RECYCLING_DRAIN_PER_S / corea_particles_per_mM()

        # One ATP in, one AMP and one pyrophosphate out: three phosphates become
        # one plus two, so the transfer closes phosphate internally and needs no
        # correction in the closure check.
        drain = ChargingDrain()
        @test Set(contributed_states(drain)) == Set([:M_atp_c, :M_amp_c, :M_ppi_c])
        c = contributions([0.0], Float64[], 0.0, drain, [3.6529])
        @test c[1] < 0 && c[2] > 0 && c[3] > 0
        @test c[2] ≈ -c[1] && c[3] ≈ -c[1]
        @test 3 * c[1] + 1 * c[2] + 2 * c[3] ≈ 0.0 atol = 1e-15

        # The anti-circularity guard, made mechanical. The demand lives in the
        # test double; the production module must not carry it, or phase 9 could
        # calibrate `k_chg` against a number this module already assumed.
        #
        # Digit separators are stripped before the search. A guard that tested
        # the bare spellings alone passed while the module header carried
        # "3,484,518" — a test that cannot fail is not evidence (spec §3), and
        # this one could not.
        src = read(joinpath(@__DIR__, "..", "src", "organisms", "coreA",
                            "nucleotide_recycling.jl"), String)
        digits_only = replace(src, "," => "", "_" => "")
        @test !occursin("553", digits_only)
        @test !occursin("3484518", digits_only)
    end

end
