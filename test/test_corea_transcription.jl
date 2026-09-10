using Test
using InferCell
using StaticArrays
using Random
using Statistics: mean, median

# Spec §11 phase 10 — transcription: seventeen genes, one constitutive
# reaction each, the corrected base mapping, five cost counters, and rate
# constants the 60 s driver rebuilds but the module never rebuilds itself.
#
# Doubles are in `corea_transcription_doubles.jl` and `hybrid_test_models.jl`,
# both included first by `runtests.jl`. `caught` comes from
# `corea_test_models.jl`.

# Spearman's rank correlation, written here rather than pulled from StatsBase,
# which is not a dependency. No ties occur in either vector.
function _spearman(x, y)
    rx = invperm(sortperm(x)); ry = invperm(sortperm(y))
    n = length(x)
    return 1 - 6 * sum((rx .- ry) .^ 2) / (n * (n^2 - 1))
end

@testset "Phase 10 — transcription" begin
    m = CoreATranscription()
    genes = transcription_genes(m)

    @testset "10.1 the per-gene extract" begin
        @test length(genes) == 17
        @test allunique(g.locus for g in genes)

        # Every gene's four base counts sum to its transcript length.
        for g in genes
            @test g.counts.A + g.counts.C + g.counts.G + g.counts.U == g.length
        end

        # The totals the spec states.
        @test sum(g.counts.A for g in genes) == 7236
        @test sum(g.counts.C for g in genes) == 2078
        @test sum(g.counts.G for g in genes) == 3094
        @test sum(g.counts.U for g in genes) == 5868

        # Four spot rows, byte for byte against the upstream sources.
        by = Dict(g.locus => g for g in genes)
        pgi = by[:JCVISYN3A_0445]
        @test (pgi.length, pgi.counts.A, pgi.counts.C, pgi.counts.G, pgi.counts.U) ==
              (1284, 521, 125, 197, 441)
        @test pgi.ptn_count == 266
        @test pgi.mean_mrna ≈ 0.4403 atol = 5e-5
        @test (by[:JCVISYN3A_0607].length, by[:JCVISYN3A_0607].ptn_count) == (1017, 1355)
        @test by[:JCVISYN3A_0607].mean_mrna ≈ 2.1781 atol = 5e-5
        @test (by[:JCVISYN3A_0779].length, by[:JCVISYN3A_0779].ptn_count) == (2238, 831)
        @test (by[:JCVISYN3A_0694].length, by[:JCVISYN3A_0694].ptn_count) == (270, 290)

        # A missing gene must abort naming the locus, not build a sixteen-gene
        # model. The reader is what the module calls, so mutate a copy of the
        # extract rather than the generator's five upstream files.
        mktempdir() do dir
            path = joinpath(dir, "truncated.tsv")
            kept = filter(l -> !occursin("JCVISYN3A_0694", l),
                          readlines(TRANSCRIPTION_EXTRACT))
            write(path, join(kept, '\n') * "\n")
            err = caught(() -> read_transcription_genes(path))
            @test err isa ArgumentError
            msg = sprint(showerror, err)
            @test occursin("JCVISYN3A_0694", msg)
            @test occursin("ptsH", msg)
        end

        # A row whose base counts do not sum to its length is refused, so a
        # hand-edit that breaks the invariant cannot load.
        mktempdir() do dir
            path = joinpath(dir, "bad_sum.tsv")
            lines = map(readlines(TRANSCRIPTION_EXTRACT)) do l
                startswith(l, "JCVISYN3A_0445") || return l
                f = split(l, '\t'); f[3] = string(parse(Int, f[3]) + 1)
                join(f, '\t')
            end
            write(path, join(lines, '\n') * "\n")
            err = caught(() -> read_transcription_genes(path))
            @test err isa ArgumentError
            @test occursin("sum to", sprint(showerror, err))
        end

        err = caught(() -> read_transcription_genes(joinpath("no", "such.tsv")))
        @test err isa ArgumentError
        @test occursin("extract_transcription_genes.jl", sprint(showerror, err))
    end

    @testset "10.2 the sub-model and its reactions" begin
        @test formalism(m) === :jump
        @test inference_mode(m) === :simulation
        @test length(states(m)) == 22
        @test !any(is_registered, states(m))
        @test states(m)[1:17] == [transcript_state(g.locus) for g in genes]
        @test states(m)[18:22] == [c.counter for c in TRANSCRIPTION_COUNTERS]

        # Genes are a fixed quantity, not a state: no locus appears among the
        # states, and the propensity is the rate constant alone.
        @test !any(s -> occursin("JCVISYN3A", string(s)) &&
                        !startswith(string(s), "mRNA_"), states(m))

        rxns = reactions(m)
        @test length(rxns) == 17
        @test all(r -> r isa Reaction, rxns)

        # Firing one gene's reaction: its transcript up by one, the energy
        # counter up by that gene's length, each monomer counter up by that
        # gene's base count, and no other gene's transcript touched.
        i = 5                       # GAPD, the largest promoter strength
        g = genes[i]
        u = zeros(Int, 22)
        p = [pp.value for pp in model_free_params(parameters(m))]
        # Zeroth order: the propensity is the rate constant and nothing else.
        @test rxns[i].rate(u, p, 0.0, Int[]) == p[17 + i]
        @test rxns[i].rate(fill(99, 22), p, 1234.0, Int[]) == p[17 + i]
        rxns[i].affect!(u, Int[])
        @test u[i] == 1
        @test sum(u[1:17]) == 1
        @test u[18] == g.length              # ATP_trsc
        @test u[19] == g.counts.A            # ATP_mRNA
        @test u[20] == g.counts.G            # GTP_mRNA
        @test u[21] == g.counts.C            # CTP_mRNA
        @test u[22] == g.counts.U            # UTP_mRNA

        # ATP is charged twice per transcript — energy and monomer.
        @test u[18] + u[19] == g.length + g.counts.A

        notes = reduction_notes(m)
        @test any(n -> occursin("replication", n), notes)
    end

    @testset "10.3 the rate constant and the corrected mapping" begin
        conc = (M_atp_c = 3.6529, M_ctp_c = 0.6874, M_gtp_c = 1.6627, M_utp_c = 2.7681)
        ks = [transcription_rate_constant(g, conc) for g in genes]

        # The band the spec's Done-when names.
        @test minimum(ks) ≈ 1.26e-3 rtol = 0.01
        @test maximum(ks) ≈ 8.29e-3 rtol = 0.01
        @test all(k -> 1.26e-3 <= k * 1.01 && k <= 8.29e-3 * 1.01, ks)

        # Under :corrected each base's count enters only its own NTP term.
        # Raising one gene's cytosine count must move the rate constant by
        # exactly the CTP term's share and leave the others alone.
        g = genes[1]
        bump_c = TranscriptionGene(g.locus, g.reaction, g.length + 1,
                                   (A = g.counts.A, C = g.counts.C + 1,
                                    G = g.counts.G, U = g.counts.U),
                                   g.first_two, g.ptn_count, g.mean_mrna)
        bump_g = TranscriptionGene(g.locus, g.reaction, g.length + 1,
                                   (A = g.counts.A, C = g.counts.C,
                                    G = g.counts.G + 1, U = g.counts.U),
                                   g.first_two, g.ptn_count, g.mean_mrna)
        # CTP is scarcer than GTP, so one more cytosine costs more than one
        # more guanine under the corrected mapping, and the reverse under the
        # published one — which is the permutation, made observable.
        dc = transcription_rate_constant(bump_c, conc) -
             transcription_rate_constant(g, conc)
        dg = transcription_rate_constant(bump_g, conc) -
             transcription_rate_constant(g, conc)
        @test dc < dg < 0
        pc = transcription_rate_constant(bump_c, conc; base_mapping = :published) -
             transcription_rate_constant(g, conc; base_mapping = :published)
        pg = transcription_rate_constant(bump_g, conc; base_mapping = :published) -
             transcription_rate_constant(g, conc; base_mapping = :published)
        @test pg < pc < 0

        # The two mappings differ, and by about a percent on the constant.
        pub = [transcription_rate_constant(g, conc; base_mapping = :published)
               for g in genes]
        @test pub != ks
        @test all(abs.(pub .- ks) ./ ks .< 0.05)

        # A longer transcript at equal promoter strength has a smaller constant.
        short = genes[findmin(g -> g.length, genes)[2]]
        long = TranscriptionGene(short.locus, short.reaction, short.length * 2,
                                 (A = short.counts.A * 2, C = short.counts.C * 2,
                                  G = short.counts.G * 2, U = short.counts.U * 2),
                                 short.first_two, short.ptn_count, short.mean_mrna)
        @test transcription_rate_constant(long, conc) <
              transcription_rate_constant(short, conc)

        # The turnover cap is reported, not assumed, and never binds.
        head = turnover_headroom(m)
        @test length(head) == 17
        @test maximum(last.(head)) ≈ 8.8516 rtol = 1e-3
        @test maximum(last.(head)) < TURNOVER_CEILING
        @test last(head[findmax(last.(head))[2]]) ==
              maximum(RNAPOL_KCAT * g.ptn_count / PROMOTER_DIVISOR for g in genes)
        @info "10.3 rate constants" k_min=minimum(ks) k_max=maximum(ks) turnover_max=maximum(last.(head)) ceiling=TURNOVER_CEILING

        # The correction is labelled; the published mapping is not.
        notes = reduction_notes(m)
        @test any(n -> occursin(":corrected", n) && occursin("1.9", n), notes)
        @test count(n -> occursin("base", n) && occursin("mapping", n), notes) == 1

        err = caught(() -> CoreATranscription(base_mapping = :permuted))
        @test err isa ArgumentError
        @test occursin(":permuted", sprint(showerror, err))
    end

    @testset "10.4 NTP concentrations from the balanced tables" begin
        ps = parameters(m)
        byname = Dict(p.name => p for p in ps)
        for (species, value, gstd, file) in
            ((:M_atp_c, 3.6529, 1.2825, "central_balanced"),
             (:M_ctp_c, 0.6874, 2.0574, "nucleotide_balanced"),
             (:M_gtp_c, 1.6627, 1.5684, "nucleotide_balanced"),
             (:M_utp_c, 2.7681, 1.3664, "nucleotide_balanced"))
            p = byname[Symbol("tx_conc_", species)]
            @test p.value ≈ value atol = 5e-5
            @test source_file(p) == file
            @test informedness(p) === :balanced
            @test p.prior.σ ≈ log(gstd) rtol = 1e-9
            # Never the first-minute setup constants.
            @test !isapprox(p.value, SETUP_CONSTANTS[species]; atol = 1e-6)
        end

        # ATP and GTP agree with the registry, whose owner integrates them.
        @test byname[:tx_conc_M_atp_c].value ≈ species_entry(:M_atp_c).initial_value atol = 5e-5
        @test byname[:tx_conc_M_gtp_c].value ≈ species_entry(:M_gtp_c).initial_value atol = 5e-5
        # CTP and UTP are chemostats the registry records no value for, which
        # is why this module carries the value and its provenance itself.
        @test species_entry(:M_ctp_c).initial_value === nothing
        @test species_entry(:M_utp_c).initial_value === nothing
    end

    @testset "10.5 one recomputation entry point, never self-called" begin
        @test length(rebuilt_params(m)) == 17
        @test rebuilt_params(m) == [rate_param(g.locus) for g in genes]

        free = model_free_params(parameters(m))
        freenames = Set(p.name for p in free)
        @test all(n -> n in freenames, rebuilt_params(m))

        # `rate_constants` reproduces the module's own law at the pools it is
        # handed, and returns one value per rebuilt name.
        p = [pp.value for pp in free]
        pools = [3.6529, 0.6874, 1.6627, 2.7681]      # TRANSCRIPTION_POOLS order
        v = rate_constants(p, 0.0, m, pools)
        @test length(v) == 17
        conc = NamedTuple{TRANSCRIPTION_POOLS}(Tuple(pools))
        for (i, g) in enumerate(genes)
            @test v[i] ≈ transcription_rate_constant(g, conc) rtol = 1e-12
        end

        # Different pools move all seventeen and nothing else about the module.
        v2 = rate_constants(p, 0.0, m, [2.0, 0.5, 1.0, 2.0])
        @test all(v2 .!= v)
        @test length(states(m)) == 22 && length(rebuilt_params(m)) == 17

        # The law reads the promoter strengths and the pools, never the slots
        # it fills: perturbing a k_tx slot changes nothing it returns.
        p3 = copy(p); p3[18:34] .*= 1000
        @test rate_constants(p3, 0.0, m, pools) ≈ v rtol = 1e-12
        # Perturbing a promoter strength does move its own constant.
        p4 = copy(p); p4[1] *= 2
        v4 = rate_constants(p4, 0.0, m, pools)
        @test v4[1] ≈ 2 * v[1] rtol = 1e-9
        @test v4[2:end] ≈ v[2:end] rtol = 1e-12
    end

    @testset "10.6 eleven edges, no inputs" begin
        es = coupling(m)
        @test length(es) == 11
        @test isempty(inputs(m))

        rc = filter(e -> edge_kind(e) === :rate_constant, es)
        @test length(rc) == 4
        @test [e.species for e in rc] == collect(TRANSCRIPTION_POOLS)
        @test all(e -> e.direction === :in, rc)
        @test all(e -> e.cadence === :piecewise_constant && e.interval == 60.0, rc)

        dc = filter(e -> edge_kind(e) === :deferred_counter, es)
        @test length(dc) == 5
        @test [(e.counter, e.species) for e in dc] ==
              [(c.counter, c.species) for c in TRANSCRIPTION_COUNTERS]
        @test all(e -> e.direction === :in, dc)
        # All five take the published policy, so none is a labelled deviation.
        @test all(e -> e.clip === :clamped_deficit_carried, dc)
        @test !any(deviates_from_published, dc)

        cl = filter(e -> edge_kind(e) === :clamped, es)
        @test length(cl) == 2
        @test Set(e.species for e in cl) == Set((:M_ctp_c, :M_utp_c))
        @test all(e -> e.origin === :ours, cl)
        @test all(deviates_from_published, cl)
        @test Dict(e.species => e.held_value for e in cl)[:M_ctp_c] ≈ 0.6874
        @test Dict(e.species => e.held_value for e in cl)[:M_utp_c] ≈ 2.7681

        # Three kinds coexist on ATP inbound (two counters, one rate constant)
        # and on CTP and UTP (counter, rate constant, clamp).
        atp = filter(e -> e.species === :M_atp_c, es)
        @test Set(edge_kind(e) for e in atp) == Set((:deferred_counter, :rate_constant))
        @test count(e -> edge_kind(e) === :deferred_counter, atp) == 2
        ctp = filter(e -> e.species === :M_ctp_c, es)
        @test Set(edge_kind(e) for e in ctp) ==
              Set((:deferred_counter, :rate_constant, :clamped))

        # …and the resolver accepts that composition rather than conflicting.
        graph = resolve_coupling([m])
        @test graph !== nothing

        # Each counter reports what it debits and what the drain produces.
        drains = counter_drains(m)
        @test length(drains) == 5
        d = Dict(x.counter => x for x in drains)
        @test d[:ATP_trsc].species === :M_atp_c
        @test Set(d[:ATP_trsc].produces) == Set((:M_adp_c, :M_pi_c))
        for c in (:ATP_mRNA, :GTP_mRNA, :CTP_mRNA, :UTP_mRNA)
            @test d[c].produces == (:M_ppi_c,)
        end
        # Transcription is a second source of pyrophosphate, which the scoping
        # note's PPA argument attributed to charging alone.
        @test count(x -> :M_ppi_c in x.produces, drains) == 4
    end

    @testset "10.7 both promoter-proxy declarations" begin
        labels = reduction_declarations([m])
        texts = [l.description for l in labels]
        @test any(t -> occursin("proxy", t) && occursin("back-calculated", t), texts)
        circ = filter(t -> occursin("circular", t), texts)
        @test length(circ) == 1
        # The second declaration names every affected parameter.
        for g in genes
            @test occursin(string(rate_param(g.locus)), circ[1])
        end

        # The proxy's values, to three decimals.
        s = Dict(promoter_strengths(m))
        @test s[:JCVISYN3A_0445] ≈ 1.478 atol = 5e-4     # PGI
        @test s[:JCVISYN3A_0607] ≈ 7.528 atol = 5e-4     # GAPD
        @test s[:JCVISYN3A_0779] ≈ 4.617 atol = 5e-4     # ptsG
        @test s[:JCVISYN3A_0203] ≈ 1.033 atol = 5e-4     # GK1

        # Two genes' promoter strengths are in the ratio of their copy numbers.
        cn = protein_copy_numbers(m)
        @test s[:JCVISYN3A_0607] / s[:JCVISYN3A_0445] ≈
              cn[:JCVISYN3A_0607] / cn[:JCVISYN3A_0445] rtol = 1e-12

        # The clamps are labelled deviations; the counters are not.
        @test count(l -> l.category === :clamp, labels) == 2

        # One number, three consumers. The cross-check against the metabolic
        # modules' enzyme concentrations needs phases 6 and 7, neither of which
        # has landed on this branch.
        @test_skip "10.7 copy numbers agree with the metabolic modules' — " *
                   "blocked on phase 6 (central glycolysis) and phase 7 " *
                   "(phosphotransferase transport), which declare the loci to " *
                   "compare against"
    end

    @testset "10.8 a full cycle, calibration and the elasticities" begin
        # Transcription needs a consumer or the counts only accumulate; the
        # decay double supplies the published per-gene law. §12, 2026-09-10.
        decay = ToyTranscriptDecay(genes)
        @test written_states(decay) == [transcript_state(g.locus) for g in genes]
        @test issubset(Set(written_states(decay)), Set(inputs(decay)))

        cycle = 6300.0
        reps = 8
        totals = zeros(17)
        nsave = 0
        maxcounter = zeros(Int, 5)
        for seed in 1:reps
            Random.seed!(seed)
            prob = build_problem([CoreATranscription(seed = seed), decay];
                                 tspan = (0.0, cycle))
            sol = solve(prob, SSAStepper(); saveat = 10.0)
            u = reduce(hcat, sol.u)                      # 23 states × times
            @test all(x -> x isa Integer && x >= 0, u)
            totals .+= vec(sum(u[1:17, :]; dims = 2))
            nsave += size(u, 2)
            maxcounter .= max.(maxcounter, u[18:22, end])
            # The counters accumulate monotonically: nothing debits them
            # without a driver.
            for r in 18:22
                @test issorted(u[r, :])
            end
        end
        avg = totals ./ nsave
        meas = [g.mean_mrna for g in genes]
        ratio = avg ./ meas

        @info "10.8 calibration" replicates=reps within2x=count(x -> 0.5 <= x <= 2.0, ratio) spearman=_spearman(avg, meas) worst=maximum(ratio) worst_gene=genes[argmax(ratio)].locus

        # §3's external-validation bound: a rank correlation across the
        # seventeen, and agreement within a factor of two for at least fifteen.
        @test _spearman(avg, meas) >= 0.7
        @test count(x -> 0.5 <= x <= 2.0, ratio) >= 15

        # The counters accrued, and ATP twice over.
        @test all(maxcounter .> 0)

        # The elasticity of the rate constants to the four pools, measured
        # rather than argued, under both mappings. Channel 4's gain: the
        # NMonoSum term is about 4.5% of the denominator and transcript length
        # is the rest, so the only ODE-to-stochastic channel in Core A′ is
        # live, bidirectional and weak.
        #
        # The spec states its bands to two and three significant figures, so
        # they are compared at the precision they are stated — half a unit in
        # the last digit — rather than as hard inequalities on rounded numbers.
        elasticity(mm, pools) = begin
            base = last.(transcription_rate_constants(mm))
            up = last.(transcription_rate_constants(
                mm, NamedTuple{TRANSCRIPTION_POOLS}(
                    Tuple(s in pools ? 1.001 * mm.conc[s] : mm.conc[s]
                          for s in TRANSCRIPTION_POOLS))))
            (log.(up) .- log.(base)) ./ log(1.001)
        end

        corrected = CoreATranscription(base_mapping = :corrected)
        published = CoreATranscription(base_mapping = :published)

        e_all = elasticity(corrected, TRANSCRIPTION_POOLS)
        e_gtp = elasticity(corrected, (:M_gtp_c,))
        p_gtp = elasticity(published, (:M_gtp_c,))
        @info "10.8 elasticity" all_pools_corrected=extrema(e_all) gtp_corrected=extrema(e_gtp) gtp_published=extrema(p_gtp) gtp_ratio=extrema(p_gtp ./ e_gtp)

        # All four pools together, under the corrected mapping: 0.044 to 0.051.
        @test minimum(e_all) ≈ 0.044 atol = 5e-4
        @test maximum(e_all) ≈ 0.051 atol = 5e-4
        @test all(0.0435 .<= e_all .<= 0.0515)

        # GTP alone, under each mapping.
        @test minimum(e_gtp) ≈ 0.0079 atol = 5e-5
        @test maximum(e_gtp) ≈ 0.0117 atol = 5e-5
        @test all(0.00785 .<= e_gtp .<= 0.01175)
        @test all(0.0156 .<= p_gtp .<= 0.0196)

        # The 1.9× the correction is worth, measured rather than argued: the
        # published mapping weights the GTP term by the uracil count, and over
        # these seventeen genes U outnumbers G by 5868 to 3094.
        @test all(p_gtp .> e_gtp)
        @test median(p_gtp ./ e_gtp) ≈ 1.9 atol = 0.2
        @test sum(g.counts.U for g in genes) / sum(g.counts.G for g in genes) ≈
              5868 / 3094
    end
end
