using Test
using InferCell
using Random
using Statistics: mean

# Spec §11 phase 12 — transcript decay: seventeen first-order reactions at one
# global constant over transcript length, returning each transcript's four
# monomers and accruing its energy cost.
#
# `TranscriptSource` is in `corea_decay_doubles.jl`; `HeldGlycolytic` in
# `nucleotide_test_models.jl`; `caught` in `corea_test_models.jl`. All are
# included first by `runtests.jl`.

# Composed jump-state index by name, off the composition itself.
_jidx(models, s) = findfirst(==(s), reduce(vcat, states.(models)))

@testset "Phase 12 — transcript decay" begin
    dec = CoreATranscriptDecay()
    tx = CoreATranscription()
    genes = transcription_genes(tx)

    @testset "12.1 one global constant over transcript length" begin
        @test length(reactions(dec)) == 17
        free = filter(p -> !p.fixed, parameters(dec))
        @test length(free) == 1
        @test only(free).name === :krnadeg
        @test only(free).value == (18 / 452) * 88
        @test only(free).provenance.informedness === :asserted

        # Each constant is krnadeg / n_g, so k·n is one number and the
        # half-life is a function of length alone.
        ks = last.(transcript_decay_constants(dec))
        @test all(i -> ks[i] * genes[i].length ≈ RNADEG_KCAT, eachindex(genes))
        hl = last.(transcript_half_lives(dec))
        @test sortperm(hl) == sortperm([g.length for g in genes])
        @test all(i -> hl[i] ≈ log(2) * genes[i].length / RNADEG_KCAT, eachindex(genes))

        # Two genes of equal length would share a half-life; the propensity
        # reads nothing per-gene but the length.
        w = fill(1, 17)
        rates = [r.rate(zeros(5), [RNADEG_KCAT], 0.0, w) for r in reactions(dec)]
        @test rates ≈ ks

        # Phase 10's standalone double runs the same law, so the two cannot
        # disagree.
        @test all(i -> ks[i] == transcript_decay_constant(genes[i]), eachindex(genes))

        # The upstream hand-tuning comment is recorded, not silently inherited.
        doc = string(@doc RNADEG_KCAT)
        @test occursin("INSTEAD OF", doc)
        @test occursin("0.00578/2", doc)          # the dead `krnadeg` it is not
    end

    @testset "12.2 decrements the transcript it does not own" begin
        @test inputs(dec) == [transcript_state(g.locus) for g in genes]
        @test written_states(dec) == inputs(dec)
        @test !any(s -> s in inputs(dec), states(dec))

        # Firing gene 5 lowers transcript 5 and no other.
        for i in (1, 5, 17)
            u = zeros(5); w = fill(3, 17)
            reactions(dec)[i].affect!(u, w)
            @test w[i] == 2
            @test all(w[j] == 3 for j in 1:17 if j != i)
        end

        # First order: zero copies, zero propensity.
        w = fill(0, 17)
        @test all(r -> r.rate(zeros(5), [RNADEG_KCAT], 0.0, w) == 0.0, reactions(dec))

        # A run from one copy each decays to zero and never below.
        src = TranscriptSource(genes; n0 = 1)
        prob = build_problem([src, dec]; tspan = (0.0, 20_000.0))
        Random.seed!(12)
        sol = solve(prob, SSAStepper(); saveat = 100.0)
        u = reduce(hcat, sol.u)
        @test all(>=(0), u[1:17, :])
        @test all(==(0), u[1:17, end])          # ~45 half-lives of the longest (443 s)

        # Remove the declaration and the write throws at its first firing,
        # naming the module and the state.
        mute = CoreATranscriptDecay(written = Symbol[])
        prob = build_problem([TranscriptSource(genes; n0 = 5), mute];
                             tspan = (0.0, 1000.0))
        err = caught(() -> solve(prob, SSAStepper()))
        @test err isa ArgumentError
        msg = sprint(showerror, err)
        @test occursin("CoreATranscriptDecay", msg)
        @test occursin("mRNA_JCVISYN3A", msg)
        @test occursin("written_states", msg)
    end

    @testset "12.3 each firing returns exactly its own base counts" begin
        # Against phase 10's extract, read through the transcription module, so
        # the two modules cannot disagree on a base count.
        @test [g.locus for g in dec.genes] == [g.locus for g in genes]
        order = [c.counter for c in DECAY_COUNTERS]
        @test order == [:ATP_mRNAdeg, :AMP_mRNAdeg, :GMP_mRNAdeg,
                        :CMP_mRNAdeg, :UMP_mRNAdeg]
        for (i, g) in enumerate(genes)
            u = zeros(Int, 5); w = fill(1, 17)
            reactions(dec)[i].affect!(u, w)
            @test u == [g.length, g.counts.A, g.counts.G, g.counts.C, g.counts.U]
            @test u[2] + u[3] + u[4] + u[5] == g.length
        end
    end

    @testset "12.4 two counters reach recycling, two are chemostat-exempt" begin
        edges = coupling(dec)
        @test length(edges) == 5
        @test all(e -> e isa DeferredCounterEdge, edges)
        @test [(e.counter, e.species, e.direction) for e in edges] ==
              [(:ATP_mRNAdeg, :M_atp_c, :in), (:AMP_mRNAdeg, :M_amp_c, :out),
               (:GMP_mRNAdeg, :M_gmp_c, :out), (:CMP_mRNAdeg, :M_ctp_c, :out),
               (:UMP_mRNAdeg, :M_utp_c, :out)]
        @test all(e -> e.clip === :clamped_deficit_carried, edges)

        # Decay alone: transcription's own CTP/UTP counters would put both
        # in the exemptions whatever decay declared.
        @test isempty(filter(s -> s in (:M_ctp_c, :M_utp_c),
                             resolve_coupling(AbstractSubModel[TranscriptSource(genes)]).chemostat_exemptions))
        r = resolve_coupling(AbstractSubModel[TranscriptSource(genes), dec])
        @test :M_ctp_c in r.chemostat_exemptions
        @test :M_utp_c in r.chemostat_exemptions
        credited = [x.species for x in r.edges
                    if x.declared_by === :CoreATranscriptDecay && x.edge.direction === :out]
        @test sort(credited) == sort([:M_amp_c, :M_gmp_c, :M_ctp_c, :M_utp_c])

        # A counter aimed at the wrong pool, a duplicate, or an unknown name is
        # refused at construction rather than silently crediting elsewhere.
        wrong = (DECAY_COUNTERS[1], merge(DECAY_COUNTERS[3], (species = :M_amp_c,)))
        err = caught(() -> CoreATranscriptDecay(counters = wrong))
        @test err isa ArgumentError && occursin("GMP_mRNAdeg", sprint(showerror, err))
        @test caught(() -> CoreATranscriptDecay(
            counters = (DECAY_COUNTERS[2], DECAY_COUNTERS[2]))) isa ArgumentError
        @test caught(() -> CoreATranscriptDecay(
            counters = (merge(DECAY_COUNTERS[2], (counter = :TMP_mRNAdeg,)),))) isa ArgumentError

        # The chemostat credit is ours and registered as such.
        notes = reduction_notes(dec)
        @test length(notes) == 1
        @test occursin("CTP and UTP chemostats", only(notes))

        # A hybrid build cannot take the CTP and UTP counters until task 13.9
        # gives a chemostatted pool an ownerless path (§12, 2026-09-10 E). Pinned,
        # so the day the framework learns it, this fails and points there.
        full = AbstractSubModel[NucleotideRecycling(), HeldGlycolytic(),
                                TranscriptSource(genes), dec]
        err = caught(() -> build_problem(full; tspan = (0.0, 60.0)))
        @test err isa ArgumentError
        @test occursin("CoreATranscriptDecay declares a DeferredCounterEdge " *
                       "debiting :M_ctp_c", sprint(showerror, err))

        # **Executed, not only declared.** Without the two chemostat counters the
        # hybrid builds, and the hook credits what decay returned into the
        # pools recycling closes. Recycling and the glycolytic double conserve
        # guanylate and adenylate exactly, so every particle of change is decay's.
        dec3 = CoreATranscriptDecay(counters = DECAY_COUNTERS[1:3])
        ms = AbstractSubModel[NucleotideRecycling(), HeldGlycolytic(),
                              TranscriptSource(genes; n0 = 5), dec3]
        d = build_problem(ms; tspan = (0.0, 600.0))
        odes = reduce(vcat, states.(ms[1:2]))
        oi(s) = findfirst(==(s), odes)
        guan(u) = (u[oi(:M_gtp_c)] + u[oi(:M_gdp_c)] + u[oi(:M_gmp_c)]) * d.factor
        aden(u) = (u[oi(:M_atp_c)] + u[oi(:M_adp_c)] + u[oi(:M_amp_c)]) * d.factor
        g0, a0 = guan(d.ode.u), aden(d.ode.u)
        Random.seed!(124)
        rec = run_handshake!(d, 600)
        n_end = rec.jump[end][1:17]
        fired = 5 .- n_end
        @test sum(fired) > 0
        # The run ends on a drain, so every accrual has been paid.
        @test all(==(0), rec.jump[end][18:20])

        g_returned = sum(fired[i] * genes[i].counts.G for i in 1:17)
        a_returned = sum(fired[i] * genes[i].counts.A for i in 1:17)
        atp_paid = sum(fired[i] * genes[i].length for i in 1:17)
        dg = guan(rec.ode[end]) - g0
        da = aden(rec.ode[end]) - a0
        @info "12.4 executed credit" fired=sum(fired) g_returned dg a_returned atp_paid da
        # Bounded at 1e-4 of the moiety rather than at one particle: the stiff
        # solve and the fractional carry leave about one particle of residue in
        # a ~80,000-particle adenylate pool, while a missing or doubled credit
        # is thousands.
        @test abs(dg - g_returned) < 1e-4 * g0
        # ATP_mRNAdeg's ADP and phosphate are not credited until task 13.10, so
        # the adenylate moiety loses what decay's energy cost debits.
        @test abs(da - (a_returned - atp_paid)) < 1e-4 * a0
    end

    @testset "12.5 the decay energy counter" begin
        drains = counter_drains(dec)
        e = only(filter(c -> c.counter === :ATP_mRNAdeg, drains))
        @test e.species === :M_atp_c
        @test e.direction === :in
        @test e.produces == (:M_adp_c, :M_pi_c)
        for (i, g) in enumerate(genes)
            u = zeros(Int, 5); w = fill(2, 17)
            reactions(dec)[i].affect!(u, w)
            reactions(dec)[i].affect!(u, w)
            @test u[1] == 2 * g.length                # one ATP per nucleotide
        end
        # The composed model reports which registry species it debits.
        r = resolve_coupling(AbstractSubModel[tx, dec])
        atp = filter(x -> x.edge isa DeferredCounterEdge &&
                          x.edge.counter === :ATP_mRNAdeg, r.edges)
        @test only(atp).species === :M_atp_c
        @test only(atp).declared_by === :CoreATranscriptDecay
    end

    @testset "12.6 and 12.7 monomer closure and the guanylate return" begin
        cycle = 6300.0
        reps = 8
        g_returned = Float64[]
        for seed in 1:reps
            Random.seed!(seed)
            ms = AbstractSubModel[CoreATranscription(seed = seed), dec]
            prob = build_problem(ms; tspan = (0.0, cycle))
            sol = solve(prob, SSAStepper(); saveat = cycle)
            u0, u1 = sol.u[1], sol.u[end]
            ix(s) = _jidx(ms, s)
            pol = (A = u1[ix(:ATP_mRNA)], C = u1[ix(:CTP_mRNA)],
                   G = u1[ix(:GTP_mRNA)], U = u1[ix(:UTP_mRNA)])
            ret = (A = u1[ix(:AMP_mRNAdeg)], C = u1[ix(:CMP_mRNAdeg)],
                   G = u1[ix(:GMP_mRNAdeg)], U = u1[ix(:UMP_mRNAdeg)])
            r = assert_monomer_closure(genes, u0[1:17], u1[1:17], pol, ret)
            @test all(==(0), values(r))
            # Energy closes the same way: one ATP_mRNAdeg per nucleotide returned.
            @test u1[ix(:ATP_mRNAdeg)] == sum(ret)
            push!(g_returned, ret.G)
        end

        # The mutation: one gene returns one guanine too many. Closure fails,
        # and names G and only G.
        bad = [g.locus === :JCVISYN3A_0607 ?
               TranscriptionGene(g.locus, g.reaction, g.length + 1,
                                 merge(g.counts, (G = g.counts.G + 1,)),
                                 g.first_two, g.ptn_count, g.mean_mrna) : g
               for g in genes]
        mdec = CoreATranscriptDecay(genes = bad)
        ms = AbstractSubModel[CoreATranscription(seed = 1), mdec]
        Random.seed!(1)
        sol = solve(build_problem(ms; tspan = (0.0, cycle)), SSAStepper();
                    saveat = cycle)
        u0, u1 = sol.u[1], sol.u[end]
        jx(s) = _jidx(ms, s)
        pol = (A = u1[jx(:ATP_mRNA)], C = u1[jx(:CTP_mRNA)],
               G = u1[jx(:GTP_mRNA)], U = u1[jx(:UTP_mRNA)])
        ret = (A = u1[jx(:AMP_mRNAdeg)], C = u1[jx(:CMP_mRNAdeg)],
               G = u1[jx(:GMP_mRNAdeg)], U = u1[jx(:UMP_mRNAdeg)])
        res = monomer_closure(genes, u0[1:17], u1[1:17], pol, ret)
        @test res.G < 0
        @test res.A == 0 && res.C == 0 && res.U == 0
        err = caught(() -> assert_monomer_closure(genes, u0[1:17], u1[1:17], pol, ret))
        msg = sprint(showerror, err)
        @test occursin("moiety G", msg)
        @test !occursin("moiety A", msg)

        # 12.7: the guanylate leak this module closes. The steady return is the
        # transcription flux's guanine, Σ k_g · T · G_g, whatever decay's
        # constant; the pool is the registry's GTP + GDP + GMP in particles.
        ks = last.(transcription_rate_constants(tx))
        expected = sum(ks[i] * cycle * genes[i].counts.G for i in eachindex(genes))
        pool = sum(species_entry(s).initial_value
                   for s in (:M_gtp_c, :M_gdp_c, :M_gmp_c)) * corea_particles_per_mM()
        @info "12.7 guanylate return" per_cycle=mean(g_returned) expected pool ratio=mean(g_returned) / pool
        @test mean(g_returned) ≈ expected rtol = 0.1
        @test expected ≈ 61_531 rtol = 1e-3
        @test pool ≈ 39_806 rtol = 1e-3
    end
end
