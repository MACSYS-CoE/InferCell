using Test
using InferCell
using Random
using Statistics: mean, median

import InferCell: states, parameters, reactions, inputs, written_states,
                  coupling, rebuilt_params, reduction_notes, formalism

# Spec §11 phase 11 — translation: one reaction per transcript plus ptsG
# translocation, with the rate constant read from the lumped charged-tRNA pool.
#
# `TranscriptSource` is in `corea_decay_doubles.jl`; `HeldGlycolytic` in
# `nucleotide_test_models.jl`; `caught` in `corea_test_models.jl`. All are
# included first by `runtests.jl`.

# Transcripts that are born and die with no coupling: transcription and decay
# with every counter and edge stripped, so a hybrid build can carry live
# transcripts before task 13.9 lets transcription's CTP and UTP clamps compose.
birth_death() = AbstractSubModel[
    CoreATranscription(counters = (), edges = CouplingEdge[], rebuilt = Symbol[]),
    CoreATranscriptDecay(counters = ())]

# Every module translation feeds: glycolysis and recycling in translated mode,
# so their fifteen enzyme slots read translation's counts, PTS for the carrier
# credits and the volume channel, charging for the tRNA transfer, and live
# transcripts. One composed type, built once per test and reused, because its
# first build compiles for minutes.
full_models(tl = CoreATranslation()) = AbstractSubModel[
    CentralGlycolysis(enzymes = :translated), PtsTransport(),
    NucleotideRecycling(enzymes = :translated), TrnaCharging(),
    birth_death()..., tl]

# Composed index of state `s` among the ODE or jump modules of `ms`, in order.
_block_idx(ms, s, f) = findfirst(==(s), reduce(vcat, [states(m) for m in ms if formalism(m) === f]))

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

    @testset "11.4 ptsG translocation, credited to the carrier PTS owns" begin
        iptsg = findfirst(g -> g.locus === TL_PTSG, genes)
        tr = reactions(tl)[18]
        # Only ptsG has it: the one translocation reads the cytosolic count.
        u = zeros(Int, n_state); u[cyto] = 3
        @test tr.rate(u, Float64[], 0.0, fill(0, 17)) ≈ 3 * 50 / 746 rtol = 1e-14
        @test TRANSLOC_KCAT / (res[iptsg] + 1) == only(q.value for q in parameters(tl)
                                                       if q.name === :tl_transloc_k)
        before = copy(u)
        tr.affect!(u, fill(0, 17))
        changed = findall(u .!= before)
        @test Set(changed) == Set([cyto, iptsg, ix(:ATP_transloc), ix(:ptsG_transloc)])
        @test u[cyto] == 2 && u[iptsg] == 1
        @test u[ix(:ATP_transloc)] == 74                  # int(746/10)
        @test u[ix(:ptsG_transloc)] == 1
        # No cytosolic protein has it: no other reaction ever lowers a count.
        for i in 1:17
            v = zeros(Int, n_state); reactions(tl)[i].affect!(v, fill(1, 17))
            @test all(>=(0), v)
        end
        # The three cytosolic carriers are credited on translation, ptsG only
        # on translocation, and the enzymes never.
        credit = Dict(c.counter => c.locus for c in TRANSLATION_COUNTERS if c.debits === nothing)
        for i in 1:17
            v = zeros(Int, n_state); reactions(tl)[i].affect!(v, fill(1, 17))
            for (c, locus) in credit
                @test v[ix(c)] == (c !== :ptsG_transloc && genes[i].locus === locus ? 1 : 0)
            end
        end

        # Composed: translocation is what raises the carrier PTS owns, and the
        # volume chain sees it through PTS's membrane flag.
        ms = full_models()
        d = build_problem(ms; tspan = (0.0, 600.0))
        ig = _block_idx(ms, :M_ptsg_c, :ode)
        igp = _block_idx(ms, :M_ptsg_P_c, :ode)
        jc = _block_idx(ms, :ptsG_transloc, :jump)
        a0 = growth_census(d).area_nm2
        # The hook writes each pool back in whole particles and carries the
        # remainder (check 0), so a pair's particles move by up to half a
        # particle on a handshake with no credit. Pool plus carried remainder is
        # the ledger that closes exactly.
        rem = d.rounding.remainders
        carriers(u) = (u[ig] + u[igp]) * d.factor + rem[ig] + rem[igp]
        Random.seed!(1104)
        credited = 0
        for _ in 1:120
            c0 = carriers(d.ode.u)
            pending = d.jump.u[jc]
            handshake_step!(d)
            # The PTS cascade conserves each carrier pair exactly, so the
            # ledger moves by the credit and by nothing else.
            @test carriers(d.ode.u) ≈ c0 + pending rtol = 1e-9
            credited += pending
        end
        @info "11.4 ptsG translocated and credited over 120 s" credited
        @test credited > 0
        g = growth_census(d)
        @test Set(s.species for s in g.states) == Set([:M_ptsg_c, :M_ptsg_P_c])
        @test g.area_nm2 > a0
    end

    @testset "11.5 the boundary, and the catalytic channel it feeds" begin
        es = coupling(tl)
        @test length(es) == 9
        rc = only(e for e in es if e isa RateConstantEdge)
        @test rc.species === :M_trna_chg_c && rc.direction === :in
        @test rc.cadence === :piecewise_constant && rc.interval == 60.0
        dc = [(e.counter, e.species, e.direction) for e in es if e isa DeferredCounterEdge]
        @test length(dc) == 8
        @test Set(dc) == Set([
            (:GTP_translat, :M_gtp_c, :in),
            (:tRNA_translat, :M_trna_chg_c, :in), (:tRNA_translat, :M_trna_c, :out),
            (:ATP_transloc, :M_atp_c, :in),
            (:ptsI_translat, :M_ptsi_c, :out), (:ptsH_translat, :M_ptsh_c, :out),
            (:Crr_translat, :M_crr_c, :out), (:ptsG_transloc, :M_ptsg_c, :out)])
        @test all(e.clip === :clamped_deficit_carried for e in es if e isa DeferredCounterEdge)
        # Translation declares no catalytic edge (the consumers own those), no
        # volume edge (PtsTransport owns ptsG's), and nothing as mass: a jump
        # module cannot continuously write an ODE state.
        @test !any(e -> e isa Union{CatalyticEdge, VolumeEdge, MassEdge,
                                    CurrencyEdge, ClampedEdge}, es)
        # GTP → GDP + Pi upstream; the products wait on task 13.10.
        gtp = only(c for c in counter_drains(tl) if c.counter === :GTP_translat)
        @test gtp.produces == (:M_gdp_c, :M_pi_c) && gtp.credits === nothing
        @test resolve_coupling(tl) isa CouplingGraph

        # Composed, every declared channel is lowered: fifteen catalytic
        # exchanges reading thirteen of translation's counts, eight counter
        # channels, one rebuild.
        ms = full_models()
        d = build_problem(ms; tspan = (0.0, 120.0))
        tl_counts = Set(states(ms[end]))
        @test length(d.catalytic) == 15
        @test length(unique(c.species for c in d.catalytic)) == 13
        @test all(c -> c.species in tl_counts, d.catalytic)
        @test Set((b.counter, b.species, b.sign) for b in d.debits) == Set(
            (c, sp, dir === :in ? -1 : 1) for (c, sp, dir) in dc)
        r = only(d.rebuilds)
        @test r.declared_by === :CoreATranslation && r.species == [:M_trna_chg_c]
        @test r.names == rebuilt_params(tl)

        # A catalytic edge moves no matter, so a conservation check refuses it
        # rather than counting it as zero.
        cats = [e for m in ms for e in coupling(m) if e isa CatalyticEdge]
        @test length(cats) == 15
        for e in cats
            @test !carries_mass(e)
            @test caught(() -> mass_contribution(e)) isa ArgumentError
        end

        # Executed: at every handshake each enzyme slot holds its live count,
        # and the counts move, so the modules run on translated protein rather
        # than on nominal stand-ins.
        Random.seed!(1105)
        c0 = [d.jump.u[c.count_idx] for c in d.catalytic]
        for _ in 1:60
            before = [d.jump.u[c.count_idx] for c in d.catalytic]
            handshake_step!(d)
            @test all(d.ode.p[c.param_idx] == counts_to_mM(before[k], d.factor)
                      for (k, c) in enumerate(d.catalytic))
        end
        @test any([d.jump.u[c.count_idx] for c in d.catalytic] .> c0)
        # The rebuild ran and wrote the law at the live charged pool.
        handshake_step!(d)
        @test r.n_refreshes >= 1
    end

    # Phase 9's composition with translation itself as the consumer: the
    # recycling module, the glycolytic double, charging, live transcripts, and
    # translation carrying only the tRNA transfer, so the pool's behaviour is
    # the transfer's alone.
    trna_only = [c for c in TRANSLATION_COUNTERS if c.counter === :tRNA_translat]
    trna_models(tlx; charging = TrnaCharging()) = AbstractSubModel[
        NucleotideRecycling(), HeldGlycolytic(), charging, birth_death()..., tlx]
    function trna_run(tlx, n; seed = 1106, charging = TrnaCharging())
        ms = trna_models(tlx; charging = charging)
        d = build_problem(ms; tspan = (0.0, float(n)))
        Random.seed!(seed)
        rec = run_handshake!(d, n)
        return ms, d, rec, _block_idx(ms, :M_trna_c, :ode), _block_idx(ms, :M_trna_chg_c, :ode)
    end

    @testset "11.6 the debit, and the collapse without it" begin
        tl1 = CoreATranslation(counters = trna_only)
        ms, d, rec, iu, ic = trna_run(tl1, 1200)
        frac = [u[ic] / (u[iu] + u[ic]) for u in rec.ode]
        # Particles, pool plus carried remainder: the transfer is exact, and
        # nothing grows the cell here, so the factor is fixed.
        ledger(u) = (u[iu] + u[ic]) * d.factor
        rem = d.rounding.remainders
        @info "11.6 steady tRNA split with translation consuming" first=frac[60] last=frac[end] mean_2nd_half=mean(frac[601:end]) min=minimum(frac[601:end]) max=maximum(frac[601:end])
        # A stochastic steady state, neither saturated nor empty. Translation's
        # demand follows the transcripts, 0 to 2 copies each, so the split
        # wanders: over a whole cycle (seed 1106) its 300 s means run 0.71 to
        # 0.90 with no trend, and a single ptsG firing debits 745 residues
        # against a pool that relaxes in about 1.5 s. So each block's mean is
        # asserted inside a band far from the collapse's 1.0, rather than two
        # windows asserted equal.
        blocks = [mean(frac[k:k+299]) for k in 1:300:1200]
        @test all(0.6 < b < 0.95 for b in blocks)
        @test all(0.1 < x < 0.9999 for x in frac)
        # tRNA is conserved through the whole run: the pair plus its carries
        # equals the initial pair.
        pool0 = 0.25 * d.factor
        @test ledger(d.ode.u) + rem[iu] + rem[ic] ≈ pool0 rtol = 1e-9

        # The witness. With no transfer at all, nothing consumes charged tRNA:
        # the pool charges completely and the charging flux collapses.
        tl0 = CoreATranslation(counters = NamedTuple[])
        _, d0, rec0, iu0, ic0 = trna_run(tl0, 600)
        f0 = rec0.ode[end][ic0] / (rec0.ode[end][iu0] + rec0.ode[end][ic0])
        @info "11.6 without the transfer" charged_fraction=f0 uncharged_mM=rec0.ode[end][iu0]
        @test f0 > 0.999
        @test rec0.ode[end][iu0] < 1e-3 * rec.ode[end][iu]

        # The defect as first committed: the uncharged pool credited and the
        # charged pool never debited. That does not collapse, it creates tRNA —
        # every residue translated adds one uncharged tRNA from nothing.
        tl2 = CoreATranslation(counters = trna_only,
                               edges = filter(e -> !(e isa DeferredCounterEdge &&
                                                     e.direction === :in), coupling(tl1)))
        _, d2, rec2, iu2, ic2 = trna_run(tl2, 600)
        grown = (rec2.ode[end][iu2] + rec2.ode[end][ic2]) * d2.factor - pool0
        @info "11.6 credit without the debit" tRNA_created_particles=grown
        @test grown > 50_000
    end

    @testset "11.7 check 7's census on the charged-tRNA counter" begin
        # At the nominal pool over twenty minutes: the counter never clips.
        tl1 = CoreATranslation(counters = trna_only)
        ms, d, rec, iu, ic = trna_run(tl1, 1200; seed = 1107)
        census = clipping_census(d)
        chg_min = minimum(u[ic] for u in rec.ode) * d.factor
        @info "11.7 census at the nominal pool" census.clipped census.drains min_charged_particles=chg_min
        @test census.drains == 1200
        @test census.clipped == 0
        # And it can fire: a pool of 0.01 mM holds about 160 charged tRNA
        # against a demand near 380 per second, and the census sees it.
        _, ds, _, _, _ = trna_run(CoreATranslation(counters = trna_only), 120;
                                  seed = 1107, charging = TrnaCharging(pool_mM = 0.01))
        cs = clipping_census(ds)
        @info "11.7 census at a 0.01 mM pool" cs.clipped cs.drains
        @test cs.clipped > 0
    end

    @testset "11.8 residue-to-energy closure, and the done-when's full cycle" begin
        cycle = 6300.0
        ms = AbstractSubModel[CoreATranscription(seed = 1108), CoreATranscriptDecay(), tl]
        Random.seed!(1108)
        sol = solve(build_problem(ms; tspan = (0.0, cycle)), SSAStepper(); saveat = 60.0)
        off = sum(length(states(m)) for m in ms[1:2])
        sl(u) = u[off .+ (1:n_state)]
        u0 = sl(sol.u[1])
        for u in sol.u
            v = sl(u)
            # Seventeen protein counts, non-negative integers, at every write.
            @test all(x -> x >= 0 && isinteger(x), v[1:cyto])
            made = proteins_made(tl, u0, v)
            @test assert_residue_energy_closure(tl, made, v[ix(:GTP_translat)] - u0[ix(:GTP_translat)]) == 0
            @test v[ix(:tRNA_translat)] - u0[ix(:tRNA_translat)] == sum(made .* res)
        end
        made = proteins_made(tl, u0, sl(sol.u[end]))
        @info "11.8 proteins made over one cycle" total=sum(made) residues=sum(made .* res)
        @test all(>(0), made)

        # The mutation: GAPD charged for half its residues. The closure fails
        # and names GAPD, and only GAPD.
        igapd = findfirst(g -> g.locus === :JCVISYN3A_0607, genes)
        bad = copy(res); bad[igapd] ÷= 2
        tlm = CoreATranslation(residues = bad)
        ms2 = AbstractSubModel[CoreATranscription(seed = 1108), CoreATranscriptDecay(), tlm]
        Random.seed!(1108)
        sol2 = solve(build_problem(ms2; tspan = (0.0, 600.0)), SSAStepper(); saveat = 600.0)
        v0, v1 = sl(sol2.u[1]), sl(sol2.u[end])
        made2 = proteins_made(tlm, v0, v1)
        @test made2[igapd] > 0
        r = residue_energy_closure(tlm, made2, v1[ix(:GTP_translat)] - v0[ix(:GTP_translat)])
        @test r == -2 * made2[igapd] * (res[igapd] - bad[igapd])
        err = caught(() -> assert_residue_energy_closure(tlm, made2,
                                                          v1[ix(:GTP_translat)] - v0[ix(:GTP_translat)]))
        msg = sprint(showerror, err)
        @test occursin("JCVISYN3A_0607 (GAPD)", msg)
        @test count("JCVISYN3A_", msg) == 1
    end

    @testset "11.9 no protein degradation: the stochastic block has 52 reactions" begin
        # 17 transcription + 17 decay + 17 translation + 1 translocation, the
        # scoping note's count. ptnDegRate is dead upstream (§12, 2026-09-23 E).
        tx, dec = CoreATranscription(), CoreATranscriptDecay()
        @test length(reactions(tx)) + length(reactions(dec)) + length(reactions(tl)) == 52
        # And no translation reaction lowers a protein count: only
        # translocation lowers one, the cytosolic ptsG, by the one it inserts.
        for (i, r) in enumerate(reactions(tl))
            u = fill(5, n_state); r.affect!(u, fill(1, 17))
            lowered = findall(u .< 5)
            @test lowered == (i == 18 ? [cyto] : Int[])
        end
    end
end
