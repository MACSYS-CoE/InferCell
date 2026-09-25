using Test
using InferCell
using Random
using StaticArrays: SVector

# Spec §11 phase 14a: the balance checks on the assembled Core A′, in the
# scoping note's order, each with a mutation showing it can fail.
#
# Every closure here is whole-cell and in particles, from
# src/organisms/coreA/validation.jl. The full ladder over six decades, the
# three rounding policies over seeds, and the mutation table are in
# `dev/scripts/corea_validation.jl`, whose result file is the numbers of record.
# This file asserts the same quantities at the pinned pair and one decade
# tighter, over the whole cycle, which spec §3 forbids shortening.

const V_SEED = 1410
const V_PINNED = (1e-10, 1e-8)

# Build, seed and record. The seed is set before the build and again before the
# run, so a rerun reproduces the jump path.
function _vrun(models, n; abstol = V_PINNED[1], reltol = V_PINNED[2],
               rounding = :fractional_carry)
    Random.seed!(V_SEED)
    d = build_problem(models; tspan = (0.0, Float64(n)), complete = true,
                      abstol, reltol, rounding)
    Random.seed!(V_SEED)
    return d, validation_run!(models, d, n)
end

_maxres(run, name) = maximum(abs, closure_residual(run, name))

# The moieties that fail with each mutant, and the one that should.
_failing(run) = sort([m.name for m in run.moieties
                      if _maxres(run, m.name) > moiety_bound(run, m.name)])

const V_CYCLE = round(Int, COREA_CYCLE_S)
const V_MUT = 600

@testset "Phase 14a: the balance checks on the assembled Core A′" begin
    ms = validation_models()
    d, run = _vrun(ms, V_CYCLE)
    _, tight = _vrun(ms, V_CYCLE; abstol = V_PINNED[1] / 10, reltol = V_PINNED[2] / 10)
    names = [m.name for m in run.moieties]

    @testset "the metered assembly is the assembly" begin
        # The meter adds two states and changes no rate law: the nine states
        # PtsTransport owns and its four contributions are bitwise the
        # unmetered module's, at the assembled model's initial state.
        inner = PtsTransport()
        metered = MeteredPtsTransport(inner)
        @test module_id(metered) === :PtsTransport
        @test states(metered)[1:9] == states(inner)
        @test coupling(metered) == coupling(inner)
        own = SVector{11}([run.ode[1][findfirst(==(s), run.ode_names)] for s in states(metered)])
        ins = SVector{4}([run.ode[1][findfirst(==(s), run.ode_names)] for s in inputs(inner)])
        pl = SVector(COREA_INITIAL_RADIUS_NM)     # the one free slot, :r_cell_nm
        @test dynamics(own, pl, 0.0, metered, ins)[1:9] ==
              dynamics(SVector{9}(own[1:9]), pl, 0.0, inner, ins)
        @test contributions(own, pl, 0.0, metered, ins) ==
              contributions(SVector{9}(own[1:9]), pl, 0.0, inner, ins)
        # The unmetered assembly composes the plain module, and the doubles'
        # default composition is exactly the metered assembly.
        @test corea_models()[2] isa PtsTransport
        @test corea_models(metered = true)[2] isa MeteredPtsTransport
        @test typeof.(ms) == typeof.(corea_models(metered = true))
    end

    @testset "14.1 check 0: the round trip on the assembly" begin
        # Fractional carry keeps every remainder within half a particle over
        # 6,300 handshakes, and the whole-cell closures stay at roundoff. They
        # are not flat: each handshake's dilution and write-back costs about an
        # ulp, so the residual accumulates with handshake count (phosphate grows
        # 9.4× from 630 to 6,300, job 17539616). That is the accumulation
        # N_restarts puts in tol_C, so the bound that tests it is tol_C, and at
        # a millionth of it this is far tighter than the tolerance principle's
        # three orders.
        @test run.max_remainder <= 0.5
        for name in (:adenylate, :guanylate, :phosphate, :carrier_ptsI)
            @test _maxres(run, name) < 1e-6 * moiety_bound(run, name)
        end

        # A rejected policy is visible at once: deterministic rounding injects
        # whole particles, orders of magnitude above the carry's roundoff, so no
        # downstream residual can be mistaken for a rounding artefact. Only the
        # deterministic policy is asserted here. Stochastic rounding is recorded
        # by the driver over five seeds, where its smallest margin over the
        # carry is 3.4e8 (dev/scripts/corea_validation_result.md).
        _, det = _vrun(ms, 630; rounding = :deterministic)
        for name in (:adenylate, :guanylate, :phosphate)
            @test _maxres(det, name) > 1e6 * max(_maxres(run, name), 1e-12)
        end
        @test det.max_remainder == 0.0     # nothing is carried, so nothing telescopes
    end

    @testset "14.2 check 1: non-negativity, naming the first violation" begin
        @test first_negative(run) === nothing
        @test assert_nonnegative(run) === nothing
        @test minimum(run.jump_min) >= 0
        # The pinned tolerances, read back off the run the bounds come from.
        @test (run.abstol, run.reltol) == V_PINNED

        # The mutation: translation's debits unclamped, so a clip takes the pool
        # it draws below zero rather than carrying the shortfall.
        # Over the full cycle: when a counter first clips depends on the seed, and
        # phase 13b's first GTP clip was at 953 s.
        _, neg = _vrun(validation_models(translation = CoreATranslation(clip = :unclamped)), V_CYCLE)
        v = first_negative(neg)
        @test v !== nothing
        err = caught(() -> assert_nonnegative(neg))
        @test err isa ArgumentError
        @test occursin(":$(v.species)", err.msg)
        @test occursin("t = $(v.t) s", err.msg)
        @test v.species in (:M_gtp_c, :M_trna_chg_c)
    end

    @testset "14.3 check 2: carbon" begin
        b = moiety_bound(run, :carbon)
        @test assert_conserved(run, :carbon; bound = b) <= b
        # Thirteen weighted terms from three modules: the composed right-hand
        # side cancels carbon to within n ulps per evaluation, the second
        # gate as amended 2026-09-25, not to one.
        g = rhs_gate(d, run, :carbon)
        @test g.gate in (:bitwise, :one_ulp, :n_ulps)
        @test g.n == 13
        ca = carbon_accounts(run)
        @test ca.glucose_in > 1e6
        # Two lactate per glucose consumed. This is the carbon closure at the
        # endpoints written as a ratio, not an independent result.
        @test abs(ca.homolactic - 2) < 1e-10
        # The independent quantity: the exporter carries the lactate, and nearly
        # all of what is formed leaves the cell.
        @test ca.exported > 0.99

        # The mutation: LDH written as making two lactate. Carbon alone fails.
        _, mut = _vrun(validation_models(glycolysis = glycolysis_mutation(:R_LDH_L, :M_lac__L_c, 2)), V_MUT)
        @test _failing(mut) == [:carbon]
        err = caught(() -> assert_conserved(mut, :carbon; bound = moiety_bound(mut, :carbon)))
        @test err isa ArgumentError && occursin("carbon moiety", err.msg)
    end

    @testset "14.4 check 3: redox at composition scope" begin
        # Phase 6's check, restated in particles with N_restarts, through the
        # same `assert_conserved` and `conservation_bound` phase 6 now calls.
        @test rhs_gate(d, run, :redox).gate === :bitwise
        b = moiety_bound(run, :redox)
        @test assert_conserved(run, :redox; bound = b) <= b

        _, mut = _vrun(validation_models(glycolysis = glycolysis_mutation(:R_GAPD, :M_nad_c, 2)), V_MUT)
        @test _failing(mut) == [:redox]
        err = caught(() -> assert_conserved(mut, :redox; bound = moiety_bound(mut, :redox)))
        @test err isa ArgumentError && occursin("redox moiety", err.msg)
    end

    @testset "14.5 check 4: adenylate and guanylate over a full cycle" begin
        # The interval is a full cycle: 144 s of it looked fine before the dead
        # end was found.
        @test last(run.t) == COREA_CYCLE_S
        @test run.n_restarts == 6300
        for name in (:adenylate, :guanylate)      # separately, so a failure names one
            b = moiety_bound(run, name)
            @test assert_conserved(run, name; bound = b) <= b
            @test rhs_gate(d, run, name).gate in (:bitwise, :one_ulp)
        end

        # The kinase removed: charging's AMP has no way back, so ATP decays. Under
        # mass action the decay is exponential rather than a constant drain, so
        # the assertion is a threshold crossing, reconciled against the scoping
        # note's 144 s for a constant drain. (The pool does reach exactly zero
        # later, once the write-back rounds it to whole particles.) Adenylate,
        # guanylate and phosphate still close.
        _, ko = _vrun(validation_models(recycling = kinase_removed()), 600)
        ia = findfirst(==(:M_atp_c), ko.ode_names)
        k = findfirst(u -> u[ia] < 0.01 * ko.ode[1][ia], ko.ode)
        @test k !== nothing
        @test 144 < ko.t[k] < 600
        for name in (:adenylate, :guanylate, :phosphate)
            @test _maxres(ko, name) <= moiety_bound(ko, name)
        end

        # Mutations, one per moiety.
        _, adk = _vrun(validation_models(recycling = MutatedRecycling(:adk1_adp_coefficient;
                                                                      enzymes = :translated)), V_MUT)
        # ADK1 writing one ADP also strands a phosphate, which is chemistry, not
        # a leak between checks.
        @test _failing(adk) == [:adenylate, :phosphate]
        _, gk = _vrun(validation_models(recycling = MutatedRecycling(:gk1_gdp_created;
                                                                     enzymes = :translated)), V_MUT)
        @test _failing(gk) == [:guanylate, :phosphate]

        # Across the boundary: GAPD's counters charge ten A for ten U, so the ODE
        # block pays ten ATP per transcript more than the sequence the closure
        # reads from the extract. Only adenylate can see it. GAPD because it is
        # the most transcribed of the seventeen, so 600 s holds events.
        _, bnd = _vrun(validation_models(transcription =
                           CoreATranscription(genes = shifted_genes(:JCVISYN3A_0607;
                                                                   A = 10, U = -10))), V_MUT)
        @test _failing(bnd) == [:adenylate]
    end

    @testset "14.6 check 4b: phosphate, through the chemostats and the transcripts" begin
        b = moiety_bound(run, :phosphate)
        @test assert_conserved(run, :phosphate; bound = b) <= b
        g = rhs_gate(d, run, :phosphate)
        @test g.gate in (:bitwise, :one_ulp, :n_ulps)
        @test g.n == 21

        # The CTP and UTP chemostats carry phosphate across the boundary, weighted
        # by what each counter moves, not by the species.
        rows = chemostat_census(d)
        @test Set(r.counter for r in rows) ==
              Set([:CTP_mRNA, :UTP_mRNA, :CMP_mRNAdeg, :UMP_mRNAdeg])
        @test all(r -> r.particles > 0, rows)
        phos = only(m for m in run.moieties if m.name === :phosphate)
        @test Dict(phos.chemostat) == Dict(:CTP_mRNA => 3.0, :UTP_mRNA => 3.0,
                                           :CMP_mRNAdeg => 1.0, :UMP_mRNAdeg => 1.0)

        # Charging's pyrophosphate needs no correction: ATP's three phosphates
        # become AMP's one and pyrophosphate's two, inside the ODE block. Asserted
        # on the composed right-hand side rather than assumed.
        ch = TrnaCharging()
        i(s) = findfirst(==(s), run.ode_names)
        u = run.samples[end].u
        uin = SVector{length(inputs(ch))}([u[i(s)] for s in inputs(ch)])
        own = SVector{2}([u[i(s)] for s in states(ch)])
        fp = model_free_params(parameters(ch))
        pc = SVector{length(fp)}([q.value for q in fp])
        c = Dict(zip(contributed_states(ch), contributions(own, pc, 0.0, ch, uin)))
        @test c[:M_atp_c] < 0
        @test c[:M_amp_c] == -c[:M_atp_c] && c[:M_ppi_c] == -c[:M_atp_c]
        # So per unit of flux the coefficients are exactly (-1, +1, +1), and the
        # phosphate they carry, 3·(-1) + 1·1 + 2·1, is exactly zero. (Weighting
        # the rates themselves would test the rounding of 3v, not the chemistry.)
        coeff = Dict(k => v / -c[:M_atp_c] for (k, v) in c)
        @test 3 * coeff[:M_atp_c] + 1 * coeff[:M_amp_c] + 2 * coeff[:M_ppi_c] === 0.0

        _, ppa = _vrun(validation_models(recycling = MutatedRecycling(:ppa_phosphate_coefficient;
                                                                      enzymes = :translated)), V_MUT)
        @test _failing(ppa) == [:phosphate]
    end

    @testset "14.7 check 5: the four carriers, up to what translation adds" begin
        for c in (:carrier_ptsI, :carrier_ptsH, :carrier_Crr, :carrier_ptsG)
            b = moiety_bound(run, c)
            @test assert_conserved(run, c; bound = b) <= b
            @test rhs_gate(d, run, c).gate === :bitwise
        end
        # Translation did add carriers, so the subtraction is doing work: each
        # carrier's two forms hold more particles at the end than at the start.
        idx(sp) = findfirst(==(sp), run.ode_names)
        held(k, a, b) = (run.ode[k][idx(a)] + run.ode[k][idx(b)]) * run.factor[k]
        for (a, b) in ((:M_ptsi_c, :M_ptsi_P_c), (:M_ptsh_c, :M_ptsh_P_c),
                       (:M_crr_c, :M_crr_P_c), (:M_ptsg_c, :M_ptsg_P_c))
            @test held(length(run.ode), a, b) > held(1, a, b) + 1
        end

        _, leak = _vrun(validation_models(pts = MeteredCarrierLeak()), V_MUT)
        # A created phospho-HPr is a created phosphate too.
        @test _failing(leak) == [:carrier_ptsH, :phosphate]
        err = caught(() -> assert_conserved(leak, :carrier_ptsH;
                                            bound = moiety_bound(leak, :carrier_ptsH)))
        @test err isa ArgumentError && occursin("carrier_ptsH", err.msg)
    end

    @testset "the tolerance principle across the checks" begin
        # Every closure is roundoff: at least three orders below its own tol_C
        # at the pinned pair and one decade tighter, and flat rather than
        # falling — the exception's signature, not the fall's (spec §3).
        for name in names
            for r in (run, tight)
                @test _maxres(r, name) < 1e-3 * moiety_bound(r, name)
            end
            a, b = _maxres(run, name), _maxres(tight, name)
            @test max(a, b) <= 100 * max(min(a, b), 1e-12)
        end
    end
end
