using Test
using InferCell
using Random
using LinearAlgebra

# Spec §11 phase 14b: the derivative, ensemble and census checks on the
# assembled Core A′, each with a mutation showing it can fail.
#
# The full-scale runs are in `dev/scripts/corea_validation_14b.jl` (check 1b's
# ensemble over a cycle, check 6, task 9.9's rerun and F5) and
# `dev/scripts/corea_census.jl` (check 7 over seeds and 200 prior draws).
# Their result files are the numbers of record. This file asserts each check's
# mechanism and its mutation at a size the suite can afford.

const B_SEED = 1410
const B_SHORT = 600

@testset "Phase 14b: the derivative, ensemble and census checks" begin

    # ------------------------------------------------------------------
    @testset "14.8 check 8: the summation theorems on the frozen variant" begin
        cq = frozen_control_problem()
        @test length(cq.multipliers) == 24          # 22 ODE reactions and two demands

        # The conserved combinations: every moiety of 14a that the frozen
        # variant keeps lies in the left null space of its stoichiometry.
        # Carbon is not among them: glucose enters and lactate leaves.
        cm = conservation_matrices(cq)
        names = cq.names[cq.internal]
        inL(pairs) = (v = zeros(length(names));
                      for (s, w) in pairs; v[findfirst(==(s), names)] = w; end;
                      norm(v - cm.L' * (cm.L * v)) / norm(v))
        for m in corea_moieties(carbon = false)
            @test inL(m.ode) < 1e-12
        end
        @test inL((:M_trna_c => 1, :M_trna_chg_c => 1)) < 1e-12
        # And one more, which no 14a check asserts (see the result file).
        @test size(cm.L, 1) == 10

        ss = steady_state(cq)
        @test ss.residual <= 1e-12 * ss.scale
        @test all(<(0), real.(ss.eigenvalues))       # stable in its conserved class
        @test all(>=(0), ss.x)

        cc = control_coefficients(cq, ss)
        # GK1's GMP has no source once expression is frozen, so it carries no
        # steady flux, and its flux has no control coefficient.
        @test cc.zero_fluxes == [:R_GK1]
        dev = assert_summation(cc)
        @test dev.ccc <= 1e-6
        @test dev.fcc <= 1e-6

        # The mutation spec §3 names: without k_chg's multiplier the identities
        # fail, and the check names the quantity.
        mut = omit_multiplier(cc, :R_charging)
        err = try
            assert_summation(mut); nothing
        catch e
            e
        end
        @test err isa ArgumentError
        @test occursin("Check 8 fails", err.msg)
        live = [k for k in eachindex(mut.multipliers) if !(mut.multipliers[k] in mut.zero_fluxes)]
        @test maximum(abs.(mut.fcc_sum[live] .- 1)) > 1e-2
        # The charging flux's own row loses exactly its elasticity to itself.
        k = findfirst(==(:R_charging), cc.multipliers)
        @test cc.fcc[k, k] > 0.01

        # Gene-level coefficients sum a gene's reactions; they satisfy the
        # concentration identity only with the non-gene multipliers added back.
        Cg, genes = group_coefficients(cc.ccc, cc.multipliers, frozen_gene_groups())
        @test length(genes) == 13
        rest = [j for (j, m) in enumerate(cc.multipliers)
                if !any(m in last(g) for g in frozen_gene_groups())]
        @test maximum(abs.(vec(sum(Cg; dims = 2)) .+ vec(sum(cc.ccc[:, rest]; dims = 2)))) <= 1e-6

        # The problem is not mutated by the analysis.
        @test cq.p == frozen_control_problem().p
    end

    # ------------------------------------------------------------------
    @testset "14.2 check 1b: the particle floor and its split" begin
        ms = corea_models()
        Random.seed!(B_SEED)
        d = build_corea(tspan = (0.0, 120.0))
        Random.seed!(B_SEED)
        vr = validation_run!(ms, d, 120; moieties = corea_moieties(carbon = false))
        rep = particle_floor(vr)
        @test length(rep) == length(vr.ode_names)            # every state, as §3 asks
        @test issorted([r.min_particles for r in rep])
        @test all(r -> r.flagged == (r.min_particles < PARTICLE_FLOOR), rep)
        # 13DPG starts at ~198 particles, so it is flagged from the first handshake.
        @test :M_13dpg_c in flagged_states(rep)

        # The mutation: a state far above the floor, scaled below it, is named.
        i = findfirst(==(:M_fdp_c), vr.ode_names)
        @test !(:M_fdp_c in flagged_states(rep))
        for u in vr.ode
            u[i] *= 1e-5
        end
        @test :M_fdp_c in flagged_states(particle_floor(vr))

        # The split of spec §3 (amended 2026-09-26): a pool with a median under
        # CONTINUUM_MEDIAN particles is excluded; the three others with the
        # smallest medians are cross-checked; an external state is in neither.
        # :d is emptied once by a clip, so its minimum is small and its median
        # is not, and it is not a small pool.
        row(s, mn, md) = (species = s, min_particles = mn, t = 0.0, median_particles = md,
                          below_one = 0.0, flagged = mn < PARTICLE_FLOOR)
        synth = [row(:ext, 0.0, 50.0), row(:a, 0.1, 0.5), row(:d, 0.0, 6000.0), row(:b, 1.0, 20.0),
                 row(:c, 2.0, 5.0), row(:e, 4.0, 100.0), row(:f, 5.0, 400.0), row(:g, 900.0, 1e4)]
        sp = langevin_pools(synth; external = [:ext])
        @test sp.excluded == [:a, :c]
        @test sp.cross_check == [:b, :e, :f]
    end

    # ------------------------------------------------------------------
    @testset "14.2 check 1b: the local-linearization Langevin double" begin
        # The covariance integral against the scalar Ornstein–Uhlenbeck formula,
        # including a mode 2000 times faster than the step.
        for (λ, h) in ((-3.0, 0.7), (-2.0e4, 0.05), (0.5, 0.3))
            Σ = noise_covariance(fill(λ, 1, 1), fill(2.0, 1, 1), h)
            @test Σ[1] ≈ 2.0 * (exp(2λ * h) - 1) / (2λ) rtol = 1e-8
        end
        # And a stiff, non-normal pair against Simpson quadrature.
        J = [-1.0 5.0; 0.0 -2.0e3]
        D = [1.0 0.3; 0.3 2.0]
        h = 0.05
        n = 200_000
        s = range(0, h; length = n + 1)
        w = [k == 1 || k == n + 1 ? 1 : (iseven(k) ? 4 : 2) for k in 1:n+1] .* (h / n / 3)
        quad = sum(w[k] * exp(J * s[k]) * D * exp(J' * s[k]) for k in eachindex(s))
        @test noise_covariance(J, D, h) ≈ quad rtol = 1e-6

        ms = corea_models()
        pn = unique(Symbol[q.name for m in ms if formalism(m) === :ode
                           for q in model_free_params(parameters(m))])
        lb = LangevinBlock(ms, pn)
        @test length(lb.channels) == 42                   # every reaction, both directions
        # The stoichiometry read off the right-hand side is the published one.
        col(c) = Dict(lb.names[i] => lb.N[i, c] for i in 1:lb.n if lb.N[i, c] != 0)
        gapd = col(findfirst(==(:kcatF_R_GAPD), lb.channel_names))
        sgn = gapd[:M_13dpg_c]
        @test Dict(k => v * sgn for (k, v) in gapd) ==
              Dict(:M_g3p_c => -1, :M_nad_c => -1, :M_pi_c => -1, :M_13dpg_c => 1, :M_nadh_c => 1)
        ppa = col(findfirst(==(:kcatF_R_PPA), lb.channel_names))
        @test abs(ppa[:M_pi_c]) == 2 && abs(ppa[:M_ppi_c]) == 1

        Random.seed!(B_SEED)
        d = build_corea(tspan = (0.0, 60.0))
        Random.seed!(B_SEED)
        fp = record_frozen_path(lb, ms, d, 60; h = 0.05)
        # With the noise off the replay is the reference: the frozen path holds
        # everything the handshakes did.
        tr, _ = replay(lb, fp, MersenneTwister(1); noise = false)
        @test maximum(maximum(abs.(tr[k] .- fp.u_ref[k])) for k in eachindex(tr)) < 1e-12
        # On a species no counter touches, `a` is the integrator's own error
        # against the pinned Rodas5P. Upper glycolysis is slow and small there.
        ig = findfirst(==(:M_fdp_c), lb.names)
        @test maximum(abs(fp.a[k][ig]) / fp.u_ref[k][ig] for k in eachindex(fp.a)) < 1e-2

        b = langevin_band(lb, fp, [:M_atp_c, :M_pep_c]; n = 4, seed = 3)
        @test all(isfinite, b.lo) && all(b.lo .<= b.hi)
        @test any(b.hi .> b.lo)                             # it is noisy

        # The mutation: a reference three times the true PEP leaves the band,
        # and the excursion names the pool.
        ipep = findfirst(==(:M_pep_c), lb.names)
        biased = [(v = copy(u); v[ipep] *= 3; v) for u in fp.u_ref]
        bb = langevin_band(lb, fp, [:M_atp_c, :M_pep_c]; n = 4, seed = 3, reference = biased)
        @test bb.excursion[2] > 0.5
        @test bb.excursion[1] == b.excursion[1]

        # The partition (§12, 2026-09-26 G). Every channel conserves the seven
        # moieties, so a replay leaves them only through what its clamps
        # inject, and never by more. With the partition the noise rarely
        # clamps (4 or 5 steps per replay in job 17558769, against 368 to 390
        # without) and the carriers stay on the reference.
        Ω = fp.Ω[end]
        for st in b.stats
            @test all(st.drift .<= st.booked .+ 1e-3 / Ω)
            @test st.clamped <= 20
            @test 0.1 < st.frozen < 0.5                     # 27.5% in job 17558769
        end
        carriers = findall(in((:carrier_ptsI, :carrier_ptsH, :carrier_Crr, :carrier_ptsG)), first(b.stats).moieties)
        @test length(carriers) == 4
        @test maximum(maximum(st.drift[carriers]) for st in b.stats) * Ω < 1.0
        # The mutation: every channel noisy, the old integrator. The noise
        # clamps on hundreds of steps and the carrier totals move by whole
        # particles, 5 to 64 of them per carrier in job 17558769.
        b0 = langevin_band(lb, fp, [:M_atp_c, :M_pep_c]; n = 4, seed = 3, min_particles = 0)
        for st in b0.stats
            @test all(st.drift .<= st.booked .+ 1e-3 / Ω)
            @test st.clamped > 200
            @test st.frozen == 0
        end
        @test minimum(maximum(st.drift[carriers]) for st in b0.stats) * Ω > 10
    end

    # ------------------------------------------------------------------
    @testset "14.8 check 7: the per-counter clip record" begin
        run7(ms) = begin
            Random.seed!(B_SEED)
            d = build_problem(ms; tspan = (0.0, Float64(B_SHORT)), complete = true)
            Random.seed!(B_SEED)
            rec = ClipRecord(d)
            for _ in 1:B_SHORT
                handshake_step!(d)
                record_clips!(rec, d)
            end
            rec, d
        end
        rec, d = run7(validation_models())
        # One record per drain, and it agrees with the driver's own census.
        c = clipping_census(d)
        @test rec.drains == c.drains == B_SHORT
        @test (sum(values(rec.clipped); init = 0) > 0) == (c.clipped > 0)
        @test isempty(clip_summary(rec)) == (c.clipped == 0)
        @test :tRNA_translat in rec.counters

        # The mutation: a tRNA pool a fifth of the asserted one buffers under
        # two seconds of demand, so translation's counter clips on it, and the
        # record names that counter.
        mrec, md = run7(validation_models(charging = TrnaCharging(pool_mM = 0.05)))
        @test clipping_census(md).clipped > 0
        rows = clip_summary(mrec)
        @test !isempty(rows)
        @test any(r -> r.counter === :tRNA_translat, rows)
        r = only(r for r in rows if r.counter === :tRNA_translat)
        @test r.drains > 0 && r.first > 0 && r.max_deficit > 0
    end

    # ------------------------------------------------------------------
    @testset "14.9 check 6: the nominal trajectory refuses a doubling time" begin
        d = build_corea(tspan = (0.0, 60.0))
        err = try
            doubling_time(d); nothing
        catch e
            e
        end
        @test err isa ErrorException && occursin("refused", err.msg)
        rc = only(reporting_constraints(d))
        @test rc.quantity === :doubling_time && rc.verdict === :refused
        @test rc.instead == [:fractional_growth, :time_to_threshold]
        for _ in 1:60
            handshake_step!(d)
        end
        g = growth_report(d)
        @test g.fractional > 1.0
    end

    # ------------------------------------------------------------------
    @testset "14b.5 F5's two comparisons, on constructed inputs" begin
        @test spearman([1, 2, 3, 4], [10, 20, 30, 40]) ≈ 1
        @test spearman([1, 2, 3, 4], [4, 3, 2, 1]) ≈ -1
        @test spearman([1, 2, 2, 3], [1, 2, 2, 3]) ≈ 1       # ties averaged
        m = collect(1.0:17.0)
        ok = transcript_comparison(1.5 .* m, m)
        @test ok.pass && ok.within == 17 && ok.rho ≈ 1
        bad = transcript_comparison(reverse(m), m)
        @test !bad.pass
        far = transcript_comparison(vcat(3 .* m[1:3], m[4:end]), m)
        @test far.within == 14 && !far.pass                  # three genes out of 2×

        L = collect(300.0:100.0:1900.0)
        fc = fold_change_report(2.2 .- L ./ 1e4, L)
        @test fc.pass && fc.slope < 0
        @test !fold_change_report(fill(1.6, 17) .- L ./ 1e5, L).pass      # median below 1.7
        @test !fold_change_report(2.0 .+ L ./ 1e4, L).pass                # slope positive
        hi = 2.0 .- L ./ 1e5; hi[1] = 3.5
        @test fold_change_report(hi, L).above == 1
    end
end
