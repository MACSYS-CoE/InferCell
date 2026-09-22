using Test
using InferCell
using StaticArrays: SA, setindex

# Spec §11 phase 9 — the lumped tRNA charging step, in the ODE block.

@testset "Core A′ tRNA charging" begin

    @testset "9.1 one mass-action reaction owning the tRNA pair" begin
        m = TrnaCharging(k_chg = 0.5, pool_mM = 0.2, charged_fraction = 0.8)

        # The state set is the registry's tRNA group, in strictly increasing
        # registry order, so a reordered or extended group cannot slip through.
        @test states(m) == species_in_group(:trna)
        @test states(m) == [:M_trna_c, :M_trna_chg_c]
        @test all(diff(species_index.(states(m))) .> 0)
        @test formalism(m) === :ode

        # Hand-checked derivative. At [trna] = 0.04 mM, [ATP] = 3.0 mM and
        # k_chg = 0.5 /(mM·s), v = 0.5 · 0.04 · 3.0 = 0.06 mM/s exactly in
        # decimal; the tRNA pair moves by ∓v and ATP, AMP, PPi by −v, +v, +v.
        u = SA[0.04, 0.16]
        p = Float64[]                       # k_chg fixed: nothing in the vector
        v = 0.5 * 0.04 * 3.0
        @test charging_flux(u, p, 0.0, m, SA[3.0]) == v
        @test dynamics(u, p, 0.0, m, SA[3.0]) == SA[-v, v]
        @test contributions(u, p, 0.0, m, SA[3.0]) == SA[-v, v, v]
        @test v ≈ 0.06 rtol = 1e-15

        # Freed, k_chg is read from the parameter vector and not the struct.
        mf = TrnaCharging(k_chg = 0.5, pool_mM = 0.2, charged_fraction = 0.8,
                          free_k_chg = true)
        @test charging_flux(u, [2.0], 0.0, mf, SA[3.0]) == 2.0 * 0.04 * 3.0

        # The initial conditions reach the composed u0 by name.
        @test [p.value for p in ic_params(parameters(m))] ≈ [0.04, 0.16] rtol = 1e-15
    end

    @testset "9.2 the pool is asserted and k_chg is derived from it" begin
        m = TrnaCharging()
        byname = Dict(p.name => p for p in parameters(m))

        # All three asserted quantities load as `:asserted`, and all three reach
        # the composed model's asserted-prior enumeration.
        for n in (:k_chg, :trna_pool_mM, :trna_charged_fraction)
            @test informedness(byname[n]) === :asserted
        end
        asserted = Set(l.subject for l in reduction_declarations([m])
                       if l.category === :asserted_prior)
        @test asserted == Set([:k_chg, :trna_pool_mM, :trna_charged_fraction])
        @test Set(p.name for p in asserted_prior_params(parameters(m))) ==
              Set([:k_chg, :trna_pool_mM, :trna_charged_fraction])

        # k_chg is computed, not typed in: recompute it here from the published
        # count, independently of the module's own constants, and check the
        # defaults against the values recorded in spec §12 (2026-09-23).
        demand_mM = (3_484_518 / 6300) / corea_particles_per_mM()
        @test demand_mM ≈ 0.027408 rtol = 1e-4
        @test byname[:k_chg].value ≈ demand_mM / (0.2 * 0.25 * 3.6529) rtol = 1e-14
        @test byname[:k_chg].value ≈ 0.15006 rtol = 1e-4
        @test charging_derivation(m).k_chg == byname[:k_chg].value

        # The initial pools: 1,009 uncharged and 4,036 charged particles.
        u0 = [p.value for p in ic_params(parameters(m))]
        @test u0 ≈ [0.05, 0.20] rtol = 1e-14
        @test u0 .* corea_particles_per_mM() ≈ [1009.02, 4036.08] rtol = 1e-5

        # Change the pool and k_chg follows: doubling the total halves it, and
        # moving the split from 0.8 to 0.9 charged doubles it.
        @test TrnaCharging(pool_mM = 0.5).k_held ≈ m.k_held / 2 rtol = 1e-14
        @test TrnaCharging(charged_fraction = 0.9).k_held ≈ 2 * m.k_held rtol = 1e-12

        # At the nominal pools the derived constant delivers the demand exactly.
        atp = species_entry(:M_atp_c).initial_value
        @test charging_flux(SA[u0...], Float64[], 0.0, m, SA[atp]) ≈ demand_mM rtol = 1e-14

        # An override is recorded in provenance, not hidden.
        o = TrnaCharging(k_chg = 0.5)
        @test occursin("override",
                       only(p for p in parameters(o) if p.name === :k_chg).provenance.identifier)

        # Only k_chg can be freed; the pool quantities are always fixed.
        @test [p.name for p in free_params(parameters(TrnaCharging(free_k_chg = true)))] == [:k_chg]
        @test isempty(free_params(parameters(m)))
        @test_throws ArgumentError TrnaCharging(charged_fraction = 1.0)
        @test_throws ArgumentError TrnaCharging(pool_mM = 0.0)
    end

    @testset "9.3 three currency edges, and they execute" begin
        m = TrnaCharging()
        edges = coupling(m)
        @test Set((edge_kind(e), e.species, e.direction) for e in edges) == Set([
            (:currency, :M_atp_c, :in),
            (:currency, :M_amp_c, :out),
            (:currency, :M_ppi_c, :out),
        ])
        @test length(edges) == 3
        @test all(e -> e isa CurrencyEdge, edges)
        # No tRNA edge: this module owns the pair, and no other module writes it
        # continuously (phase 11's debit is a deferred counter, task 11.5).
        @test !any(e -> e.species in CHARGING_STATES, edges)
        @test Set(contributed_states(m)) == Set(e.species for e in edges)

        # With phase 8: resolves, so no species and direction is both mass and
        # currency. The resolver throws on exactly that (`_check_kind_agreement`).
        @test resolve_coupling(AbstractSubModel[NucleotideRecycling(), m]) isa CouplingGraph

        # The four ODE modules together: task 8.5's deferred assertion, which
        # phase 9 is the last of phases 6 to 9 to be able to write.
        four = AbstractSubModel[CentralGlycolysis(), PtsTransport(),
                                NucleotideRecycling(), m]
        g = resolve_coupling(four)
        @test g isa CouplingGraph
        @test count(r -> r.declared_by === :TrnaCharging, g.edges) == 3

        # The contributions EXECUTE. Evaluate the composed right-hand side at the
        # initial state, and again with the uncharged pool zeroed so charging's
        # flux is zero and nothing else changes: nothing else in the composition
        # reads M_trna_c. The difference in each touched slot is the charging
        # term, and it is the flux, not zero.
        ms = AbstractSubModel[NucleotideRecycling(), HeldGlycolytic(), m]
        prob = build_problem(ms; tspan = (0.0, 1.0))
        layout = reduce(vcat, states.(ms))
        ix(s) = findfirst(==(s), layout)
        u0 = prob.u0
        uz = setindex(u0, 0.0, ix(:M_trna_c))
        d = prob.f(u0, prob.p, 0.0) - prob.f(uz, prob.p, 0.0)
        v = charging_flux(SA[u0[ix(:M_trna_c)], u0[ix(:M_trna_chg_c)]],
                          Float64[], 0.0, m, SA[u0[ix(:M_atp_c)]])
        @test v ≈ charging_demand_mM_per_s() rtol = 1e-12
        @test d[ix(:M_atp_c)] ≈ -v rtol = 1e-9
        @test d[ix(:M_amp_c)] ≈ v rtol = 1e-9
        @test d[ix(:M_ppi_c)] ≈ v rtol = 1e-9
        @test d[ix(:M_trna_c)] ≈ -v rtol = 1e-12
        @test d[ix(:M_trna_chg_c)] ≈ v rtol = 1e-12
        # And nothing it does not declare: ADP, Pi and the guanylates unmoved.
        for s in (:M_adp_c, :M_pi_c, :M_gtp_c, :M_gdp_c, :M_gmp_c)
            @test d[ix(s)] == 0.0
        end
    end

    @testset "9.4 two declarations: the lumping and the formalism" begin
        labels = [l for l in reduction_declarations([TrnaCharging()])
                  if l.category === :lumping]
        @test length(labels) == 2
        @test all(l -> l.subject === :TrnaCharging, labels)
        lump, form = labels
        # The lumping names what it replaces and why the published products won.
        @test occursin("20 per-amino-acid chains of 5 reactions", lump.description)
        @test occursin("AMP + PPi", lump.description)
        @test occursin("2 ATP -> 2 ADP + 2 Pi", lump.description)
        @test occursin("by fiat", lump.description)
        # The formalism is a separate declaration, not a clause of the lumping.
        @test occursin("deterministically", form.description)
        @test occursin("scoping note", form.description)
        @test !occursin("deterministic", lump.description)
    end

    @testset "9.5 the translation-demand double" begin
        d = TranslationDemand()
        # At the nominal charged pool it consumes the published demand, derived
        # in the double from the residue count and not read off the module.
        chg0 = CHARGING_POOL_DEFAULTS.charged_fraction * CHARGING_POOL_DEFAULTS.pool_mM
        @test dynamics(SA[0.0], Float64[], 0.0, d, SA[chg0])[1] ≈ TL_DEMAND_MM_PER_S rtol = 1e-14
        @test TL_DEMAND_PER_S ≈ 553.098 rtol = 1e-6
        # One charged consumed, one uncharged returned.
        c = contributions(SA[0.0], Float64[], 0.0, d, SA[chg0])
        @test contributed_states(d) == [:M_trna_c, :M_trna_chg_c]
        @test c[1] == -c[2] && c[1] > 0

        # Without a consumer the pool saturates and the flux collapses, which is
        # why the double is an acceptance criterion and not scaffolding.
        ms = charging_models(demand = nothing)
        sol = recycling_solve(ms; horizon = 600.0)
        chg, unc = layout_index(ms, :M_trna_chg_c), layout_index(ms, :M_trna_c)
        frac_end = sol.u[end][chg] / (sol.u[end][chg] + sol.u[end][unc])
        @test frac_end > 0.999
        # With it, the split stays near nominal.
        ms2 = charging_models()
        sol2 = recycling_solve(ms2; horizon = 600.0)
        frac2 = sol2.u[end][chg] / (sol2.u[end][chg] + sol2.u[end][unc])
        @test 0.5 < frac2 < 0.95
        @info "9.5" frac_end_no_consumer = frac_end frac_end_with_double = frac2
    end

    # ------------------------------------------------------------------
    # 9.6 — this module's own check. Spec §3's exception, gate 1 with the
    # tol_C ladder (amended 2026-09-23, §12): asserted, not argued.
    # ------------------------------------------------------------------
    @testset "9.6 tRNA conservation over a full cycle" begin
        ms = charging_models()
        tr = trna_total(ms)
        U, C = layout_index(ms, :M_trna_c), layout_index(ms, :M_trna_chg_c)
        sol = recycling_solve(ms)
        @test sol.t[end] == CYCLE_S
        total = tr(sol.u[1])
        @test total ≈ CHARGING_POOL_DEFAULTS.pool_mM rtol = 1e-15

        # The flux is live, so the check is not passed by a pool at rest: over
        # the cycle the pair turns over ~ 3.5 million residues' worth.
        CNT = layout_index(ms, :chg_residues_mM)
        @test sol.u[end][CNT] * corea_particles_per_mM() > 3.0e6

        # Gate 1's criterion: the composed right-hand side returns the pair's
        # summed derivative as bitwise 0.0 at every save point. Not `≈`.
        prob = build_problem(ms; tspan = (0.0, CYCLE_S))
        wsum(u) = (d = prob.f(u, prob.p, 0.0); d[U] + d[C])
        @test all(u -> wsum(u) === 0.0, sol.u)

        # The ladder, six decades, each rung at least three orders below its
        # own tol_C and the largest within 100x of the smallest. Measured
        # 2026-09-23: 10.2 down to 5.7 orders, spread 51x. The flat 100-ulp
        # bound is not used: the residual is ~1e-13 mM of solver roundoff set by
        # the composition's larger states, 60 to 3,097 ulps of this small sum.
        rungs = [(1e-4, 1e-2), (1e-5, 1e-3), (1e-6, 1e-4), (1e-7, 1e-5),
                 (1e-8, 1e-6), (1e-9, 1e-7), (1e-10, 1e-8)]
        residuals = Float64[]
        for (a, r) in rungs
            l = (a, r) == (ABSTOL_R, RELTOL_R) ? sol :
                recycling_solve(ms; abstol = a, reltol = r)
            push!(residuals, max_drift(l, tr))
            @test residuals[end] <
                  tolerance_bound(l, ((U, 1), (C, 1)); abstol = a, reltol = r) / 1e3
        end
        @test log10(rungs[1][1] / rungs[end][1]) >= 6
        @test maximum(residuals) < 100 * minimum(residuals)
        @info "9.6 tRNA ladder" residuals ulps = residuals ./ eps(total)

        # The mutation: charging that creates tRNA rather than transferring it
        # fails both halves, by orders of magnitude rather than by a factor.
        mms = charging_models(charging = MutatedCharging(TrnaCharging()), counter = false)
        mtr = trna_total(mms)
        mprob = build_problem(mms; tspan = (0.0, CYCLE_S))
        mU, mC = layout_index(mms, :M_trna_c), layout_index(mms, :M_trna_chg_c)
        md = mprob.f(mprob.u0, mprob.p, 0.0)
        @test abs(md[mU] + md[mC]) > 1e6 * eps(total)
        msol = recycling_solve(mms; horizon = 600.0)
        @test max_drift(msol, mtr) > 1e6 * maximum(residuals)
        @info "9.6 mutation" per_eval = md[mU] + md[mC] drift_600s = max_drift(msol, mtr)
    end

end
