# spec/phases/approximate-method.md task B4: route (ii), the mean field and the
# metabolism as one ODE (§12 2026-10-23). The toys are phase 3–5's doubles from
# hybrid_test_models.jl, whose closed forms the assertions use.

using Test
using InferCell
using ForwardDiff
using StaticArrays
using Random
using LinearAlgebra: norm
using OrdinaryDiffEq: Rodas5P

# The toy's derivative written out by hand: transcription at the rebuilt k, decay,
# translation charging `cost` ATP particles per protein, and the pool's
# Michaelis–Menten reaction at the protein's concentration.
function _toy_route_rhs(u, k_tx, m::ToyRebuiltExpression, pool::ToyPool, factor)
    mrna, prot, cost_c, atp, adp = u
    γ, ktl = m.params[4].value, m.params[5].value
    kcat, km = pool.params[1].value, pool.params[2].value
    v = kcat * (prot / factor) * atp / (km + atp)
    accrual = m.cost * ktl * mrna
    return [k_tx - γ * mrna, ktl * mrna, accrual, -v - accrual / factor, v]
end

@testset "B4: route (ii)" begin
    factor = InferCell.corea_particles_per_mM()

    @testset "catalytic, drain and rebuild on the toy" begin
        expr = ToyRebuiltExpression(; k_tx = 2.0, mrna0 = 3, protein0 = 40)
        pool = toy_slow_pool()
        models = AbstractSubModel[pool, expr]
        rt = build_mean_field_problem(models; tspan = (0.0, 200.0))
        @test rt.state_names == [:M_atp_c, :M_adp_c, :toy_mrna, :M_ptsg_c, :atp_cost]
        @test rt.n_idx == 0
        @test rt.rebuild_times == [59.0, 119.0, 179.0]
        @test Set(w.param_slot for w in mean_field_written_params(rt)) ==
              Set([:enzyme_conc, :k_tx_rb])
        @test isempty(unexecuted_edges(rt.models, rt))

        # The right-hand side, at an arbitrary state, against the hand-written one:
        # the catalytic count and the continuous drain.
        p = copy(rt.prob.p)
        ik = findfirst(==(:k_tx_rb), rt.param_names)
        p[ik] = 1.7
        u = [5.0, 0.4, 7.0, 123.0, 40.0]
        du = similar(u)
        rt.prob.f(du, u, p, 0.0)
        perm = [3, 4, 5, 1, 2]           # route order → the hand-written order
        want = _toy_route_rhs(u[perm], 1.7, expr, pool, factor)
        @test du[perm] ≈ want rtol = 1e-12

        sol = solve_mean_field(rt; saveat = [0.0, 59.0, 179.0, 200.0], abstol = 1e-12,
                               reltol = 1e-12)
        ia, id, ic = 1, 2, 5
        # ATP is conserved with ADP and the counter, in particles.
        total(u) = (u[ia] + u[id]) * factor + u[ic]
        @test all(total(x) ≈ total(sol.u[1]) for x in sol.u)

        # Window 1 runs at R₀, from the initial pools: the transcript is the
        # closed form at that constant until 59 s.
        k0 = toy_rebuilt_k(expr, [pool.params[4].value])
        γ = expr.params[4].value
        @test sol.u[2][3] ≈ k0 / γ + (3 - k0 / γ) * exp(-γ * 59) rtol = 1e-8
        # The last switch read the pools at 179 s, and nothing after it did.
        @test sol.prob.p[ik] ≈ toy_rebuilt_k(expr, [sol.u[3][ia]]) rtol = 1e-12

        # A re-solve at another θ recomputes R₀; the problem's own p is untouched.
        before = copy(rt.prob.p)
        p2 = mean_field_parameters(rt, [:k_tx_max_rb => 9.0])
        sol2 = solve_mean_field(rt; p = p2, saveat = [59.0], abstol = 1e-12, reltol = 1e-12)
        k0′ = 9.0 * pool.params[4].value / (expr.params[3].value + pool.params[4].value)
        @test sol2.u[1][3] ≈ k0′ / γ + (3 - k0′ / γ) * exp(-γ * 59) rtol = 1e-8
        @test rt.prob.p == before
        @test solve_mean_field(rt; saveat = [200.0]).u[end] ==
              solve_mean_field(rt; saveat = [200.0]).u[end]
        @test_throws ArgumentError mean_field_parameters(rt, [:k_tx_rb => 1.0])
        @test_throws ArgumentError mean_field_parameters(rt, [:no_such => 1.0])

        subjects = [l.subject for l in route_declarations(rt)]
        @test :continuous_drain in subjects && :rebuild_snapshot in subjects &&
              :unrounded_counts in subjects
    end

    @testset "growth: counts conserved under dilution, the volume by the law" begin
        n0 = 353
        expr = ToyExpression(; k_tx = 4.0, mrna0 = 2, protein0 = 50, cost = 1,
                             edges = CouplingEdge[DeferredCounterEdge(species = :M_ptsi_c,
                                                                      direction = :out,
                                                                      counter = :atp_cost)])
        pool = ToyMembranePool(; ptsi_particles = n0)
        rt = build_mean_field_problem(AbstractSubModel[pool, expr]; tspan = (0.0, 600.0))
        @test rt.n_idx == length(rt.prob.u0)
        @test rt.dilute_idxs == [1, 2, 3]            # ATP, ADP, ptsI; not lac_e
        sol = solve_mean_field(rt; saveat = 0.0:100.0:600.0, abstol = 1e-12, reltol = 1e-12)
        names = rt.state_names
        i(s) = findfirst(==(s), names)
        f0 = _toy_factor0()
        for u in sol.u
            g = toy_growth(u[rt.n_idx], n0)
            fac = 1e-3 * g.volume_litres * 6.02214076e23
            # The membrane state's particles are its credits, and N tracks them.
            @test u[i(:M_ptsi_c)] * fac ≈ n0 + u[i(:atp_cost)] rtol = 1e-9
            @test u[rt.n_idx] ≈ n0 + u[i(:atp_cost)] rtol = 1e-9
            # ATP and ADP do not react (kcat = 0), so their particles are fixed.
            @test u[i(:M_atp_c)] * fac ≈ 3.6529 * f0 rtol = 1e-9
            @test u[i(:M_lac__L_e)] == 0.5
        end
        @test toy_growth(sol.u[end][rt.n_idx], n0).fractional > 1.0
        @test :continuous_growth in [l.subject for l in route_declarations(rt)]
    end

    @testset "an overdrawn pool stops the solve" begin
        # A large cost on a small, unreplenished pool.
        expr = ToyRebuiltExpression(; k_tx = 2.0, mrna0 = 50, cost = 10_000)
        rt = build_mean_field_problem(AbstractSubModel[ToyPool(; kcat = 0.0), expr];
                                      tspan = (0.0, 600.0))
        err = try
            solve_mean_field(rt; saveat = [600.0])
            nothing
        catch e
            e
        end
        @test err isa PoolOverdraw && occursin("went below zero", sprint(showerror, err))
    end

    @testset "refusals" begin
        # Two modules declaring the slot the catalytic edge fills, as in the hybrid.
        @test_throws ArgumentError build_mean_field_problem(
            AbstractSubModel[ToyPool(), ToyPool2(), ToyExpression()])
        # A rate-constant edge with nothing rebuilt from it.
        @test_throws ArgumentError build_mean_field_problem(
            AbstractSubModel[toy_slow_pool(), ToyRebuiltExpression(; rebuilt = Symbol[])])
        # Nothing to swap.
        @test_throws ArgumentError build_mean_field_problem(AbstractSubModel[ToyPool2()])
        # A membrane protein on the jump side has no count route (ii) can grow from.
        @test_throws ArgumentError build_mean_field_problem(
            AbstractSubModel[ToyMembranePool(; membrane = Symbol[],
                                             edges = CouplingEdge[CatalyticEdge(
                                                 species = :M_ptsg_c, direction = :in,
                                                 param_slot = :enzyme_conc)]),
                             ToyGrowingExpression()])
    end

    @testset "M0" begin
        ms = m0_models()
        rt = build_mean_field_problem(ms; tspan = (0.0, M0_HORIZON_S), complete = true)
        @test count(m -> m isa MeanField, rt.models) == 1
        @test rt.n_idx > 0
        @test rt.rebuild_times == collect(59.0:60.0:599.0)
        written = Set(w.param_slot for w in mean_field_written_params(rt))
        @test written == Set(vcat([c.param_slot for c in rt.catalytic], [:r_cell_nm],
                                  rebuilt_params(last(rt.models))))
        @test isempty(intersect(written, Set(M0_TARGETS)))
        # The plain all-ODE build still refuses the catalytic edges, and says
        # where route (ii) is.
        err = try
            build_problem(mean_field_swap(ms))
            nothing
        catch e
            e
        end
        @test err isa ArgumentError && occursin("build_mean_field_problem", err.msg)

        # AV3: ForwardDiff through the solve and its callbacks against central
        # differences, on the panel's log particle counts at every save, in log θ.
        truth = draw_truth(Xoshiro(16_085), ms; names = M0_TARGETS, purpose = :calibration).values
        p0 = mean_field_parameters(rt, truth)
        base = Dict(truth)
        rows = [findfirst(==(o), rt.state_names) for o in M0_PANEL]
        function obs(x)
            p = mean_field_parameters(rt, [n => base[n] * exp(x[j]) for (j, n) in enumerate(M0_TARGETS)];
                                      p = p0)
            sol = solve_mean_field(rt; p, saveat = 60.0:60.0:600.0, abstol = 1e-14, reltol = 1e-12)
            return reduce(vcat, [log.(u[rows]) for u in sol.u])
        end
        x0 = zeros(length(M0_TARGETS))
        J = ForwardDiff.jacobian(obs, x0)
        h = 1e-4
        fd = reduce(hcat, map(eachindex(x0)) do j
            e = zeros(length(x0))
            e[j] = h
            (obs(x0 + e) - obs(x0 - e)) / 2h
        end)
        @test all(isfinite, J)
        for j in eachindex(x0)
            @test norm(J[:, j]) > 0
            @test J[:, j] ≈ fd[:, j] rtol = 1e-6
        end
    end
end
