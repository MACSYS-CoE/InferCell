# KIM-86: `predict` records time series of chosen ODE-block species, and the
# default path is unchanged. M0 at 120 s; the reference guards the default path
# (test_predict_reference.jl).

using Test
using InferCell
using InferCell: AbstractSubModel, LAC_EXPORT_METER
using Random
using Statistics: quantile

# Three draws of M0's targets near their nominal values, each with a σ_b.
function _trajectory_test_posterior(; seed = 86)
    ms = collect(AbstractSubModel, m0_models())
    nominal = nominal_parameter_values(ms)
    rng = Xoshiro(seed)
    X = [nominal[t] * exp(0.05 * randn(rng)) for _ in 1:3, t in M0_TARGETS]
    return PosteriorDraws(M0_TARGETS, X, [0.2, 0.21, 0.22], [1, 1, 1], [1, 2, 3], nothing,
                          nothing, (kind = :test,))
end

@testset "KIM-86: predict records trajectories" begin
    ms = collect(AbstractSubModel, m0_models())
    pd = _trajectory_test_posterior()
    M = LAC_EXPORT_METER
    S = [M, :M_gtp_c]

    @testset "recording leaves the run alone, and reads particles on the grid" begin
        a = predict(ms, pd; horizon = 120.0, n_cells = 2, rng = Xoshiro(7), threads = false)
        b = predict(ms, pd; horizon = 120.0, n_cells = 2, rng = Xoshiro(7), threads = false,
                    trajectory = S, every = 60.0)
        @test isequal(a.mu, b.mu) && a.seeds == b.seeds
        @test all(c -> c.trajectory === nothing, reduce(vcat, vec(a.cells)))
        c = b.cells[1, 1][1]
        @test c.trajectory.t == [0.0, 60.0, 120.0] && c.trajectory.species == S
        @test size(c.trajectory.values) == (3, 2)
        @test c.trajectory.values[end, 1] - c.trajectory.values[1, 1] ≈ c.q rtol = 1e-12
        @test b.provenance.settings.trajectory == S && b.provenance.settings.every == 60.0
        @test a.provenance.settings.trajectory == Symbol[]

        # A horizon that `every` does not divide still ends the grid at the horizon.
        e = predict(ms, pd; horizon = 90.0, draws = [1], n_cells = 1, rng = Xoshiro(7),
                    threads = false, trajectory = S)
        @test e.cells[1, 1][1].trajectory.t == [0.0, 60.0, 90.0]
    end

    @testset "refused before any cell runs" begin
        calls = Ref(0)
        counting = run -> (calls[] += 1; 1.0)
        bad(; kw...) = @test_throws ArgumentError predict(ms, pd; horizon = 120.0,
                                                          quantity = counting, kw...)
        bad(; trajectory = [:M_nope])
        @test occursin("M_nope", sprint(showerror, try
            predict(ms, pd; horizon = 120.0, trajectory = [:M_nope])
        catch err
            err
        end))
        bad(; trajectory = [transcript_state(:JCVISYN3A_0607)])
        bad(; trajectory = Symbol[])
        bad(; trajectory = S, every = 0.7)
        bad(; trajectory = S, every = 0.0)
        bad(; trajectory = S, every = -60.0)
        bad(; every = 120.0)                 # with nothing to record
        @test calls[] == 0
    end

    @testset "the mean field, and the other two entry points" begin
        mf = predict(ms, pd; horizon = 120.0, forward = :mean_field, trajectory = S,
                     threads = false)
        c = mf.cells[1, 1][1]
        @test c.trajectory.t == [0.0, 60.0, 120.0]
        @test c.trajectory.values[end, 1] - c.trajectory.values[1, 1] ≈ c.q rtol = 1e-12
        v = haldane_values(metered(ms), draw_values(pd, 1))
        fc = forward_cells(ms, v, [11]; horizon = 60.0, threads = false, trajectory = S)
        @test only(fc).trajectory.t == [0.0, 60.0]
        @test forward_cells(ms, v, [11]; horizon = 60.0, threads = false)[1].trajectory === nothing
        @test_throws ArgumentError forward_cells(ms, v, [11]; horizon = 60.0, trajectory = [:M_nope])
        ic = intervention_cells(ms, draw_values(pd, 1),
                                Knockdown(; param = :kcatF_R_ENO, dose = 0.5), [11];
                                horizon = 60.0, threads = false, trajectory = S)
        @test only(ic.baseline).trajectory.species == S
        @test only(ic.intervened).trajectory.species == S
        @test_throws ArgumentError intervention_cells(ms, draw_values(pd, 1),
                                                      Knockdown(; param = :kcatF_R_ENO, dose = 0.5),
                                                      [11]; horizon = 60.0, every = 0.7,
                                                      trajectory = S)
    end

    @testset "a cell that stalls holds no trajectory, and the others keep theirs" begin
        j = findfirst(==(:kcatF_R_ENO), pd.names)
        X = copy(pd.values)
        X[2, j] = -X[2, j]
        sp = PosteriorDraws(pd.names, X, pd.sigma_b, pd.chain, pd.iteration, nothing, nothing,
                            (kind = :test_stall,))
        s = predict(ms, sp; horizon = 60.0, n_cells = 2, rng = Xoshiro(68), threads = false,
                    trajectory = S)
        @test all(c -> c.failure === :not_advanced && c.trajectory === nothing, s.cells[2, 1])
        @test all(c -> c.trajectory !== nothing, vcat(s.cells[1, 1], s.cells[3, 1]))
    end
end

@testset "KIM-86: bands and the split over time" begin
    M = LAC_EXPORT_METER
    pd = _trajectory_test_posterior()

    @testset "the weighted quantile is type 7 at equal weights" begin
        x = [1.0, 2.0, 3.0, 10.0]
        for q in (0.0, 0.05, 0.5, 0.9, 1.0)
            @test InferCell._weighted_quantile(x, ones(4), q) ≈ quantile(x, q)
            @test InferCell._weighted_quantile(x, fill(0.3, 4), q) ≈ quantile(x, q)
        end
        @test InferCell._weighted_quantile(x, [1.0, 1.0, 1.0, 0.0], 0.5) ≈ 2.0
        @test InferCell._weighted_quantile([5.0], [2.0], 0.9) == 5.0
    end

    ms = d11_models()
    p = predict(ms, pd; horizon = 120.0, completion = :prior, n_held = 2, n_parameter_draws = 3,
                n_cells = 2, rng = Xoshiro(3), threads = false, trajectory = [M, :M_gtp_c])

    @testset "the three layers end at what predict reports" begin
        b = trajectory_bands(p, M; baseline = true)
        @test b.t[end] == 120.0 && b.n_failed == 0
        @test all(==(0.0), (b.pairs.mean[1], b.rows.mean[1], b.cells.mean[1]))
        @test b.pairs.mean[end] ≈ p.summary.uniform.mean rtol = 1e-12
        @test b.pairs.lo[0.9][end] ≈ quantile(vec(p.mu), 0.05) rtol = 1e-12
        @test b.pairs.hi[0.9][end] ≈ quantile(vec(p.mu), 0.95) rtol = 1e-12
        rowmeans = vec(sum(p.mu; dims = 2) ./ 2)
        @test b.rows.lo[0.9][end] ≈ quantile(rowmeans, 0.05) rtol = 1e-12
        qs = [c.q for cs in vec(p.cells) for c in cs]
        @test b.cells.hi[0.9][end] ≈ quantile(qs, 0.95) rtol = 1e-12
        for L in (b.pairs, b.rows, b.cells)
            @test all(L.lo[0.9] .<= L.lo[0.5] .<= L.hi[0.5] .<= L.hi[0.9])
        end
        # A species' level, unbaselined, starts where it starts.
        g = trajectory_bands(p, :M_gtp_c; levels = (0.8,))
        @test g.pairs.mean[1] > 0 && collect(keys(g.pairs.lo)) == [0.8]
    end

    @testset "the split over time ends at variance_split" begin
        s = trajectory_split(p, M; baseline = true)
        @test s.t[1] > 0 && s.t[end] == 120.0
        v = variance_split(p)
        # At two cells per pair the noise can swamp 120 s of spread, and clip the total.
        @test isapprox(s.held_share.i_first[end], v.held_share.i_first; rtol = 1e-10, nans = true)
        @test isapprox(s.held_share.h_first[end], v.held_share.h_first; rtol = 1e-10, nans = true)
        @test s.total[end] ≈ v.total rtol = 1e-10 atol = 1e-6

        # The mean field has no noise, so its total cannot clip to nothing. On M0 no
        # parameter is held, so the held share is 0. The reduced model's route would
        # add minutes of compile to the suite for the same arithmetic.
        mf = predict(collect(AbstractSubModel, m0_models()), pd; horizon = 120.0,
                     forward = :mean_field, threads = false, trajectory = [M])
        sm, vm = trajectory_split(mf, M; baseline = true), variance_split(mf)
        @test vm.total > 0
        @test sm.total[end] ≈ vm.total rtol = 1e-10
        @test sm.held_share.i_first[end] == vm.held_share.i_first == 0
        @test length(sm.t) == 2 && sm.total[1] != sm.total[2]
    end

    @testset "weights: a pair with none is absent, and a set that weighs nothing is refused" begin
        w = [1.0 0.0; 1.0 0.0; 1.0 0.0]
        b = trajectory_bands(p, M; baseline = true, weights = w)
        @test b.pairs.mean[end] ≈ sum(p.mu[:, 1]) / 3 rtol = 1e-12
        @test b.pairs.hi[0.9][end] ≈ quantile(p.mu[:, 1], 0.95) rtol = 1e-12
        @test_throws ArgumentError trajectory_bands(p, M; weights = ones(2, 2))
        @test_throws ArgumentError trajectory_bands(p, M; weights = :joint)
        @test_throws ArgumentError trajectory_bands(p, M; weights = zeros(3, 2))
        @test_throws ArgumentError trajectory_split(p, M; weights = :joint)

        never = ViabilityRule(; label = "floor", allow_clipped = true, min_viable_fraction = 0.5,
                              floors = [:M_gtp_c => 1e12])
        r = predict(ms, pd; horizon = 120.0, completion = :viability_conditioned_prior,
                    viability = never, n_viability_cells = 0, n_held = 2, n_parameter_draws = 3,
                    n_cells = 2, rng = Xoshiro(3), threads = false, trajectory = [M])
        @test_throws ArgumentError trajectory_bands(r, M)
        @test trajectory_bands(r, M; weights = ones(3, 2)).t[end] == 120.0
    end

    @testset "what cannot be read is refused" begin
        e = predict(ms, pd; horizon = 120.0, completion = :prior, n_parameter_draws = 3,
                    n_cells = 2, rng = Xoshiro(3), threads = false)
        @test_throws ArgumentError trajectory_bands(e, M)
        @test_throws ArgumentError trajectory_split(e, M)
        @test_throws ArgumentError trajectory_bands(p, :M_gtp_c_not_recorded)
    end

    @testset "M0: failed cells are skipped and counted, and row blocks keep trajectories" begin
        m0 = collect(AbstractSubModel, m0_models())
        j = findfirst(==(:kcatF_R_ENO), pd.names)
        X = copy(pd.values)
        X[2, j] = -X[2, j]
        sp = PosteriorDraws(pd.names, X, pd.sigma_b, pd.chain, pd.iteration, nothing, nothing,
                            (kind = :test_stall,))
        s = predict(m0, sp; horizon = 60.0, n_cells = 2, rng = Xoshiro(68), threads = false,
                    trajectory = [M])
        b = trajectory_bands(s, M; baseline = true)
        @test b.n_failed == 2
        @test b.pairs.hi[0.9][end] ≈ quantile(filter(isfinite, s.mu), 0.95) rtol = 1e-12

        bl = [predict(m0, pd; horizon = 60.0, n_cells = 2, draws = [i], rng = Xoshiro(i),
                      threads = false, trajectory = [M]) for i in 1:2]
        @test combine_rows(bl).cells[2, 1][1].trajectory !== nothing
    end

    @testset "the release recipe: a fixed constant joins the held set, its reverse derived" begin
        rel = d11_models()
        rel[1] = CentralGlycolysis(enzymes = :translated,
            free = [:kcatF_R_ENO, :kcatR_R_ENO, :kcatF_R_FBA, :kcatR_R_FBA,
                    :kcatF_R_LDH_L, :kcatR_R_LDH_L])
        h = held_parameters(metered(rel), pd)
        @test :kcatF_R_LDH_L in h.held && !(:kcatR_R_LDH_L in h.held)
        @test any(s -> s.name === :kcatR_R_LDH_L && s.reason === :haldane_derived, h.skipped)
    end
end
