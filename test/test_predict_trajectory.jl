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
