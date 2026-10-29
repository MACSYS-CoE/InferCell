# KIM-27: predict, end to end at smoke size. M0 for the forward descriptions,
# δ and the intervention; the reduced model (d11_models) for completion, since
# M0 holds no parameter. Every run is 60 s.

using Test
using InferCell
using InferCell: AbstractSubModel
using Random
using Statistics: mean, var

# A made-up posterior over M0's targets, near their nominal values.
function _toy_posterior(ms, n; seed = 27)
    rng = Xoshiro(seed)
    nominal = nominal_parameter_values(ms)
    X = [nominal[t] * exp(0.1 * randn(rng)) for _ in 1:n, t in M0_TARGETS]
    return PosteriorDraws(M0_TARGETS, X, 0.2 .+ 0.01 .* rand(rng, n), fill(1, n),
                          collect(1:n), nothing, nothing, (kind = :test,))
end

@testset "KIM-27: predict" begin
    H = 60.0
    ms = collect(AbstractSubModel, m0_models())
    pd = _toy_posterior(ms, 3)

    @testset "hybrid on M0: the population mean and its parts" begin
        p = predict(ms, pd; horizon = H, n_cells = 2, rng = Xoshiro(1))
        @test p isa Prediction && size(p.mu) == (3, 1) && p.draws == 1:3
        @test p.theta_I == [draw_values(pd, m) for m in 1:3]
        @test p.sigma_b == pd.sigma_b                     # carried, not simulated
        @test isempty(p.held) && p.theta_H == [Pair{Symbol, Float64}[]]
        @test size(p.seeds) == (2, 1, 3) && allunique(p.seeds)   # fresh cells per draw
        @test all(isfinite, p.mu) && all(isfinite, p.sigma2) && all(==(2), p.n_ok)
        @test p.mu[2] == mean(c.q for c in p.cells[2])
        s = p.summary.uniform
        @test s.mean ≈ mean(p.mu)
        @test s.var ≈ var(vec(p.mu))
        @test s.mc ≈ mean(p.sigma2 ./ 2)
        @test s.parameter_var == max(s.var - s.mc, 0) && s.flagged == (s.var < s.mc)
        @test s.intrinsic_var ≈ mean(p.sigma2) && s.n_pairs == 3 && s.ess ≈ 3
        @test !s.mc_undefined
        @test sum(values(p.statuses)) == 6
        @test p.effect_summary === nothing && p.delta === nothing
    end

    @testset "mean field on M0, uncorrected and δ-corrected" begin
        m = predict(ms, pd; horizon = H, forward = :mean_field, rng = Xoshiro(2))
        @test size(m.seeds) == (0, 1, 3) && all(isnan, m.sigma2)
        @test m.summary.uniform.mc == 0 && isnan(m.summary.uniform.intrinsic_var)
        c = predict(ms, pd; horizon = H, forward = :mean_field, delta = :ensemble, n_ref = 2,
                    rng = Xoshiro(2))
        @test c.delta isa DeltaCorrection && length(c.delta.seeds) == 2
        @test c.mu == m.mu                                  # the uncorrected m(θ) is kept
        @test c.summary.uniform.mean ≈ mean(c.mu) + c.delta.delta
        # A correction computed once is reused, and checked against the prediction.
        r = predict(ms, pd; horizon = H, forward = :mean_field, delta = c.delta)
        @test r.summary.uniform.mean ≈ c.summary.uniform.mean
        @test_throws ArgumentError predict(ms, pd; horizon = 2H, forward = :mean_field,
                                           delta = c.delta)
        @test_throws ArgumentError predict(ms, pd; horizon = H, forward = :mean_field,
                                           delta = c.delta, quantity = :exported)
    end

    @testset "a zero-dose knockdown through predict: ΔQ ≡ 0" begin
        p = predict(ms, pd; horizon = H, n_cells = 2, intervention = Knockdown(; dose = 0.0),
                    n_parameter_draws = 2, rng = Xoshiro(3))
        @test all(==(0.0), p.effect) && p.effect_summary.uniform.mean == 0
        @test [c.seed for c in p.cells_intervened[1]] == [c.seed for c in p.cells[1]]
    end

    @testset "completion on the reduced model, crossed" begin
        dms = collect(AbstractSubModel, d11_models())
        p = predict(dms, pd; horizon = H, completion = :prior, n_held = 2,
                    n_parameter_draws = 2, width = 2.0, rng = Xoshiro(4))
        @test size(p.mu) == (2, 2) && length(p.theta_H) == 2
        @test p.held == held_parameters(metered(dms), pd).held && length(p.held) == 15
        @test p.theta_H[1] != p.theta_H[2] && first.(p.theta_H[1]) == p.held
        @test allunique(p.seeds) && size(p.seeds) == (1, 2, 2)
        @test all(isfinite, p.mu)
        @test p.summary.uniform.n_pairs == 4
        # One cell per pair: no σ̂², so the MC term is undefined, and says so.
        @test p.summary.uniform.mc_undefined && isnan(p.summary.uniform.parameter_var)

        # A rule nothing passes: every pair is rejected, and keeps its reason.
        never = ViabilityRule(; label = "toy: never", allow_clipped = true,
                              min_viable_fraction = 1, floors = [:M_gtp_c => 1e15])
        v = predict(dms, pd; horizon = H, completion = :viability_conditioned_prior,
                    viability = never, n_held = 2, n_parameter_draws = 2, width = 2.0,
                    rng = Xoshiro(4), intervention = Knockdown(; dose = 0.3))
        @test !any(v.viable) && v.conditioning.joint.acceptance == 0
        @test v.conditioning.per_theta_I.failed_rows == [1, 2]
        @test all(d -> d.counts[:below_floor] == 1, v.decisions)
        @test v.summary.joint.n_pairs == 0 && isnan(v.summary.joint.mean)
        @test all(isempty, v.cells_intervened)              # rejected pairs are not run
        @test v.mu == p.mu                                   # same rng, same draws and cells
        # A rule everything passes weighs every pair.
        always = ViabilityRule(; label = "toy: always", allow_clipped = true,
                               min_viable_fraction = 0.01)
        a = predict(dms, pd; horizon = H, completion = :viability_conditioned_prior,
                    viability = always, n_held = 2, n_parameter_draws = 2, width = 2.0,
                    rng = Xoshiro(4))
        @test all(a.viable) && a.conditioning.joint.acceptance == 1
        @test sum(a.weights.per_theta_I) ≈ 1 && sum(a.weights.joint) ≈ 1
        @test a.summary.joint.mean ≈ p.summary.uniform.mean
    end

    @testset "refusals, before any work" begin
        bad(; kw...) = @test_throws ArgumentError predict(ms, pd; horizon = H, kw...)
        bad(; forward = :ode)
        bad(; completion = :joint)
        bad(; estimand = :median)
        bad(; quantity = :growth)
        bad(; viability = ViabilityRule(; label = "x", allow_clipped = true,
                                        min_viable_fraction = 1))
        bad(; completion = :viability_conditioned_prior)
        bad(; delta = :ensemble)                       # δ is for the mean field
        bad(; forward = :mean_field, delta = :ens)
        bad(; held = [:S_JCVISYN3A_0445])              # nothing is held without completion
        bad(; width = 2.0)
        bad(; n_parameter_draws = 4)
        bad(; draws = [1, 9])
        bad(; completion = :prior, held = [:kcatR_R_ENO])   # derived, not drawn
        bad(; theta_ref = draw_values(pd, 1))
        bad(; intervention = Knockdown(; param = :enz_R_ENO, dose = 0.5))   # driver-written
        stray = PosteriorDraws([:no_such], ones(2, 1), nothing, [1, 1], [1, 2], nothing,
                               nothing, (kind = :test,))
        @test_throws ArgumentError predict(ms, stray; horizon = H)
        # A posterior value for a reverse constant would be replaced by the
        # derived one, so a posterior naming one is refused, not quietly overridden.
        rev = PosteriorDraws(vcat(M0_TARGETS, :kcatR_R_ENO), hcat(pd.values, fill(123.0, 3)),
                             nothing, [1, 1, 1], [1, 2, 3], nothing, nothing, (kind = :test,))
        @test_throws ArgumentError predict(ms, rev; horizon = H)
        @test (name = :kcatR_R_ENO, reason = :haldane_derived) in held_parameters(ms, rev).skipped
    end
end
