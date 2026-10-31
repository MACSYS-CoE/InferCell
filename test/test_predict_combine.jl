# KIM-40: `combine_rows`, M0's row blocks put back together. M0, every run 60 s,
# with a made-up posterior near the nominal values, as test_predict.jl's. KIM-70's
# requirements: the summaries are recomputed on the concatenated rows, never
# averaged; it is tested against one call; one DeltaCorrection is shared; and
# whether two posteriors' blocks share cell seeds.

using Test
using InferCell
using InferCell: AbstractSubModel
using Random
using Serialization
using Statistics: mean, var

function _comb_posterior(ms, n; seed = 40)
    rng = Xoshiro(seed)
    nominal = nominal_parameter_values(ms)
    X = [nominal[t] * exp(0.1 * randn(rng)) for _ in 1:n, t in M0_TARGETS]
    return PosteriorDraws(M0_TARGETS, X, 0.2 .+ 0.01 .* rand(rng, n), fill(1, n),
                          collect(1:n), nothing, nothing, (kind = :test,))
end

# `p` with some fields replaced.
function _comb_with(p::Prediction; kw...)
    f = Dict{Symbol, Any}(n => getfield(p, n) for n in fieldnames(Prediction))
    for (k, v) in kw
        f[k] = v
    end
    return Prediction((f[n] for n in fieldnames(Prediction))...)
end

# Rows `r` of a one-call prediction as a block of their own, with every per-row
# result, summary and count planted wrong, so a combination that read them fails.
function _comb_slice(p::Prediction, r)
    n = length(r)
    pv = p.provenance
    bad = (; uniform = (; mean = -1.0))
    return _comb_with(p; draws = p.draws[r], theta_I = p.theta_I[r],
                      sigma_b = p.sigma_b === nothing ? nothing : p.sigma_b[r],
                      seeds = p.seeds[:, :, r], cells = p.cells[r, :],
                      cells_intervened = p.cells_intervened[r, :], mu = fill(-1.0, n, 1),
                      sigma2 = fill(-1.0, n, 1), n_ok = fill(-1, n, 1),
                      effect = fill(-1.0, n, 1), effect_sigma2 = fill(-1.0, n, 1),
                      viable = p.viable[r, :], decisions = p.decisions[r, :],
                      viability_seeds = p.viability_seeds[:, :, r],
                      viability_cells = p.viability_cells[r, :],
                      weights = (; uniform = fill(1.0, n, 1)),
                      conditioning = (; uniform = (; acceptance = -1.0, ess = -1.0)),
                      summary = bad, effect_summary = p.effect_summary === nothing ? nothing : bad,
                      statuses = Dict(:clean => -1), statuses_intervened = Dict(:clean => -1),
                      provenance = merge(pv, (; pairs = pv.pairs[r, :],
                                              settings = merge(pv.settings,
                                                               (; n_parameter_draws = n)))))
end

@testset "KIM-40: combine_rows" begin
    H = 60.0
    ms = collect(AbstractSubModel, m0_models())
    pd = _comb_posterior(ms, 4)
    kd = Knockdown(; dose = 0.3)

    @testset "the mean field: blocks equal one call over their rows" begin
        dc = delta_correction(ms, haldane_values(metered(ms), InferCell._theta_ref(pd)),
                              [101, 102, 103]; horizon = H, intervention = kd)
        kw = (; horizon = H, forward = :mean_field, delta = dc, intervention = kd)
        one = predict(ms, pd; kw..., draws = 1:4, rng = Xoshiro(10))
        # The second block's δ is a copy, as a block reads it back from a file.
        dc2 = deserialize(seekstart(let io = IOBuffer(); serialize(io, dc); io; end))
        b1 = predict(ms, pd; kw..., draws = 1:2, rng = Xoshiro(11))
        b2 = predict(ms, pd; kw..., delta = dc2, draws = 3:4, rng = Xoshiro(12))
        c = combine_rows([b1, b2])
        @test all(isfinite, c.mu) && all(isfinite, c.effect)
        @test c.draws == one.draws && isequal(c.mu, one.mu) && isequal(c.effect, one.effect)
        @test isequal(c.summary, one.summary) && isequal(c.effect_summary, one.effect_summary)
        @test c.summary.uniform.mean ≈ mean(c.mu) + dc.delta
        @test c.delta === dc && c.statuses == one.statuses
        @test c.statuses_intervened == one.statuses_intervened
        @test c.provenance.settings == one.provenance.settings
        @test c.provenance.block_rows == [2, 2]
        @test c.provenance.fork_seed == [b1.provenance.fork_seed, b2.provenance.fork_seed]
        # One δ, shared: blocks whose δ differ, or that computed their own, are refused.
        dc3 = delta_correction(ms, haldane_values(metered(ms), InferCell._theta_ref(pd)),
                               [104, 105, 106]; horizon = H, intervention = kd)
        @test_throws ArgumentError combine_rows([b1, predict(ms, pd; kw..., delta = dc3,
                                                             draws = 3:4)])
        own = predict(ms, pd; kw..., delta = :ensemble, n_ref = 2, draws = 3:4,
                      rng = Xoshiro(13))
        @test_throws ArgumentError combine_rows([b1, own])
        @test_throws ArgumentError combine_rows([b1, predict(ms, pd; horizon = H,
                                                             forward = :mean_field,
                                                             intervention = kd, draws = 3:4)])
    end

    @testset "the hybrid: one call's rows, split and recombined" begin
        p = predict(ms, pd; horizon = H, n_cells = 3, intervention = kd, draws = 1:4,
                    rng = Xoshiro(20))
        @test all(isfinite, p.mu) && all(isfinite, p.effect)
        c = combine_rows([_comb_slice(p, 1:1), _comb_slice(p, 2:4)])
        for f in (:forward, :completion, :estimand, :quantity, :horizon, :draws, :theta_I,
                  :sigma_b, :held, :theta_H, :seeds, :cells, :mu, :sigma2, :n_ok,
                  :intervention, :cells_intervened, :effect, :effect_sigma2, :viability,
                  :viable, :decisions, :viability_seeds, :viability_cells, :weights,
                  :conditioning, :delta, :summary, :effect_summary, :statuses,
                  :statuses_intervened)
            @test isequal(getfield(c, f), getfield(p, f))
        end
        @test c.provenance.pairs == p.provenance.pairs
        @test c.provenance.settings == p.provenance.settings
        @test c.provenance.block_rows == [1, 3]
    end

    @testset "the hybrid: blocks run as separate calls" begin
        kw = (; horizon = H, n_cells = 3, intervention = kd)
        b1 = predict(ms, pd; kw..., draws = [1, 3], rng = Xoshiro(21))
        b2 = predict(ms, pd; kw..., draws = [2, 4], rng = Xoshiro(22))
        c = combine_rows([b1, b2])
        @test c.draws == [1, 3, 2, 4] && c.theta_I == vcat(b1.theta_I, b2.theta_I)
        @test c.seeds == cat(b1.seeds, b2.seeds; dims = 3) && allunique(c.seeds)
        @test isequal(c.mu, vcat(b1.mu, b2.mu)) && isequal(c.sigma2, vcat(b1.sigma2, b2.sigma2))
        @test isequal(c.effect, vcat(b1.effect, b2.effect))
        @test c.cells == vcat(b1.cells, b2.cells)
        @test sum(values(c.statuses)) == 12 && sum(values(c.statuses_intervened)) == 12
        # Equal weights: the summary by hand, over all four rows.
        mu, s2, n = vec(c.mu), vec(c.sigma2), vec(c.n_ok)
        @test all(isfinite, mu) && all(isfinite, s2) && all(==(3), n)
        s = c.summary.uniform
        @test s.mean ≈ mean(mu) && s.var ≈ var(mu) && s.mc ≈ mean(s2 ./ n)
        @test s.intrinsic_var ≈ mean(s2) && s.parameter_var ≈ max(var(mu) - mean(s2 ./ n), 0)
        @test s.n_pairs == 4 && s.ess ≈ 4
        e = c.effect_summary.uniform
        @test e.mean ≈ mean(vec(c.effect)) && e.var ≈ var(vec(c.effect))
        @test e.n_failed_intervened == 0
        # The cell seeds come from the rng whatever the draws, so another posterior
        # on the same rng and block shape runs the same seeds: common cell seeds.
        pd2 = _comb_posterior(ms, 4; seed = 41)
        q = predict(ms, pd2; kw..., draws = [1, 3], rng = Xoshiro(21))
        @test q.seeds == b1.seeds && q.theta_I != b1.theta_I
        # Refused: blocks on the same rng (their seeds repeat), and repeated draws.
        @test_throws ArgumentError combine_rows([b1, predict(ms, pd; kw..., draws = [2, 4],
                                                             rng = Xoshiro(21))])
        @test_throws ArgumentError combine_rows([b1, predict(ms, pd; kw..., draws = [3, 4],
                                                             rng = Xoshiro(23))])
        # Refused: blocks that are not one prediction's rows.
        @test_throws ArgumentError combine_rows(Prediction[])
        @test_throws ArgumentError combine_rows([b1, _comb_with(b2; completion = :prior)])
        @test_throws ArgumentError combine_rows([b1, _comb_with(b2; horizon = 2H)])
        @test_throws ArgumentError combine_rows([b1, _comb_with(b2; quantity = :exported)])
        @test_throws ArgumentError combine_rows([b1, _comb_with(b2; forward = :mean_field)])
        @test_throws ArgumentError combine_rows([b1, _comb_with(b2; intervention = Knockdown(;
                                                                       dose = 0.2))])
        pv = b2.provenance
        @test_throws ArgumentError combine_rows([b1, _comb_with(b2; provenance = merge(pv,
                                                     (; settings = merge(pv.settings,
                                                                         (; n_cells = 4)))))])
        @test_throws ArgumentError combine_rows([b1, _comb_with(b2; provenance = merge(pv,
                                                     (; source = (kind = :other,))))])
        @test_throws ArgumentError combine_rows([b1, _comb_with(b2; estimand = :other)])
        extra = [pv.skipped; (; name = :not_a_parameter, reason = :inferred)]
        @test_throws ArgumentError combine_rows([b1, _comb_with(b2; provenance = merge(pv,
                                                     (; skipped = extra)))])
        @test_throws ArgumentError combine_rows([b1, _comb_with(b2; provenance = merge(pv,
                                                     (; version = v"0.0.1")))])
    end
end
