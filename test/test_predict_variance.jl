# KIM-66: the variance split on the crossed matrix, and the MC term under unequal
# weights, on made-up matrices with a known split. No cell runs here; the split's
# bias under V is measured by dev/scripts/variance_split_check_kim66.jl.

using Test
using InferCell
using Random
using Statistics: mean, std, var

# A crossed matrix: θ_I = x down the rows, θ_H = y across, both standard normal,
# μ = a x + b y + g x y, and each μ̂ the mean of n cells with noise sd s. The
# split is known: Var_I[E(μ | x)] = a², Var_H[E(μ | y)] = b², interaction g².
function _crossed(rng, M, K; a = 1.0, b = 0.7, g = 0.5, n = 4, s = 1.0)
    x, y = randn(rng, M), randn(rng, K)
    cells = [a * x[m] + b * y[k] + g * x[m] * y[k] .+ s .* randn(rng, n) for m in 1:M, k in 1:K]
    return mean.(cells), var.(cells) ./ n
end

@testset "KIM-66: the variance split and the MC term" begin
    @testset "unmasked: unbiased for the known split, cells' noise removed" begin
        rng = Xoshiro(66)
        M, K, R = 20, 10, 400
        a, b, g = 1.0, 0.7, 0.5
        est = map(1:R) do _
            mu, τ = _crossed(rng, M, K; a, b, g)
            InferCell._variance_split(mu, τ, fill(1 / (M * K), M, K); held = true)
        end
        for (name, f, truth) in (("total", e -> e.total, a^2 + b^2 + g^2),
                                 ("I first, between", e -> e.i_first.between, a^2),
                                 ("H first, between", e -> e.h_first.between, b^2),
                                 ("I first, within", e -> e.i_first.within, b^2 + g^2),
                                 ("H first, within", e -> e.h_first.within, a^2 + g^2))
            v = f.(est)
            @test abs(mean(v) - truth) <= 4 * std(v) / sqrt(R)
        end
        e = first(est)
        @test e.i_first.between + e.i_first.within ≈ e.total
        @test e.h_first.between + e.h_first.within ≈ e.total
        @test e.held_share.i_first ≈ e.i_first.within / e.total
        @test e.gap ≈ e.h_first.within - e.i_first.between
        @test e.rows == M && e.columns == K && e.pairs == M * K
    end

    # Under a mask the U-statistics' pair weights are not the weights', so the
    # estimator centres μ̂ first: a Q in particles sits ~10⁵ from zero.
    @testset "invariant to where μ sits, masked or not" begin
        rng = Xoshiro(68)
        M, K = 30, 10
        mu, τ = _crossed(rng, M, K)
        v = rand(rng, M, K) .< 0.6
        vw = viability_weights(v)
        for w in (fill(1 / (M * K), M, K), vw.per_theta_I, vw.joint)
            e, f = (InferCell._variance_split(x, τ, w; held = true) for x in (mu, mu .+ 1e5))
            @test e.total ≈ f.total rtol = 1e-8
            @test e.i_first.between ≈ f.i_first.between rtol = 1e-8
            @test e.h_first.between ≈ f.h_first.between rtol = 1e-8
        end
    end

    @testset "rows only, and a column with no μ̂" begin
        rng = Xoshiro(67)
        mu, τ = _crossed(rng, 30, 1; b = 0.0, g = 0.0)
        e = InferCell._variance_split(mu, τ, fill(1 / 30, 30, 1); held = false)
        @test e.held_varies == false && e.h_first.between == 0
        @test e.i_first.between == e.total ≈ var(vec(mu)) - mean(τ)
        mu2, τ2 = _crossed(rng, 5, 4)
        mu2[2, 3] = NaN
        e2 = InferCell._variance_split(mu2, τ2, fill(1 / 20, 5, 4); held = true)
        @test e2.pairs == 19 && isfinite(e2.total)
        τ2[1, 1] = NaN                                   # a pair with one finished cell
        e3 = InferCell._variance_split(mu2, τ2, fill(1 / 20, 5, 4); held = true)
        @test e3.mc_undefined && isnan(e3.total) && isfinite(e3.i_first.between)
        @test_throws ArgumentError InferCell._variance_split(mu2, τ2, zeros(5, 4); held = true)
    end

    # KIM-63's item 1: Σ v(1 − v)τ / (1 − Σv²) is the noise in the reliability-
    # weighted variance for any weights, where Σ vτ is right only for equal ones.
    @testset "the summary's MC term under unequal weights" begin
        rng = Xoshiro(63)
        w = [0.7 0.1 0.1 0.1]
        mu = [1.0 2.0 3.0 4.0]
        n, s = 4, [1.0 2.0 3.0 4.0]
        v = w ./ sum(w)
        truth = sum(v .* (mu .- sum(v .* mu)) .^ 2) / (1 - sum(abs2, v))
        R = 20_000
        d = map(1:R) do _
            cells = [mu[i] .+ s[i] .* randn(rng, n) for i in 1:4]
            q = InferCell._summary(reshape(mean.(cells), 1, 4), reshape(var.(cells), 1, 4),
                                   fill(n, 1, 4), w, 0.0, :hybrid)
            q.var - q.mc
        end
        @test abs(mean(d) - truth) <= 4 * std(d) / sqrt(R)
        τ = s .^ 2 ./ n
        @test sum(v .* (1 .- v) .* τ) / (1 - sum(abs2, v)) != sum(v .* τ)
    end
end
