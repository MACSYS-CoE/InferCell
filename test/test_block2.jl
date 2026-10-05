# spec/phases/16-recovery.md task 16a.4, check V3: block 2 from sufficient
# statistics, alternating with the exact transcript bridge, recovers the exact
# transcript-only posterior of one gene with frozen rates.

using Random
using Statistics
using Test
using InferCell
using Distributions: LogNormal, Normal, logpdf, params

# Gillespie for one birth-death gene, observed at every `dt`.
function _simulate_gene(rng, k, μ, m0, dt, nobs)
    m, t = m0, 0.0
    obs = Int[]
    for i in 1:nobs
        tend = i * dt
        while true
            total = k + μ * m
            t += randexp(rng) / total
            t >= tend && (t = tend; break)
            m += rand(rng) * total < k ? 1 : -1
        end
        push!(obs, m)
    end
    return obs
end

# Occupancy ∫ m dt and the birth and death counts of one bridged window.
function _window_stats(times, deltas, a, τ)
    occ, m, t = 0.0, a, 0.0
    for (ti, d) in zip(times, deltas)
        occ += m * (ti - t)
        m += d
        t = ti
    end
    occ += m * (τ - t)
    return count(==(1), deltas), count(==(-1), deltas), occ
end

# Monte Carlo SE of a mean by batch means.
function _batch_se(x; nb = 50)
    L = length(x) ÷ nb
    means = [mean(@view x[(b - 1) * L + 1:b * L]) for b in 1:nb]
    return std(means) / sqrt(nb)
end

@testset "16a.4: block 2 against the exact posterior (V3)" begin
    rng = Xoshiro(1604)
    kp, μp = LogNormal(log(0.04), log(2.0)), LogNormal(log(0.015), log(2.0))
    k_true, μ_true = 0.05, 0.02
    dt, nobs, m0 = 60.0, 60, 2
    obs = _simulate_gene(rng, k_true, μ_true, m0, dt, nobs)
    counts = vcat(m0, obs)
    cap = bridge_cap(maximum(counts))

    # The exact posterior of (ln k, ln μ): transition probabilities between
    # consecutive exact counts, and the priors carried to log scale. A coarse
    # grid over the prior finds the posterior; a fine one from −7 to +12 of its SDs
    # carries the marginals, whose CDFs are built by trapezoid. The box runs to
    # +12 SDs because k and μ rise together along the ridge that fixes the mean
    # count, so the upper tails are heavy.
    function logpost(a, b)
        P = exp(InferCell.generator(TranscriptChain(exp(a), exp(b), cap)) * dt)
        return sum(log(P[counts[i] + 1, counts[i + 1] + 1]) for i in 1:nobs) +
               logpdf(Normal(params(kp)...), a) + logpdf(Normal(params(μp)...), b)
    end
    function marginals(ua, ub)
        lp = [logpost(a, b) for a in ua, b in ub]
        w = exp.(lp .- maximum(lp))
        return vec(sum(w; dims = 2)), vec(sum(w; dims = 1))
    end
    moments(u, m) = (c = sum(u .* m) / sum(m); (c, sqrt(sum((u .- c) .^ 2 .* m) / sum(m))))
    ca, cb = marginals(range(log(0.04) - 4log(2), log(0.04) + 4log(2); length = 81),
                       range(log(0.015) - 4log(2), log(0.015) + 4log(2); length = 81))
    (mk, sk) = moments(range(log(0.04) - 4log(2), log(0.04) + 4log(2); length = 81), ca)
    (mμ, sμ) = moments(range(log(0.015) - 4log(2), log(0.015) + 4log(2); length = 81), cb)
    uk = range(mk - 7sk, mk + 12sk; length = 451)
    uμ = range(mμ - 7sμ, mμ + 12sμ; length = 451)
    fk, fμ = marginals(uk, uμ)
    function trapezoid_cdf(u, f)
        c = [0.0; cumsum([(f[i] + f[i + 1]) / 2 * step(u) for i in 1:length(u) - 1])]
        return c ./ c[end]
    end
    cdf_k, cdf_μ = trapezoid_cdf(uk, fk), trapezoid_cdf(uμ, fμ)
    function grid_quantile(cdf, u, p)
        i = findfirst(>=(p), cdf)
        return u[i - 1] + (p - cdf[i - 1]) / (cdf[i] - cdf[i - 1]) * step(u)
    end
    # The fine grid's edges hold nothing.
    @test fk[1] / maximum(fk) < 1e-8 && fk[end] / maximum(fk) < 1e-8
    @test fμ[1] / maximum(fμ) < 1e-8 && fμ[end] / maximum(fμ) < 1e-8

    # The Gibbs sampler: bridges given the rates, then block 2 given the path.
    k, μ = 0.04, 0.015
    nsweep, burn = 20_000, 1_000
    draws = zeros(nsweep, 2)
    for s in 1:(nsweep + burn)
        nb = nd = 0
        occ = 0.0
        bd = TranscriptChain(k, μ, cap)
        tables = Dict{Int, BridgeTable}()
        for i in 1:nobs
            a, b = counts[i], counts[i + 1]
            t = get!(() -> BridgeTable(bd, dt, b), tables, b)
            times, deltas = sample_bridge(rng, t, a)
            b1, d1, o1 = _window_stats(times, deltas, a, dt)
            nb += b1
            nd += d1
            occ += o1
        end
        # Transcription's exposure is k·T; decay's is μ·∫m dt.
        k = rate_update(rng, k, nb, k * nobs * dt, kp)
        μ = rate_update(rng, μ, nd, μ * occ, μp)
        s > burn && (draws[s - burn, :] = [log(k), log(μ)])
    end

    # Each exact quantile's CDF value, estimated from the chain, within 3 MC SE.
    for (j, cdf, u) in ((1, cdf_k, uk), (2, cdf_μ, uμ))
        for p in (0.05, 0.25, 0.5, 0.75, 0.95)
            q = grid_quantile(cdf, u, p)
            ind = Float64.(draws[:, j] .<= q)
            @test abs(mean(ind) - p) <= 3 * _batch_se(ind)
        end
    end
    # The posterior is informed: narrower than the prior on both.
    @test std(draws[:, 1]) < 0.9 * log(2) && std(draws[:, 2]) < 0.9 * log(2)

    @testset "the update's pieces" begin
        # The conditional carries the Jacobian: with no data it is the prior on u.
        f = rate_log_conditional(0, 0.0, kp)
        @test f(-3.0) ≈ logpdf(Normal(params(kp)...), -3.0)
        @test f(-3.0) ≈ logpdf(kp, exp(-3.0)) + (-3.0)
        # An accepted value past the linear range throws.
        @test_throws ArgumentError rate_update(Xoshiro(1), 5.0, 1000, 5.0 * 10, kp; bound = 1.0)
    end
end
