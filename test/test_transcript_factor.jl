# spec/phases/approximate-method.md task B1: the exact transcript factor, AV1 to
# AV3.

using Random
using Statistics: mean, std
using Test
using InferCell
using ForwardDiff
using Distributions: Chisq, cdf

# One gene's count after `τ`, by Gillespie, with births at `k0` until
# `τ − switch` and at `k1` after, and deaths at `μ` per transcript. No cap: the
# truncation the factor makes is what this checks.
function _ssa_switch(rng, a, k0, k1, μ, τ, switch)
    c, t = a, 0.0
    tsw = τ - switch
    while true
        k = t < tsw ? k0 : k1
        total = k + μ * c
        tnext = t + randexp(rng) / total
        # A rate change is a boundary: restart the clock there, by memorylessness.
        if t < tsw && tnext >= tsw
            t = tsw
            continue
        end
        tnext >= τ && return c
        c += rand(rng) * total < k ? 1 : -1
        t = tnext
    end
end

# Pearson χ² p-value, pooling bins until each expects at least 5, as in
# test_transcript_bridge.jl.
function _chisq_p_tf(observed::Vector{Int}, p::Vector{Float64})
    n = sum(observed)
    obs, ex = Int[], Float64[]
    o, e = 0, 0.0
    for j in eachindex(p)
        o += observed[j]
        e += n * p[j]
        if e >= 5
            push!(obs, o)
            push!(ex, e)
            o, e = 0, 0.0
        end
    end
    if !isempty(ex)
        obs[end] += o
        ex[end] += e
    end
    length(obs) >= 2 || return 1.0
    return 1 - cdf(Chisq(length(obs) - 1), sum((obs .- ex) .^ 2 ./ ex))
end

@testset "B1: the transcript factor" begin
    # Typical Core A′ rates (k ~ 0.02 /s, krnadeg/n ~ 0.009 /s, as in
    # test_transcript_bridge.jl) and an amplified case with many events.
    y0 = [2, 5]
    y = [3 1 2 4; 7 9 6 8]
    μ = [0.009, 0.05]
    f = TranscriptFactor(y0, y; kmax = [0.1, 2.0], μmin = [0.005, 0.03])
    G, M = size(y)

    @testset "the cap holds the truncation mass under 1e-12 at the prior's corner" begin
        for (g, (kx, μn)) in enumerate(zip([0.1, 2.0], [0.005, 0.03]))
            top = max(y0[g], maximum(y[g, :]))
            @test f.caps[g] >= bridge_cap(top)
            @test truncation_mass(TranscriptChain(kx, μn, f.caps[g]), 60.0, top) < 1e-12
        end
    end

    @testset "AV1: constant constants give the 60 s transition probability" begin
        k = repeat([0.02, 1.0], 1, M + 1)
        expected = sum(log(transition_probability(TranscriptChain(k[g, 1], μ[g], f.caps[g]), 60.0,
                                                  m == 1 ? y0[g] : y[g, m - 1], y[g, m]))
                       for g in 1:G, m in 1:M)
        @test transcript_loglik(f, k, μ) ≈ expected atol = 1e-10
    end

    @testset "AV1: a switch at 59 s gives P(59) · P(1) from dense matrices" begin
        # Rates differing by at least 2× across every rebuild.
        k = [0.02 0.06 0.01 0.03 0.08; 1.0 0.3 1.5 0.5 1.2]
        expected = 0.0
        for g in 1:G, m in 1:M
            a = m == 1 ? y0[g] : y[g, m - 1]
            P59 = transition_matrix(TranscriptChain(k[g, m], μ[g], f.caps[g]), 59.0)
            P1 = transition_matrix(TranscriptChain(k[g, m + 1], μ[g], f.caps[g]), 1.0)
            expected += log((P59 * P1)[a + 1, y[g, m] + 1])
        end
        @test transcript_loglik(f, k, μ) ≈ expected atol = 1e-10
        # Reading the next window's constants for the whole window is caught.
        k_late = hcat(k[:, 2:end], k[:, end])
        @test abs(transcript_loglik(f, k_late, μ) - expected) > 1e-3
    end

    @testset "AV1: the cSMC's window weight averages to the factor" begin
        # D16.3's weight, P_{R_{m−1}}(y_m | y_{m−1}; 60) · P_{R_m}(y_m | c; 1) /
        # P_{R_{m−1}}(y_m | c; 1), with c the 60 s bridge's count at 59 s, has
        # expectation Σ_c P_{R_{m−1}}(59)[a, c] · P_{R_m}(1)[c, b].
        rng = Xoshiro(8901)
        a, b, k0, k1, μ1 = 5, 9, 1.0, 0.3, 0.05
        ff = TranscriptFactor([a], fill(b, 1, 1); kmax = 1.0, μmin = μ1)
        cap = ff.caps[1]
        c0, c1 = TranscriptChain(k0, μ1, cap), TranscriptChain(k1, μ1, cap)
        lp60 = log(transition_probability(c0, 60.0, a, b))
        n = 20_000
        w = map(1:n) do _
            times, deltas = sample_bridge(rng, c0, 60.0, a, b)
            c = a + sum(d for (t, d) in zip(times, deltas) if t < 59.0; init = 0)
            exp(window_log_weight(lp60, log(transition_probability(c1, 1.0, c, b)),
                                  log(transition_probability(c0, 1.0, c, b)), 0.0))
        end
        target = exp(transcript_loglik(ff, [k0 k1], [μ1]))
        se = std(w) / sqrt(n)
        @test abs(mean(w) - target) < 3se
    end

    @testset "AV2: SSA frequencies with the rate switching at 59 s" begin
        # Two cases at 10⁵ replicates each, each χ² at p > 0.001.
        rng = Xoshiro(8902)
        for (a, k0, k1, μg) in [(2, 0.02, 0.06, 0.009), (5, 1.0, 0.3, 0.05)]
            ff = TranscriptFactor([a], fill(a, 1, 1); kmax = max(k0, k1), μmin = μg)
            p = InferCell._window_distribution(ff, 1, a, k0, k1, μg)
            @test sum(p) ≈ 1 atol = 1e-10
            counts = zeros(Int, length(p))
            for _ in 1:100_000
                c = _ssa_switch(rng, a, k0, k1, μg, 60.0, 1.0)
                c < length(p) || error("an SSA count passed the cap")
                counts[c + 1] += 1
            end
            @test _chisq_p_tf(counts, p) > 0.001
        end
    end

    @testset "AV3: ForwardDiff matches central differences" begin
        k = [0.02 0.06 0.01 0.03 0.08; 1.0 0.3 1.5 0.5 1.2]
        u0 = vcat(vec(log.(k)), log.(μ))
        nk = length(k)
        ℓ(u) = transcript_loglik(f, reshape(exp.(u[1:nk]), size(k)), exp.(u[(nk + 1):end]))
        grad = ForwardDiff.gradient(ℓ, u0)
        h = 1e-5
        fd = map(eachindex(u0)) do i
            e = zeros(length(u0))
            e[i] = h
            (ℓ(u0 + e) - ℓ(u0 - e)) / 2h
        end
        @test all(isfinite, grad)
        @test grad ≈ fd rtol = 1e-6
    end

    @testset "cells sum, and bad input is refused" begin
        k = repeat([0.02, 1.0], 1, M + 1)
        f2 = TranscriptFactor([1, 4], [2 2 1 3; 6 5 7 7]; kmax = [0.1, 2.0], μmin = [0.005, 0.03])
        @test transcript_loglik([f, f2], k, μ) ≈ transcript_loglik(f, k, μ) + transcript_loglik(f2, k, μ)
        @test_throws DimensionMismatch transcript_loglik(f, k[:, 1:M], μ)
        @test_throws DimensionMismatch transcript_loglik(f, k, μ[1:1])
        @test_throws DimensionMismatch TranscriptFactor([1], y; kmax = 1.0, μmin = 0.01)
        @test_throws ArgumentError TranscriptFactor(y0, -y; kmax = 1.0, μmin = 0.01)
        @test_throws ArgumentError TranscriptFactor(y0, y; kmax = 1.0, μmin = 0.01, switch = 60.0)
    end
end
