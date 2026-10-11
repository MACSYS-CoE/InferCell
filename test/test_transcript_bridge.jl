# spec/phases/16-recovery.md task 16a.3: exact per-gene transcript bridges and
# their transition probabilities.

using Random
using Test
using InferCell
using Distributions: Chisq, cdf

# Pearson χ² p-value of sampled counts against exact probabilities, pooling
# bins until each expects at least 5.
function _chisq_p(observed::Vector{Int}, p::Vector{Float64})
    n = sum(observed)
    obs, exp_ = Int[], Float64[]
    o, e = 0, 0.0
    for j in eachindex(p)
        o += observed[j]
        e += n * p[j]
        if e >= 5
            push!(obs, o)
            push!(exp_, e)
            o, e = 0, 0.0
        end
    end
    # The remainder joins the last bin.
    if !isempty(exp_)
        obs[end] += o
        exp_[end] += e
    end
    length(obs) >= 2 || return 1.0
    stat = sum((obs .- exp_) .^ 2 ./ exp_)
    return 1 - cdf(Chisq(length(obs) - 1), stat)
end

@testset "16a.3: transcript bridges" begin
    # Typical Core A′ rates at the 15.7 truth (k ~ 0.02 /s, krnadeg/n ~ 0.009 /s),
    # and an amplified case of the kind V5 needs.
    typical = TranscriptChain(0.02, 0.009, bridge_cap(2))
    amplified = TranscriptChain(1.0, 0.2, bridge_cap(10))

    @testset "transition probabilities match a matrix exponential to 1e-12" begin
        for bd in (typical, amplified), τ in (1.0, 20.0, 60.0)
            P = transition_matrix(bd, τ)
            E = exp(InferCell.generator(bd) * τ)
            @test maximum(abs.(P .- E)) <= 1e-12
            @test maximum(abs.(sum(P; dims = 2) .- 1)) <= 1e-12
        end
    end

    @testset "the Poisson weights stop once the tail is negligible" begin
        # The stop once ran to 10λ + 200 terms whatever λ was (the #81 review):
        # about 211 at λ near 1, where about 20 carry everything.
        for λ in (0.5, 1.0, 20.0, 600.0)
            w = InferCell._poisson_weights(λ)
            @test length(w) <= λ + 10sqrt(λ) + 25
            # The dropped tail, summed far past the stop, is below the tolerance.
            n0 = length(w)
            tail = sum(exp(-λ + n * log(λ) - sum(log, 1:n; init = 0.0)) for n in n0:(n0 + 500))
            @test tail < 1e-17
        end
    end

    @testset "the truncation mass is below 1e-12 at the dataset cap" begin
        @test truncation_mass(typical, 60.0, 2) < 1e-12
        # The dataset rule is a floor, not a guarantee: at amplified rates it
        # drops more than 1e-12, and adequate_cap grows it until it does not.
        @test truncation_mass(amplified, 20.0, 10) > 1e-12
        cap = adequate_cap(1.0, 0.2, 20.0, 10; floor = bridge_cap(10))
        @test cap > bridge_cap(10)
        @test truncation_mass(TranscriptChain(1.0, 0.2, cap), 20.0, 10) < 1e-12
        @test adequate_cap(0.02, 0.009, 60.0, 2; floor = bridge_cap(2)) == bridge_cap(2)
        # It discriminates: a cap near the counts reached drops real mass.
        @test truncation_mass(TranscriptChain(1.0, 0.2, 8), 20.0, 5) > 1e-3
    end

    @testset "bridges hit both endpoints, and their birth counts are exact" begin
        # 4 cases at 10⁵ draws each, each χ² at p > 0.01: a family-wise false
        # failure rate of about 4%.
        cases = [(typical, 60.0, 1, 1), (typical, 60.0, 0, 2),
                 (amplified, 1.0, 3, 4), (amplified, 20.0, 5, 12)]
        rng = Xoshiro(1603)
        for (bd, τ, a, b) in cases
            t = BridgeTable(bd, τ, b)
            exact = bridge_birth_distribution(bd, τ, a, b; jmax = 80)
            births = zeros(Int, length(exact))
            ok = true
            for _ in 1:100_000
                times, deltas = sample_bridge(rng, t, a)
                ok &= a + sum(deltas; init = 0) == b && issorted(times) &&
                      all(x -> 0 <= x < τ, times) && all(d -> abs(d) == 1, deltas)
                births[count(==(1), deltas) + 1] += 1
            end
            @test ok
            @test _chisq_p(births, exact) > 0.01
        end
        # An unreachable endpoint is refused.
        @test_throws ArgumentError sample_bridge(rng, TranscriptChain(0.0, 0.1, 5), 1.0, 0, 3)
    end
end
