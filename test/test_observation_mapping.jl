# Spec §11 task 15.9: observations are mapped to states by species name, and the
# ABC summary keeps the ensemble's spread (§12, 2026-09-23 G).

using Test
using InferCell
using OrdinaryDiffEq: Tsit5, solve
using JumpProcesses: SSAStepper
using Random
using Turing: logjoint
using Statistics: std

std_row(x) = std(x)

@testset "Observation mapping by name (spec §11 task 15.9)" begin
    txl = TranscriptionTranslation()
    prob = build_problem(txl; tspan = (0.0, 50.0))
    sol = solve(prob, Tsit5())
    data = observe(sol, 0.0:5.0:50.0, txl; sigma = 0.3, rng = Xoshiro(9))
    θ = (theta = [0.9, 1.1, 0.12, 0.05, 0.4],)

    lp(d) = logjoint(first(build_turing_model([txl], d, prob)), θ)
    base = lp(data)
    @test isfinite(base)

    @testset "permuting the observed species leaves the log-density unchanged" begin
        perm = ObservedData(data.times, data.observations[[2, 1], :], data.species[[2, 1]])
        @test lp(perm) == base
        # The positional mapping this replaces would score the swap differently.
        swapped = ObservedData(data.times, data.observations[[2, 1], :], data.species)
        @test lp(swapped) != base
    end

    @testset "a subset of states builds and scores" begin
        sub = ObservedData(data.times, data.observations[[2], :], [:protein])
        @test isfinite(lp(sub))
        @test lp(sub) != base
    end

    @testset "an unknown species throws, naming it" begin
        bad = ObservedData(data.times, data.observations, [:mRNA, :not_a_state])
        err = try
            build_turing_model([txl], bad, prob)
            nothing
        catch e
            e
        end
        @test err isa ArgumentError && occursin("not_a_state", err.msg)
        @test_throws DimensionMismatch ObservedData(data.times, data.observations, [:mRNA])
    end
end

@testset "The ABC summary keeps the spread (spec §11 task 15.9)" begin
    model = StochasticGeneExpression()
    prob = build_problem(model; tspan = (0.0, 50.0))
    times = [10.0, 25.0, 50.0]
    Random.seed!(15)
    trajs = [solve(prob, SSAStepper(); saveat = times) for _ in 1:40]
    sp = states(model)

    s = compute_summary_stats(trajs, sp; times)
    @test length(s) == 2 * length(sp) * length(times)
    # Per time: the means, then the standard deviations across replicates.
    vals = [trajs[j](times[2])[i] for i in eachindex(sp), j in eachindex(trajs)]
    blk = s[(2 * length(sp) + 1):(4 * length(sp))]
    @test blk[1:length(sp)] ≈ vec(sum(vals; dims = 2)) ./ length(trajs)
    @test blk[(length(sp) + 1):end] ≈ [std_row(vals[i, :]) for i in eachindex(sp)]
    @test all(>(0), blk[(length(sp) + 1):end])
    # The mean-only form is still available, and is the first half of each block.
    m = compute_summary_stats(trajs, sp; times, spread = false)
    @test m[1:length(sp)] == s[1:length(sp)]

    # By name: a reversed subset summarises the same numbers, reordered.
    r = compute_summary_stats(trajs, reverse(sp); times, state_names = sp)
    @test r[1:length(sp)] == reverse(s[1:length(sp)])
    @test_throws ArgumentError compute_summary_stats(trajs, [:nope]; times, state_names = sp)

    # observe keeps the spread, and it matches the summary.
    d = observe(trajs, times, model)
    @test d.spread !== nothing && size(d.spread) == size(d.observations)
    @test vec(vcat(d.observations, d.spread)) == s
    # One trajectory has no spread to keep.
    @test observe(trajs[1:1], times, model).spread === nothing
end
