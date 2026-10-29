# KIM-33: predict's provenance and rng. Two calls with the same rng give the
# same result, the global generator is left as it was, and every pair records
# what it ran and from what.

using Test
using InferCell
using InferCell: AbstractSubModel
using Random

@testset "KIM-33: provenance and rng" begin
    H = 60.0
    ms = collect(AbstractSubModel, m0_models())
    nominal = nominal_parameter_values(ms)
    rng0 = Xoshiro(33)
    X = [nominal[t] * exp(0.1 * randn(rng0)) for _ in 1:4, t in M0_TARGETS]
    pd = PosteriorDraws(M0_TARGETS, X, fill(0.2, 4), [1, 1, 2, 2], [1, 2, 1, 2], nothing,
                        nothing, (kind = :test, path = "made up"))
    go(rng) = predict(ms, pd; horizon = H, n_cells = 2, n_parameter_draws = 2, rng)

    @testset "the same rng gives the same result" begin
        # From two different global states, so the result rests on `rng` alone.
        Random.seed!(1)
        a = go(Xoshiro(7))
        Random.seed!(2)
        b = go(Xoshiro(7))
        @test a.draws == b.draws && a.seeds == b.seeds
        @test a.mu == b.mu && a.sigma2 == b.sigma2
        @test [c.q for cs in a.cells for c in cs] == [c.q for cs in b.cells for c in cs]
        @test a.provenance.fork_seed == b.provenance.fork_seed
        # Another rng draws other posterior draws or other cells.
        c = go(Xoshiro(8))
        @test c.seeds != a.seeds
    end

    @testset "the global generator is left as it was" begin
        Random.seed!(330)
        go(Xoshiro(1))
        after = rand()
        Random.seed!(330)
        @test after == rand()
        # With the global generator as `rng`, the caller's stream moves by exactly
        # the fork's one UInt64, as infer's does.
        Random.seed!(331)
        p = go(Random.default_rng())
        after = rand()
        Random.seed!(331)
        fork = rand(UInt64)
        @test after == rand() && p.provenance.fork_seed == fork
        # And consecutive calls on it differ.
        Random.seed!(332)
        @test go(Random.default_rng()).seeds != go(Random.default_rng()).seeds
    end

    @testset "what each pair ran, and from what" begin
        p = predict(ms, pd; horizon = H, n_cells = 2, draws = [3, 2], rng = Xoshiro(9),
                    intervention = Knockdown(; dose = 0.2))
        pr = p.provenance
        @test pr.version == pkgversion(InferCell)
        @test pr.source == pd.source
        @test pr.settings.forward === :hybrid && pr.settings.n_cells == 2
        @test pr.settings.intervention == Knockdown(; dose = 0.2) && pr.settings.delta === :none
        @test size(pr.pairs) == (2, 1)
        @test [x.draw for x in vec(pr.pairs)] == [3, 2]
        @test pr.pairs[1].chain == 2 && pr.pairs[1].iteration == 1
        @test pr.pairs[2].seeds == p.seeds[:, 1, 2]
        @test all(x -> x.theta_H_source === :none && x.viable && isempty(x.failed), pr.pairs)
        @test !isempty(pr.reductions) && all(l -> l isa ReductionLabel, pr.reductions)
        @test Set(s.name for s in pr.skipped if s.reason === :inferred) == Set(M0_TARGETS)
        m = predict(ms, pd; horizon = H, forward = :mean_field, n_parameter_draws = 1,
                    rng = Xoshiro(9))
        # The mean field's reductions include route (ii)'s own.
        @test any(l -> l.subject === :continuous_drain, m.provenance.reductions)
        @test m.provenance.pairs[1].seeds == Int[]
    end
end
