# KIM-34: posterior_predictive, the toy models' predictive, on predict's ingest
# and rng path. It reads every chain, maps columns by name, and draws from a
# forked stream. The integration-gated test_txl.jl and test_joint_inference.jl
# are its original call sites.

using Test
using InferCell
using Random
using MCMCChains: Chains

@testset "KIM-34: posterior_predictive on the shared ingest and rng" begin
    txl = TranscriptionTranslation()
    names = [:k_tx, :k_tl, :gamma_mRNA, :gamma_protein, :sigma_obs]
    # Two chains of five draws; k_tx is 1 in chain 1 and 3 in chain 2.
    X = zeros(5, 5, 2)
    X[:, :, 1] .= [1.0 2.0 0.5 0.1 0.3]
    X[:, :, 2] .= [3.0 2.0 0.5 0.1 0.3]
    ch = Chains(X, names)
    kw = (; tspan = (0.0, 10.0), saveat = 0.0:1.0:10.0)

    ppc = posterior_predictive(txl, ch; n_samples = 10, rng = Xoshiro(1), kw...)
    @test ppc isa PosteriorPredictive && length(ppc.solutions) == 10
    @test ppc.species == [:mRNA, :protein] && ppc.times == collect(0.0:1.0:10.0)
    # Every chain is read, not only the first.
    @test sort(unique(s.prob.p[1] for s in ppc.solutions)) == [1.0, 3.0]

    k(p) = [s.prob.p[1] for s in p.solutions]
    a = posterior_predictive(txl, ch; n_samples = 4, rng = Xoshiro(2), kw...)
    b = posterior_predictive(txl, ch; n_samples = 4, rng = Xoshiro(2), kw...)
    @test k(a) == k(b) && length(a.solutions) == 4

    Random.seed!(34)
    posterior_predictive(txl, ch; n_samples = 3, rng = Xoshiro(3), kw...)
    after = rand()
    Random.seed!(34)
    @test after == rand()

    # A vector of models and a PosteriorDraws take the same path.
    @test length(posterior_predictive([txl], posterior_draws(ch); n_samples = 2, kw...).solutions) == 2
    # The ABC path resamples by weight.
    abc = ABCPosterior(permutedims(X[:, 1:4, 1]), [0.0, 0.0, 1.0, 0.0, 0.0], names[1:4], 0.1, 2)
    r = posterior_predictive(txl, abc; n_samples = 3, kw...)
    @test length(r.solutions) == 3

    @test_throws ArgumentError posterior_predictive(txl, Chains(X[:, 1:3, :], names[1:3]);
                                                    n_samples = 2, kw...)
    @test_throws ErrorException posterior_predictive(txl, ch; n_samples = 2)
end
