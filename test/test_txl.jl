using Test
using InferCell
using OrdinaryDiffEq: Tsit5, solve
using Random
using Statistics: quantile

@testset "TranscriptionTranslation Integration" begin
    Random.seed!(42)

    # --- Setup ---
    txl = TranscriptionTranslation()
    prob = build_problem(txl; tspan=(0.0, 50.0))
    sol = solve(prob, Tsit5())
    data = observe(sol, 0.0:2.5:50.0, txl; sigma=0.3)

    # --- Parameter recovery via NUTS ---
    chain = infer(txl, data; n_samples=1000)

    # True parameter values (defaults from TranscriptionTranslation constructor)
    true_vals = Dict(
        :k_tx          => 1.0,
        :k_tl          => 2.0,
        :gamma_mRNA    => 0.5,
        :gamma_protein => 0.1,
        :sigma_obs     => 0.3,
    )

    @testset "parameter recovery (90% CI)" begin
        for (name, true_val) in true_vals
            samples = vec(chain[string(name)].data)
            lo, hi = quantile(samples, [0.05, 0.95])
            @test lo <= true_val <= hi
        end
    end

    # --- Posterior predictive ---
    @testset "posterior predictive" begin
        ppc = posterior_predictive(txl, chain;
            tspan=(0.0, 50.0), saveat=0.0:2.5:50.0, n_samples=50)

        @test ppc isa PosteriorPredictive
        @test length(ppc.solutions) == 50
        @test ppc.species == [:mRNA, :protein]
    end
end
