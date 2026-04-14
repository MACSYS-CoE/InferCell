using Test
using InferCell
using OrdinaryDiffEq: Tsit5, solve
using Random
using Statistics: quantile

@testset "Joint TX/TL + Metabolism Inference (Integration)" begin
    Random.seed!(123)

    # --- Setup ---
    txl = TranscriptionTranslation()
    metab = LightMetabolism()
    models = [txl, metab]

    prob = build_problem(models; tspan=(0.0, 50.0))
    sol = solve(prob, Tsit5())
    data = observe(sol, 0.0:2.5:50.0, models; sigma=0.3)

    @test size(data.observations) == (5, 21)
    @test data.species == [:mRNA, :protein, :ATP, :NTP, :AA]

    # --- Identifiability check ---
    ident = check_identifiability(models, prob, collect(0.0:2.5:50.0))
    @test ident.full_rank == true

    # --- Joint NUTS inference ---
    chain = infer(models, data; n_samples=1000)

    true_vals = Dict(
        :k_tx          => 1.0,
        :k_tl          => 2.0,
        :gamma_mRNA    => 0.5,
        :gamma_protein => 0.1,
        :k_atp         => 1.0,
        :k_ntp         => 0.2,
        :k_aa          => 0.4,
        :sigma_obs     => 0.3,
    )

    @testset "parameter recovery (95% CI)" begin
        for (name, true_val) in true_vals
            samples = vec(chain[string(name)].data)
            lo, hi = quantile(samples, [0.025, 0.975])
            @test lo <= true_val <= hi
        end
    end

    # --- Posterior predictive ---
    @testset "posterior predictive" begin
        ppc = posterior_predictive(models, chain;
            tspan=(0.0, 50.0), saveat=0.0:2.5:50.0, n_samples=50)

        @test ppc isa PosteriorPredictive
        @test length(ppc.solutions) == 50
        @test ppc.species == [:mRNA, :protein, :ATP, :NTP, :AA]
    end
end
