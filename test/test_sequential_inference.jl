using Test
using InferCell
using OrdinaryDiffEq: Tsit5, solve
using MCMCChains: Chains
using Random
using Statistics: quantile, var

@testset "Sequential Inference (Integration)" begin
    Random.seed!(789)

    # True parameters (shared between ODE and SSA models)
    true_params = (k_tx=1.0, k_tl=2.0, gamma_mRNA=0.5, gamma_protein=0.1)

    # --- ODE block: TranscriptionTranslation ---
    txl = TranscriptionTranslation(; true_params..., sigma_obs=0.3)
    ode_prob = build_problem(txl; tspan=(0.0, 50.0))
    ode_sol = solve(ode_prob, Tsit5())
    ode_times = 0.0:2.5:50.0
    ode_data = observe(ode_sol, ode_times, txl; sigma=0.3)

    # --- SSA block: StochasticGeneExpression ---
    sge = StochasticGeneExpression(; true_params...)
    ssa_tspan = (0.0, 100.0)
    ssa_times = collect(0.0:10.0:100.0)
    ssa_prob = build_problem(sge; tspan=ssa_tspan)
    n_data_replicates = 200
    trajectories = [solve(ssa_prob, SSAStepper(); saveat=ssa_times) for _ in 1:n_data_replicates]
    ssa_data = observe(trajectories, ssa_times, sge)

    # --- Sequential inference ---
    result = sequential_infer(txl, ode_data, sge, ssa_data;
                              ode_kwargs=(; n_samples=500),
                              ssa_kwargs=(; n_particles=200, n_populations=6,
                                            n_replicates=50, tspan=ssa_tspan),
                              boundary_method=:kde, verbose=true)

    @test result.ode_chain isa Chains
    @test result.ssa_posterior isa ABCPosterior
    @test length(result.boundary.priors) == 4
    @test all(result.boundary.conditioned)

    @testset "ODE parameter recovery (90% CI)" begin
        ode_true = Dict(:k_tx => 1.0, :k_tl => 2.0, :gamma_mRNA => 0.5,
                        :gamma_protein => 0.1, :sigma_obs => 0.3)
        for (name, true_val) in ode_true
            samples = vec(result.ode_chain[string(name)].data)
            lo, hi = quantile(samples, [0.05, 0.95])
            @test lo <= true_val <= hi
        end
    end

    @testset "SSA parameter recovery (wide tolerance)" begin
        for (j, (name, true_val)) in enumerate(pairs(true_params))
            posterior_mean = sum(result.ssa_posterior.weights .* result.ssa_posterior.particles[j, :])
            @test 0.1 * true_val < posterior_mean < 10.0 * true_val
        end
    end

    @testset "conditioned posteriors tighter than unconditioned" begin
        # The control run differs from the conditioned one only in its priors,
        # so seed immediately before it: the comparison is then reproducible
        # rather than a single draw from two uncoupled RNG streams.
        Random.seed!(789)
        uncond = infer(sge, ssa_data;
                       n_samples=200, n_populations=6,
                       n_replicates=50, tspan=ssa_tspan)

        # Weighted variance of row `j` of an ABC posterior.
        function weighted_var(post, j)
            mean = sum(post.weights .* post.particles[j, :])
            return sum(post.weights .* (post.particles[j, :] .- mean).^2)
        end

        # This is the whole claim of the boundary protocol: conditioning the
        # SSA block on the ODE block's posterior narrows the shared
        # parameters. Assert the direction, with no slack. The previous
        # `cond_var < uncond_var * 5.0` passed when conditioning made the
        # posterior five times *wider*, so it could not fail for the reason
        # the testset exists.
        #
        # NOT VERIFIED. This testset has never executed to completion — not in
        # the 48 GitHub runs, and not on Slurm, where job 16617632 timed out at
        # six hours. The assertion states the claim the protocol is supposed to
        # make; whether the claim holds is unknown. Do not read a green suite
        # here as evidence until it has actually run.
        for j in 1:4
            @test weighted_var(result.ssa_posterior, j) < weighted_var(uncond, j)
        end
    end
end
