# Standalone integration test for sequential inference
# Runs directly (not via Pkg.test) so we get live progress output

using InferCell
using OrdinaryDiffEq: Tsit5, solve
using MCMCChains: Chains
using Random
using Statistics: quantile
using Dates
using Test

println("=== Sequential Inference Integration Test ===")
println("Start: $(Dates.now())")
flush(stdout)

Random.seed!(789)

# True parameters
true_params = (k_tx=1.0, k_tl=2.0, gamma_mRNA=0.5, gamma_protein=0.1)

# --- ODE block ---
println("\n1. Setting up ODE block...")
flush(stdout)
txl = TranscriptionTranslation(; true_params..., sigma_obs=0.3)
ode_prob = build_problem(txl; tspan=(0.0, 50.0))
ode_sol = solve(ode_prob, Tsit5())
ode_data = observe(ode_sol, 0.0:2.5:50.0, txl; sigma=0.3)
println("   ODE data: $(size(ode_data.observations))")
flush(stdout)

# --- SSA block ---
println("\n2. Setting up SSA block...")
flush(stdout)
sge = StochasticGeneExpression(; true_params...)
ssa_tspan = (0.0, 100.0)
ssa_times = collect(0.0:10.0:100.0)
ssa_prob = build_problem(sge; tspan=ssa_tspan)
trajectories = [solve(ssa_prob, SSAStepper(); saveat=ssa_times) for _ in 1:200]
ssa_data = observe(trajectories, ssa_times, sge)
println("   SSA data: $(size(ssa_data.observations))")
flush(stdout)

# --- Sequential inference ---
println("\n3. Running sequential_infer...")
flush(stdout)
t_start = time()
result = sequential_infer(txl, ode_data, sge, ssa_data;
                          ode_kwargs=(; n_samples=500),
                          ssa_kwargs=(; n_particles=200, n_populations=6,
                                        n_replicates=50, tspan=ssa_tspan),
                          boundary_method=:kde, verbose=true)
elapsed = time() - t_start
println("   Sequential inference complete in $(round(elapsed; digits=1))s")
flush(stdout)

# --- Check results ---
println("\n4. Checking results...")
flush(stdout)

@testset "Sequential inference results" begin
    @test result.ode_chain isa Chains
    @test result.ssa_posterior isa ABCPosterior
    @test length(result.boundary.priors) == 4
    @test all(result.boundary.conditioned)
end

@testset "ODE parameter recovery (90% CI)" begin
    ode_true = Dict(:k_tx => 1.0, :k_tl => 2.0, :gamma_mRNA => 0.5,
                    :gamma_protein => 0.1, :sigma_obs => 0.3)
    for (name, true_val) in ode_true
        samples = vec(result.ode_chain[string(name)].data)
        lo, hi = quantile(samples, [0.05, 0.95])
        recovered = lo <= true_val <= hi
        println("   $name: true=$true_val, CI=[$lo, $hi] $(recovered ? "PASS" : "FAIL")")
        @test recovered
    end
end

@testset "SSA parameter recovery (wide tolerance)" begin
    for (j, (name, true_val)) in enumerate(pairs(true_params))
        posterior_mean = sum(result.ssa_posterior.weights .* result.ssa_posterior.particles[j, :])
        println("   $name: true=$true_val, posterior_mean=$(round(posterior_mean; digits=3))")
        @test 0.1 * true_val < posterior_mean < 10.0 * true_val
    end
end

# --- Comparison: conditioned vs unconditioned ---
println("\n5. Running unconditioned ABC-SMC for comparison...")
flush(stdout)
t_uncond = time()
uncond = infer(sge, ssa_data;
               n_samples=200, n_populations=6,
               n_replicates=50, tspan=ssa_tspan, verbose=true)
elapsed_uncond = time() - t_uncond
println("   Unconditioned ABC-SMC complete in $(round(elapsed_uncond; digits=1))s")
flush(stdout)

@testset "conditioned posteriors tighter than unconditioned" begin
    for j in 1:4
        cond_mean = sum(result.ssa_posterior.weights .* result.ssa_posterior.particles[j, :])
        cond_var = sum(result.ssa_posterior.weights .* (result.ssa_posterior.particles[j, :] .- cond_mean).^2)

        uncond_mean = sum(uncond.weights .* uncond.particles[j, :])
        uncond_var = sum(uncond.weights .* (uncond.particles[j, :] .- uncond_mean).^2)

        ratio = uncond_var > 0 ? round(uncond_var / cond_var; digits=1) : Inf
        println("   Param $j: cond_var=$(round(cond_var;digits=4)), uncond_var=$(round(uncond_var;digits=4)), ratio=$(ratio)x")
        @test cond_var < uncond_var * 5.0
    end
end

println("\n=== All tests complete at $(Dates.now()) ===")
