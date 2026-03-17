using OrdinaryDiffEq
using SciMLSensitivity
using SciMLBase: ReturnCode
using Turing
using Distributions
using StaticArrays
using Random

Random.seed!(42)

# --- Step 1: ODE dynamics (out-of-place) ---
function txl_dynamics(u, p, t)
    mRNA, protein = u
    k_tx, k_tl, γ_mRNA, γ_protein = p
    SA[k_tx - γ_mRNA * mRNA, k_tl * mRNA - γ_protein * protein]
end

# --- Step 2: Synthetic data ---
true_params = SA[1.0, 2.0, 0.5, 0.1]  # k_tx, k_tl, γ_mRNA, γ_protein
σ_true = 0.3
u0 = SA[0.0, 0.0]
tspan = (0.0, 50.0)
times = 0.0:2.5:50.0

prob = ODEProblem(txl_dynamics, u0, tspan, true_params)
sol_true = solve(prob, Tsit5(); saveat=times)
data = Array(sol_true) .+ σ_true .* randn(size(Array(sol_true)))

println("Synthetic data generated: $(size(data, 2)) timepoints, $(size(data, 1)) species")

# --- Step 3: Turing model ---
@model function txl_inference(data, times, prob)
    k_tx ~ LogNormal(0, 1)
    k_tl ~ LogNormal(0, 1)
    γ_mRNA ~ LogNormal(0, 1)
    γ_protein ~ LogNormal(0, 1)
    σ ~ truncated(Normal(0, 1); lower=0)

    p = SA[k_tx, k_tl, γ_mRNA, γ_protein]
    sol = solve(remake(prob, p=p), Tsit5(); saveat=times,
                sensealg=ForwardDiffSensitivity())

    if sol.retcode !== ReturnCode.Success
        Turing.@addlogprob! -Inf
        return
    end

    for i in eachindex(times)
        data[:, i] ~ MvNormal(sol[:, i], σ)
    end
end

# --- Step 4: Sample ---
println("Starting NUTS sampling (1000 draws)...")
t_start = time()
model = txl_inference(data, collect(times), prob)
chain = sample(model, NUTS(), 1000)
elapsed = time() - t_start

# --- Step 5: Validate ---
println("\n", chain)
println("\nElapsed time: $(round(elapsed; digits=1)) seconds")

# Check 90% credible intervals
true_vals = Dict(:k_tx => 1.0, :k_tl => 2.0, :γ_mRNA => 0.5, :γ_protein => 0.1, :σ => 0.3)
all_pass = true
println("\n--- Parameter recovery (90% CI) ---")
for (name, true_val) in true_vals
    samples = chain[name] |> vec
    lo, hi = quantile(samples, [0.05, 0.95])
    recovered = lo ≤ true_val ≤ hi
    status = recovered ? "PASS" : "FAIL"
    global all_pass = all_pass && recovered
    println("  $name: true=$true_val, 90% CI=[$( round(lo;digits=3)), $(round(hi;digits=3))] → $status")
end

perf_pass = elapsed < 600
println("\nPerformance: $(round(elapsed; digits=1))s $(perf_pass ? "< 10min → PASS" : ">= 10min → FAIL")")
println("\n$(all_pass && perf_pass ? "✓ ALL CRITERIA PASSED" : "✗ SOME CRITERIA FAILED")")
