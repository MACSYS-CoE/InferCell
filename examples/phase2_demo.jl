# InferCell Phase 2 Demo: Unified inference API for ODE and SSA models
#
# Demonstrates that the same `infer()` call automatically dispatches to
# the correct inference backend based on the model's formalism:
#   1. ODE model (differentiable)  → NUTS via Turing.jl
#   2. SSA model (simulation-based) → ABC-SMC (likelihood-free)
#
# These are two independent parameter recovery experiments, not a comparison.

using InferCell
using JumpProcesses
using CairoMakie
using Random
using Statistics: mean, quantile

Random.seed!(42)

tspan = (0.0, 50.0)
obs_times = collect(0.0:2.5:50.0)

# ============================================================
# Part 1: ODE model → infer() dispatches to NUTS
# ============================================================
println("=" ^ 60)
println("Part 1: ODE model (TranscriptionTranslation)")
println("  formalism = :ode → inference_mode = :differentiable → NUTS")
println("=" ^ 60); flush(stdout)

ode_true = (k_tx=1.0, k_tl=2.0, gamma_mRNA=0.5, gamma_protein=0.1)
ode_model = TranscriptionTranslation(; ode_true..., sigma_obs=0.3)
ode_prob = build_problem(ode_model; tspan=tspan)
ode_sol = solve(ode_prob, Tsit5(); saveat=obs_times)
ode_data = observe(ode_sol, obs_times, ode_model; sigma=0.3)

println("Starting inference (1000 draws)..."); flush(stdout)
t1 = time()
chain = infer(ode_model, ode_data; n_samples=1000)
println("Inference complete in $(round(time()-t1; digits=1))s"); flush(stdout)

println("\n--- ODE Parameter Recovery ---")
for (name, tv) in pairs(ode_true)
    samples = vec(chain[string(name)].data)
    lo, hi = quantile(samples, [0.05, 0.95])
    status = lo <= tv <= hi ? "PASS" : "FAIL"
    println("  $name: true=$tv, median=$(round(quantile(samples, 0.5); digits=3)), 90% CI=[$(round(lo;digits=3)), $(round(hi;digits=3))] → $status")
end
flush(stdout)

# ============================================================
# Part 2: SSA model → infer() dispatches to ABC-SMC
# ============================================================
println("\n" * "=" ^ 60)
println("Part 2: SSA model (StochasticGeneExpression)")
println("  formalism = :jump → inference_mode = :simulation → ABC-SMC")
println("=" ^ 60)

ssa_true = (k_tx=1.0, k_tl=2.0, gamma_mRNA=0.5, gamma_protein=0.1)
ssa_model = StochasticGeneExpression(; ssa_true...)
ssa_prob = build_problem(ssa_model; tspan=tspan)

n_data_replicates = 200
println("Generating observed data ($n_data_replicates SSA replicates)..."); flush(stdout)
ssa_trajectories = [solve(ssa_prob, SSAStepper(); saveat=obs_times) for _ in 1:n_data_replicates]
ssa_data = observe(ssa_trajectories, obs_times, ssa_model)

println("Starting inference (300 particles, 5 populations)...")
flush(stdout)
t2 = time()
abc_result = infer(ssa_model, ssa_data;
                   n_samples=300, n_populations=5, alpha=0.5,
                   n_replicates=50, tspan=tspan, verbose=true)
println("Inference complete in $(round(time()-t2; digits=1))s")

println("\n--- SSA Parameter Recovery ---")
for (j, (name, tv)) in enumerate(pairs(ssa_true))
    weighted_mean = sum(abc_result.weights .* abc_result.particles[j, :])
    sorted_idx = sortperm(abc_result.particles[j, :])
    cum_w = cumsum(abc_result.weights[sorted_idx])
    lo_idx = sorted_idx[findfirst(x -> x >= 0.05, cum_w)]
    hi_idx = sorted_idx[findfirst(x -> x >= 0.95, cum_w)]
    lo = abc_result.particles[j, lo_idx]
    hi = abc_result.particles[j, hi_idx]
    status = lo <= tv <= hi ? "PASS" : "FAIL"
    println("  $name: true=$tv, wmean=$(round(weighted_mean; digits=3)), 90% CI=[$(round(lo;digits=3)), $(round(hi;digits=3))] → $status")
end

# ============================================================
# Figure: Independent parameter recovery for each backend
# ============================================================
println("\nGenerating figure...")
mkpath("examples/figures")

param_names = [:k_tx, :k_tl, :gamma_mRNA, :gamma_protein]
param_labels = ["k_tx", "k_tl", "γ_mRNA", "γ_protein"]

fig = Figure(size=(1000, 600))

# Row 1: ODE / NUTS posteriors
for (col, (pname, lab)) in enumerate(zip(param_names, param_labels))
    tv = ode_true[pname]
    ax = Axis(fig[1, col]; xlabel=lab, ylabel=col == 1 ? "Density" : "")
    samples = vec(chain[string(pname)].data)
    hist!(ax, samples; bins=30, normalization=:pdf, color=(:steelblue, 0.6))
    vlines!(ax, [tv]; color=:black, linewidth=2, linestyle=:dash)
end

# Row 2: SSA / ABC-SMC posteriors
for (col, (pname, lab)) in enumerate(zip(param_names, param_labels))
    tv = ssa_true[pname]
    j = col
    ax = Axis(fig[2, col]; xlabel=lab, ylabel=col == 1 ? "Density" : "")
    abc_samples = abc_result.particles[j, :]
    hist!(ax, abc_samples; bins=30, normalization=:pdf, color=(:coral, 0.6))
    vlines!(ax, [tv]; color=:black, linewidth=2, linestyle=:dash)
end

# Row labels
Label(fig[1, 0], "ODE → NUTS"; fontsize=14, font=:bold, rotation=π/2, tellheight=false)
Label(fig[2, 0], "SSA → ABC-SMC"; fontsize=14, font=:bold, rotation=π/2, tellheight=false)

Label(fig[0, :],
      "infer(): Automatic backend dispatch by model formalism";
      fontsize=16, font=:bold)

# Legend (shared)
legend_elems = [PolyElement(color=(:steelblue, 0.6)), PolyElement(color=(:coral, 0.6)),
                LineElement(color=:black, linewidth=2, linestyle=:dash)]
Legend(fig[3, :], legend_elems, ["NUTS posterior", "ABC-SMC posterior", "True value"];
       orientation=:horizontal, tellwidth=false)

save("examples/figures/phase2_parameter_recovery.png", fig; px_per_unit=2)
println("Saved examples/figures/phase2_parameter_recovery.png")

println("\nPhase 2 demo complete.")
