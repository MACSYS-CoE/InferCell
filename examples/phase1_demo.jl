# InferCell Phase 1 Demo: Transcription-Translation Inference
#
# End-to-end workflow:
#   1. Define sub-model
#   2. Forward simulation
#   3. Generate synthetic data
#   4. Bayesian parameter inference (NUTS)
#   5. Posterior predictive checks
#   6. Diagnostic plots

using InferCell
using CairoMakie
using Random
using Statistics: quantile

Random.seed!(42)

# --- 1. Define sub-model ---
txl = TranscriptionTranslation()

# --- 2. Forward simulation ---
prob = build_problem(txl; tspan=(0.0, 50.0))
true_sol = solve(prob, Tsit5())

println("Forward simulation complete.")
println("  mRNA steady state:    $(round(Array(true_sol)[1, end]; digits=2))")
println("  Protein steady state: $(round(Array(true_sol)[2, end]; digits=2))")

# --- 3. Synthetic data ---
times = 0.0:2.5:50.0
data = observe(true_sol, times, txl; sigma=0.3)

println("Synthetic data: $(size(data.observations, 2)) timepoints, $(size(data.observations, 1)) species")

# --- 4. Structural identifiability check ---
ident = check_identifiability(txl, prob, collect(times))
println("\nIdentifiability check:")
println("  Sensitivity matrix rank: $(ident.rank) / $(ident.n_params) parameters")
println("  Full rank: $(ident.full_rank)")

# --- 5. Bayesian inference ---
println("\nStarting NUTS sampling (1000 draws)...")
t_start = time()
chain = infer(txl, data; n_samples=1000)
elapsed = time() - t_start
println("Sampling complete in $(round(elapsed; digits=1)) seconds.")

# Parameter recovery check
true_vals = Dict(
    "k_tx"          => 1.0,
    "k_tl"          => 2.0,
    "gamma_mRNA"    => 0.5,
    "gamma_protein" => 0.1,
    "sigma_obs"     => 0.3,
)

println("\n--- Parameter Recovery (90% CI) ---")
all_pass = true
for (name, tv) in true_vals
    samples = vec(chain[name].data)
    lo, hi = quantile(samples, [0.05, 0.95])
    recovered = lo <= tv <= hi
    global all_pass = all_pass && recovered
    status = recovered ? "PASS" : "FAIL"
    println("  $name: true=$tv, 90% CI=[$(round(lo;digits=3)), $(round(hi;digits=3))] → $status")
end
println(all_pass ? "\n✓ ALL PARAMETERS RECOVERED" : "\n✗ SOME PARAMETERS NOT RECOVERED")

# --- 6. Posterior predictive ---
ppc = posterior_predictive(txl, chain;
    tspan=(0.0, 50.0), saveat=range(0, 50; length=200), n_samples=100)

# --- 7. Diagnostic plots ---
mkpath("examples/figures")

# -- Plot 1: Posterior distributions --
param_names = ["k_tx", "k_tl", "gamma_mRNA", "gamma_protein", "sigma_obs"]
param_true  = [1.0,     2.0,    0.5,          0.1,             0.3]
param_labels = ["k_tx", "k_tl", "γ_mRNA", "γ_protein", "σ_obs"]

fig1 = Figure(size=(1000, 700))
for (i, (pname, tv, lab)) in enumerate(zip(param_names, param_true, param_labels))
    row, col = fldmod1(i, 3)
    ax = Axis(fig1[row, col]; xlabel=lab, ylabel="Density")
    samples = vec(chain[pname].data)
    hist!(ax, samples; bins=40, normalization=:pdf, color=(:steelblue, 0.6))
    vlines!(ax, [tv]; color=:red, linewidth=2)
    lo, hi = quantile(samples, [0.05, 0.95])
    vspan!(ax, lo, hi; color=(:orange, 0.15))
end
Label(fig1[0, :], "Posterior Distributions (red = true value, orange = 90% CI)";
      fontsize=16, font=:bold)
save("examples/figures/posteriors.png", fig1; px_per_unit=2)
println("\nSaved examples/figures/posteriors.png")

# -- Plot 2: Posterior predictive trajectories --
fig2 = Figure(size=(900, 400))
ax_m = Axis(fig2[1, 1]; xlabel="Time", ylabel="Concentration", title="mRNA")
ax_p = Axis(fig2[1, 2]; xlabel="Time", ylabel="Concentration", title="Protein")

for sol in ppc.solutions
    lines!(ax_m, ppc.times, Array(sol)[1, :]; color=(:steelblue, 0.08))
    lines!(ax_p, ppc.times, Array(sol)[2, :]; color=(:steelblue, 0.08))
end

# True trajectory
t_fine = ppc.times
true_fine = solve(prob, Tsit5(); saveat=t_fine)
lines!(ax_m, t_fine, Array(true_fine)[1, :]; color=:red, linewidth=2, label="True")
lines!(ax_p, t_fine, Array(true_fine)[2, :]; color=:red, linewidth=2, label="True")

# Noisy data
scatter!(ax_m, collect(times), data.observations[1, :]; color=:black, markersize=6, label="Data")
scatter!(ax_p, collect(times), data.observations[2, :]; color=:black, markersize=6, label="Data")

axislegend(ax_m; position=:rb)
axislegend(ax_p; position=:rb)
Label(fig2[0, :], "Posterior Predictive (blue = samples, red = truth)";
      fontsize=16, font=:bold)
save("examples/figures/posterior_predictive.png", fig2; px_per_unit=2)
println("Saved examples/figures/posterior_predictive.png")

# -- Plot 3: Trace plots --
fig3 = Figure(size=(1000, 700))
for (i, (pname, tv, lab)) in enumerate(zip(param_names, param_true, param_labels))
    row, col = fldmod1(i, 3)
    ax = Axis(fig3[row, col]; xlabel="Iteration", ylabel=lab)
    samples = vec(chain[pname].data)
    lines!(ax, 1:length(samples), samples; color=(:steelblue, 0.6), linewidth=0.5)
    hlines!(ax, [tv]; color=:red, linewidth=1.5)
end
Label(fig3[0, :], "Trace Plots (red = true value)"; fontsize=16, font=:bold)
save("examples/figures/traces.png", fig3; px_per_unit=2)
println("Saved examples/figures/traces.png")

println("\nPhase 1 demo complete.")
