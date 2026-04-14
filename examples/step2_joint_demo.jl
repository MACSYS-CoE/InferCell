# InferCell Step 2 Demo: Joint ODE Inference (TX/TL + Metabolism)
#
# End-to-end workflow:
#   1. Compose sub-models
#   2. Forward simulation
#   3. Generate synthetic data
#   4. Structural identifiability check
#   5. Joint Bayesian inference (NUTS)
#   6. Posterior predictive checks
#   7. Diagnostic plots

using InferCell
using CairoMakie
using Random
using Statistics: quantile

Random.seed!(42)

# --- 1. Compose sub-models ---
txl = TranscriptionTranslation()
metab = LightMetabolism()
models = [txl, metab]

println("Composed system: $(length(reduce(vcat, states.(models)))) states, " *
        "$(length(unique_params(reduce(vcat, model_free_params.(parameters.(models)))))) ODE params")

# --- 2. Forward simulation ---
prob = build_problem(models; tspan=(0.0, 50.0))
true_sol = solve(prob, Tsit5())

println("\nForward simulation complete.")
state_names = reduce(vcat, states.(models))
final = Array(true_sol)[:, end]
for (s, v) in zip(state_names, final)
    println("  $s steady state: $(round(v; digits=2))")
end

# --- 3. Synthetic data ---
times = 0.0:2.5:50.0
data = observe(true_sol, times, models; sigma=0.3)

println("\nSynthetic data: $(size(data.observations, 2)) timepoints, $(size(data.observations, 1)) species")

# --- 4. Structural identifiability check ---
ident = check_identifiability(models, prob, collect(times))
println("\nIdentifiability check:")
println("  Sensitivity matrix rank: $(ident.rank) / $(ident.n_params) parameters")
println("  Full rank: $(ident.full_rank)")

# --- 5. Joint Bayesian inference ---
println("\nStarting joint NUTS sampling (1000 draws, 8 parameters)...")
t_start = time()
chain = infer(models, data; n_samples=1000)
elapsed = time() - t_start
println("Sampling complete in $(round(elapsed; digits=1)) seconds.")

# Parameter recovery check
true_vals = Dict(
    "k_tx"          => 1.0,
    "k_tl"          => 2.0,
    "gamma_mRNA"    => 0.5,
    "gamma_protein" => 0.1,
    "k_atp"         => 1.0,
    "k_ntp"         => 0.2,
    "k_aa"          => 0.4,
    "sigma_obs"     => 0.3,
)

println("\n--- Parameter Recovery (90% CI) ---")
all_pass = true
for (name, tv) in sort(collect(true_vals); by=first)
    samples = vec(chain[name].data)
    lo, hi = quantile(samples, [0.05, 0.95])
    recovered = lo <= tv <= hi
    global all_pass = all_pass && recovered
    status = recovered ? "PASS" : "FAIL"
    println("  $name: true=$tv, 90% CI=[$(round(lo;digits=3)), $(round(hi;digits=3))] → $status")
end
println(all_pass ? "\nALL PARAMETERS RECOVERED" : "\nSOME PARAMETERS NOT RECOVERED")
println("Wall-clock: $(round(elapsed; digits=1)) seconds")

# --- 6. Posterior predictive ---
ppc = posterior_predictive(models, chain;
    tspan=(0.0, 50.0), saveat=range(0, 50; length=200), n_samples=100)

# --- 7. Diagnostic plots ---
mkpath("examples/figures")

# -- Plot 1: Posterior distributions (4x2 grid) --
param_names  = ["gamma_mRNA", "gamma_protein", "k_aa", "k_atp",
                "k_ntp", "k_tl", "k_tx", "sigma_obs"]
param_true   = [0.5, 0.1, 0.4, 1.0, 0.2, 2.0, 1.0, 0.3]
param_labels = ["gamma_mRNA", "gamma_protein", "k_aa", "k_atp",
                "k_ntp", "k_tl", "k_tx", "sigma_obs"]

fig1 = Figure(size=(1000, 800))
for (i, (pname, tv, lab)) in enumerate(zip(param_names, param_true, param_labels))
    row, col = fldmod1(i, 2)
    ax = Axis(fig1[row, col]; xlabel=lab, ylabel="Density")
    samples = vec(chain[pname].data)
    hist!(ax, samples; bins=40, normalization=:pdf, color=(:steelblue, 0.6))
    vlines!(ax, [tv]; color=:red, linewidth=2)
    lo, hi = quantile(samples, [0.05, 0.95])
    vspan!(ax, lo, hi; color=(:orange, 0.15))
end
Label(fig1[0, :], "Joint Posteriors (red = true, orange = 90% CI)";
      fontsize=16, font=:bold)
save("examples/figures/step2_posteriors.png", fig1; px_per_unit=2)
println("\nSaved examples/figures/step2_posteriors.png")

# -- Plot 2: Posterior predictive trajectories (5 panels) --
species_names = ["mRNA", "protein", "ATP", "NTP", "AA"]

fig2 = Figure(size=(1200, 800))
for (i, sp) in enumerate(species_names)
    row, col = fldmod1(i, 3)
    ax = Axis(fig2[row, col]; xlabel="Time", ylabel="Concentration", title=sp)

    for sol in ppc.solutions
        lines!(ax, ppc.times, Array(sol)[i, :]; color=(:steelblue, 0.08))
    end

    # True trajectory
    true_fine = solve(prob, Tsit5(); saveat=ppc.times)
    lines!(ax, ppc.times, Array(true_fine)[i, :]; color=:red, linewidth=2, label="True")

    # Noisy data
    scatter!(ax, collect(times), data.observations[i, :]; color=:black, markersize=6, label="Data")

    if i == 1
        axislegend(ax; position=:rb)
    end
end
Label(fig2[0, :], "Posterior Predictive (blue = samples, red = truth)";
      fontsize=16, font=:bold)
save("examples/figures/step2_posterior_predictive.png", fig2; px_per_unit=2)
println("Saved examples/figures/step2_posterior_predictive.png")

# -- Plot 3: Trace plots (4x2 grid) --
fig3 = Figure(size=(1000, 800))
for (i, (pname, tv, lab)) in enumerate(zip(param_names, param_true, param_labels))
    row, col = fldmod1(i, 2)
    ax = Axis(fig3[row, col]; xlabel="Iteration", ylabel=lab)
    samples = vec(chain[pname].data)
    lines!(ax, 1:length(samples), samples; color=(:steelblue, 0.6), linewidth=0.5)
    hlines!(ax, [tv]; color=:red, linewidth=1.5)
end
Label(fig3[0, :], "Trace Plots (red = true value)"; fontsize=16, font=:bold)
save("examples/figures/step2_traces.png", fig3; px_per_unit=2)
println("Saved examples/figures/step2_traces.png")

println("\nStep 2 joint inference demo complete.")
