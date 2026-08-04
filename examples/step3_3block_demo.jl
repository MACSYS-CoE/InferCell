# InferCell Step 3: Three-Block Boundary Protocol
# [TX/TL + Metabolism] (joint ODE/NUTS) → cut posterior → [Gene Expr] (SSA/ABC-SMC)

using InferCell
using CairoMakie
using Random
using Statistics: quantile, mean, std
using MCMCChains: Chains
using Distributions: LogNormal, pdf
using Serialization

Random.seed!(42)

# --- 1. Define models ---
true_rates = (k_tx=1.0, k_tl=2.0, gamma_mRNA=0.5, gamma_protein=0.1)

txl = TranscriptionTranslation(; true_rates..., sigma_obs=0.3)
metab = LightMetabolism()  # uses default k_atp=1.0, k_ntp=0.2, k_aa=0.4, shares k_tx, k_tl
sge = StochasticGeneExpression(; true_rates...)

ode_models = [txl, metab]
println("ODE block: TX/TL + Metabolism ($(length(reduce(vcat, states.(ode_models)))) states)")
println("SSA block: StochasticGeneExpression ($(length(states(sge))) states)")
flush(stdout)

# --- 2. Generate data ---
# ODE data (5 species: mRNA, protein, ATP, NTP, AA)
ode_prob = build_problem(ode_models; tspan=(0.0, 50.0))
ode_sol = solve(ode_prob, Tsit5())
ode_times = 0.0:2.5:50.0
ode_data = observe(ode_sol, ode_times, ode_models; sigma=0.3)
println("ODE data: $(size(ode_data.observations)) — species: $(ode_data.species)")
flush(stdout)

# SSA data (2 species: mRNA, protein — discrete)
ssa_tspan = (0.0, 100.0)
ssa_times = collect(0.0:10.0:100.0)
ssa_prob = build_problem(sge; tspan=ssa_tspan)
trajectories = [solve(ssa_prob, SSAStepper(); saveat=ssa_times) for _ in 1:200]
ssa_data = observe(trajectories, ssa_times, sge)
println("SSA data: $(size(ssa_data.observations)) — species: $(ssa_data.species)")
flush(stdout)

# --- 3. Joint ODE inference (NUTS on 5-state, 8-param system) ---
println("\n--- Stage 1: Joint ODE inference (NUTS, 1000 draws, 8 params) ---")
flush(stdout)
t_ode = time()
ode_chain = infer(ode_models, ode_data; n_samples=1000)
elapsed_ode = time() - t_ode
println("ODE inference complete in $(round(elapsed_ode; digits=1)) seconds")
flush(stdout)

# ODE recovery
ode_true = Dict("k_tx" => 1.0, "k_tl" => 2.0, "gamma_mRNA" => 0.5,
                "gamma_protein" => 0.1, "k_atp" => 1.0, "k_ntp" => 0.2,
                "k_aa" => 0.4, "sigma_obs" => 0.3)
println("\nODE Parameter Recovery (90% CI):")
for (name, tv) in sort(collect(ode_true); by=first)
    samples = vec(ode_chain[name].data)
    lo, hi = quantile(samples, [0.05, 0.95])
    status = lo <= tv <= hi ? "PASS" : "FAIL"
    println("  $name: true=$tv, 90% CI=[$(round(lo;digits=3)), $(round(hi;digits=3))] -> $status")
end
flush(stdout)

# --- 4. Boundary conditioning ---
println("\n--- Stage 2: Boundary conditioning (KDE) ---")
bc = boundary_condition(ode_chain, sge; method=:kde)
println("Conditioned $(sum(bc.conditioned)) / $(length(bc.param_names)) SSA parameters:")
for (name, cond) in zip(bc.param_names, bc.conditioned)
    println("  $name: $(cond ? "KDE from ODE posterior" : "original prior")")
end
flush(stdout)

# --- 5. Conditioned ABC-SMC ---
println("\n--- Stage 3: Conditioned ABC-SMC (200 particles, 6 populations) ---")
flush(stdout)
t_cond = time()
cond_result = InferCell._infer_abc([sge], ssa_data;
    priors=bc.priors, param_names=bc.param_names,
    n_particles=200, n_populations=6,
    n_replicates=50, tspan=ssa_tspan, verbose=true)
elapsed_cond = time() - t_cond
println("Conditioned ABC-SMC complete in $(round(elapsed_cond; digits=1)) seconds")
flush(stdout)

# --- 6. Unconditioned ABC-SMC ---
println("\n--- Baseline: Unconditioned ABC-SMC (same settings) ---")
flush(stdout)
t_uncond = time()
uncond_result = infer(sge, ssa_data;
    n_samples=200, n_populations=6,
    n_replicates=50, tspan=ssa_tspan, verbose=true)
elapsed_uncond = time() - t_uncond
println("Unconditioned ABC-SMC complete in $(round(elapsed_uncond; digits=1)) seconds")
flush(stdout)

# --- Save all inference results ---
mkpath("examples/results")
results = Dict(
    "ode_chain" => ode_chain,
    "cond_result" => cond_result,
    "uncond_result" => uncond_result,
    "boundary" => bc,
    "ode_data" => ode_data,
    "ssa_data" => ssa_data,
    "ode_true" => ode_true,
    "true_rates" => true_rates,
    "elapsed_ode" => elapsed_ode,
    "elapsed_cond" => elapsed_cond,
    "elapsed_uncond" => elapsed_uncond,
)
serialize("examples/results/step3_3block_results.jls", results)
println("\nSaved inference results to examples/results/step3_3block_results.jls")
flush(stdout)

# --- 7. Results ---
param_names_str = ["k_tx", "k_tl", "gamma_mRNA", "gamma_protein"]
param_true_vals = [1.0, 2.0, 0.5, 0.1]

println("\n--- Parameter Recovery Comparison ---")
println("Parameter       | True  | Cond. mean (var)        | Uncond. mean (var)")
println(repeat("-", 75))
for (j, (name, tv)) in enumerate(zip(param_names_str, param_true_vals))
    cmean = sum(cond_result.weights .* cond_result.particles[j, :])
    cvar = sum(cond_result.weights .* (cond_result.particles[j, :] .- cmean).^2)
    umean = sum(uncond_result.weights .* uncond_result.particles[j, :])
    uvar = sum(uncond_result.weights .* (uncond_result.particles[j, :] .- umean).^2)
    ratio = uvar > 0 ? round(uvar / cvar; digits=1) : Inf
    println("  $(rpad(name, 15))| $tv   | $(round(cmean;digits=3)) ($(round(cvar;digits=4))) | $(round(umean;digits=3)) ($(round(uvar;digits=4)))  var_ratio=$(ratio)x")
end
println("\nWall-clock: ODE=$(round(elapsed_ode;digits=1))s, " *
        "Cond.ABC=$(round(elapsed_cond;digits=1))s, " *
        "Uncond.ABC=$(round(elapsed_uncond;digits=1))s, " *
        "Total=$(round(elapsed_ode+elapsed_cond+elapsed_uncond;digits=1))s")
flush(stdout)

# --- 8. Figures ---
mkpath("examples/figures")

# Plot 1: All 8 ODE posteriors
all_names  = ["gamma_mRNA", "gamma_protein", "k_aa", "k_atp", "k_ntp", "k_tl", "k_tx", "sigma_obs"]
all_true   = [0.5, 0.1, 0.4, 1.0, 0.2, 2.0, 1.0, 0.3]

fig1 = Figure(size=(1000, 800))
for (i, (pname, tv)) in enumerate(zip(all_names, all_true))
    row, col = fldmod1(i, 2)
    ax = Axis(fig1[row, col]; xlabel=pname, ylabel="Density")
    samples = vec(ode_chain[pname].data)
    hist!(ax, samples; bins=40, normalization=:pdf, color=(:steelblue, 0.6))
    vlines!(ax, [tv]; color=:red, linewidth=2)
    lo, hi = quantile(samples, [0.05, 0.95])
    vspan!(ax, lo, hi; color=(:orange, 0.15))
end
Label(fig1[0, :], "Joint ODE Posteriors: TX/TL + Metabolism (red = true, orange = 90% CI)";
      fontsize=14, font=:bold)
save("examples/figures/step3_3block_ode_posteriors.png", fig1; px_per_unit=2)
println("\nSaved examples/figures/step3_3block_ode_posteriors.png")

# Plot 2: Boundary priors vs originals (4 shared params)
fig2 = Figure(size=(900, 600))
x_ranges = [(0.0, 4.0), (0.0, 8.0), (0.0, 3.0), (0.0, 1.0)]
for (i, (name, tv, xr)) in enumerate(zip(param_names_str, param_true_vals, x_ranges))
    row, col = fldmod1(i, 2)
    ax = Axis(fig2[row, col]; xlabel=name, ylabel="Density")
    xs = range(xr[1], xr[2]; length=200)
    original = LogNormal(0, 1)
    lines!(ax, xs, [pdf(original, x) for x in xs]; color=:gray, linewidth=2, linestyle=:dash, label="Original prior")
    kde_prior = bc.priors[i]
    lines!(ax, xs, [pdf(kde_prior, x) for x in xs]; color=:steelblue, linewidth=2, label="KDE prior")
    vlines!(ax, [tv]; color=:red, linewidth=2, label="True value")
    i == 1 && axislegend(ax; position=:rt)
end
Label(fig2[0, :], "Boundary Conditioning: Original vs KDE-Fitted Priors (from joint ODE posterior)";
      fontsize=14, font=:bold)
save("examples/figures/step3_3block_boundary_priors.png", fig2; px_per_unit=2)
println("Saved examples/figures/step3_3block_boundary_priors.png")

# Plot 3: Conditioned vs unconditioned SSA posteriors
fig3 = Figure(size=(900, 600))
for (i, (name, tv)) in enumerate(zip(param_names_str, param_true_vals))
    row, col = fldmod1(i, 2)
    ax = Axis(fig3[row, col]; xlabel=name, ylabel="Density")
    n_resample = 2000
    uncond_idx = [InferCell._weighted_sample(uncond_result.weights, Random.default_rng()) for _ in 1:n_resample]
    uncond_samples = [uncond_result.particles[i, k] for k in uncond_idx]
    cond_idx = [InferCell._weighted_sample(cond_result.weights, Random.default_rng()) for _ in 1:n_resample]
    cond_samples = [cond_result.particles[i, k] for k in cond_idx]
    hist!(ax, uncond_samples; bins=40, normalization=:pdf, color=(:gray, 0.4), label="Unconditioned")
    hist!(ax, cond_samples; bins=40, normalization=:pdf, color=(:steelblue, 0.5), label="Conditioned")
    vlines!(ax, [tv]; color=:red, linewidth=2, label="True value")
    i == 1 && axislegend(ax; position=:rt)
end
Label(fig3[0, :], "SSA Posteriors: Conditioned on Joint ODE (blue) vs Unconditioned (gray)";
      fontsize=14, font=:bold)
save("examples/figures/step3_3block_ssa_comparison.png", fig3; px_per_unit=2)
println("Saved examples/figures/step3_3block_ssa_comparison.png")

println("\nStep 3 three-block demo complete.")
