using OrdinaryDiffEq
using SciMLSensitivity
using SciMLBase: ReturnCode
using Turing
using Distributions
using StaticArrays
using Random
using CairoMakie

Random.seed!(42)

# --- ODE + synthetic data (same as spike) ---
function txl_dynamics(u, p, t)
    mRNA, protein = u
    k_tx, k_tl, γ_mRNA, γ_protein = p
    SA[k_tx - γ_mRNA * mRNA, k_tl * mRNA - γ_protein * protein]
end

true_params = SA[1.0, 2.0, 0.5, 0.1]
σ_true = 0.3
u0 = SA[0.0, 0.0]
tspan = (0.0, 50.0)
times = collect(0.0:2.5:50.0)

prob = ODEProblem(txl_dynamics, u0, tspan, true_params)
sol_true = solve(prob, Tsit5(); saveat=times)
data = Array(sol_true) .+ σ_true .* randn(size(Array(sol_true)))

# --- Turing model + sample ---
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

println("Sampling...")
model = txl_inference(data, times, prob)
chain = sample(model, NUTS(), 1000)

mkpath("spike/figures")

# ── Plot 1: Posterior distributions ──
param_names = [:k_tx, :k_tl, :γ_mRNA, :γ_protein, :σ]
true_vals   = [1.0, 2.0, 0.5, 0.1, 0.3]
labels      = ["k_tx", "k_tl", "γ_mRNA", "γ_protein", "σ"]

fig1 = Figure(size=(1000, 700))
for (i, (pname, tv, lab)) in enumerate(zip(param_names, true_vals, labels))
    row, col = fldmod1(i, 3)
    ax = Axis(fig1[row, col]; xlabel=lab, ylabel="Density")
    samples = vec(chain[pname])
    hist!(ax, samples; bins=40, normalization=:pdf, color=(:steelblue, 0.6))
    vlines!(ax, [tv]; color=:red, linewidth=2, label="true")
    lo, hi = quantile(samples, [0.05, 0.95])
    vspan!(ax, lo, hi; color=(:orange, 0.15))
end
Label(fig1[0, :], "Posterior Distributions (red = true value, orange band = 90% CI)";
      fontsize=16, font=:bold)
save("spike/figures/posteriors.png", fig1; px_per_unit=2)
println("Saved spike/figures/posteriors.png")

# ── Plot 2: Posterior predictive trajectories ──
fig2 = Figure(size=(900, 400))
ax_m = Axis(fig2[1, 1]; xlabel="Time", ylabel="Concentration", title="mRNA")
ax_p = Axis(fig2[1, 2]; xlabel="Time", ylabel="Concentration", title="Protein")

# Draw 100 posterior predictive trajectories
n_draws = 100
idxs = rand(1:length(chain), n_draws)
t_fine = range(0, 50; length=200)
for idx in idxs
    p_draw = SA[chain[:k_tx].data[idx],
                chain[:k_tl].data[idx],
                chain[:γ_mRNA].data[idx],
                chain[:γ_protein].data[idx]]
    sol_draw = solve(remake(prob, p=p_draw), Tsit5(); saveat=t_fine)
    lines!(ax_m, t_fine, sol_draw[1, :]; color=(:steelblue, 0.08))
    lines!(ax_p, t_fine, sol_draw[2, :]; color=(:steelblue, 0.08))
end

# True trajectory
lines!(ax_m, t_fine, solve(prob, Tsit5(); saveat=t_fine)[1, :];
       color=:red, linewidth=2, label="True")
lines!(ax_p, t_fine, solve(prob, Tsit5(); saveat=t_fine)[2, :];
       color=:red, linewidth=2, label="True")

# Noisy data
scatter!(ax_m, times, data[1, :]; color=:black, markersize=6, label="Data")
scatter!(ax_p, times, data[2, :]; color=:black, markersize=6, label="Data")

axislegend(ax_m; position=:rb)
axislegend(ax_p; position=:rb)
Label(fig2[0, :], "Posterior Predictive Trajectories (blue = samples, red = truth)";
      fontsize=16, font=:bold)
save("spike/figures/posterior_predictive.png", fig2; px_per_unit=2)
println("Saved spike/figures/posterior_predictive.png")

# ── Plot 3: Trace plots ──
fig3 = Figure(size=(1000, 700))
for (i, (pname, tv, lab)) in enumerate(zip(param_names, true_vals, labels))
    row, col = fldmod1(i, 3)
    ax = Axis(fig3[row, col]; xlabel="Iteration", ylabel=lab)
    samples = vec(chain[pname])
    lines!(ax, 1:length(samples), samples; color=(:steelblue, 0.6), linewidth=0.5)
    hlines!(ax, [tv]; color=:red, linewidth=1.5)
end
Label(fig3[0, :], "Trace Plots (red = true value)"; fontsize=16, font=:bold)
save("spike/figures/traces.png", fig3; px_per_unit=2)
println("Saved spike/figures/traces.png")

println("\nDone — all plots in spike/figures/")
