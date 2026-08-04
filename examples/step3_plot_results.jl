# Step 3 publication-quality figures
# Loads serialized results — no inference needed

using CairoMakie
using LaTeXStrings
using Serialization
using Statistics: quantile, mean, std
using StatsBase: fit, Histogram
using Random
using InferCell
using PairPlots

set_theme!(theme_latexfonts())

# --- Load results ---
results = deserialize("examples/results/step3_3block_results.jls")
ode_chain    = results["ode_chain"]
cond_result  = results["cond_result"]
uncond_result = results["uncond_result"]
bc           = results["boundary"]
ode_true     = results["ode_true"]
true_rates   = results["true_rates"]

mkpath("examples/figures")
cols = Makie.wong_colors()

# ============================================================
# Figure 1: Corner plot — joint ODE posteriors (8 parameters)
# ============================================================

param_names = [:k_tx, :k_tl, :gamma_mRNA, :gamma_protein,
               :k_atp, :k_ntp, :k_aa, :sigma_obs]

param_labels = Dict{Symbol,Any}(
    :k_tx          => L"k_{tx}",
    :k_tl          => L"k_{tl}",
    :gamma_mRNA    => L"\gamma_{mRNA}",
    :gamma_protein => L"\gamma_{prot}",
    :k_atp         => L"k_{atp}",
    :k_ntp         => L"k_{ntp}",
    :k_aa          => L"k_{aa}",
    :sigma_obs     => L"\sigma",
)

truth_table = (
    k_tx = 1.0, k_tl = 2.0, gamma_mRNA = 0.5, gamma_protein = 0.1,
    k_atp = 1.0, k_ntp = 0.2, k_aa = 0.4, sigma_obs = 0.3,
)

# Build a table (NamedTuple of vectors) for PairPlots
table = (; (k => vec(ode_chain[string(k)].data) for k in param_names)...)

fig1 = pairplot(
    table => (
        PairPlots.Contourf(sigmas=1:3),
        PairPlots.Scatter(filtersigma=3, markersize=1, color=(cols[1], 0.1)),
        PairPlots.MarginDensity(linewidth=1.5),
    ),
    PairPlots.Truth(truth_table, color=cols[2], linewidth=1.5),
    labels=param_labels,
    figure=(; size=(700, 700), figure_padding=(4, 8, 4, 4)),
)

save("examples/figures/step3_corner.png", fig1; px_per_unit=4)
println("Saved examples/figures/step3_corner.png")

# ============================================================
# Figure 2: SSA posteriors — conditioned vs unconditioned
# ============================================================

ssa_names  = ["k_tx", "k_tl", "gamma_mRNA", "gamma_protein"]
ssa_labels = [L"k_{tx}", L"k_{tl}", L"\gamma_{mRNA}", L"\gamma_{prot}"]
ssa_true   = [1.0, 2.0, 0.5, 0.1]

# Resample from weighted particles
n_resample = 5000
rng = Random.MersenneTwister(42)

function weighted_resample(result, param_idx, n, rng)
    idx = [InferCell._weighted_sample(result.weights, rng) for _ in 1:n]
    return [result.particles[param_idx, k] for k in idx]
end

fig2 = Figure(size=(700, 280), figure_padding=(2, 8, 2, 2))

for (i, (name, label, tv)) in enumerate(zip(ssa_names, ssa_labels, ssa_true))
    ax = Axis(fig2[1, i];
        xlabel=label, ylabel = i == 1 ? L"Density" : "",
        topspinevisible=false, rightspinevisible=false,
        xgridvisible=false, ygridvisible=false,
        xlabelsize=14, ylabelsize=14,
        xticklabelsize=12, yticklabelsize=12,
        yticklabelsvisible = i == 1,
    )

    uncond_s = weighted_resample(uncond_result, i, n_resample, rng)
    cond_s   = weighted_resample(cond_result, i, n_resample, rng)

    # Determine shared bin edges from the union of both samples
    lo = min(minimum(uncond_s), minimum(cond_s))
    hi = max(maximum(uncond_s), maximum(cond_s))
    edges = range(lo, hi; length=51)
    bw = edges[2] - edges[1]

    # Unconditioned step histogram
    h_u = fit(Histogram, uncond_s, edges)
    dens_u = h_u.weights ./ (sum(h_u.weights) * bw)
    stairs!(ax, collect(edges), vcat(dens_u, dens_u[end]);
            color=(cols[3], 0.7), linewidth=2.0, step=:post,
            label="Unconditioned")

    # Conditioned step histogram
    h_c = fit(Histogram, cond_s, edges)
    dens_c = h_c.weights ./ (sum(h_c.weights) * bw)
    stairs!(ax, collect(edges), vcat(dens_c, dens_c[end]);
            color=(cols[1], 0.9), linewidth=1.5, step=:post,
            label="Conditioned")

    # True value
    vlines!(ax, [tv]; color=:black, linewidth=1.2, linestyle=:dash,
            label="True value")

    if i == 1
        Legend(fig2[1, 1], ax;
            tellwidth=false, tellheight=false,
            halign=:right, valign=:top,
            margin=(8, 8, 8, 8),
            padding=(4, 4, 3, 3),
            labelsize=10,
            framevisible=false,
        )
    end
end

save("examples/figures/step3_ssa_posteriors.png", fig2; px_per_unit=4)
println("Saved examples/figures/step3_ssa_posteriors.png")

println("Done.")
