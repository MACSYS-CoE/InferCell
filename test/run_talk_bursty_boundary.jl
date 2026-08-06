# ISAB talk figure: information flow across a module boundary in the
# bursty-regulator cell.
#
# Three distinct modules, not one subsystem described twice:
#   1. LightMetabolism            (ODE)  ATP/NTP/AA, enzyme-gated production
#   2. TranscriptionTranslation   (ODE)  enzyme gene; its protein feeds (1)
#   3. BurstyGeneExpression       (SSA)  telegraph promoter, discrete counts
#
# (1) and (2) are state-coupled both ways. (3) is parameter-coupled to (2)
# through k_tl, gamma_mRNA, gamma_protein -- the shared ribosomes and shared
# degradation machinery. 3 shared of 11 unique.
#
# Runs ABC-SMC on the SSA block twice: once under its own priors, once under
# KDE priors fitted to the NUTS posterior of the ODE block. The expected result
# is that the three SHARED parameters sharpen onto truth while k_on, k_off and
# k_tx_burst stay as broad as their priors -- nothing in the bulk metabolic data
# knows about promoter switching. That contrast is the figure.
#
# Exports CSV for dev/talks/ISAB/figures/plot_bursty_boundary.py.

using InferCell
using OrdinaryDiffEq: Tsit5, solve
using JumpProcesses: SSAStepper
using Random
using Statistics: mean, std
using Dates

println("=== ISAB boundary figure: bursty-regulator cell ===")
println("Start: $(Dates.now())")
flush(stdout)

Random.seed!(2026)

const OUTDIR = get(ENV, "INFERCELL_OUTDIR", "examples/results")
mkpath(OUTDIR)

# ABC settings. The proven configuration in run_integ_iterative.jl is 200
# particles; 400 here because this is a 6-dimensional posterior being shown at
# projector size and the contours need to be smooth.
const N_PARTICLES  = parse(Int, get(ENV, "INFERCELL_N_PARTICLES", "400"))
const N_POPULATIONS = parse(Int, get(ENV, "INFERCELL_N_POPULATIONS", "6"))
const N_REPLICATES = parse(Int, get(ENV, "INFERCELL_N_REPLICATES", "50"))
const N_NUTS       = parse(Int, get(ENV, "INFERCELL_N_NUTS", "1000"))

println("Settings: particles=$N_PARTICLES populations=$N_POPULATIONS " *
        "replicates=$N_REPLICATES nuts=$N_NUTS")
flush(stdout)

# Ground truth. Same values as run_integ_iterative.jl so the two runs are
# comparable; k_tl, gamma_mRNA and gamma_protein appear in both tuples and must
# agree -- that is what makes them shared.
true_ode = (k_tx=1.0, k_tl=2.0, gamma_mRNA=0.5, gamma_protein=0.1)
true_ssa = (k_on=0.2, k_off=0.5, k_tx_burst=10.0, k_tl=2.0,
            gamma_mRNA=0.5, gamma_protein=0.1)

# --- Differentiable block: metabolism + TX/TL, coupled ---
println("\n1. Differentiable block (metabolism + TX/TL)...")
flush(stdout)
txl = TranscriptionTranslation(; true_ode..., sigma_obs=0.3)
metab = LightMetabolism()  # mRNA_source=:mRNA, enzyme_source=:protein
ode_prob = build_problem([txl, metab]; tspan=(0.0, 50.0))
ode_sol = solve(ode_prob, Tsit5())
ode_data = observe(ode_sol, 0.0:2.5:50.0, [txl, metab]; sigma=0.3)
println("   ODE data: $(size(ode_data.observations))")
flush(stdout)

# --- Simulation block: bursty regulator ---
println("\n2. Simulation block (bursty regulator)...")
flush(stdout)
bge = BurstyGeneExpression(; true_ssa...)
ssa_tspan = (0.0, 100.0)
ssa_times = collect(0.0:10.0:100.0)
ssa_prob = build_problem(bge; tspan=ssa_tspan)
trajectories = [solve(ssa_prob, SSAStepper(); saveat=ssa_times) for _ in 1:200]
ssa_data = observe(trajectories, ssa_times, bge)
println("   SSA data: $(size(ssa_data.observations))")
flush(stdout)

ssa_kwargs = (; n_particles=N_PARTICLES, n_populations=N_POPULATIONS,
                n_replicates=N_REPLICATES, tspan=ssa_tspan)

# --- Conditioned: NUTS on the ODE block, then ABC under boundary priors ---
println("\n3. sequential_infer: NUTS -> boundary -> ABC-SMC...")
flush(stdout)
t0 = time()
cond = sequential_infer([txl, metab], ode_data, [bge], ssa_data;
                        ode_kwargs=(; n_samples=N_NUTS),
                        ssa_kwargs=ssa_kwargs,
                        boundary_method=:kde, verbose=true)
println("   Conditioned run done in $(round(time() - t0; digits=1))s")
flush(stdout)

# --- Unconditioned baseline: same block, same data, own priors ---
println("\n4. Unconditioned ABC-SMC baseline...")
flush(stdout)
t1 = time()
uncond = infer(bge, ssa_data;
               n_samples=N_PARTICLES, n_populations=N_POPULATIONS,
               n_replicates=N_REPLICATES, tspan=ssa_tspan, verbose=true)
println("   Unconditioned run done in $(round(time() - t1; digits=1))s")
flush(stdout)

# Both posteriors must order their rows identically or the figure overlays the
# wrong marginals. Assert rather than assume: sequential_infer derives its order
# from boundary_condition and the baseline derives its own inside _infer_abc.
@assert cond.ssa_posterior.param_names == uncond.param_names """
Parameter orders differ:
  conditioned:   $(cond.ssa_posterior.param_names)
  unconditioned: $(uncond.param_names)
"""

ssa_names = cond.ssa_posterior.param_names
conditioned_mask = cond.boundary.conditioned
println("\n5. Boundary summary")
println("   SSA parameters: $(ssa_names)")
println("   Conditioned:    $(conditioned_mask)  " *
        "($(sum(conditioned_mask)) of $(length(ssa_names)))")
flush(stdout)

# --- Effect size, printed so the log alone answers "did it work?" ---
function weighted_moments(particles, weights, j)
    m = sum(weights .* particles[j, :])
    v = sum(weights .* (particles[j, :] .- m) .^ 2)
    return m, sqrt(v)
end

truth_lookup = Dict(pairs(true_ssa))
println("\n6. Per-parameter effect (shared parameters should tighten):")
println(rpad("param", 16), rpad("shared", 8), rpad("truth", 10),
        rpad("uncond mean±sd", 22), rpad("cond mean±sd", 22), "sd ratio")
for (j, name) in enumerate(ssa_names)
    um, us = weighted_moments(uncond.particles, uncond.weights, j)
    cm, cs = weighted_moments(cond.ssa_posterior.particles,
                              cond.ssa_posterior.weights, j)
    ratio = cs > 0 ? us / cs : Inf
    println(rpad(string(name), 16),
            rpad(conditioned_mask[j] ? "yes" : "no", 8),
            rpad(string(round(truth_lookup[name]; digits=3)), 10),
            rpad("$(round(um; digits=3)) ± $(round(us; digits=3))", 22),
            rpad("$(round(cm; digits=3)) ± $(round(cs; digits=3))", 22),
            round(ratio; digits=2))
end
flush(stdout)

# --- Export ---
# Written by hand rather than via DelimitedFiles/CSV: neither is a dependency of
# this project, and the files are small.
function write_matrix(path, rows::Vector{Symbol}, m::AbstractMatrix)
    open(path, "w") do io
        println(io, join(rows, ","))
        for k in axes(m, 2)
            println(io, join((m[j, k] for j in axes(m, 1)), ","))
        end
    end
end

function write_vector(path, header::String, v::AbstractVector)
    open(path, "w") do io
        println(io, header)
        for x in v
            println(io, x)
        end
    end
end

println("\n7. Exporting to $OUTDIR ...")
write_matrix(joinpath(OUTDIR, "isab_cond_particles.csv"), ssa_names,
             cond.ssa_posterior.particles)
write_matrix(joinpath(OUTDIR, "isab_uncond_particles.csv"), ssa_names,
             uncond.particles)
write_vector(joinpath(OUTDIR, "isab_cond_weights.csv"), "weight",
             cond.ssa_posterior.weights)
write_vector(joinpath(OUTDIR, "isab_uncond_weights.csv"), "weight",
             uncond.weights)

open(joinpath(OUTDIR, "isab_ssa_truth.csv"), "w") do io
    println(io, "param,truth,conditioned")
    for (j, name) in enumerate(ssa_names)
        println(io, "$name,$(truth_lookup[name]),$(conditioned_mask[j])")
    end
end

# The ODE chain backs the backup identifiability slide.
ode_names = Symbol.(String.(names(cond.ode_chain, :parameters)))
open(joinpath(OUTDIR, "isab_ode_chain.csv"), "w") do io
    println(io, join(ode_names, ","))
    samples = [vec(cond.ode_chain[string(n)].data) for n in ode_names]
    for i in eachindex(samples[1])
        println(io, join((s[i] for s in samples), ","))
    end
end

open(joinpath(OUTDIR, "isab_ode_truth.csv"), "w") do io
    println(io, "param,truth")
    for (name, val) in pairs(true_ode)
        println(io, "$name,$val")
    end
end

println("   Wrote isab_{cond,uncond}_{particles,weights}.csv, " *
        "isab_ssa_truth.csv, isab_ode_chain.csv, isab_ode_truth.csv")
println("\n=== Complete at $(Dates.now()) ===")
