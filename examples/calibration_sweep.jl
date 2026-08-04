# Multi-seed calibration sweep for the iterative boundary protocol (issue #11).
#
# For each seed in this shard's slice, generate synthetic ODE+SSA data, run both
# the single-pass (`sequential_infer`) and iterative (`iterative_infer`) Level-2
# protocols, then extract CI bounds at nominal levels [0.80, 0.90, 0.95, 0.99]
# for the 4 ODE parameters. Writes one CSV per shard; the plot script
# concatenates all shards.
#
# Shard selection: reads SLURM_ARRAY_TASK_ID from the environment (defaults to 0
# for local runs). Pass `--smoke` to run a fast 2-seed smoke test with tiny
# sample counts.

using InferCell
using OrdinaryDiffEq: Tsit5, solve
using MCMCChains: Chains
using Random
using Statistics: quantile
using Dates
using Printf

# --- Shard / smoke configuration ---

const SMOKE = "--smoke" in ARGS
const SHARD_ID = parse(Int, get(ENV, "SLURM_ARRAY_TASK_ID", "0"))
const SEEDS_PER_SHARD = 5
const ALL_SEEDS = 1001:1020

const seeds_for_shard = if SMOKE
    1001:1002
else
    lo = first(ALL_SEEDS) + SHARD_ID * SEEDS_PER_SHARD
    hi = lo + SEEDS_PER_SHARD - 1
    lo:hi
end

const N_SAMPLES = SMOKE ? 50 : 300
const SSA_KW = if SMOKE
    (; n_particles=50, n_populations=3, n_replicates=20)
else
    (; n_particles=200, n_populations=6, n_replicates=50)
end

const NOMINAL_LEVELS = [0.80, 0.90, 0.95, 0.99]
const ODE_PARAMS = [:k_tx, :k_tl, :gamma_mRNA, :gamma_protein]
const ODE_TRUTH = Dict(:k_tx => 1.0, :k_tl => 2.0,
                      :gamma_mRNA => 0.5, :gamma_protein => 0.1)

const RESULTS_DIR = joinpath(@__DIR__, "results")
mkpath(RESULTS_DIR)
const OUT_PATH = joinpath(RESULTS_DIR,
    SMOKE ? "step4_calibration_shard_smoke.csv" :
            "step4_calibration_shard_$(SHARD_ID).csv")

println("=== Calibration sweep ===")
println("Shard ID:       $SHARD_ID  (smoke=$SMOKE)")
println("Seeds:          $(collect(seeds_for_shard))")
println("ODE n_samples:  $N_SAMPLES")
println("SSA kwargs:     $SSA_KW")
println("Output:         $OUT_PATH")
println("Start:          $(Dates.now())")
flush(stdout)

# --- Helpers ---

function ci_bounds(samples::AbstractVector{<:Real}, level::Float64)
    α = (1 - level) / 2
    return quantile(samples, [α, 1 - α])
end

# Build the same closed-loop biology setup as test/run_integ_iterative.jl.
# Truth values are fixed; only the noise realisation varies by seed.
function build_data(seed::Int)
    Random.seed!(seed)

    true_ode = (k_tx=1.0, k_tl=2.0, gamma_mRNA=0.5, gamma_protein=0.1)
    true_ssa = (k_on=0.2, k_off=0.5, k_tx_burst=10.0, k_tl=2.0,
                gamma_mRNA=0.5, gamma_protein=0.1)

    txl = TranscriptionTranslation(; true_ode..., sigma_obs=0.3)
    metab = LightMetabolism()
    ode_prob = build_problem([txl, metab]; tspan=(0.0, 50.0))
    ode_sol = solve(ode_prob, Tsit5())
    ode_data = observe(ode_sol, 0.0:2.5:50.0, [txl, metab]; sigma=0.3)

    bge = BurstyGeneExpression(; true_ssa...)
    ssa_tspan = (0.0, 100.0)
    ssa_times = collect(0.0:10.0:100.0)
    ssa_prob = build_problem(bge; tspan=ssa_tspan)
    trajectories = [solve(ssa_prob, SSAStepper(); saveat=ssa_times) for _ in 1:200]
    ssa_data = observe(trajectories, ssa_times, bge)

    return (; txl, metab, bge, ode_data, ssa_data, ssa_tspan)
end

function write_rows!(io, seed, protocol, ode_chain, n_iters, kl_trace, runtime_sec)
    kl_str = isempty(kl_trace) ? "" : join(kl_trace, ";")
    n_iters_str = ismissing(n_iters) ? "" : string(n_iters)
    for param in ODE_PARAMS
        samples = vec(ode_chain[string(param)].data)
        truth = ODE_TRUTH[param]
        for level in NOMINAL_LEVELS
            lo, hi = ci_bounds(samples, level)
            covered = lo <= truth <= hi
            @printf(io, "%d,%s,%s,%.4f,%.6f,%.6f,%.6f,%d,%s,%.2f,\"%s\"\n",
                    seed, string(protocol), string(param), level,
                    lo, hi, truth, Int(covered),
                    n_iters_str, runtime_sec, kl_str)
        end
    end
    flush(io)
end

# --- Sweep ---
# Write incrementally so a timeout still preserves completed seeds.

open(OUT_PATH, "w") do io
    println(io, "seed,protocol,param,level,lo,hi,truth,covered,n_iters,runtime_sec,kl_trace")
    flush(io)

    for seed in seeds_for_shard
        println("\n--- Seed $seed ---")
        flush(stdout)

        setup = build_data(seed)
        ode_kwargs = (; n_samples=N_SAMPLES)
        ssa_kwargs = merge(SSA_KW, (; tspan=setup.ssa_tspan))

        # Single-pass
        Random.seed!(seed * 31 + 1)
        t0 = time()
        single = sequential_infer([setup.txl, setup.metab], setup.ode_data,
                                  [setup.bge], setup.ssa_data;
                                  ode_kwargs=ode_kwargs, ssa_kwargs=ssa_kwargs,
                                  boundary_method=:kde, verbose=false)
        single_runtime = time() - t0
        println("  single-pass: $(round(single_runtime; digits=1))s")
        flush(stdout)
        write_rows!(io, seed, :single, single.ode_chain,
                    missing, Float64[], single_runtime)

        # Iterative
        Random.seed!(seed * 31 + 2)
        t1 = time()
        iter = iterative_infer([setup.txl, setup.metab], setup.ode_data,
                               [setup.bge], setup.ssa_data;
                               ode_kwargs=ode_kwargs, ssa_kwargs=ssa_kwargs,
                               boundary_method=:kde,
                               max_iters=4, kl_tol=0.05,
                               keep_history=false, verbose=false,
                               rng=Xoshiro(seed * 31 + 2))
        iter_runtime = time() - t1
        println("  iterative:   $(round(iter_runtime; digits=1))s, n_iters=$(iter.n_iters), kl_trace=$(iter.kl_trace)")
        flush(stdout)
        write_rows!(io, seed, :iter, iter.ode_chain,
                    iter.n_iters, iter.kl_trace, iter_runtime)
    end
end

println("\nWrote results to $OUT_PATH")
println("End: $(Dates.now())")
flush(stdout)
