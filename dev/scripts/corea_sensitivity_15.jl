# Spec §11 tasks 15.6 and 15.8 (§12 2026-09-29, phase 15 planning): the
# ensemble sensitivity of every candidate observable to D11's six targets and
# the polymerase direction, on the sampler model.
#
# Each seed runs 15 configurations, nominal and ln θ ± 0.2 on each of seven
# columns, under common random numbers. Each run is a full cycle, recorded
# every 60 s by `emit_observables!`. The candidate observables are every
# non-protein ODE state (as particles), every transcript and the volume;
# protein counts are ruled out of the likelihood (§4 D8). The merge builds the
# Jacobian, its split-half noise floor and the influence audit.
#
# Usage: julia --project dev/scripts/corea_sensitivity_15.jl <task> <ntasks> <nseeds>
# Regenerate: sbatch dev/scripts/corea_sensitivity_15.slurm, then
#             julia --project dev/scripts/corea_sensitivity_15_merge.jl

using InferCell
using Random
using Serialization
using Dates

const TASK = parse(Int, get(ARGS, 1, get(ENV, "SLURM_ARRAY_TASK_ID", "0")))
const NTASK = parse(Int, get(ARGS, 2, "1"))
const NSEEDS = parse(Int, get(ARGS, 3, "100"))
const SMOKE = get(ENV, "SMOKE", "0") == "1"
const CYCLE = SMOKE ? 120.0 : COREA_CYCLE_S
const DELTA = 0.2
const SEEDS = 20_000 .+ (1:NSEEDS)
const DIR = joinpath(@__DIR__, "sensitivity_15")
const COMMIT = strip(read(`git rev-parse --short HEAD`, String))
mkpath(DIR)

const COLUMNS = vcat(D11_TARGETS, [POLYMERASE_DIRECTION])
const CONFIGS = vcat([(:nominal, 0)], [(c, s) for c in COLUMNS for s in (+1, -1)])

function run_config(seed, column, sign)
    Random.seed!(seed)
    ms = d11_models(; smoothing = COREA_SAMPLER.smoothing)
    d = build_problem(ms; tspan = (0.0, CYCLE), complete = true,
                      rounding = COREA_SAMPLER.rounding)
    column === :nominal || set_parameters!(d, ms, perturbed_values(ms, column, exp(sign * DELTA)))
    Random.seed!(seed)
    r = emit_observables!(d, ms; every = 60)
    proteins = Set(protein_count_states(ms))
    keep = [i for (i, s) in enumerate(r.species) if !(s in proteins)]
    return r.species[keep], r.times, vcat(r.counts[keep, :], permutedims(r.volume))
end

out = Dict{Tuple{Int, Symbol, Int}, Matrix{Float64}}()
species, times = Symbol[], Float64[]
for (j, seed) in enumerate(SEEDS)
    (j - 1) % NTASK == TASK || continue
    for (column, sign) in CONFIGS
        println("$(now()) seed $seed $column $sign"); flush(stdout)
        sp, ts, x = run_config(seed, column, sign)
        global species, times = vcat(sp, [:volume_litres]), ts
        out[(seed, column, sign)] = x
    end
end
serialize(joinpath(DIR, "task_$(TASK).jls"),
          (commit = COMMIT, job = get(ENV, "SLURM_ARRAY_JOB_ID", get(ENV, "SLURM_JOB_ID", "none")),
           smoke = SMOKE, delta = DELTA, species = species, times = times,
           columns = COLUMNS, runs = out))
println("Done")
