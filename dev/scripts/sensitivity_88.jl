# Issue #88, tasks A1 and A2: the ensemble sensitivity that the linearised
# information bound and R are built from, on M0 (A1) or on Core A′ (A2, a
# rerun of 15.8 that keeps its raw outputs).
#
# Each seed runs nominal and ln θ ± 0.2 on each target, under common random
# numbers, as `corea_sensitivity_15.jl` does. Every run keeps the full species ×
# times matrix (non-protein states as particles, plus volume), so the nominal
# runs are also the per-cell ensemble R is computed from.
#
# - `m0`: `build_m0()` as the V6/V7 sampler and `m0_dataset` build it, M0's
#   five targets, 600 s.
# - `corea`: 15.8's sampler model (`d11_models` with `COREA_SAMPLER`), D11's six
#   (15.8's polymerase column is dropped), one 6300 s cycle.
#
# Task files are written to OUTDIR, which must lie outside the worktrees, so a
# worktree cleanup cannot delete them as it deleted 15.8's.
#
# Usage: julia --project dev/scripts/sensitivity_88.jl <m0|corea> <task> <ntasks> <nseeds> <outdir>
# Regenerate: dev/scripts/sensitivity_88.sh <m0|corea> <nseeds> <ntasks> <outdir>

using InferCell
using Random
using Serialization
using Dates

const MODEL = Symbol(ARGS[1])
const TASK = parse(Int, ARGS[2])
const NTASK = parse(Int, ARGS[3])
const NSEEDS = parse(Int, ARGS[4])
const DIR = ARGS[5]
const SMOKE = get(ENV, "SMOKE", "0") == "1"
const DELTA = 0.2
const SEEDS = 20_000 .+ (1:NSEEDS)
const COMMIT = strip(read(`git rev-parse --short HEAD`, String))

MODEL in (:m0, :corea) || error("model must be m0 or corea, not $MODEL")
any(w -> occursin(w, abspath(DIR)), ("/.worktrees/", "/.claude/worktrees/")) &&
    error("OUTDIR $DIR is inside a worktree")

const HORIZON = SMOKE ? 120.0 : MODEL === :m0 ? M0_HORIZON_S : COREA_CYCLE_S
const COLUMNS = MODEL === :m0 ? M0_TARGETS : D11_TARGETS
const CONFIGS = vcat([(:nominal, 0)], [(c, s) for c in COLUMNS for s in (+1, -1)])

models() = MODEL === :m0 ? m0_models() : d11_models(; smoothing = COREA_SAMPLER.smoothing)
# `build_m0` is this call on `m0_models()` with the default rounding.
const BUILD_KW = MODEL === :m0 ? (;) : (; rounding = COREA_SAMPLER.rounding)
build(ms) = build_problem(ms; tspan = (0.0, HORIZON), complete = true, BUILD_KW...)

function run_config(seed, column, sign)
    Random.seed!(seed)
    ms = models()
    d = build(ms)
    column === :nominal || set_parameters!(d, ms, perturbed_values(ms, column, exp(sign * DELTA)))
    Random.seed!(seed)
    r = emit_observables!(d, ms; every = 60)
    proteins = Set(protein_count_states(ms))
    keep = [i for (i, s) in enumerate(r.species) if !(s in proteins)]
    return vcat(r.species[keep], [:volume_litres]), r.times,
           vcat(r.counts[keep, :], permutedims(r.volume))
end

function main()
    mkpath(DIR)
    out = Dict{Tuple{Int, Symbol, Int}, Matrix{Float64}}()
    species, times = Symbol[], Float64[]
    for (j, seed) in enumerate(SEEDS)
        (j - 1) % NTASK == TASK || continue
        for (column, sign) in CONFIGS
            println("$(now()) seed $seed $column $sign"); flush(stdout)
            species, times, out[(seed, column, sign)] = run_config(seed, column, sign)
        end
    end
    serialize(joinpath(DIR, "task_$(TASK).jls"),
              (model = MODEL, commit = COMMIT,
               job = get(ENV, "SLURM_ARRAY_JOB_ID", get(ENV, "SLURM_JOB_ID", "none")),
               smoke = SMOKE, delta = DELTA, horizon = HORIZON, species, times,
               columns = COLUMNS, runs = out))
    println("Done")
end

main()
