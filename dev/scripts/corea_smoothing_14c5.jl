# Spec §11 task 14c.5, closed as §12 (2026-09-29, the later entry) specifies:
# the sampler model against the published one over 5,000 seeds.
#
# Each seed runs a full cycle under the published model (clamped drain,
# fractional carry) and the sampler model (COREA_SAMPLER: drain smoothed at one
# particle, pools carried continuously). Two things are recorded per seed:
# - the coupled half, in the run: every paired difference at a 60 s save point
#   before the first handshake at which any jump state differs, relative to
#   max(|published|, floor). Only the largest and the count are written.
# - the decoupled half's input: every candidate observable at the end of the
#   cycle, under both models. The merge gates them with `ensemble_agreement`
#   and builds the null from the published runs alone.
#
# Candidate observables are 14c.5's 50: every ODE state as a particle count
# with its carried remainder, every transcript count, and the cell volume.
# The seeds are 14c's 100 (job 17625498) and 4,900 more, all rerun here at one
# commit.
#
# Usage: julia --project dev/scripts/corea_smoothing_14c5.jl <task> <ntasks>
#        Merge with corea_smoothing_14c5_merge.jl.
# Regenerate: sbatch dev/scripts/corea_smoothing_14c5.slurm

using InferCell
using Printf
using Random
using Dates

const TASK = parse(Int, get(ARGS, 1, get(ENV, "SLURM_ARRAY_TASK_ID", "0")))
const NTASK = parse(Int, get(ARGS, 2, "1"))
const SMOKE = get(ENV, "SMOKE", "0") == "1"
const CYCLE = SMOKE ? 120 : round(Int, COREA_CYCLE_S)
const SEEDS = SMOKE ? [14_800, 15_001] :
              vcat([1310, 1410], 14_800 .+ (0:7), 15_000 .+ (1:4990))
const EVERY = 60
const THRESHOLD = 0.01
const ODE_FLOOR = 500.0
const TX_FLOOR = 1.0
const DIR = joinpath(@__DIR__, "smoothing_14c5")
const HEADER = "# commit $(strip(read(`git rev-parse --short HEAD`, String))), job " *
               "$(get(ENV, "SLURM_ARRAY_JOB_ID", get(ENV, "SLURM_JOB_ID", "none"))), " *
               "task $TASK of $NTASK, $(now()), Julia $VERSION$(SMOKE ? ", SMOKE" : "")"
mkpath(DIR)
@assert length(unique(SEEDS)) == length(SEEDS) == (SMOKE ? 2 : 5000)

_names(ms, f) = InferCell._block_names(ms, f)
const TX = Set(transcript_state(g.locus) for g in read_transcription_genes())

floor_of(o) = o === :volume_litres ? 0.0 :
              startswith(String(o), "mRNA_") ? TX_FLOOR : ODE_FLOOR

function observe(d, ode_names, jump_names)
    vals = Pair{Symbol, Float64}[]
    for (i, s) in enumerate(ode_names)
        push!(vals, s => d.ode.u[i] * d.factor + d.rounding.remainders[i])
    end
    for (j, s) in enumerate(jump_names)
        s in TX && push!(vals, s => float(d.jump.u[j]))
    end
    push!(vals, :volume_litres => d.volume_litres)
    return vals
end

# One full cycle. The published run (`reference === nothing`) returns its jump
# state at every handshake and its observables at every save point. The sampler
# run compares against them: a mismatch at handshake k means the save at k is
# already decoupled, as in 14c's merge.
function run_cycle(seed, variant; reference = nothing)
    Random.seed!(seed)
    kw = variant === :sampler ? COREA_SAMPLER : (;)
    ms = corea_models(; smoothing = get(kw, :smoothing, nothing))
    d = build_corea(; tspan = (0.0, Float64(CYCLE)), kw...)
    ode_names, jump_names = _names(ms, :ode), _names(ms, :jump)
    Random.seed!(seed)
    jumps = Vector{Vector{Int}}(undef, CYCLE)
    saves = Dict{Int, Vector{Pair{Symbol, Float64}}}(0 => observe(d, ode_names, jump_names))
    t_dec = -1
    worst = (r = 0.0, t = 0, o = :none)
    ncomp = 0
    compare(t) = for ((o, v), (o2, p)) in zip(saves[t], reference.saves[t])
        @assert o === o2
        r = (v - p) / max(abs(p), floor_of(o))
        ncomp += 1
        abs(r) > abs(worst.r) && (worst = (r = r, t = t, o = o))
    end
    reference === nothing || compare(0)
    wall = @elapsed for k in 1:CYCLE
        handshake_step!(d)
        jumps[k] = collect(Int, d.jump.u)
        if reference !== nothing && t_dec < 0 && jumps[k] != reference.jumps[k]
            t_dec = k
        end
        if k % EVERY == 0 || k == CYCLE
            saves[k] = observe(d, ode_names, jump_names)
            reference !== nothing && t_dec < 0 && compare(k)
        end
    end
    c = clipping_census(d)
    return (jumps = jumps, saves = saves, t_dec = t_dec, worst = worst, ncomp = ncomp,
            wall = wall, drains = c.drains, clipped = c.clipped)
end

open(joinpath(DIR, "task_$(TASK).tsv"), "w") do io
    println(io, HEADER)
    println(io, join(("seed", "variant", "observable", "value"), '\t'))
    for (j, seed) in enumerate(SEEDS)
        (j - 1) % NTASK == TASK || continue
        println("$(now()) seed $seed"); flush(stdout)
        pub = run_cycle(seed, :published)
        smp = run_cycle(seed, :sampler; reference = pub)
        for (v, r) in ((:published, pub), (:sampler, smp))
            for (o, x) in r.saves[CYCLE]
                println(io, join((seed, v, o, repr(x)), '\t'))
            end
            println(io, "# census\t$seed\t$v\t$(r.drains)\t$(r.clipped)\t",
                    @sprintf("%.1f", r.wall))
        end
        println(io, "# couple\t$seed\t$(smp.t_dec)\t$(smp.ncomp)\t$(repr(smp.worst.r))\t",
                "$(smp.worst.t)\t$(smp.worst.o)")
        flush(io)
    end
end
println("Done")
