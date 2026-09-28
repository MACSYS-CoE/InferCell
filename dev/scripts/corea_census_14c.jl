# Spec §11 tasks 14c.2, 14c.3 and 14c.6: check 7's census, rerun three ways.
# All three run the published clamped drain; they ask where it clips, not what
# smoothing does (that is 14c.5).
#
# MODE (first argument):
# - d11:   200 draws of D11's six targets from their priors under 14c.1's
#          Haldane rule, every other parameter at its published value. Scored
#          on T3 as K5's prior half (14c.2).
# - broad: 14b's broad census — every informed kinetic constant of glycolysis
#          and recycling — drawn under 14c.1 instead of independently. Recorded,
#          not gated (14c.3).
# - gtp:   14b's published-parameter seeds with the GTP branch's capacity, PGK3
#          and PYK3, scaled up. Scaling both catalytic constants of a reaction
#          by s is scaling its capacity E·kcat by s at every handshake, whatever
#          count the catalytic channel writes into E, and it keeps the reaction's
#          equilibrium constant. Recorded, not gated (14c.6).
#
# Draw i runs at seed DRAW_SEED0 + i, 14b's seeds, and draw 0 is the nominal
# values through the freed build: its census must equal the published build's
# at seed 14800, which is the evidence that `set_parameters!` changes nothing
# but which vector holds a value.
#
# Every task writes one TSV. `corea_census_14c_merge.jl` reads all three modes
# and writes the result file. A failed cycle is recorded with its message.
#
# Usage: julia --project dev/scripts/corea_census_14c.jl <mode> <task> <ntasks>
# Regenerate: sbatch dev/scripts/corea_census_14c.slurm <mode>, for each mode,
#             then julia --project dev/scripts/corea_census_14c_merge.jl

using InferCell
using Printf
using Random
using Dates

const MODE = Symbol(ARGS[1])
MODE in (:d11, :broad, :gtp) || error("mode must be d11, broad or gtp, not $MODE")
const TASK = parse(Int, get(ARGS, 2, get(ENV, "SLURM_ARRAY_TASK_ID", "0")))
const NTASK = parse(Int, get(ARGS, 3, "1"))
const SMOKE = get(ENV, "SMOKE", "0") == "1"
const CYCLE = SMOKE ? 120 : round(Int, COREA_CYCLE_S)
const N_DRAWS = SMOKE ? 2 : 200
const DRAW_SEED0 = 14_800
const PUB_SEEDS = SMOKE ? [DRAW_SEED0] : vcat([1310, 1410], DRAW_SEED0 .+ (0:7))
const GTP_SCALES = SMOKE ? [1, 10] : [1, 2, 5, 10, 100]
const GTP_CONSTANTS = [:kcatF_R_PGK3, :kcatR_R_PGK3, :kcatF_R_PYK3, :kcatR_R_PYK3]
const OUT = joinpath(@__DIR__, "census_14c", String(MODE), "task_$(TASK).tsv")

_informed(m) = Symbol[q.name for q in parameters(m)
                      if q.role === :rate && informedness(q) === :balanced]

# The composition each mode runs, and the names it draws or sets.
function mode_models()
    if MODE === :d11
        return d11_models(), D11_TARGETS
    elseif MODE === :broad
        g0 = CentralGlycolysis(enzymes = :translated)
        r0 = NucleotideRecycling(enzymes = :translated)
        drawn = vcat(_informed(g0), _informed(r0))
        # Every reverse constant is freed, informed or not, since a draw derives it.
        freed(m) = unique(vcat(_informed(m), [r.reverse for r in haldane_relations(m)]))
        ms = corea_models()
        ms[1] = CentralGlycolysis(enzymes = :translated, free = freed(g0))
        ms[3] = NucleotideRecycling(enzymes = :translated, free = freed(r0))
        return ms, drawn
    else
        ms = corea_models()
        ms[3] = NucleotideRecycling(enzymes = :translated, free = GTP_CONSTANTS)
        return ms, GTP_CONSTANTS
    end
end

"One full cycle with a clip record. Returns the record, the census and the wall time."
function census_cycle(ms, values; seed)
    Random.seed!(seed)
    d = build_problem(ms; tspan = (0.0, Float64(CYCLE)), complete = true)
    set_parameters!(d, ms, values)
    Random.seed!(seed)
    rec = ClipRecord(d)
    wall = @elapsed for _ in 1:CYCLE
        handshake_step!(d)
        record_clips!(rec, d)
    end
    return rec, clipping_census(d), wall
end

function row(io, kind, id, seed, status, rec, census, wall, msg = "")
    clips = rec === nothing ? "" :
        join((@sprintf("%s:%d:%.1f:%.3g", r.counter, r.drains, r.first, r.max_deficit)
              for r in clip_summary(rec)), ";")
    println(io, join((kind, id, seed, status,
                      census === nothing ? "" : census.drains,
                      census === nothing ? "" : census.clipped,
                      census === nothing ? "" : @sprintf("%.6g", census.fraction),
                      @sprintf("%.1f", wall), clips, replace(msg, r"\s+" => " ")), '\t'))
    flush(io)
end

mkpath(dirname(OUT))
open(OUT, "w") do io
    println(io, "# commit $(strip(read(`git rev-parse --short HEAD`, String))), job ",
            get(ENV, "SLURM_JOB_ID", "none"), ", mode $MODE, task $TASK of $NTASK, ",
            "$(now()), Julia $VERSION", SMOKE ? ", SMOKE" : "")
    println(io, join(("kind", "id", "seed", "status", "drains", "clipped", "fraction",
                      "wall_s", "clips(counter:drains:first_s:max_deficit)", "message"), '\t'))
    ms, names = mode_models()
    nominal = nominal_parameter_values(ms)
    jobs = MODE === :gtp ?
        [(:scale, s, seed) for s in GTP_SCALES for seed in PUB_SEEDS] :
        vcat([(:published, 1, DRAW_SEED0)],
             [(:draw, i, DRAW_SEED0 + i) for i in 0:N_DRAWS])
    TASK == 0 && println("Mode $MODE: $(length(names)) names; $(length(jobs)) cycles over $NTASK tasks")
    for (j, (kind, i, seed)) in enumerate(jobs)
        (j - 1) % NTASK == TASK || continue
        println("$(now()) $kind $i seed $seed"); flush(stdout)
        try
            values = if kind === :scale
                [n => i * nominal[n] for n in names]
            elseif kind === :draw && i > 0
                v = draw_parameters(MersenneTwister(seed), ms, names)
                assert_haldane(ms, v)
                v
            else
                Pair{Symbol, Float64}[]          # published, or draw 0
            end
            rec, census, wall = census_cycle(kind === :published ? corea_models() : ms,
                                             values; seed)
            row(io, kind, i, seed, "ok", rec, census, wall)
        catch e
            row(io, kind, i, seed, "failed", nothing, nothing, 0.0,
                first(sprint(showerror, e), 300))
        end
    end
end
println("Done: $OUT")
