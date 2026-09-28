# Spec §11 task 14.8, check 7: the clipping census on the assembled Core A′,
# scored for K5 (spec §8). It is run and scored, not gated (§12, 2026-09-25 C).
#
# Two parts, one Slurm array:
# - published parameters over PUB_SEEDS, one full cycle each;
# - N_DRAWS prior draws, one full cycle each. A draw varies every ODE-block
#   kinetic constant with an informed (`:balanced`) prior, in glycolysis and in
#   recycling. Asserted and uninformed priors stay at their values, and no slot
#   the driver writes is drawn (§11 phase 14b, decided 2026-09-26).
#
# Draw 0 is the nominal values through the freed build. Its census must equal
# the published build's at the same seed: that is the evidence that freeing the
# constants changes nothing but which vector holds them.
#
# Every task writes one TSV of rows. `corea_census_merge.jl` reads them all and
# writes the result file. A draw whose solve fails is recorded as failed with
# its message, never dropped.
#
# Usage: julia --project dev/scripts/corea_census.jl <task> <ntasks>
# Regenerate: sbatch dev/scripts/corea_census.slurm, then
#             julia --project dev/scripts/corea_census_merge.jl

using InferCell
using Printf
using Random
using Dates

const TASK = parse(Int, get(ARGS, 1, get(ENV, "SLURM_ARRAY_TASK_ID", "0")))
const NTASK = parse(Int, get(ARGS, 2, get(ENV, "NTASK", "1")))
const SMOKE = get(ENV, "SMOKE", "0") == "1"
const CYCLE = SMOKE ? 120 : round(Int, COREA_CYCLE_S)
const N_DRAWS = SMOKE ? 2 : 200
const DRAW_SEED0 = 14_800             # draw i runs at seed DRAW_SEED0 + i
const PUB_SEEDS = SMOKE ? [DRAW_SEED0] :
    vcat([1310, 1410], DRAW_SEED0 .+ (0:7))   # 13b's and 14a's seeds, and draw 0's
const OUT = joinpath(@__DIR__, "census", "task_$(TASK).tsv")

# The modules the draws vary: the assembled composition with every informed
# kinetic constant of the two balanced modules freed at its value.
_informed(m) = Symbol[q.name for q in parameters(m)
                      if q.role === :rate && informedness(q) === :balanced]
function census_models()
    g0 = CentralGlycolysis(enzymes = :translated)
    r0 = NucleotideRecycling(enzymes = :translated)
    ms = corea_models()
    ms[1] = CentralGlycolysis(enzymes = :translated, free = _informed(g0))
    ms[3] = NucleotideRecycling(enzymes = :translated, free = _informed(r0))
    return ms
end

ode_param_names(ms) = unique(Symbol[q.name for m in ms if formalism(m) === :ode
                                    for q in model_free_params(parameters(m))])

"One full cycle with a clip record. Returns the record, the census and the wall time."
function census_cycle(ms; seed, draw = nothing)
    Random.seed!(seed)
    d = build_problem(ms; tspan = (0.0, Float64(CYCLE)), complete = true)
    if draw !== nothing
        names = ode_param_names(ms)
        for (s, v) in draw
            d.ode.p[findfirst(==(s), names)] = v
        end
    end
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
            get(ENV, "SLURM_JOB_ID", "none"), ", task $TASK of $NTASK, $(now()), ",
            "Julia $VERSION", SMOKE ? ", SMOKE" : "")
    println(io, join(("kind", "id", "seed", "status", "drains", "clipped", "fraction",
                      "wall_s", "clips(counter:drains:first_s:max_deficit)", "message"), '\t'))
    jobs = vcat([(:published, i, s) for (i, s) in enumerate(PUB_SEEDS)],
                [(:draw, i, DRAW_SEED0 + i) for i in 0:N_DRAWS])   # draw 0 is nominal
    ms_pub = corea_models()
    ms_cen = census_models()
    names = ode_param_names(ms_cen)
    drawn = [q for m in ms_cen if formalism(m) === :ode for q in model_free_params(parameters(m))
             if q.role === :rate && informedness(q) === :balanced]
    TASK == 0 && println("Drawing $(length(drawn)) informed kinetic constants; ",
                         "$(length(jobs)) cycles over $NTASK tasks")
    for (j, (kind, i, seed)) in enumerate(jobs)
        (j - 1) % NTASK == TASK || continue
        println("$(now()) $kind $i seed $seed"); flush(stdout)
        draw = nothing
        if kind === :draw
            rng = MersenneTwister(seed)
            vals = [rand(rng, q.prior) for q in drawn]
            draw = i == 0 ? [q.name => q.value for q in drawn] :
                            [q.name => v for (q, v) in zip(drawn, vals)]
        end
        try
            rec, census, wall = census_cycle(kind === :draw ? ms_cen : ms_pub; seed, draw)
            row(io, kind, i, seed, "ok", rec, census, wall)
        catch e
            row(io, kind, i, seed, "failed", nothing, nothing, 0.0, first(sprint(showerror, e), 300))
        end
    end
end
println("Done: $OUT")
