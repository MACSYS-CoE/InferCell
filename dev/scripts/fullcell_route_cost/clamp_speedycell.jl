# KIM-195: SpeedyCell's own dNTP pools and replication cost counters, minute by minute,
# the reference for the `clamp` stage of route_cost_meanfield.jl (which form of the
# counter payment the mean field should use).
#
#   julia --project=dev/scripts/fullcell_timing_cycle \
#       dev/scripts/fullcell_route_cost/clamp_speedycell.jl <outdir>
#
# Runs in KIM-183's environment. Procid 1 (the mean field's initial state), seeds 1..8,
# over the fit window (54 recorded minutes); SMOKE=1: 2 seeds, 6 minutes. Writes
# <outdir>/clamp/speedycell.csv: per seed and minute, the 4 dNTP pools and their
# counters at the minute's end, and each pool's minimum over the minute's hooks.

using SpeedyCell, SpeedyCell.Cell

const SMOKE = get(ENV, "SMOKE", "0") == "1"
const DNTP = [("M_datp_c", "dATP_DNArep"), ("M_dttp_c", "dTTP_DNArep"), ("M_dctp_c", "dCTP_DNArep"),
              ("M_dgtp_c", "dGTP_DNArep")]

out = joinpath(abspath(ARGS[1]), "clamp"); mkpath(out)
m = CellModel()
sp = [split(l, ',')[2] for l in Iterators.drop(eachline(joinpath(SpeedyCell.Cell.MODEL_DIR, "cme", "species.csv")), 1)]
ix(s) = something(findfirst(==(s), sp))
pools = [ix(p) for (p, _) in DNTP]; ctrs = [ix(c) for (_, c) in DNTP]
seeds = SMOKE ? (1:2) : (1:8); minutes = SMOKE ? 6 : 54
t = @elapsed runs = [run_cell(m, 1; seed = s, minutes, track = [[i] for i in pools]) for s in seeds]
open(joinpath(out, "speedycell.csv"), "w") do io
    println(io, "seed,minute,", join(first.(DNTP), ','), ',', join(last.(DNTP), ','), ',',
            join(["min_hook_" * p for (p, _) in DNTP], ','))
    for (s, r) in zip(seeds, runs), k in 1:size(r.record, 2)
        mh = k > 1 ? r.track_min[:, k - 1] : fill(-1, length(pools))
        println(io, s, ',', k - 1, ',', join(Int.(r.record[pools, k]), ','), ',', join(Int.(r.record[ctrs, k]), ','),
                ',', join(mh, ','))
    end
end
open(joinpath(out, "speedycell_meta.txt"), "w") do io
    println(io, "host ", gethostname(), " job ", get(ENV, "SLURM_JOB_ID", "none"), " seeds ", seeds,
            " minutes ", minutes, " wall ", round(t; digits = 1), " s")
    println(io, "syn3a-ensemble ", readchomp(`git -C /fred/oz022/tkimpson/syn3a-ensemble rev-parse HEAD`))
end
println("wrote ", out)
