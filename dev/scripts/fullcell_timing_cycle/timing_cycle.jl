# KIM-183: wall time of one full syn3A cell cycle in SpeedyCell, the syn3a-ensemble copy
# (/fred/oz022/tkimpson/syn3a-ensemble/code/SpeedyCell, developed by path in this
# directory's environment; nothing is written to that tree).
#
#   MODE=serial    one thread: a warm-up cell, then PROCIDS one after another
#   MODE=ensemble  one process, Threads.nthreads() threads: a warm-up cell, then PROCIDS
#                  under Threads.@threads :dynamic, timed as a whole
#
# Environment: MODE, PROCIDS ("a-b"), SEED (default 1), MINUTES (recorded minutes,
# default 121: blocks 0 … 120, so 120 simulated minutes, syn3a-ensemble's T_end), OUT
# (output prefix; writes OUT_cells.csv, OUT_minutes.csv, OUT_meta.txt).
#
# Every cell runs with ode_failure=:record, so a dead cell's failed ODE step ends that
# cell and is recorded instead of throwing; any other error is caught and recorded too.
# No cell is replaced or rerun. The cycle readouts per cell, from the per-minute record
# (code/cell_cycle/events.py's definitions):
#   t_SA2   first minute with CellSA >= 2 × its minute-0 value (the cycle end)
#   t_Vcap  first minute with CellV >= 670 (the volume ceiling, 6.70e-17 L)
#   t_short first minute with fewer than 59 hooks (collapse)
# and the wall time to minute k is the run's setup (total less the minutes' sum) plus
# the first k minutes' wall times.

using SpeedyCell.Cell
using InteractiveUtils, LinearAlgebra, Printf, SHA, Statistics, Dates

const SYN3A = "/fred/oz022/tkimpson/syn3a-ensemble"
const CELLV_CAP = 670
const HEALTHY_HOOKS = 59

function int_range(name, default)
    f = parse.(Int, split(get(ENV, name, default), '-'))
    length(f) <= 2 && first(f) <= last(f) || error("$name must be a or a-b with a <= b")
    first(f):last(f)
end

first_minute(f, n) = (k = findfirst(f, 0:(n - 1)); k === nothing ? -1 : k - 1)

# Wall time from run_cell's entry to the end of simulated minute k (k <= minutes run).
function wall_to(r, k)
    setup = r.timings.total - sum(r.minute_wall_s; init=0.0)
    k < 0 || k > length(r.minute_wall_s) ? NaN : setup + sum(view(r.minute_wall_s, 1:k); init=0.0)
end

function provenance(io)
    println(io, "date ", now())
    println(io, "host ", gethostname())
    cpu = try
        only(filter(l -> startswith(l, "Model name"), readlines(`lscpu`)))
    catch
        "Model name: " * Sys.cpu_info()[1].model
    end
    println(io, "lscpu ", cpu)
    println(io, "slurm_job ", get(ENV, "SLURM_JOB_ID", "none"), " cpus_per_task ", get(ENV, "SLURM_CPUS_PER_TASK", "?"))
    println(io, "threads ", Threads.nthreads(), " gcthreads ", Threads.ngcthreads(),
            " sweep_threads ", Base.JLOptions().nsweepthreads, " blas_threads ", BLAS.get_num_threads())
    println(io, "infercell_git ", readchomp(`git -C $(@__DIR__) rev-parse HEAD`))
    println(io, "syn3a_ensemble_git ", readchomp(`git -C $SYN3A rev-parse HEAD`))
    st = readchomp(`git -C $SYN3A status --porcelain -- code/SpeedyCell data/model`)
    println(io, "syn3a_ensemble_code_data_clean ", isempty(st))
    println(io, "speedycell_path ", pathof(SpeedyCell))
    for (root, _, files) in walkdir(joinpath(SYN3A, "data", "model")), f in sort(files)
        p = joinpath(root, f)
        println(io, "sha256 ", bytes2hex(open(sha256, p)), "  ", relpath(p, SYN3A))
    end
    versioninfo(io)
end

# One cell, never throwing: (run or nothing, outer wall seconds, error string).
function one_cell(m, p, seed, minutes)
    t0 = time()
    try
        r = run_cell(m, p; seed=seed, minutes=minutes, ode_failure=:record)
        return r, time() - t0, ""
    catch e
        return nothing, time() - t0, first(split(sprint(showerror, e), '\n'))
    end
end

function write_outputs(out, m, jobs, results, minutes)
    species = m.cme.species
    sa, v = findfirst(==("CellSA"), species), findfirst(==("CellV"), species)
    open(out * "_cells.csv", "w") do io
        println(io, "procid,seed,status,error,wall_outer_s,wall_total_s,setup_s,minutes_run,nevents,failure,",
                "t_SA2,t_Vcap,t_short,wall_to_105_s,wall_to_SA2_s,wall_to_120_s,cme_s,ode_s,ode_setup_s,",
                "hook_other_s,rebuild_s,io_s,thread,t_start,t_end")
        for ((p, s), (r, w, err, th, ts, te)) in zip(jobs, results)
            if r === nothing
                println(io, join((p, s, "error", replace(err, ',' => ';'), @sprintf("%.3f", w)), ","),
                        repeat(",", 17), ",", th, ",", ts, ",", te)
                continue
            end
            n = size(r.record, 2)
            x0 = r.record[sa, 1]
            tsa = first_minute(k -> r.record[sa, k + 1] >= 2 * x0, n)
            tv = first_minute(k -> r.record[v, k + 1] >= CELLV_CAP, n)
            tsh = (k = findfirst(<(HEALTHY_HOOKS), r.hooks_per_minute); k === nothing ? -1 : k)
            status = r.failure !== nothing ? "ode_failure" : "ok"
            tm = r.timings
            f = x -> isnan(x) ? "" : @sprintf("%.3f", x)
            println(io, join((p, s, status, "", f(w), f(tm.total), f(wall_to(r, 0)), n - 1, Int(r.nevents),
                              r.failure === nothing ? "" : "$(r.failure[1]):$(r.failure[2])",
                              tsa, tv, tsh, f(wall_to(r, 105)), f(tsa < 0 ? NaN : wall_to(r, tsa)),
                              f(wall_to(r, 120)), f(tm.cme), f(tm.ode), f(tm.ode_setup), f(tm.hook_other),
                              f(tm.rebuild), f(tm.io), th, @sprintf("%.3f", ts), @sprintf("%.3f", te)), ","))
        end
    end
    open(out * "_minutes.csv", "w") do io
        println(io, "procid,seed,minute,CellSA,CellV,hooks,wall_s,events")
        for ((p, s), (r, _...)) in zip(jobs, results)
            r === nothing && continue
            for k in 0:(size(r.record, 2) - 1)
                h, w, e = k == 0 ? (0, 0.0, 0) : (r.hooks_per_minute[k], r.minute_wall_s[k], Int(r.minute_events[k]))
                println(io, p, ',', s, ',', k, ',', r.record[sa, k + 1], ',', r.record[v, k + 1], ',', h, ',',
                        @sprintf("%.4f", w), ',', e)
            end
        end
    end
end

function main()
    mode = get(ENV, "MODE", "serial")
    mode in ("serial", "ensemble") || error("MODE must be serial or ensemble")
    procids = int_range("PROCIDS", "1-20")
    seed = parse(Int, get(ENV, "SEED", "1"))
    minutes = parse(Int, get(ENV, "MINUTES", "121"))
    minutes >= 3 || error("MINUTES must be at least 3")
    out = get(ENV, "OUT", "") == "" ? error("set OUT") : ENV["OUT"]
    mkpath(dirname(out))
    meta = open(out * "_meta.txt", "w")
    for io in (stdout, meta)
        println(io, "mode $mode procids $procids seed $seed minutes $minutes out $out")
        provenance(io)
    end
    flush(stdout)

    t_model = @elapsed m = CellModel()
    # Warm-up: the timed path (same kwargs), 3 recorded minutes so restart_clock! runs,
    # on a procid outside every timed set (1000). Compilation lands here, before the clock.
    t_warm = @elapsed one_cell(m, 1000, seed, 3)
    for io in (stdout, meta)
        @printf(io, "CellModel() %.2f s; warm-up (procid 1000, 3 minutes, includes compilation) %.2f s\n", t_model, t_warm)
    end
    flush(stdout)

    jobs = [(p, seed) for p in procids]
    results = Vector{Any}(undef, length(jobs))
    gc0 = Base.gc_time_ns()
    t0 = time()
    if mode == "serial"
        for (j, (p, s)) in enumerate(jobs)
            ts = time()
            r, w, err = one_cell(m, p, s, minutes)
            results[j] = (r, w, err, Threads.threadid(), ts - t0, time() - t0)
            @printf("procid %d: %.2f s, %s\n", p, w, r === nothing ? "error: $err" :
                    r.failure === nothing ? "ok" : "ode failure at minute $(r.failure[1])")
            flush(stdout)
        end
    else
        Threads.@threads :dynamic for j in eachindex(jobs)
            p, s = jobs[j]
            ts = time()
            r, w, err = one_cell(m, p, s, minutes)
            results[j] = (r, w, err, Threads.threadid(), ts - t0, time() - t0)
        end
    end
    wall = time() - t0
    gc_s = (Base.gc_time_ns() - gc0) / 1e9
    write_outputs(out, m, jobs, results, minutes)

    walls = [x[2] for x in results]
    nerr = count(x -> x[1] === nothing, results)
    nfail = count(x -> x[1] !== nothing && x[1].failure !== nothing, results)
    for io in (stdout, meta)
        @printf(io, "%s: %d cells (%d errors, %d ODE failures) in %.1f s wall; GC %.1f s (%.1f%% of wall)\n",
                mode, length(jobs), nerr, nfail, wall, gc_s, 100 * gc_s / wall)
        @printf(io, "per-cell outer wall: median %.2f s, min %.2f, max %.2f\n", median(walls), minimum(walls), maximum(walls))
        mode == "ensemble" && @printf(io, "cells per node-hour: %.0f (%d cells / %.1f s, %d threads)\n",
                                      length(jobs) * 3600 / wall, length(jobs), wall, Threads.nthreads())
    end
    close(meta)
end

main()
