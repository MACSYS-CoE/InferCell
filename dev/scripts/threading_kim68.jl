# KIM-68: why `predict`'s threaded cells scale poorly (M5 prep). KIM-38 measured an
# efficiency of 0.48 at 16 threads and 0.16 at 64, on a 128-cell reduced-model
# matrix at 6,300 s. The PR #123 review read it as GC contention, not BLAS. This
# measures it on the same matrix: B7's held-out `v6` fit, `completion = :prior`,
# M 4 × K 8 × 4 cells, `rng = Xoshiro(18_438_301)`. These are KIM-38's seeds, from
# M4's block, reused so that the before and after compare like for like. Nothing
# here reads an evidence posterior. The instrumentation wraps
# `InferCell._hybrid_cell` here; `src` is not touched.
#
# Usage (dev/scripts/threading_kim68.sh submits every stage):
#   julia --project [--threads=N] [--gcthreads=…] dev/scripts/threading_kim68.jl <stage> <outdir> [<label>]
# Stages:
#   serial    8 of the matrix's cells, one at a time, one thread: per-cell
#             `@timed` (time, gctime, bytes)
#   threaded  2 cells serially first (this node's baseline), then all 128 cells
#             under each loop in LOOPS, joined by "+" (`chunked`: the frozen `_each`,
#             `Threads.@threads` over the jobs; `greedy`: the new `_each`). Per
#             cell: `@timed`, the thread at its start and stop, and its start and
#             stop times. Per pass: `Base.gc_num()` and `Base.gc_time_ns()`
#             differences, and MaxRSS. With PREDICT=1, also one full `predict` call
#             on the matrix, KIM-38's measurement. <label> names the variant (GC
#             flags or BLAS threads, set by the launcher).
#   procs     one of 4 Slurm tasks on one node (SLURM_PROCID), each running a
#             quarter of the cells with the new `_each` on its own threads, in its
#             own process and heap
#   allocs    `Profile.Allocs` on one 60-s cell: the top allocation sites by
#             bytes, by InferCell frame and by type
#   report    reads every stage's output and writes summary.md
#
# Efficiency = cells × this node's serial median / (threads × wall).

include(joinpath(@__DIR__, "pilot_common.jl"))
using InferCell: AbstractSubModel
using Random
using Statistics: median
using Printf
using LinearAlgebra: BLAS
using Profile

const STAGE = Symbol(ARGS[1])
const OUT = ARGS[2]
const LABEL = length(ARGS) >= 3 ? ARGS[3] : "default"
const H = SMOKE ? 120.0 : 6300.0
const M, K, NC = SMOKE ? (2, 2, 2) : (4, 8, 4)
const SEED = 18_438_301
const LOOPS = split(get(ENV, "LOOPS", "greedy"), "+")      # "+": sbatch --export splits on commas

# The matrix's parameter sets and seeds, from a 60-s `predict` with KIM-38's rng.
# Its draws, held values and seeds come before the horizon is read.
function matrix(; threads)
    pd = posterior_draws(V6)
    ms = metered(collect(AbstractSubModel, d11_models()))
    p0 = predict(ms, pd; completion = :prior, horizon = 60.0, n_held = K, n_parameter_draws = M,
                 n_cells = NC, rng = Xoshiro(SEED), threads)
    jobs = [(m, k, c) for m in 1:M for k in 1:K for c in 1:NC]     # `_run_matrix`'s order
    values = [haldane_values(ms, vcat(p0.theta_I[m], p0.theta_H[k])) for m in 1:M, k in 1:K]
    return (; pd, ms, p0, jobs, values,
            f = InferCell._forward_setup(ms, :hybrid; horizon = H))
end

function timed_cell(x, j, t0)
    m, k, c = x.jobs[j]
    a, ta = time_ns() - t0, Threads.threadid()
    st = @timed InferCell._hybrid_cell(x.f, x.values[m, k], x.p0.seeds[c, k, m], :lactate_out)
    b, tb = time_ns() - t0, Threads.threadid()
    return (; j, seed = x.p0.seeds[c, k, m], q = st.value.q, failure = st.value.failure,
            time = st.time, gctime = st.gctime, bytes = st.bytes,
            compile = get(st, :compile_time, NaN), thread_start = ta, thread_stop = tb,
            start = a / 1e9, stop = b / 1e9)
end

gcstats(d) = (; total_time = d.total_time / 1e9, pause = d.pause, full_sweep = d.full_sweep,
              allocd = d.allocd)

# One pass over jobs `js` with `loop`, each cell timed.
function pass(x, js, loop)
    out = Vector{Any}(undef, length(js))
    GC.gc()
    g0, n0, t0 = Base.gc_num(), Base.gc_time_ns(), time_ns()
    run1 = i -> (out[i] = timed_cell(x, js[i], t0))
    InferCell._restoring_global_rng() do
        if loop == "serial"
            foreach(run1, eachindex(js))
        elseif loop == "chunked"                       # the frozen `_each`, threaded
            Threads.@threads for i in eachindex(js)
                run1(i)
            end
        elseif loop == "greedy"                        # the new `_each`
            InferCell._each(run1, length(js), true)
        else
            error("unknown loop $loop")
        end
    end
    wall = (time_ns() - t0) / 1e9
    cells = collect(NamedTuple, out)
    return (; loop, wall, cells, gc = gcstats(Base.GC_Diff(Base.gc_num(), g0)),
            gc_time = (Base.gc_time_ns() - n0) / 1e9, maxrss_gib = Sys.maxrss() / 2^30)
end

env() = (; threads = Threads.nthreads(), gcthreads = Threads.ngcthreads(),
         blas_threads = BLAS.get_num_threads(), openblas = get(ENV, "OPENBLAS_NUM_THREADS", ""),
         node = gethostname(), label = LABEL, smoke = SMOKE, horizon = H)

# The 60-s `predict` in `matrix` has compiled the cell's path before any pass.
function serial()
    x = matrix(; threads = false)
    say("matrix of $(length(x.jobs)) cells; timing 8 at $H s, one at a time")
    r = pass(x, 1:min(8, length(x.jobs)), "serial")
    say("serial: median $(round(median(c.time for c in r.cells); digits = 2)) s per cell")
    save(OUT, "serial", (; env = env(), pass = r); stage = "serial")
end

function threaded()
    x = matrix(; threads = true)
    base = pass(x, 1:2, "serial")
    say("$(LABEL):$(Threads.nthreads()) threads, $(Threads.ngcthreads()) GC threads, ",
        "BLAS $(BLAS.get_num_threads()); baseline $(round(median(c.time for c in base.cells); digits = 2)) s")
    passes = Any[]
    for loop in LOOPS
        r = pass(x, eachindex(x.jobs), loop)
        say("  $loop: $(length(r.cells)) cells in $(round(r.wall; digits = 1)) s, GC ",
            "$(round(r.gc_time; digits = 1)) s, maxrss $(round(r.maxrss_gib; digits = 1)) GiB")
        push!(passes, r)
    end
    pr = nothing
    if get(ENV, "PREDICT", "0") == "1"
        GC.gc()
        t = @elapsed p = predict(x.ms, x.pd; completion = :prior, horizon = H, n_held = K,
                                 n_parameter_draws = M, n_cells = NC, rng = Xoshiro(SEED),
                                 threads = true)
        p.seeds == x.p0.seeds && p.theta_H == x.p0.theta_H && p.draws == x.p0.draws ||
            error("the 60-s matrix's draws differ from the full call's")
        qq = [p.cells[m, k][c].q for (m, k, c) in x.jobs]
        same = all(r -> isequal([c.q for c in sort(r.cells; by = c -> c.j)], qq), passes)
        pr = (; wall = t, cells = length(qq), statuses = plain(p).statuses, same_q_as_passes = same,
              maxrss_gib = Sys.maxrss() / 2^30)
        say("  predict: $(pr.cells) cells in $(round(t; digits = 1)) s; passes' Q equal: $same")
    end
    save(OUT, "threaded_$(Threads.nthreads())_$(LABEL)", (; env = env(), base, passes, predict = pr);
         stage = "threaded")
end

function procs()
    id, n = parse(Int, ENV["SLURM_PROCID"]), parse(Int, ENV["SLURM_NTASKS"])
    x = matrix(; threads = true)
    base = pass(x, 1:2, "serial")
    # The tasks start seconds apart; a file barrier makes their cells overlap.
    touch(joinpath(OUT, "procs_ready_$id"))
    while count(startswith("procs_ready_"), readdir(OUT)) < n
        sleep(0.5)
    end
    js = collect((id + 1):n:length(x.jobs))
    r = pass(x, js, "greedy")
    say("proc $id of $n: $(length(js)) cells in $(round(r.wall; digits = 1)) s on $(Threads.nthreads()) threads")
    save(OUT, "procs_$(n)x$(Threads.nthreads())_$id", (; env = env(), id, n, base, pass = r);
         stage = "procs")
end

function allocs()
    x = matrix(; threads = false)
    f60 = InferCell._forward_setup(x.ms, :hybrid; horizon = 60.0)
    cell() = InferCell._restoring_global_rng(() -> InferCell._hybrid_cell(
        f60, x.values[1, 1], x.p0.seeds[1, 1, 1], :lactate_out))
    cell()
    st = @timed cell()
    Profile.Allocs.clear()
    rate = 0.01
    Profile.Allocs.@profile sample_rate = rate cell()
    res = Profile.Allocs.fetch()
    # The innermost frame in InferCell's own source.
    ours(s) = (f = string(s.file); occursin("InferCell", f) && occursin("/src/", f))
    function site(a)
        i = findfirst(ours, a.stacktrace)
        i === nothing && return "(no InferCell frame)"
        s = a.stacktrace[i]
        return "$(s.func) at $(basename(string(s.file))):$(s.line)"
    end
    by_site, by_type = Dict{String, Vector{Int}}(), Dict{String, Vector{Int}}()
    for a in res.allocs
        for (d, key) in ((by_site, site(a)), (by_type, string(a.type)))
            v = get!(d, key, [0, 0])
            v[1] += a.size
            v[2] += 1
        end
    end
    top(d) = sort!([(; key = k, bytes = round(Int, v[1] / rate), count = round(Int, v[2] / rate))
                    for (k, v) in d]; by = r -> -r.bytes)[1:min(30, length(d))]
    r = (; env = env(), bytes = st.bytes, time = st.time, gctime = st.gctime, sample_rate = rate,
         sampled = length(res.allocs), sites = top(by_site), types = top(by_type))
    say("60-s cell: $(round(st.bytes / 2^20; digits = 1)) MiB in $(round(st.time; digits = 2)) s")
    for s in r.sites[1:min(10, end)]
        say("  ", lpad(round(s.bytes / 2^20; digits = 1), 8), " MiB  ", s.key)
    end
    save(OUT, "allocs", r; stage = "allocs")
end

function report()
    files = filter(f -> endswith(f, ".jls"), readdir(OUT))
    rows = String[]
    serial_med(base) = median(c.time for c in base.cells)
    # Shared by a pass and by the procs: per-cell slowdown, GC's share of the cells'
    # time and of the wall (the process's own GC clock), the pauses, and the idle
    # share, 1 - Σ cell time / (threads × wall).
    function cols(ps, T, s)
        cells = [c for p in ps for c in p.cells]
        gc, pauses, wall = sum(p.gc_time for p in ps), sum(p.gc.pause for p in ps), sum(p.wall for p in ps)
        return @sprintf("%.1f | %.2f | %.2f | %d | %.2f | %.2f", median(c.time for c in cells) / s,
                        sum(c.gctime for c in cells) / sum(c.time for c in cells), gc / wall,
                        pauses, gc / pauses, 1 - sum(c.time for c in cells) / (T * maximum(p.wall for p in ps)))
    end
    for f in sort(files)
        r = deserialize(joinpath(OUT, f))
        x = r.result
        if startswith(f, "threaded_")
            s = serial_med(x.base)
            for p in x.passes
                n = length(p.cells)
                eff = n * s / (x.env.threads * p.wall)
                push!(rows, @sprintf("| %s | %s | %d | %d | %s | %d | %.0f | %.2f | %.2f | %s | %.1f |",
                                     x.env.label, x.env.node, x.env.threads, x.env.gcthreads,
                                     p.loop, n, p.wall, eff, s, cols([p], x.env.threads, s), p.maxrss_gib))
            end
            x.predict === nothing || push!(rows, @sprintf("| %s | %s | %d | %d | predict | %d | %.0f | %.2f | %.2f | | | | | | | %.1f |",
                                                         x.env.label, x.env.node, x.env.threads,
                                                         x.env.gcthreads, x.predict.cells,
                                                         x.predict.wall,
                                                         x.predict.cells * s / (x.env.threads * x.predict.wall),
                                                         s, x.predict.maxrss_gib))
        end
    end
    procs = [deserialize(joinpath(OUT, f)).result for f in files if startswith(f, "procs_")]
    if !isempty(procs)
        wall = maximum(p.pass.wall for p in procs)
        n = sum(length(p.pass.cells) for p in procs)
        s = median(serial_med(p.base) for p in procs)
        T = sum(p.env.threads for p in procs)
        push!(rows, @sprintf("| procs %dx%d | %s | %d | %d | greedy | %d | %.0f | %.2f | %.2f | %s | %.1f |",
                             length(procs), procs[1].env.threads, procs[1].env.node, T,
                             procs[1].env.gcthreads, n, wall, n * s / (T * wall), s,
                             cols([p.pass for p in procs], T, s),
                             maximum(p.pass.maxrss_gib for p in procs)))
    end
    open(joinpath(OUT, "summary.md"), "w") do io
        println(io, "# KIM-68 threading diagnosis: summary, at $COMMIT\n")
        println(io, "Efficiency = cells × the node's serial median / (threads × wall). GC's share ",
                "of the wall is the process's GC clock over its wall (summed over processes ",
                "for procs). Idle share = 1 − Σ cell time / (threads × wall).\n")
        println(io, "| variant | node | threads | GC threads | loop | cells | wall (s) | efficiency | serial median (s) | per-cell slowdown | GC share of cell time | GC share of wall | GC pauses | mean pause (s) | idle share | MaxRSS (GiB) |")
        println(io, "|---|---|---|---|---|---|---|---|---|---|---|---|---|---|---|---|")
        foreach(r -> println(io, r), rows)
        if "serial.jls" in files
            x = load(OUT, "serial")
            println(io, "\nSerial ($(x.env.node)): ", join((@sprintf("%.2f s", c.time) for c in x.pass.cells), ", "),
                    @sprintf("; GC share %.3f; %.2f GiB per cell", sum(c.gctime for c in x.pass.cells) / sum(c.time for c in x.pass.cells),
                             median(c.bytes for c in x.pass.cells) / 2^30))
        end
        if "allocs.jls" in files
            a = load(OUT, "allocs")
            println(io, @sprintf("\nAllocations, one 60-s cell (%s): %.1f MiB in %.2f s (sample rate %.2f, %d sampled)\n",
                                 a.env.node, a.bytes / 2^20, a.time, a.sample_rate, a.sampled))
            println(io, "| MiB (est.) | count (est.) | first InferCell frame |\n|---|---|---|")
            foreach(s -> println(io, @sprintf("| %.1f | %d | `%s` |", s.bytes / 2^20, s.count, s.key)), a.sites)
            println(io, "\n| MiB (est.) | count (est.) | type |\n|---|---|---|")
            foreach(s -> println(io, @sprintf("| %.1f | %d | `%s` |", s.bytes / 2^20, s.count, first(s.key, 100))), a.types[1:min(15, end)])
        end
    end
    say("wrote $(joinpath(OUT, "summary.md"))")
end

function main()
    say("KIM-68 threading, stage $STAGE ($LABEL), at $COMMIT, job $JOB, ",
        "$(Threads.nthreads()) thread(s)", SMOKE ? ", SMOKE" : "")
    STAGE === :serial ? serial() : STAGE === :threaded ? threaded() :
    STAGE === :procs ? procs() : STAGE === :allocs ? allocs() :
    STAGE === :report ? report() : error("unknown stage $STAGE")
end

main()
