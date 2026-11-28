# KIM-56: an override run at the default values reproduces the default run bitwise
# (the Spec's full-cell correctness check).
#
# For each procid, three runs at the same seed and options:
#   published   run_cell(m, p; ..., params=nothing)            the published path
#   params      run_cell(m, p; ..., params=Params(m))          lane 10a's parameterised path
#   vector      run_cell(m, p; ..., params=params_from_vector(m, ones(838)))
# and every pair must agree under `same_output` (procid, seed, minutes, jitter, events,
# min count, stop, failure, the per-minute record, spans, hooks per minute) and
# `output_digest`. Two controls on the first procid must each give a different digest
# from its published run, so the check can fail: an ODE constant, PGI.kcatF × 2, over
# the full window, and a CME constant, translation_rate_restart.py:riboKcat × 2, over 3
# minutes. They run before the cells, so a dead control fails in minutes. Options are
# KIM-180's run of record (dev/scripts/fullcell_doubling_kim180):
# seed 1, minutes=181, ode_failure=:record, no stop. Exits non-zero on any mismatch.
#
# Environment: PROCIDS ("a-b", default 1-6), SEED (default 1), MINUTES (default 181),
# OUT (output directory; writes check.csv and meta.txt).

import SpeedyCell
using SpeedyCell.Cell
using InteractiveUtils, Printf, Dates, SHA
include(joinpath(@__DIR__, "adapter.jl"))


function int_range(name, default)
    f = parse.(Int, split(get(ENV, name, default), '-'))
    length(f) <= 2 && first(f) <= last(f) || error("$name must be a or a-b with a <= b")
    first(f):last(f)
end

function provenance(io)
    println(io, "date ", now(), "  host ", gethostname(), "  slurm_job ", get(ENV, "SLURM_JOB_ID", "none"))
    println(io, "infercell_git ", readchomp(`git -C $(@__DIR__) rev-parse HEAD`))
    println(io, "syn3a_ensemble_git ", readchomp(`git -C $SYN3A rev-parse HEAD`))
    st = readchomp(`git -C $SYN3A status --porcelain -- code/SpeedyCell data/model`)
    println(io, "syn3a_ensemble_code_data_clean ", isempty(st))
    println(io, "speedycell_path ", pathof(SpeedyCell))
    versioninfo(io)
end

function main()
    procids = int_range("PROCIDS", "1-6")
    seed = parse(Int, get(ENV, "SEED", "1"))
    minutes = parse(Int, get(ENV, "MINUTES", "181"))
    out = get(ENV, "OUT", "")
    isempty(out) && error("set OUT")
    mkpath(out)
    open(io -> provenance(io), joinpath(out, "meta.txt"), "w")
    provenance(stdout)

    m = CellModel()
    names = param_names()
    n = length(names)
    paths = (published=nothing, params=Params(m), vector=params_from_vector(m, ones(n); names=names))
    sim(p, P; minutes=minutes) = run_cell(m, p; seed=seed, minutes=minutes, params=P, ode_failure=:record)

    sim(first(procids), nothing; minutes=3)   # warm-up: compiles the published path
    sim(first(procids), paths.params; minutes=3)   # and the parameterised one
    p = first(procids)
    control_differs = true
    for (name, mins) in (("PGI.kcatF", minutes), ("translation_rate_restart.py:riboKcat", 3))
        i = findfirst(==(name), names)
        i === nothing && error("no named constant $name")
        v = ones(n); v[i] = 2.0
        dc = output_digest(sim(p, params_from_vector(m, v; names=names); minutes=mins))
        dp = output_digest(sim(p, nothing; minutes=mins))
        control_differs &= dc != dp
        @printf("control: procid %d with %s x 2 over %d min, digest %016x against %016x: %s\n", p, name, mins,
                dc, dp, dc != dp ? "differs" : "SAME (the check cannot fail)")
        flush(stdout)
    end
    control_differs || (println("CHECK FAIL"); exit(1))
    mismatches = 0
    rows = String[]
    for p in procids
        t = time()
        r = map(P -> sim(p, P), paths)
        d = map(output_digest, r)
        ok_pp = same_output(r.published, r.params) && d.published == d.params
        ok_pv = same_output(r.published, r.vector) && d.published == d.vector
        mismatches += !ok_pp + !ok_pv
        collapsed = any(<(59), r.published.hooks_per_minute)
        @printf("procid %d: digest %016x  params %s  vector %s  events %d  failure %s  collapsed %s  (%.0f s)\n",
                p, d.published, ok_pp ? "same" : "DIFFERENT", ok_pv ? "same" : "DIFFERENT",
                r.published.nevents, something(r.published.failure, "none"), collapsed, time() - t)
        push!(rows, @sprintf("%d,%d,%d,%016x,%016x,%016x,%s,%s,%d,%s,%s", p, seed, minutes, d.published,
                             d.params, d.vector, ok_pp, ok_pv, r.published.nevents,
                             r.published.failure === nothing ? "" : string(r.published.failure[1]), collapsed))
        flush(stdout)
    end
    open(joinpath(out, "check.csv"), "w") do io
        println(io, "procid,seed,minutes,digest_published,digest_params,digest_vector,same_params,same_vector,nevents,failure_minute,collapsed")
        foreach(l -> println(io, l), rows)
    end
    println("$(length(procids)) cells x 3 paths, $mismatches mismatches; controls ",
            control_differs ? "differ" : "do not differ")
    ok = mismatches == 0 && control_differs
    println(ok ? "CHECK PASS" : "CHECK FAIL")
    ok || exit(1)
end

main()
