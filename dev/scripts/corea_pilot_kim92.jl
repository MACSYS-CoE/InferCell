# KIM-92: the reduced-model pilot, the kill check (Linear M8). The gates, the design
# and the scoring targets are frozen in spec/prereg/kim-92.md at tag prereg/kim-92,
# before any of this ran; this script implements them and changes none.
#
# Usage: julia --project dev/scripts/corea_pilot_kim92.jl <stage> <out> [<arg>]
# Submit through dev/scripts/corea_pilot_kim92.sh, from a worktree at a pushed commit.
# Stages:
#   dataset <out>   the pilot's dataset, corea_dataset(17_200_001; ncells = 200,
#                   fit_to = COREA_FIT_S, targets = COREA_TARGETS), to dataset.jls. The
#                   full 6300 s record of the same cells is saved beside it
#                   (dataset_full.jls) for KIM-95, after checking the window is its
#                   prefix, so the second half is never simulated twice;
#   profile <out>   item 2: the log density and its 21-coordinate gradient at the
#                   pilot's setup (default tolerances, δ = S = 0), median of 3, with
#                   time, allocation and peak RSS, and one solve's step counts. Writes
#                   profile.jls and profile.md, with what a NUTS chain of the frozen
#                   design would cost at tree depths 3 to 7.
# SMOKE=1: 3 cells, a 300 s run cut at 150 s (two saves), never a result.

include(joinpath(@__DIR__, "pilot_common.jl"))
using InferCell: _mixed_solve
using ForwardDiff
using Printf
using Random
using Statistics: median

const SEED = 17_200_001                         # prereg §1
const NCELLS = SMOKE ? 3 : 200
const HORIZON = SMOKE ? 300.0 : COREA_CYCLE_S
const FIT_TO = SMOKE ? 150.0 : COREA_FIT_S
const NREPEAT = 3                               # item 2: median of 3
# The frozen chain design (prereg §1), for the projection only.
const NADAPT, NDRAW, NCHAIN = 200, 200, 4

dataset_file(out) = joinpath(out, "dataset.jls")

pilot_setup(ds) = mixed_setup(ds; models = d11_models(), targets = COREA_TARGETS)

function truth_x(ds, s)
    tv = Dict(ds.truth.values)
    return vcat([log(tv[n]) for n in s.names], log(ds.sigma))
end

# ---------------------------------------------------------------- dataset

function run_dataset(out)
    mkpath(out)
    kw = (; ncells = NCELLS, horizon = HORIZON, targets = COREA_TARGETS)
    t_win = @elapsed ds = corea_dataset(SEED; kw..., fit_to = FIT_TO)
    t_full = @elapsed full = corea_dataset(SEED; kw...)
    # The window is the full record's first saves (KIM-91's prefix test, here at size).
    k = 1:length(ds.times)
    ds.times == full.times[k] && ds.latent == full.latent[:, :, k] &&
        ds.transcripts == full.transcripts[:, :, k] && ds.observed == full.observed[:, :, k] &&
        ds.bulk == full.bulk[:, k] && ds.truth.values == full.truth.values &&
        ds.sigma == full.sigma && ds.seeds == full.seeds ||
        error("the window is not a prefix of the full record")
    save(out, "dataset", ds; stage = "dataset")
    save(out, "dataset_full", full; stage = "dataset")
    say(@sprintf("dataset: seed %d, %d cells, %d saves to %.0f s (window) in %.0f s", SEED,
                 NCELLS, length(ds.times), ds.times[end], t_win),
        @sprintf("; full record %d saves to %.0f s in %.0f s; %d bulk rows", length(full.times),
                 full.times[end], t_full, length(ds.bulk)))
end

# ---------------------------------------------------------------- profile

# `f()` run NREPEAT times: seconds, GiB allocated, GC seconds, and peak RSS after.
function timed(f)
    map(1:NREPEAT) do _
        st = @timed f()
        (; seconds = st.time, gib = st.bytes / 2^30, gc = st.gctime, maxrss_gib = Sys.maxrss() / 2^30,
         value = st.value)
    end
end

function run_profile(out)
    ds = load(out, "dataset")
    s = pilot_setup(ds)
    ℓ(z) = mixed_logdensity(s, z)
    xs = [:truth => truth_x(ds, s), :prior_median => mixed_prior_median(s)]
    # Compile on a two-cell, two-save case, so no timed call pays for it.
    t_compile = @elapsed begin
        tiny = pilot_setup(corea_dataset(SEED; ncells = 2, horizon = 120.0, targets = COREA_TARGETS))
        ForwardDiff.gradient(z -> mixed_logdensity(tiny, z), mixed_prior_median(tiny))
    end
    say(@sprintf("profile: compiled in %.0f s; %d coordinates, %d bulk rows, %d saves to %.0f s",
                 t_compile, length(xs[1][2]), length(s.logz), length(s.saves), s.saves[end]))
    rows = []
    for (label, x) in xs
        d = timed(() -> ℓ(x))
        say(@sprintf("density at %s: %s s; log density %.6f", label,
                     join((@sprintf("%.1f", r.seconds) for r in d), ", "), d[1].value))
        g = timed(() -> ForwardDiff.gradient(ℓ, x))
        say(@sprintf("gradient at %s: %s s; |g| %.4g", label,
                     join((@sprintf("%.1f", r.seconds) for r in g), ", "),
                     sqrt(sum(abs2, g[1].value))))
        sol = _mixed_solve(s, exp.(x[1:length(s.names)])).sol
        st = sol.stats
        push!(rows, (; label, density = [Base.structdiff(r, (; value = 0)) for r in d],
                     gradient = [Base.structdiff(r, (; value = 0)) for r in g],
                     logdensity = d[1].value, grad = g[1].value,
                     nf = st.nf, njacs = st.njacs, naccept = st.naccept, nreject = st.nreject))
        say(@sprintf("solve at %s: %d RHS calls, %d Jacobians, %d accepted and %d rejected steps",
                     label, st.nf, st.njacs, st.naccept, st.nreject))
    end
    save(out, "profile", (; rows, compile_seconds = t_compile, nrows = length(s.logz),
                          saves = s.saves, abstol = s.abstol, reltol = s.reltol); stage = "profile")
    write_profile_md(out, rows, s)
end

function write_profile_md(out, rows, s)
    io = IOBuffer()
    p(xs...) = println(io, xs...)
    p("# KIM-92 item 2: the pilot's density and gradient\n")
    p("Commit `", COMMIT, "`, job ", JOB, ", on ", RAN_ON.host, " (", RAN_ON.cpu, "), ",
      Threads.nthreads(), " thread(s). SMOKE = ", SMOKE, ". The setup is the pilot's: ",
      NCELLS, " cells, ", length(s.saves), " saves to ", s.saves[end], " s, ", length(s.logz),
      " bulk rows, 21 coordinates, route (ii) at abstol ", s.abstol, " and reltol ", s.reltol,
      ", δ = S = 0. Each call ran ", NREPEAT, " times after compiling on a two-cell case. ",
      "Peak RSS is the process's, after the call.\n")
    p("| x | call | seconds (each) | median s | GiB allocated (median) | GC share | peak RSS GiB |")
    p("|---|---|---|---|---|---|---|")
    for r in rows, (what, v) in (("log density", r.density), ("gradient", r.gradient))
        sec = [q.seconds for q in v]
        p(@sprintf("| %s | %s | %s | %.1f | %.1f | %.3f | %.2f |", r.label, what,
                   join((@sprintf("%.1f", t) for t in sec), ", "), median(sec),
                   median([q.gib for q in v]), median([q.gc / q.seconds for q in v]),
                   maximum(q.maxrss_gib for q in v)))
    end
    p("\nOne route (ii) solve at each x:\n")
    p("| x | RHS calls | Jacobians | accepted steps | rejected steps | RHS calls per step |")
    p("|---|---|---|---|---|---|")
    for r in rows
        p(@sprintf("| %s | %d | %d | %d | %d | %.1f |", r.label, r.nf, r.njacs, r.naccept,
                   r.nreject, r.nf / max(r.naccept + r.nreject, 1)))
    end
    g = median([q.seconds for q in only(r for r in rows if r.label === :truth).gradient])
    p("\n## What the frozen design would cost (projection, not a measurement)\n")
    p("One NUTS iteration at tree depth d takes up to 2^d − 1 gradients. A chain is ",
      NADAPT, " + ", NDRAW, " iterations; the design is ", NCHAIN,
      " chains (prereg §1). At the median gradient at the truth, ", @sprintf("%.0f", g), " s:\n")
    p("| tree depth | gradients per iteration | days per chain | CPU-h, ", NCHAIN, " chains |")
    p("|---|---|---|---|")
    for d in 3:7
        n = 2^d - 1
        sec = g * n * (NADAPT + NDRAW)
        p(@sprintf("| %d | %d | %.1f | %.0f |", d, n, sec / 86_400, NCHAIN * sec / 3600))
    end
    p("\nmilan's limit is 7 days per job. The depth the chains reach is not known until they run.")
    write(joinpath(out, "profile.md"), String(take!(io)))
    print(read(joinpath(out, "profile.md"), String))
end

function main(args)
    stage, out = args[1], args[2]
    if stage == "dataset"
        run_dataset(out)
    elseif stage == "profile"
        run_profile(out)
    else
        error("unknown stage $stage")
    end
end

main(ARGS)
