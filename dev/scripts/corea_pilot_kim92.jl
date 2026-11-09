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
#                   design would cost at tree depths 3 to 7;
#   prep <out>      θ_ref (the bounded MAP at δ = S = 0), the 800-cell reference
#                   ensemble at θ_ref to the window's last save, δ and S (C = 200);
#   chain <out> c   chain c of 4: NUTS(NADAPT, 0.8), NDRAW kept draws, from a start
#                   jittered by 0.25 prior SD about θ_ref, seed mixed_chain_seeds(ds, 4)[c];
#   check <out> k   the convergence gate on each chain's first k kept draws (k = NDRAW,
#                   or 2 NDRAW after the extension), to gate_<k>.md;
#   extend <out> c  the one extension (prereg §1): chain c continued NDRAW more draws
#                   from its saved state, no new adaptation, seed arm 1's c-th;
#   ppc <out>       the PPC (prereg §2.2): 200 posterior draws evenly spaced through the
#                   pooled kept draws, 10 replicates of log z_rep each, threaded;
#   merge <out>     the gates, recovery, contraction and cost, to pilot_summary.md and
#                   the TSVs the figures read.
# SMOKE=1: 3 cells, a 300 s run cut at 150 s (two saves), 8 reference cells, 5 + 5
# NUTS steps at tree depth ≤ 3, 10 PPC draws × 2: a plumbing pass, never a result. NADAPT and NDRAW can
# be set for a run, but the frozen design is 200 + 200 (a change is a Decision log entry).

include(joinpath(@__DIR__, "pilot_common.jl"))
using InferCell: _mixed_solve
using ForwardDiff
using Printf
using Random
using LinearAlgebra
using Statistics: median, std, cor, quantile
using Turing: sample, NUTS
using MCMCChains
using MCMCChains: MCMCDiagnosticTools
const MDT = MCMCDiagnosticTools

envint(k, default) = (v = strip(get(ENV, k, "")); isempty(v) ? default : parse(Int, v))

# prereg §1
const SEED = 17_200_001
const NCELLS = SMOKE ? 3 : 200
const HORIZON = SMOKE ? 300.0 : COREA_CYCLE_S
const FIT_TO = SMOKE ? 150.0 : COREA_FIT_S
const NREF = SMOKE ? 8 : 800
const MAP_ITERS = SMOKE ? 5 : 200
const NADAPT = envint("NADAPT", SMOKE ? 5 : 200)
const NDRAW = envint("NDRAW", SMOKE ? 5 : 200)
# Turing's default, 10, for the run. A smoke caps it at 3: plumbing only, since one
# deep iteration took over 3 h at the smoke's size (jobs 20453041–44).
const MAX_DEPTH = SMOKE ? 3 : 10
nuts() = NUTS(NADAPT, TARGET_ACCEPT; max_depth = MAX_DEPTH)
const NCHAIN = 4
const TARGET_ACCEPT = 0.8
const JITTER_SD = 0.25
const ARM_CHAIN, ARM_EXTEND = 0, 1
# prereg §2 and §3
const RHAT_MAX, ESS_MIN = 1.02, 100
const PPC_DRAWS, PPC_REPS = SMOKE ? (10, 2) : (200, 10)
const PPC_LO, PPC_HI, PPC_MAX_OUT_SHARE = 0.005, 0.995, 0.05
const CI_LO, CI_HI = 0.05, 0.95
const REACHED = 0.5
const NREPEAT = 3                               # item 2: median of 3

pilot_setup(ds; kw...) = mixed_setup(ds; models = d11_models(), targets = COREA_TARGETS, kw...)

function truth_x(ds, s)
    tv = Dict(ds.truth.values)
    return vcat([log(tv[n]) for n in s.names], log(ds.sigma))
end

has(out, name) = isfile(joinpath(out, "$name.jls"))

# The chains' setup: the dataset with the fit's δ and S.
function fit_setup(out)
    ds = load(out, "dataset")
    p = load(out, "prep")
    return (; ds, p, s = with_reference(pilot_setup(ds), p.delta, p.S))
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

# ---------------------------------------------------------------- prep

function run_prep(out)
    ds = load(out, "dataset")
    s0 = pilot_setup(ds)
    t_map = @elapsed ref = reference_theta(s0; maxiters = MAP_ITERS)
    say(@sprintf("prep: θ_ref by %s in %.0f s, log density %.4f", ref.method, t_map, ref.logdensity))
    values = mixed_theta_values(s0, ref.x)
    m_ref = mean_field_log_bulk(s0, _mixed_solve(s0, exp.(ref.x[1:length(s0.names)])).sol)
    seeds = reference_seeds(ds, NREF)
    t_ens = @elapsed cells = reference_ensemble(values, seeds; models = d11_models(),
                                                horizon = ds.times[end])
    size(cells, 3) == length(ds.times) || error("the ensemble's saves are not the window's")
    d = delta_and_S(cells, m_ref; C = NCELLS)
    save(out, "prep", (; ref, values, m_ref, seeds, cells, delta = d.delta, S = d.S,
                       seconds = (map = t_map, ensemble = t_ens)); stage = "prep")
    say(@sprintf("prep: %d reference cells to %.0f s in %.0f s; tr S %.4g; max |δ| %.4g",
                 length(seeds), ds.times[end], t_ens, tr(d.S), maximum(abs, d.delta)))
end

# ---------------------------------------------------------------- chains

# A line every `every` iterations, so a days-long chain shows its progress.
function progress_callback(label, total; every = SMOKE ? 1 : 5)
    t0 = time()
    return function (rng, model, sampler, transition, state, i; kwargs...)
        (i % every == 0 || i == total) || return nothing
        st = transition.stat
        say(@sprintf("%s: iteration %d of %d, %.1f h; depth %d, %d steps, step size %.3g, lp %.2f",
                     label, i, total, (time() - t0) / 3600, st.tree_depth, st.n_steps,
                     st.step_size, transition.lp))
        return nothing
    end
end

function run_chain(out, c)
    f = fit_setup(out)
    seed = mixed_chain_seeds(f.ds, NCHAIN; arm = ARM_CHAIN)[c]
    rng = Xoshiro(seed)
    st = mixed_jittered_start(rng, f.s, f.p.ref.x; scale = JITTER_SD)
    say(@sprintf("chain %d: seed %d, start after %d tries, log density %.4f", c, seed, st.tries,
                 mixed_logdensity(f.s, st.x)))
    # The adaptation draws are kept in the chain (`discard_adapt = false`) so that the
    # callback, which AbstractMCMC calls only for kept draws, reports them too: the
    # adaptation is most of a chain's days. `chain_draws` drops them. With them kept,
    # the N Turing is given counts them, and the kept draws equal a default run's
    # bitwise (checked on a toy model, 2026-11-10).
    t = @elapsed chain = sample(rng, mixed_model(f.s), nuts(), NADAPT + NDRAW;
                                initial_params = st.x, progress = false, save_state = true,
                                discard_adapt = false,
                                callback = progress_callback("chain $c", NADAPT + NDRAW))
    size(chain, 1) == NADAPT + NDRAW ||
        error("chain $c holds $(size(chain, 1)) draws, not $NADAPT + $NDRAW")
    save(out, "chain_$c", (; c, seed, start = st.x, start_tries = st.tries, chain, seconds = t,
                           nadapt = NADAPT, ndraw = NDRAW); stage = "chain")
    say(@sprintf("chain %d: %d + %d iterations in %.1f h", c, NADAPT, NDRAW, t / 3600))
end

function run_extend(out, c)
    f = fit_setup(out)
    base = load(out, "chain_$c")
    seed = mixed_chain_seeds(f.ds, NCHAIN; arm = ARM_EXTEND)[c]
    rng = Xoshiro(seed)
    # No new adaptation: `resume_from` continues from the saved state, with its step
    # size and mass matrix, and Turing runs it with nadapts = 0. (Turing 0.39 ignores
    # `initial_state` unless `resume_from` is given, and would start afresh.)
    t = @elapsed chain = sample(rng, mixed_model(f.s), nuts(), base.ndraw;
                                resume_from = base.chain, progress = false, save_state = true,
                                callback = progress_callback("extend $c", base.ndraw))
    ε0, ε1 = base.chain[:step_size][end], unique(chain[:step_size])
    ε1 == [ε0] || error("the extension's step size $ε1 is not the chain's adapted $ε0")
    save(out, "extend_$c", (; c, seed, chain, seconds = t, ndraw = base.ndraw); stage = "extend")
    say(@sprintf("extend %d: %d more iterations in %.1f h", c, base.ndraw, t / 3600))
end

# Draws × parameters × chains on the log scale, mixed_parameter_names' order: each
# chain's kept draws, after its `skip` adaptation draws, followed by its extension's
# when it has run.
function chain_draws(ch, P; skip = 0)
    Array(ch[(skip + 1):end, [Symbol("x[$i]") for i in 1:P], 1])
end

function posterior_draws(out, P; extended = all(c -> has(out, "extend_$c"), 1:NCHAIN))
    X = map(1:NCHAIN) do c
        r = load(out, "chain_$c")
        size(r.chain, 1) == r.nadapt + r.ndraw ||
            error("chain $c holds $(size(r.chain, 1)) draws, not $(r.nadapt) + $(r.ndraw)")
        a = chain_draws(r.chain, P; skip = r.nadapt)
        extended ? vcat(a, chain_draws(load(out, "extend_$c").chain, P)) : a
    end
    return cat(X...; dims = 3), extended
end

pooled(X) = reshape(permutedims(X, (1, 3, 2)), :, size(X, 2))

# Rank-normalised split R̂ and bulk ESS per parameter (B7's `diagnostics`).
function diagnostics(X)
    A = permutedims(X, (1, 3, 2))                 # draws × chains × parameters
    return MDT.rhat(A), MDT.ess(A; kind = :bulk)
end

function gate(X, pnames)
    rh, es = diagnostics(X)
    fail = [i for i in eachindex(pnames) if !(rh[i] <= RHAT_MAX && es[i] >= ESS_MIN)]
    return (; rh, es, fail, pass = isempty(fail))
end

function run_check(out, k)
    f = fit_setup(out)
    pnames = mixed_parameter_names(f.s)
    X, ext = posterior_draws(out, length(pnames))
    size(X, 1) >= k || error("only $(size(X, 1)) kept draws per chain, not $k")
    g = gate(X[1:k, :, :], pnames)
    io = IOBuffer()
    p(xs...) = (println(io, xs...); println(xs...))
    p("# KIM-92: the convergence gate on ", k, " kept draws per chain\n")
    p("Commit `", COMMIT, "`, job ", JOB, ". ", NCHAIN, " chains; extension ",
      ext ? "included" : "not run", ". Gate (prereg §2.1): R̂ ≤ ", RHAT_MAX, " and bulk ESS ≥ ",
      ESS_MIN, " on every coordinate.\n")
    p("| coordinate | R̂ | bulk ESS | passes |")
    p("|---|---|---|---|")
    for i in eachindex(pnames)
        p(@sprintf("| %s | %.4f | %.1f | %s |", pnames[i], g.rh[i], g.es[i],
                   i in g.fail ? "no" : "yes"))
    end
    p("\n**The gate ", g.pass ? "passes" : "FAILS", "** on ", k, " draws per chain: ",
      length(pnames) - length(g.fail), " of ", length(pnames), " coordinates pass.",
      g.pass || ext ? "" : " The prereg's one extension (stage `extend`) runs next.")
    write(joinpath(out, "gate_$k.md"), String(take!(io)))
end

# ---------------------------------------------------------------- PPC

function run_ppc(out)
    f = fit_setup(out)
    s = f.s
    pnames = mixed_parameter_names(s)
    X, ext = posterior_draws(out, length(pnames))
    flat = pooled(X)
    idx = round.(Int, range(1, size(flat, 1); length = PPC_DRAWS))
    n = length(s.logz)
    Z = fill(NaN, n, PPC_DRAWS * PPC_REPS)
    failed = zeros(Bool, PPC_DRAWS)
    t = @elapsed Threads.@threads for j in eachindex(idx)
        x = flat[idx[j], :]
        θ = exp.(x[1:length(s.names)])
        σ = exp(x[end])
        r = _mixed_solve(s, θ)
        if r === nothing
            failed[j] = true
            continue
        end
        m = mean_field_log_bulk(s, r.sol)
        C = cholesky(Symmetric(s.S + σ^2 * I)).L
        rng = Xoshiro(SEED + 92_000_000 + j)           # one stream per draw
        for k in 1:PPC_REPS
            Z[:, (j - 1) * PPC_REPS + k] = m .+ s.delta .+ C * randn(rng, n)
        end
    end
    any(failed) && error("route (ii) failed at $(count(failed)) PPC draws")
    lo = [quantile(Z[i, :], PPC_LO) for i in 1:n]
    hi = [quantile(Z[i, :], PPC_HI) for i in 1:n]
    outside = (s.logz .< lo) .| (s.logz .> hi)
    save(out, "ppc", (; idx, lo, hi, outside, logz = s.logz, extended = ext, seconds = t,
                      draws = PPC_DRAWS, reps = PPC_REPS); stage = "ppc")
    say(@sprintf("ppc: %d draws × %d in %.0f s on %d threads; %d of %d rows outside",
                 PPC_DRAWS, PPC_REPS, t, Threads.nthreads(), count(outside), n))
end

# ---------------------------------------------------------------- merge

# CPU-h by sacct for the jobs in jobs.tsv, by name.
function sacct_cpuh(out)
    file = joinpath(out, "jobs.tsv")
    isfile(file) || return Dict{String, Float64}()
    cost = Dict{String, Float64}()
    for l in readlines(file)
        name, id = split(l, '\t')
        txt = try
            read(`sacct -j $id -X -n -P --format=CPUTimeRAW`, String)
        catch
            continue
        end
        v = sum((parse(Float64, strip(x)) for x in split(strip(txt), '\n') if !isempty(strip(x)));
                init = 0.0) / 3600
        cost[String(name)] = get(cost, String(name), 0.0) + v
    end
    return cost
end

function run_merge(out)
    f = fit_setup(out)
    s, ds = f.s, f.ds
    pnames = mixed_parameter_names(s)
    P = length(pnames)
    X, ext = posterior_draws(out, P)
    first_gate = gate(X[1:NDRAW, :, :], pnames)
    g = ext ? gate(X, pnames) : first_gate
    flat = pooled(X)
    ppc = load(out, "ppc")
    ppc.extended == ext || error("the PPC ran on a different set of draws from the merge's")
    n_out = count(ppc.outside)
    ppc_pass = n_out <= floor(Int, PPC_MAX_OUT_SHARE * length(ppc.outside))
    tx = truth_x(ds, s)
    psd = mixed_prior_sd(s)
    μ, sd = vec(mean(flat; dims = 1)), vec(std(flat; dims = 1))
    lo90 = [quantile(flat[:, q], CI_LO) for q in 1:P]
    hi90 = [quantile(flat[:, q], CI_HI) for q in 1:P]
    inside = lo90 .<= tx .<= hi90
    contraction = 1 .- sd .^ 2 ./ psd .^ 2
    R = cor(flat)
    cost = sacct_cpuh(out)

    io = IOBuffer()
    p(xs...) = (println(io, xs...); println(xs...))
    p("# KIM-92: the reduced-model pilot\n")
    p("Merge at `", COMMIT, "`, job ", JOB, ". SMOKE = ", SMOKE, ". Frozen at `prereg/kim-92`. ",
      "Dataset seed ", SEED, ", ", NCELLS, " cells, ", length(s.saves), " saves to ",
      s.saves[end], " s, ", length(s.logz), " bulk rows; ", NCHAIN, " chains of ", NADAPT,
      " + ", NDRAW, ext ? ", extended once by $(NDRAW)" : "", ".\n")
    p("## The gates\n")
    p("| gate | criterion | result | verdict |")
    p("|---|---|---|---|")
    p(@sprintf("| convergence, %d draws per chain | R̂ ≤ %.2f and bulk ESS ≥ %d, every coordinate | worst R̂ %.4f, smallest ESS %.1f, %d failing | %s |",
               NDRAW, RHAT_MAX, ESS_MIN, maximum(first_gate.rh), minimum(first_gate.es),
               length(first_gate.fail), first_gate.pass ? "pass" : "fail"))
    ext && p(@sprintf("| convergence after the extension, %d per chain | as above | worst R̂ %.4f, smallest ESS %.1f, %d failing | %s |",
                      size(X, 1), maximum(g.rh), minimum(g.es), length(g.fail),
                      g.pass ? "pass" : "fail"))
    p(@sprintf("| PPC | at most %d of %d rows outside the central 99%% band | %d outside | %s |",
               floor(Int, PPC_MAX_OUT_SHARE * length(ppc.outside)), length(ppc.outside), n_out,
               ppc_pass ? "pass" : "fail"))
    p("\n**The kill check ", g.pass && ppc_pass ? "passes" : "FAILS", "** on the frozen gates. ",
      "Cost has no cap; KIM-93 decides.\n")
    p("## Per coordinate (recovery and contraction are reported, not gated)\n")
    p("| coordinate | R̂ | bulk ESS | prior sd | posterior mean ± sd | truth | 90% interval | truth inside | (mean − truth)/sd | contraction |")
    p("|---|---|---|---|---|---|---|---|---|---|")
    for q in 1:P
        p(@sprintf("| %s | %.4f | %.1f | %.3f | %.4f ± %.4f | %.4f | [%.4f, %.4f] | %s | %+.2f | %.3f |",
                   pnames[q], g.rh[q], g.es[q], psd[q], μ[q], sd[q], tx[q], lo90[q], hi90[q],
                   inside[q] ? "yes" : "no", (μ[q] - tx[q]) / sd[q], contraction[q]))
    end
    p(@sprintf("\nTruth inside the 90%% interval: %d of %d. Reached (contraction ≥ %.1f): %d of %d.",
               count(inside), P, REACHED, count(>=(REACHED), contraction), P))
    strong = [(pnames[i], pnames[j], R[i, j]) for i in 1:P for j in (i + 1):P if abs(R[i, j]) >= 0.5]
    sort!(strong; by = t -> -abs(t[3]))
    p("\nPosterior correlations with |r| ≥ 0.5: ",
      isempty(strong) ? "none." : join((@sprintf("%s–%s %.2f", a, b, r) for (a, b, r) in strong), "; "))
    p("\n## Cost\n")
    if isempty(cost)
        p("No sacct records (jobs.tsv missing or sacct unavailable).")
    else
        p("| jobs | CPU-h (sacct) |")
        p("|---|---|")
        for k in sort(collect(keys(cost)))
            p(@sprintf("| %s | %.2f |", k, cost[k]))
        end
        p(@sprintf("| all | %.2f |", sum(values(cost))))
    end
    p(@sprintf("\nStage seconds: MAP %.0f, ensemble %.0f, chains %s, PPC %.0f.", f.p.seconds.map,
               f.p.seconds.ensemble,
               join((@sprintf("%.0f", load(out, "chain_$c").seconds) for c in 1:NCHAIN), ", "),
               ppc.seconds))
    write(joinpath(out, "pilot_summary.md"), String(take!(io)))

    open(joinpath(out, "pilot_coordinates.tsv"), "w") do o
        println(o, join(["coordinate", "rhat", "ess", "prior_mean", "prior_sd", "post_mean",
                         "post_sd", "truth", "lo90", "hi90", "contraction"], '\t'))
        pm = mixed_prior_median(s)
        for q in 1:P
            println(o, join(vcat(String(pnames[q]),
                                 [@sprintf("%.10g", v) for v in (g.rh[q], g.es[q], pm[q], psd[q],
                                                                  μ[q], sd[q], tx[q], lo90[q],
                                                                  hi90[q], contraction[q])]), '\t'))
        end
    end
    open(joinpath(out, "pilot_draws.tsv"), "w") do o
        println(o, join(String.(pnames), '\t'))
        foreach(i -> println(o, join((@sprintf("%.17g", v) for v in flat[i, :]), '\t')), axes(flat, 1))
    end
    open(joinpath(out, "pilot_ppc.tsv"), "w") do o
        println(o, "row\tlogz\tlo\thi\toutside")
        for i in eachindex(ppc.logz)
            println(o, i, '\t', ppc.logz[i], '\t', ppc.lo[i], '\t', ppc.hi[i], '\t', Int(ppc.outside[i]))
        end
    end
end

function main(args)
    stage, out = args[1], args[2]
    if stage == "dataset"
        run_dataset(out)
    elseif stage == "profile"
        run_profile(out)
    elseif stage == "prep"
        run_prep(out)
    elseif stage == "chain"
        run_chain(out, parse(Int, args[3]))
    elseif stage == "check"
        run_check(out, parse(Int, args[3]))
    elseif stage == "extend"
        run_extend(out, parse(Int, args[3]))
    elseif stage == "ppc"
        run_ppc(out)
    elseif stage == "merge"
        run_merge(out)
    else
        error("unknown stage $stage")
    end
end

main(ARGS)
