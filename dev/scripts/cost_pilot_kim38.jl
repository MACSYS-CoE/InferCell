# KIM-38: cost the evidence package and gate the run (M4). A pilot on θ held out
# of the evidence run: B7's `v6` fit (m0_v6, not one of the 50 evidence
# datasets), with every rng seed from M4's reserved block 18_438_xxx. Nothing
# here reads an evidence posterior or a coverage.
#
# Usage (dev/scripts/cost_pilot_kim38.sh submits the first four):
#   julia --project dev/scripts/cost_pilot_kim38.jl <stage> <outdir> [<kim35_outdir>]
# Stages:
#   cells      warm per-cell wall time, one thread: the M0 hybrid at 600 s, at
#              baseline and under a 25% ptsG knockdown (20 cells each); the M0
#              mean field (10 solves); the d11 hybrid at 6,300 s and at 600 s
#              (10 cells each); the d11 mean field at 6,300 s (5 solves)
#   threads    a d11 crossed matrix at 6,300 s (M = 4 × K = 8 × 4 cells = 128
#              cells), threaded; run at 16 and at 64 threads for the efficiency
#   sizing     an M0 `predict` on v6 at 600 s, 20 draws × 10 cells, hybrid and
#              mean field: Var_θ[μ] and σ² for the sizing rules, for both Q
#              candidates
#   project    reads the three above and KIM-35's arms, and projects the package
#
# Rules, declared before the run (M4 plan, approved 2026-10-30):
# - Sizing: the MC term is at most 10% of the spread over draws, so
#   n_cells ≥ 9 σ²/Var_θ; the truth's SE and δ's SE are each at most 10% of the
#   predictive SD, so n_truth, n_ref ≥ 100 σ²/Var_θ. n_cells is at least 4.
# - The package (completion plan §3.5; Linear KIM-40/41/42):
#   M0 crossed: 50 datasets × 2 posteriors × (hybrid: 100 draws × n_cells cells;
#   mean field: 100 solves; δ: n_ref cells per dataset × posterior), plus the
#   truth μ(θ*) with n_truth cells per dataset; the intervention doubles the
#   hybrid, δ and truth cells and the solves. The reduced-model budget:
#   100 × K = 20 × n_cells (d11 sizing from KIM-35) at 3 widths, with V's own
#   cells (doubling it) shown as an option, and an intervened arm as an option.
# - Wall-clock: the M0 parts as single-core array tasks at a stated concurrency;
#   the budget as one threaded 64-core node per width, at the measured efficiency.
# - The gate (KIM-38): if the projection runs past ~Nov 20, (1) drop the
#   intervention panel, (2) next, run 2 widths, not 3. Each cut goes into §12.

using InferCell
using InferCell: AbstractSubModel
using Random
using Serialization
using Statistics: mean, median, quantile, var
using Printf
using Dates

const STAGE = Symbol(ARGS[1])
const OUT = ARGS[2]
const COMMIT = strip(read(`git rev-parse --short HEAD`, String))
const JOB = get(ENV, "SLURM_JOB_ID", "none")
const V6 = ["/fred/oz022/tkimpson/infercell-89-b7-archive-2026-10-27/run/v6/chain_$c.jls" for c in 1:2]
const SMOKE = get(ENV, "SMOKE", "0") == "1"      # a plumbing pass at tiny sizes, not a result
const H_M0 = SMOKE ? 60.0 : 600.0
const H_CYCLE = SMOKE ? 120.0 : 6300.0
const SCALE = SMOKE ? 0.15 : 1.0                  # run counts, scaled down in a smoke pass
n_of(n) = max(2, round(Int, SCALE * n))

say(x...) = (println(now(), " ", x...); flush(stdout))

save(name, x) = serialize(joinpath(OUT, "$name.jls"),
                          (; result = x, commit = COMMIT, job = JOB, stage = STAGE,
                           threads = Threads.nthreads(), when = now()))

# Wall time of each run, one at a time, after one untimed run to compile.
function time_runs(label, ms, values, n; forward, horizon, seed0)
    f = InferCell._forward_setup(metered(ms), forward; horizon)
    run1(s) = forward === :hybrid ? InferCell._hybrid_cell(f, values, s, :lactate_out) :
                                   InferCell._mean_field_cell(f, values, :lactate_out)
    InferCell._restoring_global_rng(() -> run1(seed0))
    ts, st = Float64[], Dict(s => 0 for s in CELL_STATUSES)
    for i in 1:n
        local c
        t = @elapsed c = InferCell._restoring_global_rng(() -> run1(seed0 + i))
        push!(ts, t)
        st[cell_status(c)] += 1
        say("  $label $i/$n: $(round(t; digits = 2)) s, $(c.failure), q = $(c.q)")
    end
    r = (; label, forward, horizon, n, median = median(ts), p90 = quantile(ts, 0.9),
         max = maximum(ts), mean = mean(ts), times = ts, statuses = st,
         maxrss_gib = Sys.maxrss() / 2^30)
    say("$label: median $(round(r.median; digits = 2)) s, p90 $(round(r.p90; digits = 2)) s")
    return r
end

function cells()
    pd = posterior_draws(V6)
    ref = InferCell._theta_ref(pd)
    mm = metered(m0_models())
    dm = metered(d11_models())
    m0v = haldane_values(mm, ref)
    m0k = intervene(Knockdown(; dose = 0.25), mm, ref)
    d11v = haldane_values(dm, ref)
    rows = [time_runs("m0_hybrid_600", mm, m0v, n_of(20); forward = :hybrid, horizon = H_M0,
                      seed0 = 18_438_100),
            time_runs("m0_hybrid_600_knockdown", mm, m0k, n_of(20); forward = :hybrid, horizon = H_M0,
                      seed0 = 18_438_100),
            time_runs("m0_mean_field_600", mm, m0v, n_of(10); forward = :mean_field, horizon = H_M0,
                      seed0 = 0),
            time_runs("d11_hybrid_6300", dm, d11v, n_of(10); forward = :hybrid, horizon = H_CYCLE,
                      seed0 = 18_438_200),
            time_runs("d11_hybrid_600", dm, d11v, n_of(10); forward = :hybrid, horizon = H_M0,
                      seed0 = 18_438_200),
            time_runs("d11_mean_field_6300", dm, d11v, n_of(5); forward = :mean_field,
                      horizon = H_CYCLE, seed0 = 0)]
    save("cells", rows)
end

function threads()
    pd = posterior_draws(V6)
    dms = d11_models()
    kw = (; completion = :prior, threads = true)
    # Compiles the threaded path, at a short horizon.
    predict(dms, pd; kw..., horizon = 60.0, n_held = 2, n_parameter_draws = 2, n_cells = 2,
            rng = Xoshiro(18_438_300))
    say("compiled; running the 128-cell matrix on $(Threads.nthreads()) threads")
    local p
    t = @elapsed p = predict(dms, pd; kw..., horizon = H_CYCLE, n_held = SMOKE ? 2 : 8,
                             n_parameter_draws = SMOKE ? 2 : 4, n_cells = SMOKE ? 2 : 4,
                             rng = Xoshiro(18_438_301))
    r = (; threads = Threads.nthreads(), wall = t, cells = sum(length, p.cells),
         statuses = p.statuses, maxrss_gib = Sys.maxrss() / 2^30)
    say("threads $(r.threads): $(r.cells) cells in $(round(t; digits = 1)) s, ",
        "maxrss $(round(r.maxrss_gib; digits = 1)) GiB, statuses $(r.statuses)")
    save("threads_$(Threads.nthreads())", r)
end

# μ̂ and σ̂² of a quantity read from each pair's finished cells.
function moments(p::Prediction, q::Symbol)
    xs = [[getproperty(c.accounts, q) for c in cs if c.failure === :none] for cs in p.cells]
    mu = [isempty(x) ? NaN : mean(x) for x in xs]
    s2 = [length(x) >= 2 ? var(x) : NaN for x in xs]
    return mu, s2, [length(x) for x in xs]
end

function sizing_numbers(p, q)
    mu, s2, n = moments(p, q)
    ok = isfinite.(mu)
    v = var(mu[ok])
    mc = p.forward === :hybrid ? mean(s2[ok] ./ n[ok]) : 0.0
    σ2 = p.forward === :hybrid ? mean(s2[ok]) : NaN
    pv = v - mc
    return (; q, forward = p.forward, mean = mean(mu[ok]), var = v, mc, var_theta = pv,
            sigma2 = σ2, n_cells_rule = isfinite(σ2) && pv > 0 ? ceil(Int, 9σ2 / pv) : missing,
            n_ref_rule = isfinite(σ2) && pv > 0 ? ceil(Int, 100σ2 / pv) : missing)
end

function sizing()
    pd = posterior_draws(V6)
    ms = m0_models()
    hy = predict(ms, pd; horizon = H_M0, n_parameter_draws = n_of(20), n_cells = n_of(10),
                 rng = Xoshiro(18_438_400))
    say("hybrid done: statuses $(hy.statuses)")
    mf = predict(ms, pd; horizon = H_M0, forward = :mean_field, n_parameter_draws = n_of(20),
                 rng = Xoshiro(18_438_400))
    say("mean field done")
    rows = [sizing_numbers(p, q) for p in (hy, mf) for q in (:lactate_out, :exported)]
    for r in rows
        say(r)
    end
    save("sizing", (; rows, hybrid = hy, mean_field = mf))
end

load(name, dir = OUT) = deserialize(joinpath(dir, "$name.jls"))

function project()
    kim35 = length(ARGS) >= 3 ? ARGS[3] : error("project needs KIM-35's outdir")
    c = Dict(r.label => r for r in load("cells").result)
    th = Dict(r.result.threads => r.result for r in (load("threads_16"), load("threads_64")))
    sz = load("sizing").result.rows
    t1 = c["d11_hybrid_6300"].median
    eff = Dict(k => v.cells * t1 / (k * v.wall) for (k, v) in th)
    say("threaded efficiency (per-cell median $(round(t1; digits = 1)) s at 1 thread): ", eff)
    tm0, tm0k = c["m0_hybrid_600"].median, c["m0_hybrid_600_knockdown"].median
    tmf = c["m0_mean_field_600"].median
    D, P, M = 50, 2, 100
    io = open(joinpath(OUT, "projection.md"), "w")
    out(x...) = (println(io, x...); println(x...))
    out("# KIM-38 projection at $COMMIT, job $JOB, $(now())\n")
    out("Per run (median, one thread): M0 hybrid 600 s $(round(tm0; digits = 2)) s, ",
        "knocked down $(round(tm0k; digits = 2)) s, M0 mean field $(round(tmf; digits = 2)) s, ",
        "d11 hybrid 6,300 s $(round(t1; digits = 1)) s, d11 hybrid 600 s ",
        "$(round(c["d11_hybrid_600"].median; digits = 2)) s. Threaded efficiency: ",
        join(["$(k) threads $(round(eff[k]; digits = 2))" for k in sort(collect(keys(eff)))], ", "),
        ".\n")
    for q in (:lactate_out, :exported)
        h = only(r for r in sz if r.q === q && r.forward === :hybrid)
        n = max(4, coalesce(h.n_cells_rule, 4))
        nref = max(2, coalesce(h.n_ref_rule, 100))
        ntruth = nref
        m0 = D * P * M * n * tm0 + D * P * M * tmf + D * P * nref * tm0 + D * ntruth * tm0 + D * tmf
        iv = D * P * M * n * tm0k + D * P * M * tmf + D * P * nref * tm0k + D * ntruth * tm0k + D * tmf
        out("## Q = :$q\n")
        out(@sprintf("M0 sizing (v6, 20 draws × 10 cells): Var_θ %.4g, σ² %.4g, so n_cells %d, n_ref = n_truth %d.\n",
                     h.var_theta, h.sigma2, n, nref))
        out("| part | CPU-h | wall at 100 tasks (h) |\n|---|---|---|")
        out(@sprintf("| M0 crossed (KIM-40) | %.1f | %.1f |", m0 / 3600, m0 / 3600 / 100))
        out(@sprintf("| intervention (KIM-41) | %.1f | %.1f |", iv / 3600, iv / 3600 / 100))
        for H in (6300.0, 600.0)                   # KIM-35's stage names
            a = load("arm_H_$(Int(H))", kim35).result
            b = load("arm_I_$(Int(H))", kim35).result
            hq, iq = a.numbers[q], b.numbers[q]
            vr = max(hq.var_theta, 0) + max(iq.var_theta, 0)
            σr = mean((hq.sigma2, iq.sigma2))
            nr = vr > 0 ? max(4, ceil(Int, 9σr / vr)) : 4
            tr = H == 6300.0 ? t1 : c["d11_hybrid_600"].median
            e64 = eff[64]
            for (label, f) in (("", 1), (" with V's own cells", 2))
                bud = 3 * 100 * 20 * nr * f * tr
                out(@sprintf("| reduced budget at %d s, n_cells %d%s (KIM-42) | %.1f | %.1f per width on one 64-core node |",
                             Int(H), nr, label, bud / 3600, bud / 3 / (64 * e64) / 3600))
            end
        end
        out("")
    end
    close(io)
end

function main()
    mkpath(OUT)
    say("KIM-38 cost pilot, stage $STAGE, at $COMMIT, job $JOB, $(Threads.nthreads()) thread(s)")
    STAGE === :cells ? cells() : STAGE === :threads ? threads() : STAGE === :sizing ? sizing() :
    STAGE === :project ? project() : error("unknown stage $STAGE")
    say("done")
end

main()
