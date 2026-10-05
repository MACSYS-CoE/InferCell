# spec/phases/16-recovery.md task 16a.7, the variant choice: per-window update
# rates and cost per sweep of the conditional SMC path update, for particle Gibbs,
# ancestor sampling truncated at lag 1 and 5, and full ancestor sampling, at
# N ∈ {5, 10, 20, 50}, on M0 and on one full-scale cell (D16.3); and, after KP
# fired (§12 2026-10-02), annealed windows at N ∈ {2, 5} on M0 (16a.7a).
#
# The cells:
# - m0: cells 1 to 5 of M0 replicate 16085, the replicate block 1's V4 used, at
#   the replicate's truth and σ. Each cell's path is re-recorded from its seed
#   and checked against the dataset's latent.
# - full: Core A′ cell 30000 + c at 15.7's truth (`Xoshiro(1507)`), its path
#   recorded; the variant runs use c = 1, cell 30001.
#   The 15.7 observations are not tracked, so the panel is observed afresh,
#   lognormal at 15.7's σ = 0.1 with a one-particle floor, from `Xoshiro(1508)`
#   on cell 30001 alone, and from `Xoshiro(1508 + 1000(c − 1))` on the others. It
#   is drawn as 15.7's was, but it is not 15.7's draw.
#
# Every chain starts from the cell's true path at the true θ. That path is an
# exact draw from p(X | data, θ), so the kernel is stationary from the first
# sweep and no burn-in is dropped from the update rates.
#
# A window's update rate is the fraction of sweeps in which its events changed.
# Cost is wall-clock per sweep, and per cycle of plain simulation timed in the
# same job, so cycles per sweep is machine-independent. Each task's first sweep
# compiles, so its time is dropped from the cost.
#
# Stages: run <scale> <variant> <N> <task> <sweeps>, and merge.
# Variants: pg (lag 0), l1, l5, full (lag T), and a20, a50, a100 (16a.7a: particle
# Gibbs with each window annealed over K stages).
# Regenerate: dev/scripts/csmc_variants_16a7.sh, from a worktree at a pushed commit.
# CSMC_SMOKE=1 shrinks both cells to three windows.

using InferCell
using Random
using Serialization
using Statistics
using Printf

const SMOKE = get(ENV, "CSMC_SMOKE", "0") == "1"
const DIR = joinpath(@__DIR__, "csmc_variants_16a7", SMOKE ? "smoke" : "")
mkpath(DIR)
const M0_REPLICATE = 16085
const M0_CELLS = 5
const FULL_TRUTH_SEED = 1507
const FULL_CELL = 30_001
const FULL_SIGMA = 0.1
const FULL_NOISE_SEED = 1508
const THRESHOLD = 0.10            # D16.3's update-rate threshold, ⚠️ DRAFT
# The annealed variants, aK, are particle Gibbs with each window annealed over K
# geometric stages from β1 = 1e-4 (§12 2026-10-02, task 16a.7a).
const VARIANTS = ["pg", "l1", "l5", "full", "a20", "a50", "a100"]
const NS = [5, 10, 20, 50, 2]
const BETA1 = 1e-4
# CSMC_FLOOR=F re-observes the panel with its floor at F particles instead of
# one, keeping each observation's own noise draw: y_F = max(x, F) · y / max(x, 1).
# The likelihood's floor moves with it. The floor a cited single-cell detection
# limit implies is about 1.2e5 molecules (0.2 amol of NAD⁺; Lin et al., Anal.
# Chem. 2011), which is option 1 of the KP decision (handoff 2026-10-02).
const OBS_FLOOR = parse(Float64, get(ENV, "CSMC_FLOOR", "1.0"))

schedule(variant) = startswith(variant, "a") ?
    geometric_schedule(parse(Int, variant[2:end]); β1 = BETA1) : nothing
lag(variant, T) = startswith(variant, "a") || variant == "pg" ? 0 : variant == "l1" ? 1 : variant == "l5" ? 5 :
                  variant == "full" ? T : error("unknown variant $variant")

# Record a run from `seed` window by window: returns its path, its transcripts in
# the map's order and its panel latent, both windows as columns.
function record(d, seed, tm, rows, T)
    Random.seed!(seed)
    record_path!(d)
    tx = zeros(Int, length(tm.states), T)
    lat = zeros(length(rows), T)
    for m in 1:T
        for _ in 1:60
            handshake_step!(d)
        end
        tx[:, m] = [d.jump.u[s] for s in tm.states]
        lat[:, m] = [d.ode.u[i] * d.factor + d.rounding.remainders[i] for i in rows]
    end
    return recorded_path(d), tx, lat
end

# One cell: its base driver at t = 0 and the true θ, the sweep's spec, its data
# and its true path by window.
function cell(scale, c)
    if scale == "m0"
        horizon = SMOKE ? 180.0 : M0_HORIZON_S
        ms = m0_models()
        ds = m0_dataset(M0_REPLICATE; ncells = M0_CELLS, horizon)
        fresh = () -> (d = build_m0(tspan = (0.0, horizon)); set_parameters!(d, ms, ds.truth.values); d)
        seed, σ = ds.seeds[c], ds.sigma
    else
        horizon = SMOKE ? 180.0 : COREA_CYCLE_S
        ms = d11_models()
        truth = draw_truth(Xoshiro(FULL_TRUTH_SEED), ms; purpose = :recovery)
        fresh = () -> (d = build_problem(ms; tspan = (0.0, horizon), complete = true);
                       set_parameters!(d, ms, truth.values); d)
        seed, σ = FULL_CELL + c - 1, FULL_SIGMA
    end
    T = round(Int, horizon) ÷ 60
    Random.seed!(seed)
    base = snapshot(fresh())
    tm = transcript_map(base, ms)
    rows = [findfirst(==(s), InferCell._block_names(ms, :ode)) for s in M0_PANEL]
    path, tx, lat = record(restore(base), seed, tm, rows, T)
    if scale == "m0"
        lrows = [findfirst(==(s), ds.species) for s in M0_PANEL]
        orows = [findfirst(==(s), ds.observed_species) for s in M0_PANEL]
        lat == ds.latent[c, lrows, :] || error("cell $seed: the recorded run is not the dataset's")
        met = ds.observed[c, orows, :]
    else
        noise = Xoshiro(FULL_NOISE_SEED + 1000 * (c - 1))      # cell 30001 keeps 1508
        met = max.(lat, 1.0) .* exp.(FULL_SIGMA .* randn(noise, size(lat)))
    end
    OBS_FLOOR == 1 || (met = max.(lat, OBS_FLOOR) .* met ./ max.(lat, 1.0))
    spec = CSMCSpec(tm, rows; sigma = σ, floor = OBS_FLOOR, cap_floor = bridge_cap(maximum(tx)))
    return (; base, spec, data = CellData(tx, met), ref = window_events(path, T), T, seed, σ, ms)
end

# Plain simulation of the whole cycle, timed warm: the cost unit.
function cycle_seconds(base, seed, T)
    run() = (d = restore(base); Random.seed!(seed); for _ in 1:(60T); handshake_step!(d); end)
    run()
    return @elapsed run()
end

function run_task(scale, variant, N, task, sweeps)
    c = scale == "m0" ? task : 1
    x = cell(scale, c)
    L = lag(variant, x.T)
    out = joinpath(DIR, "$(scale)_$(variant)_N$(N)_t$(task).jls")
    cyc = cycle_seconds(x.base, x.seed, x.T)
    @printf("%s %s N = %d task %d: cell %d, T = %d, lag %d, σ = %.3g; a plain cycle takes %.2f s\n",
            scale, variant, N, task, x.seed, x.T, L, x.σ, cyc)
    rng = Xoshiro(16_074_000_000 + 100_000 * (scale == "full") + 10_000 * findfirst(==(variant), VARIANTS) +
                  100 * findfirst(==(N), NS) + task)
    ref = x.ref
    changed = falses(x.T, 0)
    secs = Float64[]
    for s in 1:sweeps
        t = @elapsed begin
            ref, ch = csmc_sweep(rng, x.base, x.spec, x.data, ref; N, lag = L,
                                 schedule = schedule(variant))
        end
        changed = hcat(changed, ch)
        push!(secs, t)
        # Saved after every sweep, so a job that runs out of time still counts.
        serialize(out, (; scale, variant, N, task, seed = x.seed, T = x.T, lag = L,
                        sigma = x.σ, cycle = cyc, changed, secs))
        @printf("  sweep %d: %.1f s (%.1f cycles); %d of %d windows changed\n",
                s, t, t / cyc, count(ch), x.T)
        flush(stdout)
    end
end

function merge()
    files = filter(endswith(".jls"), readdir(DIR; join = true))
    rs = [deserialize(f) for f in files]
    for scale in ["m0", "full"]
        println("\n## $scale\n")
        println("| Variant | N | Sweeps | Cycles per sweep | Min window rate | Windows below $(THRESHOLD) | Mean rate, windows 1–10 | Mean rate, rest |")
        println("|---|---|---|---|---|---|---|---|")
        for v in VARIANTS, N in NS
            g = [r for r in rs if r.scale == scale && r.variant == v && r.N == N]
            isempty(g) && continue
            ch = reduce(hcat, [r.changed for r in g])
            S = size(ch, 2)
            S == 0 && continue
            rate = vec(mean(ch; dims = 2))
            # Each task's first sweep compiles; it is dropped from the cost.
            cps = [s / r.cycle for r in g for s in r.secs[2:end]]
            cost = isempty(cps) ? "—" : @sprintf("%.1f", mean(cps))
            early = mean(rate[1:min(10, end)])
            rest = length(rate) > 10 ? @sprintf("%.2f", mean(rate[11:end])) : "—"
            @printf("| %s | %d | %d | %s | %.2f (window %d) | %d of %d | %.2f | %s |\n",
                    v, N, S, cost, minimum(rate), argmin(rate), count(<(THRESHOLD), rate),
                    length(rate), early, rest)
        end
    end
    println("\nPer-window rates are in each task's file under $(DIR).")
end

if abspath(PROGRAM_FILE) == @__FILE__
    if ARGS[1] == "run"
        run_task(ARGS[2], ARGS[3], parse(Int, ARGS[4]), parse(Int, ARGS[5]), parse(Int, ARGS[6]))
    elseif ARGS[1] == "merge"
        merge()
    else
        error("stage must be run or merge")
    end
end
