# spec/phases/16-recovery.md task 16a.7, check V5: the conditional SMC path update
# against brute force, at fixed θ, for block 3 alone.
#
# Two cases:
# - m0: M0 cut to ptsG alone (translation needs it), over three windows.
# - toy: the phase 4 toy with a small, fast-drained pool. Its transcription
#   constant is rebuilt from ATP, so the constants differ between particles and
#   jump across each rebuild (V5's amplified case, sub-spec §12 2026-09-30). The
#   achieved spreads are measured and reported, not assumed.
#
# The design:
# 1. Rejection: simulate the published model, and keep runs whose transcripts
#    match the data at every window. Each kept run is an independent exact draw
#    from p(X | transcripts).
# 2. Kernel: each chain starts from a different exact draw and runs K sweeps of
#    csmc_sweep, and only its final path is kept. If the kernel leaves the target
#    invariant, those final paths are independent exact draws too. So invariance
#    is tested on independent samples, without assuming mixing.
#    - Metabolite likelihood off: the target is p(X | transcripts). Starts are
#      rejection draws from the first half of the kept set, compared with the
#      second half by a two-sample χ² per functional.
#    - Metabolite likelihood on: the target is p(X | transcripts, metabolites).
#      The reference is the kept set weighted by the metabolite likelihood, with
#      σ inflated from 0.1 until the weighted ESS is at least 1,000. Starts are
#      resampled from it, and expectations are compared within 3 SE.
# The functionals are each window's transcript births and the gene's translation
# count.
#
# Stages: reject <case> <task>, prep <case>, kernel <case> <mode> <task>, merge <case>.
# Regenerate: dev/scripts/csmc_v5_16a7.sh

using InferCell
using Random
using Serialization
using Statistics
using Printf
using Distributions: Chisq, cdf, Normal, logpdf

const STAGE, CASE = ARGS[1], ARGS[2]
const DIR = joinpath(@__DIR__, "csmc_v5_16a7", CASE)
mkpath(DIR)
const T = 3
const N_PART = 5
const K_SWEEPS = 5
# Sizes, overridable for a smoke run (V5_SMOKE=1 shrinks every stage).
const SMOKE = get(ENV, "V5_SMOKE", "0") == "1"
const REJECT_TASKS = SMOKE ? 1 : 100
const RUNS_PER_TASK = SMOKE ? 300 : 1000      # m0 keeps ~9%, the toy ~23%
const KERNEL_TASKS = SMOKE ? 1 : 50
const CHAINS_PER_TASK = SMOKE ? 2 : 40
const SIGMA0 = 0.1

# --- the two cases ---------------------------------------------------------

if CASE == "m0"
    const GENES = [:JCVISYN3A_0779]
    const MODELS = m0_models(genes = GENES)
    build() = build_m0(genes = GENES, tspan = (0.0, 60.0 * T))
    const TRUTH = draw_truth(Xoshiro(1607), MODELS;
                             names = [promoter_param(:JCVISYN3A_0779), :krnadeg, :kcatF_R_ENO, :kcatF_R_FBA],
                             purpose = :calibration).values
    const PANEL = M0_PANEL
    tmap(d) = transcript_map(d, MODELS)
    translation_label = :CoreATranslation_1
elseif CASE == "toy"
    include(joinpath(@__DIR__, "..", "..", "test", "corea_test_models.jl"))
    include(joinpath(@__DIR__, "..", "..", "test", "contribution_test_models.jl"))
    include(joinpath(@__DIR__, "..", "..", "test", "jump_test_models.jl"))
    include(joinpath(@__DIR__, "..", "..", "test", "hybrid_test_models.jl"))
    # A small pool, drained hard by translation, so the rebuilt constant falls
    # by large factors across rebuilds and differs between particles.
    const F = corea_particles_per_mM()
    const ATP0 = 3000 / F
    const MODELS = [ToyPool(kcat = 0.0, atp0 = ATP0),
                    ToyRebuiltExpression(k_tx_max = 2.0, km_tx = ATP0, gamma_m = 0.5,
                                         k_tl = 1.0, cost = 20,
                                         k_tx = 2.0 * ATP0 / (ATP0 + ATP0))]
    build() = build_problem(MODELS; tspan = (0.0, 60.0 * T))
    const TRUTH = Pair{Symbol, Float64}[]
    const PANEL = [:M_atp_c]
    function tmap(d)
        names = InferCell._block_names(MODELS, :jump)
        lab(s) = findfirst(==(s), d.events.labels)
        TranscriptMap([findfirst(==(:toy_mrna), names)], [lab(:ToyRebuiltExpression_1)],
                      [lab(:ToyRebuiltExpression_2)])
    end
    translation_label = :ToyRebuiltExpression_3
else
    error("case must be m0 or toy")
end

function fresh()
    d = build()
    isempty(TRUTH) || set_parameters!(d, MODELS, TRUTH)
    return d
end
const BASE = snapshot(fresh())
const TM = tmap(BASE)
const ROWS = [findfirst(==(s), InferCell._block_names(MODELS, :ode)) for s in PANEL]
const TL = findfirst(==(translation_label), BASE.events.labels)
# The rebuilt transcription constant's slot, for the amplified case's conditions.
const KSLOT = CASE == "toy" ? only(only(BASE.rebuilds).fill_idxs) : 0

# Simulate one run from `seed`; return its windows, its transcripts at each 60m,
# and its panel's latent at each 60m.
function simulate(seed)
    d = restore(BASE)
    Random.seed!(seed)
    record_path!(d)
    tx = zeros(Int, length(TM.states), T)
    lat = zeros(length(ROWS), T)
    k = zeros(T)
    for m in 1:T
        for _ in 1:60
            handshake_step!(d)
        end
        tx[:, m] = [d.jump.u[s] for s in TM.states]
        lat[:, m] = [d.ode.u[i] * d.factor + d.rounding.remainders[i] for i in ROWS]
        KSLOT > 0 && (k[m] = d.jump.p[KSLOT])     # rebuilt at 60m, in force from 60m − 1
    end
    return window_events(recorded_path(d), T), tx, lat, k
end

# The data: one run at seed 16070, metabolites observed lognormally at σ0.
function the_data()
    ws, tx, lat, _ = simulate(16_070)
    rng = Xoshiro(16_071)
    met = max.(lat, 1.0) .* exp.(SIGMA0 .* randn(rng, size(lat)))
    return ws, tx, met
end
const TRUE_WS, DATA_TX, DATA_MET = the_data()

functionals(ws) = vcat([count(==(TM.births[1]), w.reactions) for w in ws],
                       [count(==(TL), w.reactions) for w in ws])
const FNAMES = vcat(["births w$m" for m in 1:T], ["translation w$m" for m in 1:T])

metloglik(lat, σ) = sum(logpdf(Normal(log(max(x, 1.0)), σ), log(y)) for (x, y) in zip(lat, DATA_MET))

if STAGE == "reject"
    task = parse(Int, ARGS[3])
    kept = []
    for r in 1:RUNS_PER_TASK
        seed = 1_607_000_000 + task * RUNS_PER_TASK + r
        ws, tx, lat, k = simulate(seed)
        tx == DATA_TX && push!(kept, (seed = seed, ws = ws, lat = lat, f = functionals(ws), k = k))
    end
    serialize(joinpath(DIR, "reject_$(task).jls"), kept)
    @printf("reject %s task %d: %d of %d runs match\n", CASE, task, length(kept), RUNS_PER_TASK)
elseif STAGE == "prep"
    kept = reduce(vcat, [deserialize(joinpath(DIR, "reject_$(t).jls")) for t in 0:(REJECT_TASKS - 1)])
    n = length(kept)
    half = n ÷ 2
    # The smallest σ, doubling from σ0, whose weighted ESS over the whole kept set
    # is at least 1,000.
    σ, ess = SIGMA0, 0.0
    while true
        global σ, ess
        lw = [metloglik(k.lat, σ) for k in kept]
        w = exp.(lw .- maximum(lw))
        ess = sum(w)^2 / sum(w .^ 2)
        ess >= (SMOKE ? 5 : 1000) && break
        σ *= 2
        σ > 100 && error("no σ up to 100 gives an ESS of 1,000")
    end
    serialize(joinpath(DIR, "prep.jls"), (n = n, half = half, sigma = σ, ess = ess))
    if KSLOT > 0
        # The amplified case's two conditions, measured on runs that match the
        # data, which are what the particles sample. The constant before the
        # first rebuild is the build's.
        k0 = BASE.jump.p[KSLOT]
        ks = reduce(hcat, [vcat(k0, kk.k) for kk in kept])          # (T + 1) × runs
        q(x, p) = quantile(x, p)
        for m in 1:T
            @printf("  constant from %d s: 90th/10th percentile across matching runs %.2f; max/min %.2f\n",
                    60m - 1, q(ks[m + 1, :], 0.9) / q(ks[m + 1, :], 0.1),
                    maximum(ks[m + 1, :]) / max(minimum(ks[m + 1, :]), 1e-300))
            jumps = ks[m, :] ./ max.(ks[m + 1, :], 1e-300)
            @printf("    each run's drop across that rebuild: median %.2f×, 10th percentile %.2f×\n",
                    median(jumps), q(jumps, 0.1))
        end
        @printf("  rates: the first constant is %.3g /s\n", k0)
    end
    @printf("prep %s: %d kept runs (acceptance %.4f); metabolite-on σ = %.3g (ESS %.0f)\n",
            CASE, n, n / (REJECT_TASKS * RUNS_PER_TASK), σ, ess)
elseif STAGE == "kernel"
    mode, task = ARGS[3], parse(Int, ARGS[4])
    kept = reduce(vcat, [deserialize(joinpath(DIR, "reject_$(t).jls")) for t in 0:(REJECT_TASKS - 1)])
    pr = deserialize(joinpath(DIR, "prep.jls"))
    rng = Xoshiro(16_072_000 + 1000 * (mode == "on") + task)
    if mode == "off"
        spec = CSMCSpec(TM, Int[]; sigma = 1.0, cap_floor = bridge_cap(maximum(DATA_TX)))
        pool = kept[1:pr.half]
        starts = [pool[mod1(task * CHAINS_PER_TASK + c, length(pool))].ws for c in 1:CHAINS_PER_TASK]
    else
        spec = CSMCSpec(TM, ROWS; sigma = pr.sigma, cap_floor = bridge_cap(maximum(DATA_TX)))
        lw = [metloglik(k.lat, pr.sigma) for k in kept]
        w = exp.(lw .- maximum(lw))
        w ./= sum(w)
        cw = cumsum(w)
        starts = [kept[min(searchsortedfirst(cw, rand(rng)), length(kept))].ws for _ in 1:CHAINS_PER_TASK]
    end
    data = CellData(DATA_TX, DATA_MET)
    out = map(starts) do ref
        changed = falses(T)
        for _ in 1:K_SWEEPS
            ref, ch = csmc_sweep(rng, BASE, spec, data, ref; N = N_PART)
            changed .|= ch
        end
        (f = functionals(ref), changed = collect(changed))
    end
    serialize(joinpath(DIR, "kernel_$(mode)_$(task).jls"), out)
    @printf("kernel %s %s task %d: %d chains\n", CASE, mode, task, length(out))
elseif STAGE == "merge"
    kept = reduce(vcat, [deserialize(joinpath(DIR, "reject_$(t).jls")) for t in 0:(REJECT_TASKS - 1)])
    pr = deserialize(joinpath(DIR, "prep.jls"))
    @printf("V5, case %s: %d kept runs from %d; N = %d particles, %d sweeps per chain from an exact start\n",
            CASE, pr.n, REJECT_TASKS * RUNS_PER_TASK, N_PART, K_SWEEPS)
    # Metabolite off: two-sample χ² per functional against the second half.
    off = reduce(vcat, [deserialize(joinpath(DIR, "kernel_off_$(t).jls")) for t in 0:(KERNEL_TASKS - 1)])
    refset = kept[(pr.half + 1):end]
    ok = true
    @printf("metabolite likelihood off: %d kernel draws against %d rejection draws; %d χ² tests at p > 0.01\n",
            length(off), length(refset), length(FNAMES))
    for (j, name) in enumerate(FNAMES)
        a = [o.f[j] for o in off]
        b = [k.f[j] for k in refset]
        vals = sort(unique(vcat(a, b)))
        ca = [count(==(v), a) for v in vals]
        cb = [count(==(v), b) for v in vals]
        # Pool sparse categories until each expects at least 5 in both samples.
        na, nb = length(a), length(b)
        tab = Tuple{Int, Int}[]
        oa, ob = 0, 0
        for (x, y) in zip(ca, cb)
            oa += x
            ob += y
            if (oa + ob) * min(na, nb) / (na + nb) >= 5
                push!(tab, (oa, ob))
                oa, ob = 0, 0
            end
        end
        if !isempty(tab)
            tab[end] = (tab[end][1] + oa, tab[end][2] + ob)
        end
        stat = 0.0
        for (x, y) in tab
            e1 = (x + y) * na / (na + nb)
            e2 = (x + y) * nb / (na + nb)
            stat += (x - e1)^2 / e1 + (y - e2)^2 / e2
        end
        p = length(tab) >= 2 ? 1 - cdf(Chisq(length(tab) - 1), stat) : 1.0
        global ok &= p > 0.01
        @printf("  %s: kernel mean %.3f, rejection mean %.3f; χ² %.2f on %d df, p = %.3f%s\n",
                name, mean(a), mean(b), stat, length(tab) - 1, p, p > 0.01 ? "" : "  FAIL")
    end
    # Metabolite on: weighted reference expectations within 3 SE.
    on = reduce(vcat, [deserialize(joinpath(DIR, "kernel_on_$(t).jls")) for t in 0:(KERNEL_TASKS - 1)])
    lw = [metloglik(k.lat, pr.sigma) for k in kept]
    w = exp.(lw .- maximum(lw))
    w ./= sum(w)
    ess = 1 / sum(w .^ 2)
    @printf("metabolite likelihood on, σ = %.3g, reference ESS %.0f: %d kernel draws\n", pr.sigma, ess, length(on))
    for (j, name) in enumerate(FNAMES)
        x = [k.f[j] for k in kept]
        m_ref = sum(w .* x)
        se_ref = sqrt(sum(w .* (x .- m_ref) .^ 2) / ess)
        a = [o.f[j] for o in on]
        m_k, se_k = mean(a), std(a) / sqrt(length(a))
        z = (m_k - m_ref) / sqrt(se_ref^2 + se_k^2)
        global ok &= abs(z) <= 3
        @printf("  %s: kernel %.3f ± %.3f, reference %.3f ± %.3f; z = %.2f%s\n",
                name, m_k, se_k, m_ref, se_ref, z, abs(z) <= 3 ? "" : "  FAIL")
    end
    upd = mean(reduce(hcat, [o.changed for o in off]); dims = 2)
    @printf("windows changed within %d sweeps (metabolite off): %s\n", K_SWEEPS,
            join([@sprintf("%.2f", x) for x in upd], ", "))
    println(ok ? "V5 ($CASE) passes" : "V5 ($CASE) FAILS")
end
