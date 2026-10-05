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
# V5_FRESH=1: a second, independent rejection reference, from seeds offset by
# 50,000,000, in its own directory (task 16a.7c's watch item on M0's window-2
# births, which every comparison so far made against one reference set). The
# `refcheck` stage tests earlier kernel draws, which do not depend on the
# reference, against it.
const FRESH = get(ENV, "V5_FRESH", "0") == "1"
const DIR = joinpath(@__DIR__, "csmc_v5_16a7", FRESH ? CASE * "_fresh" : CASE)
const SEED_OFFSET = FRESH ? 50_000_000 : 0
mkpath(DIR)
const T = 3
const N_PART = 5
const K_SWEEPS = 5
# V5_ANNEAL=K runs the kernel with each window annealed over K geometric stages
# from β1 = 1e-4 (task 16a.7a, §12 2026-10-02). Its kernel files carry the K, and
# it reuses the unannealed run's rejection draws, which no kernel changes.
const ANNEAL = parse(Int, get(ENV, "V5_ANNEAL", "0"))
const SCHEDULE = ANNEAL > 0 ? geometric_schedule(ANNEAL; β1 = 1e-4) : nothing
const KTAG = ANNEAL > 0 ? "_a$(ANNEAL)" : ""
# Sizes, overridable for a smoke run (V5_SMOKE=1 shrinks every stage).
const SMOKE = get(ENV, "V5_SMOKE", "0") == "1"
const REJECT_TASKS = SMOKE ? 1 : 100
# m0 keeps about 9% of runs, the toy about 3%, so the toy runs more to keep its
# start pool larger than the 2,000 chains.
const RUNS_PER_TASK = SMOKE ? 300 : (CASE == "m0" ? 1000 : 5000)
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
    # ATP decays first order, at a rate set by the live protein count through the
    # catalytic edge, so it never empties. Transcription's constant, nearly
    # linear in ATP since km_tx is far above it, falls about 2.5× a minute, and
    # differs between particles because their protein counts do. A drain by
    # translation cost emptied the pool in the first smoke run (job 17853650).
    const F = corea_particles_per_mM()
    const ATP0 = 2000 / F
    const KM_TX = 100 * ATP0
    const KMAX = 1.0 * (KM_TX + ATP0) / ATP0          # k = 1 /s at the start
    const MODELS = [ToyPool(kcat = 0.015 * F / 100, km = 1.0, atp0 = ATP0),
                    ToyRebuiltExpression(k_tx_max = KMAX, km_tx = KM_TX, gamma_m = 0.5,
                                         k_tl = 1.0, cost = 0, protein0 = 10,
                                         k_tx = KMAX * ATP0 / (KM_TX + ATP0))]
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

# --- V5's bulk "on" half (task 16a.7c, §12 2026-10-05) ----------------------
#
# A population of three cells, each with the data's transcripts, so each cell's
# exact draws are the kept runs and draws of different cells are independent.
# The data's three true cells are kept runs drawn from a fixed seed, and the
# bulk observation is their mean latent, lognormal at σ0. The reference is
# random triples of kept runs (the second half), weighted by the bulk likelihood
# of their mean, with σ_b doubled from σ0 until the weighted ESS is at least
# 1,000. The kernel scans the three cells in turn, K_SWEEPS times, each
# conditioned on the other two, from a triple resampled from the weighted first
# half. Each cell's functionals are averaged over the three cells, which are
# exchangeable. Expectations must agree within 3 SE.
const BULK_C = 3
const BULK_TRIPLES = SMOKE ? 2_000 : 200_000

bulkloglik(lats, z, σ) = sum(logpdf(Normal(log(max(x, 1.0)), σ), log(y))
                             for (x, y) in zip(sum(lats) ./ BULK_C, z))

# The bulk reference's weighted mean of functional `j` over `triples` of runs,
# and its standard error by a bootstrap over the runs themselves: the triples
# reuse a finite pool of runs, so they are not independent, and an ESS-based SE
# understates the error. (The first bulk merge, jobs 18052539 and 18052542, used
# the ESS-based SE; the fresh-reference check of 2026-10-05 found it too small.
# The merge now uses this, on a quarter of the prep's triples.)
function bulk_reference(kept, triples, z, σ, j; B = 100, rng = Xoshiro(16_093))
    fx(map) = begin
        l = [bulkloglik([kept[map[i]].lat for i in t], z, σ) for t in triples]
        w = exp.(l .- maximum(l))
        sum(w .* [mean(kept[map[i]].f[j] for i in t) for t in triples]) / sum(w)
    end
    m = fx(1:length(kept))
    boots = [fx(rand(rng, 1:length(kept), length(kept))) for _ in 1:B]
    return m, std(boots)
end

function bulk_stage(stage, args)
    kept = reduce(vcat, [deserialize(joinpath(DIR, "reject_$(t).jls")) for t in 0:(REJECT_TASKS - 1)])
    n = length(kept)
    half = n ÷ 2
    if stage == "bulkprep"
        rng = Xoshiro(16_090)
        truth = [kept[rand(rng, 1:n)] for _ in 1:BULK_C]
        z = sum(t.lat for t in truth) ./ BULK_C
        z = max.(z, 1.0) .* exp.(SIGMA0 .* randn(rng, size(z)))
        second = (half + 1):n
        triples = [rand(rng, second, BULK_C) for _ in 1:BULK_TRIPLES]
        σ, ess = SIGMA0, 0.0
        while true
            lw = [bulkloglik([kept[i].lat for i in t], z, σ) for t in triples]
            w = exp.(lw .- maximum(lw))
            ess = sum(w)^2 / sum(w .^ 2)
            ess >= (SMOKE ? 5 : 1000) && break
            σ *= 2
            σ > 100 && error("no σ_b up to 100 gives an ESS of 1,000")
        end
        serialize(joinpath(DIR, "bulkprep.jls"), (; z, sigma = σ, ess, triples))
        @printf("bulk prep %s: %d triples from %d second-half runs; σ_b = %.3g (ESS %.0f)\n",
                CASE, length(triples), n - half, σ, ess)
    elseif stage == "bulkkernel"
        task = parse(Int, args[3])
        pr = deserialize(joinpath(DIR, "bulkprep.jls"))
        rng = Xoshiro(16_091_000 + task)
        spec = CSMCSpec(TM, ROWS; sigma = pr.sigma, cap_floor = bridge_cap(maximum(DATA_TX)))
        data = CellData(DATA_TX, pr.z)
        # Exact starts: triples from the first half, resampled by weight.
        first_half = 1:half
        starts_pool = [rand(rng, first_half, BULK_C) for _ in 1:(10 * CHAINS_PER_TASK * 50)]
        lw = [bulkloglik([kept[i].lat for i in t], pr.z, pr.sigma) for t in starts_pool]
        w = exp.(lw .- maximum(lw))
        cw = cumsum(w ./ sum(w))
        out = map(1:CHAINS_PER_TASK) do _
            t = starts_pool[min(searchsortedfirst(cw, rand(rng)), length(cw))]
            refs = [kept[i].ws for i in t]
            lats = [kept[i].lat for i in t]
            for _ in 1:K_SWEEPS, c in 1:BULK_C
                others = sum(lats[j] for j in 1:BULK_C if j != c)
                bulk = bulk_window_observations(pr.z, others, BULK_C)
                refs[c], _, lats[c] = csmc_sweep(rng, BASE, spec, data, refs[c]; N = N_PART, bulk)
            end
            (f = mean(functionals.(refs)),)
        end
        serialize(joinpath(DIR, "bulkkernel_$(task).jls"), out)
        @printf("bulk kernel %s task %d: %d chains\n", CASE, task, length(out))
    else
        pr = deserialize(joinpath(DIR, "bulkprep.jls"))
        on = reduce(vcat, [deserialize(joinpath(DIR, "bulkkernel_$(t).jls")) for t in 0:(KERNEL_TASKS - 1)])
        lw = [bulkloglik([kept[i].lat for i in t], pr.z, pr.sigma) for t in pr.triples]
        w = exp.(lw .- maximum(lw))
        w ./= sum(w)
        ess = 1 / sum(w .^ 2)
        triples = pr.triples[1:(length(pr.triples) ÷ 4)]
        @printf("V5 bulk, case %s: %d cells; σ_b = %.3g, reference ESS %.0f over %d triples; %d kernel draws, N = %d, %d scans; reference SE by a bootstrap over runs on %d triples\n",
                CASE, BULK_C, pr.sigma, ess, length(pr.triples), length(on), N_PART, K_SWEEPS, length(triples))
        ok = true
        for (j, name) in enumerate(FNAMES)
            m_ref, se_ref = bulk_reference(kept, triples, pr.z, pr.sigma, j)
            a = [o.f[j] for o in on]
            m_k, se_k = mean(a), std(a) / sqrt(length(a))
            den = sqrt(se_ref^2 + se_k^2)
            z = den > 0 ? (m_k - m_ref) / den : (m_k == m_ref ? 0.0 : Inf)
            ok &= abs(z) <= 3
            @printf("  %s (mean over cells): kernel %.3f ± %.3f, reference %.3f ± %.3f; z = %.2f%s\n",
                    name, m_k, se_k, m_ref, se_ref, z, abs(z) <= 3 ? "" : "  FAIL")
        end
        println(ok ? "V5 bulk ($CASE) passes" : "V5 bulk ($CASE) FAILS")
    end
end

# The fresh reference against the kernel draws in `olddir` (an earlier run's
# csmc_v5_16a7/<case>): the off half's χ², the per-cell on half's and the bulk
# half's z, each at the σ that run's prep chose. Every functional is reported,
# with window 2's births the one in question.
function chi2_p(a, b)
    vals = sort(unique(vcat(a, b)))
    ca = [count(==(v), a) for v in vals]
    cb = [count(==(v), b) for v in vals]
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
    isempty(tab) || (tab[end] = (tab[end][1] + oa, tab[end][2] + ob))
    stat = 0.0
    for (x, y) in tab
        e1, e2 = (x + y) * na / (na + nb), (x + y) * nb / (na + nb)
        stat += (x - e1)^2 / e1 + (y - e2)^2 / e2
    end
    return length(tab) >= 2 ? 1 - cdf(Chisq(length(tab) - 1), stat) : 1.0
end

function refcheck(olddir)
    kept = reduce(vcat, [deserialize(joinpath(DIR, "reject_$(t).jls")) for t in 0:(REJECT_TASKS - 1)])
    @printf("V5 reference check, case %s: %d fresh kept runs from %d (seed offset %d)\n",
            CASE, length(kept), REJECT_TASKS * RUNS_PER_TASK, SEED_OFFSET)
    # V5_POOL=<dir>: pool another rejection set (an earlier run's) into the
    # reference, which is still exact, and halves its variance.
    pool = get(ENV, "V5_POOL", "")
    if !isempty(pool)
        extra = reduce(vcat, [deserialize(joinpath(pool, "reject_$(t).jls")) for t in 0:(REJECT_TASKS - 1)])
        kept = vcat(kept, extra)
        @printf("pooled with %d runs from %s: %d in all\n", length(extra), pool, length(kept))
    end
    zfun(a, x, w) = begin
        ess = 1 / sum(w .^ 2)
        m = sum(w .* x)
        se = sqrt(sum(w .* (x .- m) .^ 2) / ess)
        (mean(a) - m) / sqrt(se^2 + var(a) / length(a)), m
    end
    for tag in ("", "_a5")
        f0 = joinpath(olddir, "kernel_off$(tag)_0.jls")
        isfile(f0) || continue
        off = reduce(vcat, [deserialize(joinpath(olddir, "kernel_off$(tag)_$(t).jls")) for t in 0:(KERNEL_TASKS - 1)])
        on = reduce(vcat, [deserialize(joinpath(olddir, "kernel_on$(tag)_$(t).jls")) for t in 0:(KERNEL_TASKS - 1)])
        σ = deserialize(joinpath(olddir, "prep.jls")).sigma
        lw = [metloglik(k.lat, σ) for k in kept]
        w = exp.(lw .- maximum(lw))
        w ./= sum(w)
        @printf("kernel%s: off, %d draws against %d fresh runs; on at σ = %.3g, fresh ESS %.0f\n",
                tag == "" ? "" : " (annealed, K = 5)", length(off), length(kept), σ, 1 / sum(w .^ 2))
        for (j, name) in enumerate(FNAMES)
            a = [o.f[j] for o in off]
            b = [k.f[j] for k in kept]
            z, m = zfun([o.f[j] for o in on], [k.f[j] for k in kept], w)
            @printf("  %s: off kernel %.4f, fresh %.4f, χ² p = %.3f; on kernel %.4f, fresh %.4f, z = %.2f\n",
                    name, mean(a), mean(b), chi2_p(a, b), mean(o.f[j] for o in on), m, z)
        end
    end
    bp = joinpath(olddir, "bulkprep.jls")
    if isfile(bp)
        pr = deserialize(bp)
        bk = reduce(vcat, [deserialize(joinpath(olddir, "bulkkernel_$(t).jls")) for t in 0:(KERNEL_TASKS - 1)])
        rng = Xoshiro(16_092)
        triples = [rand(rng, 1:length(kept), BULK_C) for _ in 1:(BULK_TRIPLES ÷ 4)]
        @printf("bulk: %d kernel draws; %d triples over %d runs, σ_b = %.3g; reference SE by a bootstrap over runs\n",
                length(bk), length(triples), length(kept), pr.sigma)
        for (j, name) in enumerate(FNAMES)
            m, se = bulk_reference(kept, triples, pr.z, pr.sigma, j)
            a = [o.f[j] for o in bk]
            z = (mean(a) - m) / sqrt(se^2 + var(a) / length(a))
            @printf("  %s (mean over cells): kernel %.4f ± %.4f, reference %.4f ± %.4f, z = %.2f\n",
                    name, mean(a), std(a) / sqrt(length(a)), m, se, z)
        end
    end
end

if STAGE == "reject"
    task = parse(Int, ARGS[3])
    kept = []
    for r in 1:RUNS_PER_TASK
        seed = 1_607_000_000 + SEED_OFFSET + task * RUNS_PER_TASK + r
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
            any(<=(0), ks[m + 1, :]) && @printf("  %d matching runs have a zero constant from %d s\n",
                                                 count(<=(0), ks[m + 1, :]), 60m - 1)
            @printf("  constant from %d s: 90th/10th percentile across matching runs %.2f; max/min %.2f\n",
                    60m - 1, q(ks[m + 1, :], 0.9) / q(ks[m + 1, :], 0.1),
                    maximum(ks[m + 1, :]) / minimum(ks[m + 1, :]))
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
    rng = Xoshiro(16_072_000 + 1_000_000 * ANNEAL + 1000 * (mode == "on") + task)
    if mode == "off"
        spec = CSMCSpec(TM, Int[]; sigma = 1.0, cap_floor = bridge_cap(maximum(DATA_TX)))
        pool = kept[1:pr.half]
        length(pool) >= KERNEL_TASKS * CHAINS_PER_TASK ||
            @warn "the start pool ($(length(pool))) is smaller than the chain count; starts repeat"
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
            ref, ch = csmc_sweep(rng, BASE, spec, data, ref; N = N_PART, schedule = SCHEDULE)
            changed .|= ch
        end
        (f = functionals(ref), changed = collect(changed))
    end
    serialize(joinpath(DIR, "kernel_$(mode)$(KTAG)_$(task).jls"), out)
    @printf("kernel %s %s task %d: %d chains\n", CASE, mode, task, length(out))
elseif STAGE == "merge"
    kept = reduce(vcat, [deserialize(joinpath(DIR, "reject_$(t).jls")) for t in 0:(REJECT_TASKS - 1)])
    pr = deserialize(joinpath(DIR, "prep.jls"))
    @printf("V5, case %s: %d kept runs from %d; N = %d particles, %d sweeps per chain from an exact start\n",
            CASE, pr.n, REJECT_TASKS * RUNS_PER_TASK, N_PART, K_SWEEPS)
    ANNEAL > 0 && @printf("each window annealed over %d stages from β1 = 1e-4\n", ANNEAL)
    # Metabolite off: two-sample χ² per functional against the second half.
    off = reduce(vcat, [deserialize(joinpath(DIR, "kernel_off$(KTAG)_$(t).jls")) for t in 0:(KERNEL_TASKS - 1)])
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
    on = reduce(vcat, [deserialize(joinpath(DIR, "kernel_on$(KTAG)_$(t).jls")) for t in 0:(KERNEL_TASKS - 1)])
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
        den = sqrt(se_ref^2 + se_k^2)
        # A functional with no spread on either side passes only if the means agree.
        z = den > 0 ? (m_k - m_ref) / den : (m_k == m_ref ? 0.0 : Inf)
        global ok &= abs(z) <= 3
        @printf("  %s: kernel %.3f ± %.3f, reference %.3f ± %.3f; z = %.2f%s\n",
                name, m_k, se_k, m_ref, se_ref, z, abs(z) <= 3 ? "" : "  FAIL")
    end
    upd = mean(reduce(hcat, [o.changed for o in off]); dims = 2)
    @printf("windows changed within %d sweeps (metabolite off): %s\n", K_SWEEPS,
            join([@sprintf("%.2f", x) for x in upd], ", "))
    println(ok ? "V5 ($CASE) passes" : "V5 ($CASE) FAILS")
elseif STAGE in ("bulkprep", "bulkkernel", "bulkmerge")
    bulk_stage(STAGE, ARGS)
elseif STAGE == "refcheck"
    refcheck(ARGS[3])
end
