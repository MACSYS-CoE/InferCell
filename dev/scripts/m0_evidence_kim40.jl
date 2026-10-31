# KIM-40 and KIM-41: M0's evidence run (M5). The crossed comparison (KIM-40) and
# the ptsG panel (KIM-41) run in the same calls, on the M4 freeze's design
# (`93abbce`), at M0's scale (spec §12 2026-10-31, KIM-69, `5ec56d3`), and are
# scored as KIM-40's §12 entry declares before the run (below).
#
# Usage (dev/scripts/m0_evidence_kim40.sh submits every stage):
#   julia --project dev/scripts/m0_evidence_kim40.jl <stage> <outdir> [<index>]
# Stages:
#   prep       checks the 50 V7 replicates and the 50 B7 fits (names, truths,
#              draw counts) and writes the manifest: each posterior's draws and
#              blocks, and the SHA-256 of the 200 input files (B7's checked
#              against its archive's SHA256SUMS), which `score` checks again. It
#              reads draws and truths, never Q.
#   work i     one single-core task, i = 0…1,149:
#              - 0…99: δ for (dataset r, posterior p);
#              - 100…149: the truth for dataset r;
#              - 150…1,149: one hybrid block (r, p, b) of 10 draws.
#   finish i   i = 0…99, one (r, p): its blocks combined (`combine_rows`), the
#              mean field on the same 100 draws with the shared δ, and the checks.
#   score      coverage and sharpness over the datasets, for Q (KIM-40) and ΔQ
#              (KIM-41).
# Outputs, never overwritten: work/ (whole objects, which read back only at this
# commit), work_plain/ and finish/ (plain values, readable at any commit),
# score.md, score.jls and score_datasets.tsv. A file is written under tmp/ and
# moved into place, so a killed task leaves no output.
#
# The design (frozen: the M4 freeze and KIM-69):
# - Q = :lactate_out over M0_HORIZON_S = 600 s; the estimand is the population
#   mean μ(θ). completion = :none, with no viability rule (the freeze: V applies
#   on the reduced model only);
# - every call carries the KIM-41 knockdown, Knockdown(; dose = 0.3) on ptsG:
#   the hybrid's, the mean field's, δ's and the truth's;
# - per dataset r = 1…50 and posterior p, V7 replicate r (exact) or B7 fit r
#   (mixed, its two chains stacked): 100 draws, evenly spaced,
#   round.(Int, range(1, N; length = 100)) with N = 556 (V7, after its burn-in of
#   100) or 400 (B7); 10 blocks of 10 consecutive draws;
# - the hybrid: 887 cells per draw. The mean field: one solve per draw, on the
#   hybrid's draws (paired draws), uncorrected (`mu`) and δ-corrected;
# - δ: one DeltaCorrection per (r, p), 9,849 cells at θ_ref, the posterior mean of
#   log θ (`predict`'s default), computed as `predict`'s own `delta = :ensemble`
#   would, and passed to the mean field's call;
# - the truth μ(θ*) and Δμ(θ*): `predict` on θ* alone, 9,849 cells, so it runs
#   the predictions' own code path. θ* is V7 replicate r's stored truth, which
#   `prep` checks against B7 fit r's to 1e-12 relative (V7 stores it on the log
#   scale, B7 on the natural one).
#
# Seeds (KIM-40's, from M5's 19_500_000–19_599_999; KIM-41 shares them; KIM-42
# keeps 19_542_xxx):
# - block b of dataset r: rng Xoshiro(19_540_000 + 10(r − 1) + (b − 1)), for V7
#   and B7 alike. With explicit `draws` the cell seeds come from the rng alone, so
#   matched V7 and B7 blocks run the same cells (common cell seeds), and a
#   dataset's blocks run disjoint ones;
# - δ for dataset r: cell seeds `_cell_seeds(Xoshiro(19_541_000 + r), 9_849)`,
#   for V7 and B7 alike;
# - the mean field for r: rng Xoshiro(19_541_200 + r). It draws no cell seeds;
# - the truth for r: rng Xoshiro(19_541_500 + r);
# - the smoke pass: 19_541_900 to 19_541_930.
# `finish` checks that a dataset's block, δ and truth seeds are disjoint, and
# `score` that V7's and B7's cell seeds are equal.
#
# Scoring (declared before the run, in KIM-40's §12 entry). Per dataset, the
# predictive sample of each arm is its 100 values, as they are:
# - the arms: {V7, B7} × {hybrid μ̂, mean field m(θ), mean field m(θ) + δ}. For ΔQ
#   they are the hybrid's paired effect, m(g_a θ) − m(θ), and that plus
#   δ_intervened − δ. A draw with no finished cell leaves its sample, and is
#   counted;
# - the truth: μ̂(θ*) and Δμ̂(θ*), from the truth's cells, as they are;
# - **primary.** Coverage: the fraction of datasets whose truth lies in the
#   central 50% and 90% intervals of the sample (`quantile`, type 7), with a
#   binomial SE. Sharpness: the 90% interval's width, as its median over
#   datasets, and each dataset's width over V7-hybrid's (median and IQR);
# - **secondary.** The truth's PIT, the fraction of the sample below it, in 10
#   bins. The standardised error (mean − truth)/SD. Each arm's shift from
#   V7-hybrid, (mean − mean_ref)/SD_ref, the table's effects;
# - **beside them.** The hybrid's MC share mc/var, and the truth's and δ's SEs
#   over the hybrid's √parameter_var, against the 10% rules (KIM-69). For ΔQ
#   each is the SE of its cells' paired effect: for δ, that of δ_intervened − δ,
#   the correction ΔQ's arm adds. The intrinsic spread E[σ̂²] (the hybrid's
#   `intrinsic_var`) and the truth's σ̂², apart from the estimand. The failed and
#   starved cells.
#
# SMOKE=1 is a plumbing pass, never a result. It runs dataset v6 only (B7's
# held-out fit, not an evidence dataset), with its chain 1 standing in for V7 and
# its chain 2 for B7: 60 s, 4 draws in 2 blocks of 2, 3 cells per draw, and
# n_ref = n_truth = 4. Its `prep` still checks the 100 evidence files.

include(joinpath(@__DIR__, "pilot_common.jl"))
using InferCell: AbstractSubModel
using Random
using Statistics: median, quantile, std

const STAGE = Symbol(ARGS[1])
const OUT = ARGS[2]
const IDX = length(ARGS) >= 3 ? parse(Int, ARGS[3]) : -1

const V7_DIR = "/fred/oz022/tkimpson/infercell-16a9b-archive-2026-10-17/v7/sbc_v7_16a9b"
const B7_DIR = "/fred/oz022/tkimpson/infercell-89-b7-archive-2026-10-27/run"
const EVIDENCE = 1:50
const POSTERIORS = (:v7, :b7)
const DATASETS = SMOKE ? [0] : collect(EVIDENCE)     # 0 is v6, in the smoke pass
const H = SMOKE ? 60.0 : M0_HORIZON_S
const N_DRAWS = SMOKE ? 4 : 100
const N_BLOCKS = SMOKE ? 2 : 10
const N_CELLS = SMOKE ? 3 : 887
const N_REF = SMOKE ? 4 : 9_849
const N_TRUTH = N_REF
const KD = Knockdown(; dose = 0.3)
const LEVELS = (0.5, 0.9)
const ARMS = (:hybrid, :mean_field, :mean_field_delta)
const N_D, N_P = length(DATASETS), length(POSTERIORS)
const N_WORK = N_D * N_P + N_D + N_D * N_P * N_BLOCKS

block_seed(r, b) = SMOKE ? 19_541_900 + (b - 1) : 19_540_000 + N_BLOCKS * (r - 1) + (b - 1)
delta_seed(r) = SMOKE ? 19_541_910 : 19_541_000 + r
mf_seed(r) = SMOKE ? 19_541_920 : 19_541_200 + r
truth_seed(r) = SMOKE ? 19_541_930 : 19_541_500 + r

function posterior(r, p)
    r == 0 && return posterior_draws(V6[p === :v7 ? 1 : 2])
    return p === :v7 ? posterior_draws(joinpath(V7_DIR, "replicate$r.jls")) :
                       posterior_draws([joinpath(B7_DIR, "r$r", "chain_$c.jls") for c in 1:2])
end

# What `posterior` opens: V7's replicates, and B7's chains with the prep.jls beside
# them, which holds B7's truth. The smoke pass also reads v6's.
input_files() = vcat([joinpath(V7_DIR, "replicate$r.jls") for r in EVIDENCE],
                     [joinpath(B7_DIR, "r$r", f) for r in EVIDENCE
                      for f in ("chain_1.jls", "chain_2.jls", "prep.jls")],
                     SMOKE ? [V6; joinpath(dirname(V6[1]), "prep.jls")] : String[])
# Their SHA-256 sums, as `sha256sum` prints them: 64 hex digits, two spaces, the path.
input_sums() = [(; sha256 = l[1:64], path = l[67:end])
                for l in eachline(IOBuffer(read(`sha256sum $(input_files())`, String)))]

all_draws(pd) = round.(Int, range(1, length(pd); length = N_DRAWS))
block_draws(pd, b) = all_draws(pd)[((b - 1) * (N_DRAWS ÷ N_BLOCKS) + 1):(b * (N_DRAWS ÷ N_BLOCKS))]

# Work index i (0-based) → its task. δ first, then the truths, then the blocks,
# so the longest single tasks start first.
function task(i)
    0 <= i < N_WORK || error("work index $i is outside 0…$(N_WORK - 1)")
    if i < N_D * N_P
        return (; kind = :delta, r = DATASETS[i ÷ N_P + 1], p = POSTERIORS[i % N_P + 1], b = 0)
    end
    i -= N_D * N_P
    i < N_D && return (; kind = :truth, r = DATASETS[i + 1], p = :none, b = 0)
    i -= N_D
    rp, b = divrem(i, N_BLOCKS)
    return (; kind = :hybrid, r = DATASETS[rp ÷ N_P + 1], p = POSTERIORS[rp % N_P + 1], b = b + 1)
end
work_index(kind, r, p, b) = findfirst(i -> task(i) == (; kind, r, p, b), 0:(N_WORK - 1)) - 1
work_name(i) = (t = task(i); string(lpad(i, 4, '0'), "_", t.kind, "_r", t.r,
                                     t.p === :none ? "" : "_$(t.p)", t.b == 0 ? "" : "_b$(t.b)"))
finish_task(i) = (0 <= i < N_D * N_P || error("finish index $i is outside 0…$(N_D * N_P - 1)");
                  (; r = DATASETS[i ÷ N_P + 1], p = POSTERIORS[i % N_P + 1]))
finish_name(i) = (t = finish_task(i); string(lpad(i, 3, '0'), "_r", t.r, "_", t.p))

# Written under tmp/ and moved into place; an existing output is never replaced.
function put(dir, name, x)
    path = joinpath(OUT, dir, "$name.jls")
    isfile(path) && error("$path exists; this run never overwrites an output")
    mkpath(joinpath(OUT, "tmp"))
    mkpath(joinpath(OUT, dir))
    save(joinpath(OUT, "tmp"), name, x; stage = string(STAGE))
    mv(joinpath(OUT, "tmp", "$name.jls"), path)
    say("wrote $path")
end

# A saved output, refused unless this commit wrote it: a run is one commit.
function read_out(dir, name)
    x = deserialize(joinpath(OUT, dir, "$name.jls"))
    x.commit == COMMIT || error("$dir/$name was written at $(x.commit), not $COMMIT")
    return x.result
end

function prep()
    # The inputs, pinned. B7's archive has its own SHA256SUMS, so its files are
    # checked against it; V7's archive has none.
    sums = input_sums()
    b7 = dirname(B7_DIR)
    archived = Dict(joinpath(b7, l[67:end]) => l[1:64] for l in eachline(joinpath(b7, "SHA256SUMS")))
    for s in sums
        startswith(s.path, b7 * "/") && get(archived, s.path, "") != s.sha256 &&
            error("$(s.path): SHA-256 $(s.sha256), not its archive's $(get(archived, s.path, "(none)"))")
    end
    say("$(length(sums)) input files hashed; B7's $(count(s -> startswith(s.path, b7 * "/"), sums)) ",
        "match their archive's SHA256SUMS")
    rows = map(EVIDENCE) do r
        v, b = posterior(r, :v7), posterior(r, :b7)
        v.names == b.names == M0_TARGETS || error("r$r: names $(v.names), $(b.names)")
        # V7 stores the truth on the log scale and B7 on the natural one, so the two
        # may differ by roundoff.
        keys(v.truth) == keys(b.truth) || error("r$r: the truths name different parameters")
        dev = maximum(abs(v.truth[k] / b.truth[k] - 1) for k in keys(v.truth))
        dev < 1e-12 || error("r$r: V7's and B7's truths differ, by $dev relative")
        length(v) == 556 && length(b) == 400 ||
            error("r$r: $(length(v)) V7 and $(length(b)) B7 draws, not 556 and 400")
        (; r, n_v7 = length(v), n_b7 = length(b), truth = v.truth, truth_rel_dev = dev)
    end
    say("the 100 evidence files: names, truths and draw counts as expected; the truths ",
        "agree to $(maximum(x.truth_rel_dev for x in rows)) relative")
    used = map(DATASETS) do r
        pds = Dict(p => posterior(r, p) for p in POSTERIORS)
        (; r, draws = Dict(p => all_draws(pds[p]) for p in POSTERIORS),
         blocks = Dict(p => [block_draws(pds[p], b) for b in 1:N_BLOCKS] for p in POSTERIORS),
         truth = pds[:v7].truth)
    end
    for u in used, p in POSTERIORS
        allunique(u.draws[p]) && reduce(vcat, u.blocks[p]) == u.draws[p] ||
            error("r$(u.r) $p: the draws are not $N_BLOCKS disjoint blocks")
    end
    tasks = [(; i, name = work_name(i), task(i)...) for i in 0:(N_WORK - 1)]
    put(".", "manifest", (; evidence = rows, used, tasks, H, N_DRAWS, N_BLOCKS, N_CELLS, N_REF,
                          N_TRUTH, dose = KD.dose, inputs = sums))
    # For `sha256sum -c inputs.sha256`.
    write(joinpath(OUT, "inputs.sha256"), join(("$(s.sha256)  $(s.path)\n" for s in sums)))
    open(joinpath(OUT, "manifest.tsv"), "w") do io
        println(io, "index\tname\tkind\tdataset\tposterior\tblock\tdraws\trng")
        for t in tasks
            d = t.kind === :hybrid ? only(u for u in used if u.r == t.r).blocks[t.p][t.b] : Int[]
            rng = t.kind === :hybrid ? block_seed(t.r, t.b) :
                  t.kind === :delta ? delta_seed(t.r) : truth_seed(t.r)
            println(io, join((t.i, t.name, t.kind, t.r, t.p, t.b, join(d, ","), rng), '\t'))
        end
    end
end

function work(i)
    t = task(i)
    name = work_name(i)
    isfile(joinpath(OUT, "work", "$name.jls")) && error("work/$name exists")
    say("work $name: $t")
    ms = collect(AbstractSubModel, m0_models())
    local x, p
    sec = @elapsed if t.kind === :delta
        pd = posterior(t.r, t.p)
        seeds = InferCell._cell_seeds(Xoshiro(delta_seed(t.r)), N_REF)
        x = delta_correction(ms, haldane_values(metered(ms), InferCell._theta_ref(pd)), seeds;
                             horizon = H, intervention = KD, threads = false)
        # The SE of δ_intervened − δ, ΔQ's correction: its cells are paired, as the
        # truth's are.
        e, es2, en = InferCell._moments(InferCell._paired_delta(x.hybrid, x.hybrid_intervened))
        p = (; t, x.hybrid_mean, x.hybrid_se, x.n_ok, x.delta, x.delta_intervened,
             effect = e, effect_sigma2 = es2, effect_n = en, effect_se = sqrt(es2 / en),
             mean_field_q = x.mean_field.q, mean_field_intervened_q = x.mean_field_intervened.q,
             x.seeds, q = [c.q for c in x.hybrid], q_intervened = [c.q for c in x.hybrid_intervened],
             failed = count(c -> c.failure !== :none, x.hybrid),
             failed_intervened = count(c -> c.failure !== :none, x.hybrid_intervened))
    elseif t.kind === :truth
        truth = posterior(t.r, :v7).truth
        θ = reshape([truth[n] for n in M0_TARGETS], 1, :)
        pd = PosteriorDraws(M0_TARGETS, θ, nothing, [1], [1], nothing, truth,
                            (kind = :truth, dataset = t.r))
        x = predict(ms, pd; horizon = H, n_cells = N_TRUTH, intervention = KD,
                    rng = Xoshiro(truth_seed(t.r)), threads = false)
        mu, s2, n = InferCell._moments([c.q for c in x.cells[1]])
        e, es2, en = InferCell._moments(InferCell._paired_delta(x.cells[1], x.cells_intervened[1]))
        p = (; t, prediction = plain(x), mu, sigma2 = s2, n, se = sqrt(s2 / n), effect = e,
             effect_sigma2 = es2, effect_n = en, effect_se = sqrt(es2 / en),
             q = [c.q for c in x.cells[1]], q_intervened = [c.q for c in x.cells_intervened[1]])
    else
        pd = posterior(t.r, t.p)
        x = predict(ms, pd; horizon = H, n_cells = N_CELLS, draws = block_draws(pd, t.b),
                    intervention = KD, rng = Xoshiro(block_seed(t.r, t.b)), threads = false)
        p = (; t, prediction = plain(x))
    end
    say("work $name: $(round(sec / 3600; digits = 2)) h, maxrss $(round(Sys.maxrss() / 2^30; digits = 2)) GiB, ",
        "object $(round(Base.summarysize(x) / 2^20; digits = 1)) MiB")
    put("work_plain", name, merge(p, (; seconds = sec, maxrss_gib = Sys.maxrss() / 2^30,
                                      object_mib = Base.summarysize(x) / 2^20)))
    put("work", name, x)
end

function finish(i)
    (; r, p) = finish_task(i)
    name = finish_name(i)
    say("finish $name")
    ms = collect(AbstractSubModel, m0_models())
    pd = posterior(r, p)
    blocks = [read_out("work", work_name(work_index(:hybrid, r, p, b))) for b in 1:N_BLOCKS]
    hy = combine_rows(blocks)
    dc = read_out("work", work_name(work_index(:delta, r, p, 0)))
    tr = read_out("work", work_name(work_index(:truth, r, :none, 0)))
    mf = predict(ms, pd; horizon = H, forward = :mean_field, draws = all_draws(pd), delta = dc,
                 intervention = KD, rng = Xoshiro(mf_seed(r)), threads = false)
    hy.draws == mf.draws == all_draws(pd) || error("$name: the draws are not paired")
    hy.provenance.block_rows == fill(N_DRAWS ÷ N_BLOCKS, N_BLOCKS) ||
        error("$name: blocks of $(hy.provenance.block_rows) rows")
    isdisjoint(hy.seeds, dc.seeds) && isdisjoint(hy.seeds, tr.seeds) &&
        isdisjoint(dc.seeds, tr.seeds) || error("$name: the block, δ and truth seeds meet")
    put("finish", name, (; r, p, hybrid = plain(hy), mean_field = plain(mf),
                         delta = (; dc.hybrid_mean, dc.hybrid_se, dc.n_ok, dc.delta,
                                  dc.delta_intervened)))
end

interval(v, L) = (quantile(v, (1 - L) / 2), quantile(v, (1 + L) / 2))

# One arm's sample against the truth.
function score_arm(vals, truth)
    v = filter(isfinite, vals)
    n = length(v)
    ints = Dict(L => interval(v, L) for L in LEVELS)
    m, s = sum(v) / n, std(v)
    return (; n, n_no_q = length(vals) - n, mean = m, sd = s,
            cover = Dict(L => ints[L][1] <= truth <= ints[L][2] for L in LEVELS),
            width = Dict(L => ints[L][2] - ints[L][1] for L in LEVELS), intervals = ints,
            pit = count(<(truth), v) / n, z = (m - truth) / s)
end

function samples(f, target)
    h, m, d = f.hybrid, f.mean_field, f.delta
    target === :q && return (; hybrid = vec(h.mu), mean_field = vec(m.mu),
                             mean_field_delta = vec(m.mu) .+ d.delta)
    return (; hybrid = vec(h.effect), mean_field = vec(m.effect),
            mean_field_delta = vec(m.effect) .+ (d.delta_intervened - d.delta))
end

function score()
    pinned, current = read_out(".", "manifest").inputs, input_sums()
    pinned == current || error("the inputs changed since prep: ",
                               join(unique(x.path for x in symdiff(pinned, current)), ", "))
    fin = Dict((x.r, x.p) => x for x in (read_out("finish", finish_name(i)) for i in 0:(N_D * N_P - 1)))
    tru = Dict(r => read_out("work_plain", work_name(work_index(:truth, r, :none, 0))) for r in DATASETS)
    dpl = Dict((r, p) => read_out("work_plain", work_name(work_index(:delta, r, p, 0)))
               for r in DATASETS, p in POSTERIORS)
    for r in DATASETS
        fin[(r, :v7)].hybrid.seeds == fin[(r, :b7)].hybrid.seeds ||
            error("r$r: V7's and B7's cell seeds differ; the run promised common seeds")
    end
    rows = NamedTuple[]
    for target in (:q, :dq), r in DATASETS, p in POSTERIORS
        f, t = fin[(r, p)], tru[r]
        truth, tse = target === :q ? (t.mu, t.se) : (t.effect, t.effect_se)
        s = target === :q ? f.hybrid.summary.uniform : f.hybrid.effect_summary.uniform
        psd = sqrt(s.parameter_var)
        smp = samples(f, target)
        for a in ARMS
            push!(rows, (; target, r, p, arm = a, truth, truth_se = tse,
                         score_arm(getproperty(smp, a), truth)...,
                         mc_share = s.mc / s.var, truth_se_ratio = tse / psd,
                         delta_se_ratio = (target === :q ? dpl[(r, p)].hybrid_se :
                                           dpl[(r, p)].effect_se) / psd,
                         intrinsic_var = s.intrinsic_var,
                         truth_sigma2 = target === :q ? t.sigma2 : t.effect_sigma2,
                         failed = a === :hybrid ? f.hybrid.statuses : f.mean_field.statuses,
                         failed_intervened = a === :hybrid ? f.hybrid.statuses_intervened :
                                             f.mean_field.statuses_intervened))
        end
    end
    # Each arm's shift and width against V7-hybrid, on the same dataset.
    ref = Dict((x.target, x.r) => x for x in rows if x.p === :v7 && x.arm === :hybrid)
    rows = [merge(x, (; shift = (x.mean - ref[(x.target, x.r)].mean) / ref[(x.target, x.r)].sd,
                      width_ratio = x.width[0.9] / ref[(x.target, x.r)].width[0.9])) for x in rows]
    put(".", "score", (; rows, levels = LEVELS, datasets = DATASETS))
    write_score(rows)
end

pct(x) = string(round(100x; digits = 1), "%")
r3(x) = round(x; sigdigits = 3)
iqr(v) = "[$(r3(quantile(v, 0.25))), $(r3(quantile(v, 0.75)))]"

function write_score(rows)
    io = open(joinpath(OUT, "score.md"), "w")
    out(x...) = (println(io, x...); println(x...))
    out("# M0's crossed comparison (KIM-40) and ptsG panel (KIM-41): the scores\n")
    out("At $COMMIT, job $JOB, $(now()). $(SMOKE ? "**SMOKE: plumbing, not a result.** " : "")",
        "$(N_D) datasets, $(N_DRAWS) draws per posterior, $(N_CELLS) cells per draw, ",
        "n_ref = n_truth = $(N_REF), $(H) s, knockdown dose $(KD.dose).\n")
    label(p, a) = string(uppercase(String(p)), " ",
                         a === :hybrid ? "hybrid" : a === :mean_field ? "mean field" : "mean field + δ")
    for target in (:q, :dq)
        out("## ", target === :q ? "Q = lactate_out (KIM-40)" : "ΔQ under the knockdown (KIM-41)", "\n")
        out("| arm | cover 50% | cover 90% | median 90% width | width / V7 hybrid, median [IQR] | ",
            "shift from V7 hybrid (SD), median | z, mean ± SD | draws without Q |")
        out("|---|---|---|---|---|---|---|---|")
        for p in POSTERIORS, a in ARMS
            x = [y for y in rows if y.target === target && y.p === p && y.arm === a]
            n = length(x)
            cov(L) = (c = count(y -> y.cover[L], x) / n; "$(pct(c)) ± $(pct(sqrt(c * (1 - c) / n)))")
            w = [y.width[0.9] for y in x]
            wr = [y.width_ratio for y in x]
            z = [y.z for y in x]
            out("| $(label(p, a)) | $(cov(0.5)) | $(cov(0.9)) | $(r3(median(w))) | ",
                "$(r3(median(wr))) $(iqr(wr)) | $(r3(median([y.shift for y in x]))) | ",
                "$(r3(sum(z) / n)) ± $(r3(n > 1 ? std(z) : NaN)) | $(sum(y -> y.n_no_q, x)) |")
        end
        out("\nThe truth's PIT (fraction of the sample below it), in 10 bins:\n")
        out("| arm | ", join(["$(b - 1)/10" for b in 1:10], " | "), " |")
        out("|---|", repeat("---|", 10))
        for p in POSTERIORS, a in ARMS
            x = [y.pit for y in rows if y.target === target && y.p === p && y.arm === a]
            bins = [count(v -> min(floor(Int, 10v), 9) == b - 1, x) for b in 1:10]
            out("| $(label(p, a)) | ", join(bins, " | "), " |")
        end
        out("\nAgainst the 10% rules, and the intrinsic spread (the hybrid's, per posterior):\n")
        out("| posterior | MC share, median (max; > 10%) | truth SE / SD, median (max; > 10%) | ",
            "δ SE / SD, median (max) | √E[σ̂²], median | √σ̂²(θ*), median |")
        out("|---|---|---|---|---|---|")
        for p in POSTERIORS
            x = [y for y in rows if y.target === target && y.p === p && y.arm === :hybrid]
            f = [y.mc_share for y in x]
            tr = [y.truth_se_ratio for y in x]
            dr = [y.delta_se_ratio for y in x]
            out("| $(uppercase(String(p))) | $(r3(median(f))) ($(r3(maximum(f))); $(count(>(0.1), f))) | ",
                "$(r3(median(tr))) ($(r3(maximum(tr))); $(count(>(0.1), tr))) | ",
                "$(r3(median(dr))) ($(r3(maximum(dr)))) | $(r3(median(sqrt.([y.intrinsic_var for y in x])))) | ",
                "$(r3(median(sqrt.([y.truth_sigma2 for y in x])))) |")
        end
        out("")
    end
    out("## Cells by status, summed over the datasets\n")
    out("| arm | baseline | intervened |")
    out("|---|---|---|")
    for p in POSTERIORS, a in (:hybrid, :mean_field)
        x = [y for y in rows if y.target === :q && y.p === p && y.arm === a]
        sm(k) = join(["$s $(sum(y -> Dict(getproperty(y, k))[s], x))" for s in CELL_STATUSES], ", ")
        out("| $(label(p, a)) | $(sm(:failed)) | $(sm(:failed_intervened)) |")
    end
    close(io)
    open(joinpath(OUT, "score_datasets.tsv"), "w") do io
        println(io, "target\tdataset\tposterior\tarm\ttruth\ttruth_se\tmean\tsd\tlo50\thi50\tlo90\thi90\t",
                "cover50\tcover90\tpit\tz\tshift\twidth_ratio\tmc_share\tn_no_q")
        for y in rows
            println(io, join((y.target, y.r, y.p, y.arm, y.truth, y.truth_se, y.mean, y.sd,
                              y.intervals[0.5]..., y.intervals[0.9]..., y.cover[0.5], y.cover[0.9],
                              y.pit, y.z, y.shift, y.width_ratio, y.mc_share, y.n_no_q), '\t'))
        end
    end
end

function main()
    mkpath(OUT)
    say("M0 evidence (KIM-40/41), stage $STAGE$(IDX >= 0 ? " $IDX" : ""), at $COMMIT, job $JOB, ",
        "$(Threads.nthreads()) thread(s)$(SMOKE ? ", SMOKE" : "")")
    STAGE === :prep ? prep() : STAGE === :work ? work(IDX) : STAGE === :finish ? finish(IDX) :
    STAGE === :score ? score() : error("unknown stage $STAGE")
    say("done")
end

main()
