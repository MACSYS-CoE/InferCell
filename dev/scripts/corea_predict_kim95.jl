# KIM-95: scored predictions on the reduced model (Linear M9). The targets and scores
# are frozen in spec/prereg/kim-92.md §4 at tag prereg/kim-92; Tom's calls of
# 2026-11-18, made before any prediction existed, fix their form:
#   - pools and transcripts are scored as population means per save, the mean over
#     cells with no σ_b, like Q's frozen estimand. The truth is the mean over
#     N_TRUTH fresh cells at θ*;
#   - the sizes follow M0's declared rules (spec §12 2026-10-31, KIM-69): the MC term
#     at most 10% of the predictive spread, and the truth's SE at most 10% of the
#     predictive SD. The `sizing` stage reads the predictive only, never the truth;
#   - no viability rule on either arm; failed and starved cells are counted per arm;
#   - the scores are in corea_predict_kim95_score.jl.
#
# Usage: julia --project dev/scripts/corea_predict_kim95.jl <stage> <out> [<arm> <chunk>]
# Submit through dev/scripts/corea_predict_kim95.sh, from a worktree at a pushed commit.
# Stages:
#   prep <out>      the pilot's 1600 draws (KIM-92's archive), checked to the bit against
#                   its pilot_draws.tsv; θ* from its dataset; PRIOR_POOL draws from the
#                   21-coordinate prior on PRIOR_SEED; to prep.jls and prep.md;
#   sizing <out>    SIZING_DRAWS draws × SIZING_CELLS cells per arm (posterior, prior),
#                   both arms of do(ptsG), on seeds of their own: Var_θ and σ² per
#                   target, the sizes the rules give, and the cost per cell, to
#                   sizing.jls, sizing.tsv and sizing.md. It reads no truth.
#                   SIZING_ARMS (e.g. "posterior") limits the arms, and SIZING_TAG=resolve
#                   writes sizing_resolve.* on seeds of its own instead (the resolving
#                   pass, when the first did not resolve Var_θ);
#   run <out> <arm> <k>  chunk k of arm `posterior`, `prior` or `truth`, at the frozen
#                   sizes (refused until FROZEN, except under SMOKE), to
#                   rows_<arm>_<k>.jls: per pair, μ̂, σ̂² and n for every target;
#   collect <out>   the chunks joined, every cell seed checked unique and disjoint from
#                   the dataset's and the pilot's reference ensemble's, to plain TSVs
#                   (targets.tsv, predictive_<arm>_{mu,s2,n}.tsv, truth.tsv, status.tsv)
#                   and collect.md. Item 3's scores read only these.
# Targets (prereg §4): Q = lactate_out over the cycle, its population mean, at baseline
# (`Q`) and under Knockdown(dose = 0.3) on ptsG (`Q_kd`), and the paired effect (`dQ`);
# then M0_PANEL's 17 pools and the 17 genes' transcripts at each save after the fit
# window, `<species>@<t>`. The run goes to the full cycle and is cut at COREA_FIT_S
# (KIM-94: never a separate 3150 s run).
# SMOKE=1: a 120 s run cut at 60 s, 2 draws × 2 cells, 3 truth cells: plumbing only.

include(joinpath(@__DIR__, "pilot_common.jl"))
using Printf
using Random
using Statistics: median, quantile

envint(k, default) = (v = strip(get(ENV, k, "")); isempty(v) ? default : parse(Int, v))

const PILOT = "/fred/oz022/tkimpson/infercell-kim92-archive-2026-11-17"
const PILOT_FILES = [joinpath(PILOT, f) for c in 1:4 for f in ("chain_$c.jls", "extend_$c.jls")]
const HORIZON = SMOKE ? 120.0 : COREA_CYCLE_S
const AFTER = SMOKE ? 60.0 : COREA_FIT_S
const EVERY = 60.0
const DOSE = 0.3                                   # KIM-39 / M4 freeze, on ptsG
knockdown() = Knockdown(; dose = DOSE)

# Seeds: rng seeds for `predict`'s stream, a block of their own (17_95x_xxx). The cell
# seeds `predict` draws from them are checked against the dataset's in `collect`.
const PRIOR_SEED = 17_950_001
const PRIOR_POOL = 1600
const SIZING_TAG = get(ENV, "SIZING_TAG", "")
const SIZING_SEED = SIZING_TAG == "resolve" ? Dict(:posterior => 17_950_200, :prior => 17_950_201) :
                    isempty(SIZING_TAG) ? Dict(:posterior => 17_950_100, :prior => 17_950_101) :
                    error("SIZING_TAG is \"\" or \"resolve\", not $(SIZING_TAG)")
const SIZING_ARMS = Symbol.(split(get(ENV, "SIZING_ARMS", "posterior,prior"), ','))
const SIZING_NAME = isempty(SIZING_TAG) ? "sizing" : "sizing_$(SIZING_TAG)"
const ARM_SEED = Dict(:posterior => 17_951_000, :prior => 17_952_000, :truth => 17_953_000)
const SIZING_DRAWS = SMOKE ? 2 : envint("SIZING_DRAWS", 40)
const SIZING_CELLS = SMOKE ? 2 : envint("SIZING_CELLS", 8)

# The frozen sizes (the freeze, 2026-11-18, Tom's call on Linear KIM-95, before any truth
# was read). M0's rules (KIM-69), sized on Q_kd, the one Q target whose posterior Var_θ
# the resolving pass resolved (job 21126961, 40 draws × 100 cells, sizing_resolve.tsv):
# Var_θ/σ² = 0.009251, so n_cells = ⌈9 σ²/Var_θ⌉ = 973 and n_truth = ⌈100 σ²/Var_θ⌉ =
# 10,810 (M0's formula). Q's and ΔQ's posterior Var_θ stayed unresolved (Var_θ/σ² below
# about 0.0045), so their predictive spread is mostly the MC term, reported per target.
# The prior: 200 cells a draw, above every target's rule in the first sizing (job
# 21122610, max 202 over 1,802 per-save targets; 1 on Q). Chunks hold a ~6 h job each on
# 32 threads at 75.5 (posterior) and 62 (prior) CPU-s a cell. `run` refuses unless FROZEN.
const FROZEN = true
const N_DRAWS = SMOKE ? 2 : 100                    # per arm, evenly spaced over its pool
const N_CELLS = Dict(:posterior => SMOKE ? 2 : 973, :prior => SMOKE ? 2 : 200)   # per draw
const N_TRUTH = SMOKE ? 3 : 10_810                 # cells at θ*
const CHUNKS = Dict(:posterior => SMOKE ? 2 : 20, :prior => SMOKE ? 2 : 4, :truth => SMOKE ? 2 : 4)

# The rules (KIM-69): MC share ≤ F of the predictive spread, truth SE ≤ F of its SD.
const F_RULE = 0.1

const ARMS = (:posterior, :prior, :truth)
const MODELS = d11_models()
const LOCI = [g.locus for g in only(m for m in MODELS if m isa CoreATranscription).genes]
const SPECIES = vcat(M0_PANEL, [transcript_state(l) for l in LOCI])
family(sp) = sp in M0_PANEL ? "pool" : "transcript"

evenly(n_pool, n) = n == n_pool ? collect(1:n) : unique(round.(Int, range(1, n_pool; length = n)))

# ---------------------------------------------------------------- prep

function read_pilot_tsv()
    lines = readlines(joinpath(PILOT, "pilot_draws.tsv"))
    X = permutedims(reduce(hcat, [parse.(Float64, split(l, '\t')) for l in lines[2:end]]))
    return split(lines[1], '\t'), X
end

function run_prep(out)
    mkpath(out)
    pd = posterior_draws(PILOT_FILES)
    head, X = read_pilot_tsv()
    head == String.(InferCell._log_names(COREA_TARGETS)) && pd.names == COREA_TARGETS &&
        length(pd) == 1600 == size(X, 1) && pd.values == exp.(X[:, 1:20]) &&
        pd.sigma_b == exp.(X[:, 21]) && isempty(held_parameters(MODELS, pd).held) ||
        error("the pilot's draws are not its pilot_draws.tsv")
    ds = load(PILOT, "dataset")
    tv = Dict(ds.truth.values)
    pd.truth == merge(Dict(n => tv[n] for n in COREA_TARGETS), Dict(:sigma_b => ds.sigma)) ||
        error("the draws' truth is not the dataset's")
    # The prior: the pilot's own, independent normals in log space (mixed_setup's).
    s = mixed_setup(ds; models = MODELS, targets = COREA_TARGETS)
    mu, sd = mixed_prior_median(s), mixed_prior_sd(s)
    rng = Xoshiro(PRIOR_SEED)
    Z = mu' .+ sd' .* randn(rng, PRIOR_POOL, length(mu))
    prior = PosteriorDraws(COREA_TARGETS, exp.(Z[:, 1:20]), exp.(Z[:, 21]), ones(Int, PRIOR_POOL),
                           collect(1:PRIOR_POOL), nothing, pd.truth,
                           (; kind = :prior, seed = PRIOR_SEED, route = :mixed))
    theta = [tv[n] for n in COREA_TARGETS]
    truth = PosteriorDraws(COREA_TARGETS, reshape(theta, 1, :), [ds.sigma], [0], [0], nothing,
                           pd.truth, (; kind = :truth, replicate_seed = ds.replicate_seed))
    ref = load(PILOT, "prep")
    save(out, "prep", (; posterior = pd, prior, truth, dataset_seeds = collect(ds.seeds),
                       reference_seeds = collect(ref.seeds), prior_mu = mu, prior_sd = sd);
         stage = "prep")
    open(joinpath(out, "prep.md"), "w") do io
        println(io, "# KIM-95 prep\n\ncommit `", COMMIT, "`, job ", JOB, SMOKE ? ", SMOKE" : "", "\n")
        println(io, "- posterior: ", length(pd), " draws from ", PILOT, ", == pilot_draws.tsv")
        println(io, "- prior: ", PRIOR_POOL, " draws on seed ", PRIOR_SEED)
        println(io, "- θ*: dataset seed ", ds.replicate_seed, ", σ_b ", ds.sigma)
        println(io, "- dataset cells ", length(ds.seeds), ", reference cells ", length(ref.seeds))
    end
    say("prep: ", length(pd), " posterior draws == pilot_draws.tsv; ", PRIOR_POOL, " prior draws")
end

# ---------------------------------------------------------------- one prediction

predict_arm(pd, draws, n_cells, seed) =
    predict(MODELS, pd; horizon = HORIZON, forward = :hybrid, completion = :none,
            intervention = knockdown(), n_cells, draws, trajectory = SPECIES, every = EVERY,
            rng = Xoshiro(seed))

# Per pair: μ̂, σ̂² and n for each target, at baseline over its finished cells; Q_kd
# over the intervened arm's; dQ over the cells finished in both.
function pair_rows(p::Prediction)
    M = size(p.cells, 1)
    tv = [trajectory_values(p, sp; after = AFTER) for sp in SPECIES]
    t = tv[1].t
    all(x -> x.t == t, tv) || error("the species' saves differ")
    ids = vcat(["Q", "Q_kd", "dQ"], ["$sp@$(Int(round(ti)))" for sp in SPECIES for ti in t])
    T = length(ids)
    mu, s2, n = fill(NaN, M, T), fill(NaN, M, T), zeros(Int, M, T)
    sv(x) = length(x) >= 2 ? var(x) : NaN
    for m in 1:M
        b, k = p.cells[m, 1], p.cells_intervened[m, 1]
        [c.seed for c in b] == [c.seed for c in k] || error("the arms' seeds differ on pair $m")
        qb = [c.q for c in b if c.failure === :none]
        qk = [c.q for c in k if c.failure === :none]
        both = [k[i].q - b[i].q for i in eachindex(b) if b[i].failure === :none && k[i].failure === :none]
        for (j, x) in enumerate((qb, qk, both))
            mu[m, j], s2[m, j], n[m, j] = isempty(x) ? NaN : mean(x), sv(x), length(x)
        end
        # dQ is predict's paired effect; check the reading against it.
        isempty(both) || mean(both) ≈ p.effect[m, 1] ||
            error("pair $m: the paired effect is not predict's")
        j = 3
        for x in tv
            X = x.values[m, 1]                       # finished cells × saves
            for c in 1:length(t)
                j += 1
                col = X[:, c]
                mu[m, j], s2[m, j], n[m, j] = isempty(col) ? NaN : mean(col), sv(col), length(col)
            end
        end
    end
    return (; ids, t, mu, s2, n, draws = p.draws, seeds = vec(p.seeds),
            statuses = p.statuses, statuses_intervened = p.statuses_intervened)
end

# ---------------------------------------------------------------- sizing

# For one target over the pairs: σ² (the mean σ̂²), Var_θ net of the MC term, and the
# rules' sizes: n_cells ≥ (1 − F)/F · σ²/Var_θ, then n_truth ≥ σ²/(F² · pred_var).
function rule(mu, s2, n)
    ok = findall(i -> isfinite(mu[i]) && isfinite(s2[i]) && n[i] >= 2, eachindex(mu))
    length(ok) >= 2 || return (; sigma2 = NaN, var_theta = NaN, n_cells = -1, n_truth = -1)
    σ2 = mean(s2[ok])
    vθ = var(mu[ok]) - mean(s2[ok] ./ n[ok])
    vθ > 0 || return (; sigma2 = σ2, var_theta = vθ, n_cells = -1, n_truth = -1)
    nc = max(1, ceil(Int, (1 - F_RULE) / F_RULE * σ2 / vθ))
    nt = ceil(Int, σ2 / (F_RULE^2 * (vθ + σ2 / nc)))
    return (; sigma2 = σ2, var_theta = vθ, n_cells = nc, n_truth = nt)
end

function run_sizing(out)
    pp = load(out, "prep")
    res = Dict{Symbol, Any}()
    for arm in SIZING_ARMS
        pd = getproperty(pp, arm)
        draws = evenly(length(pd), SIZING_DRAWS)
        t = @elapsed p = predict_arm(pd, draws, SIZING_CELLS, SIZING_SEED[arm])
        r = pair_rows(p)
        ncell = 2 * length(draws) * SIZING_CELLS
        res[arm] = (; r..., seconds = t, cells = ncell, threads = Threads.nthreads(),
                    rules = [rule(r.mu[:, j], r.s2[:, j], r.n[:, j]) for j in eachindex(r.ids)])
        say(@sprintf("sizing %s: %d draws × %d cells × 2 arms in %.0f s on %d threads (%.1f CPU-s a cell)",
                     arm, length(draws), SIZING_CELLS, t, Threads.nthreads(),
                     t * Threads.nthreads() / ncell))
    end
    save(out, SIZING_NAME, res; stage = SIZING_NAME)
    write_sizing(out, res)
end

function write_sizing(out, res)
    arms = [a for a in (:posterior, :prior) if haskey(res, a)]
    ids = res[arms[1]].ids
    open(joinpath(out, SIZING_NAME * ".tsv"), "w") do io
        println(io, join(["target", "arm", "sigma2", "var_theta", "n_cells", "n_truth"], '\t'))
        for arm in arms, (j, id) in enumerate(ids)
            r = res[arm].rules[j]
            println(io, join([id, arm, @sprintf("%.17g", r.sigma2), @sprintf("%.17g", r.var_theta),
                              r.n_cells, r.n_truth], '\t'))
        end
    end
    fam(id) = id in ("Q", "Q_kd", "dQ") ? id : family(Symbol(first(split(id, '@'))))
    open(joinpath(out, SIZING_NAME * ".md"), "w") do io
        println(io, "# KIM-95 ", SIZING_NAME, " (reads the predictive only, never the truth)\n")
        println(io, "commit `", COMMIT, "`, job ", JOB, SMOKE ? ", SMOKE" : "", ". ",
                SIZING_DRAWS, " draws × ", SIZING_CELLS, " cells per arm (", join(arms, ", "),
                "), both arms of do(ptsG), rng seeds ", join([SIZING_SEED[a] for a in arms], ", "), ".\n")
        for arm in arms
            a = res[arm]
            println(io, @sprintf("- %s: %.0f s on %d threads, %.1f CPU-s a cell (%d cells)",
                                 arm, a.seconds, a.threads, a.seconds * a.threads / a.cells, a.cells))
            println(io, "  statuses: ", _sorted(a.statuses), "; intervened: ", _sorted(a.statuses_intervened))
        end
        println(io, "\nThe rules' sizes: n_cells = ⌈9 σ²/Var_θ⌉, then n_truth = ⌈100 σ²/(Var_θ + σ²/n_cells)⌉.",
                " −1: Var_θ not resolved (≤ 0) or too few cells.\n")
        println(io, "| target or family | arm | n targets | unresolved | n_cells median | n_cells 90% | n_cells max | n_truth median | n_truth 90% | n_truth max |")
        println(io, "|---|---|---|---|---|---|---|---|---|---|")
        for f in ("Q", "Q_kd", "dQ", "pool", "transcript"), arm in arms
            js = findall(id -> fam(id) == f, ids)
            rs = res[arm].rules[js]
            good = filter(r -> r.n_cells > 0, rs)
            q(v, p) = isempty(v) ? "—" : string(ceil(Int, quantile(v, p)))
            nc, nt = [r.n_cells for r in good], [r.n_truth for r in good]
            println(io, "| ", f, " | ", arm, " | ", length(js), " | ", length(rs) - length(good), " | ",
                    q(nc, 0.5), " | ", q(nc, 0.9), " | ", q(nc, 1.0), " | ", q(nt, 0.5), " | ",
                    q(nt, 0.9), " | ", q(nt, 1.0), " |")
        end
        # Whether Var_θ is resolved at all: var(μ̂) over its MC term, whose sampling SD
        # under Var_θ = 0 is about √(2/(M − 1)) for M draws.
        println(io, "\nResolution: var(μ̂) across draws over its MC term mean(σ̂²/n); 1 is no θ-spread,",
                " ± √(2/(M−1)) its rough sampling SD at Var_θ = 0. Var_θ/σ² beside it.\n")
        println(io, "| target or family | arm | ratio (median over targets) | ± | Var_θ/σ² (median) |")
        println(io, "|---|---|---|---|---|")
        for f in ("Q", "Q_kd", "dQ", "pool", "transcript"), arm in arms
            js = findall(id -> fam(id) == f, ids)
            a = res[arm]
            rat, vr, Ms = Float64[], Float64[], Int[]
            for j in js
                ok = findall(i -> isfinite(a.mu[i, j]) && isfinite(a.s2[i, j]) && a.n[i, j] >= 2, axes(a.mu, 1))
                length(ok) >= 3 || continue
                mc = mean(a.s2[ok, j] ./ a.n[ok, j])
                push!(rat, var(a.mu[ok, j]) / mc)
                push!(vr, (var(a.mu[ok, j]) - mc) / mean(a.s2[ok, j]))
                push!(Ms, length(ok))
            end
            isempty(rat) && continue
            println(io, "| ", f, " | ", arm, " | ", @sprintf("%.3f", median(rat)), " | ",
                    @sprintf("%.3f", sqrt(2 / (median(Ms) - 1))), " | ", @sprintf("%.4g", median(vr)), " |")
        end
    end
end

# ---------------------------------------------------------------- run

function chunk_draws(arm, n_pool, k)
    K = CHUNKS[arm]
    1 <= k <= K || error("chunk $k of $K")
    all = evenly(n_pool, N_DRAWS)
    return all[k:K:end]
end

function run_chunk(out, arm, k)
    FROZEN || SMOKE || error("the sizes are not frozen: commit and push the freeze first")
    pp = load(out, "prep")
    pd = getproperty(pp, arm)
    if arm === :truth
        K = CHUNKS[:truth]
        n_cells = N_TRUTH ÷ K + (k <= N_TRUTH % K ? 1 : 0)
        draws = [1]
    else
        n_cells, draws = N_CELLS[arm], chunk_draws(arm, length(pd), k)
    end
    t = @elapsed p = predict_arm(pd, draws, n_cells, ARM_SEED[arm] + k)
    r = pair_rows(p)
    save(out, "rows_$(arm)_$k", (; r..., arm, k, seconds = t, threads = Threads.nthreads(),
                                 cells = 2 * length(draws) * n_cells); stage = "run")
    say(@sprintf("run %s chunk %d: %d draws × %d cells × 2 arms in %.0f s on %d threads",
                 arm, k, length(draws), n_cells, t, Threads.nthreads()))
end

# ---------------------------------------------------------------- collect

f17(x) = @sprintf("%.17g", x)

# Chunks' (n, mean, sample variance) per target pooled into the whole's, as within- plus
# between-chunk sums of squares. A chunk with no finished cell (mean NaN) or one (s² NaN)
# adds only what it has.
function pool(nk, mk, sk)
    N = reduce(+, nk)
    m = reduce(+, [ifelse.(n .> 0, n .* μ, 0.0) for (n, μ) in zip(nk, mk)]) ./ N
    W = reduce(+, [ifelse.(n .> 1, (n .- 1) .* s, 0.0) for (n, s) in zip(nk, sk)])
    B = reduce(+, [ifelse.(n .> 0, n .* (μ .- m) .^ 2, 0.0) for (n, μ) in zip(nk, mk)])
    return m, (W .+ B) ./ (N .- 1), N
end

function write_wide(path, head, rows, X)
    open(path, "w") do io
        println(io, join(vcat(head, rows.ids), '\t'))
        for m in 1:size(X, 1)
            println(io, join(vcat(string.(collect(rows.lead[m])), f17.(X[m, :])), '\t'))
        end
    end
end

function run_collect(out)
    pp = load(out, "prep")
    parts = Dict(arm => [load(out, "rows_$(arm)_$k") for k in 1:CHUNKS[arm]] for arm in ARMS)
    ids = parts[:posterior][1].ids
    all(r -> r.ids == ids, (r for arm in ARMS for r in parts[arm])) || error("the targets differ")
    seeds = reduce(vcat, [r.seeds for arm in ARMS for r in parts[arm]])
    allunique(seeds) || error("a cell seed repeats across chunks or arms")
    isdisjoint(seeds, pp.dataset_seeds) || error("a cell seed is one of the dataset's")
    isdisjoint(seeds, pp.reference_seeds) || error("a cell seed is one of the reference ensemble's")
    t = parts[:posterior][1].t
    open(joinpath(out, "targets.tsv"), "w") do io
        println(io, "target\tfamily\tspecies\tt")
        for id in ids
            if id in ("Q", "Q_kd", "dQ")
                println(io, id, "\tQ\t", id, "\t", f17(HORIZON))
            else
                sp, ti = split(id, '@')
                println(io, id, '\t', family(Symbol(sp)), '\t', sp, '\t', ti)
            end
        end
    end
    for arm in (:posterior, :prior)
        pd = getproperty(pp, arm)
        rs = parts[arm]
        draws = reduce(vcat, [r.draws for r in rs])
        o = sortperm(draws)
        draws == unique(draws) || error("$arm: a draw repeats across chunks")
        lead = [(d, pd.chain[d], pd.iteration[d]) for d in draws[o]]
        for (name, f) in (("mu", r -> r.mu), ("s2", r -> r.s2), ("n", r -> r.n))
            X = reduce(vcat, [f(r) for r in rs])[o, :]
            write_wide(joinpath(out, "predictive_$(arm)_$name.tsv"), ["draw", "chain", "iteration"],
                       (; ids, lead), X)
        end
    end
    # The truth: the chunks' cells pooled.
    rs = parts[:truth]
    m, v, N = pool([vec(r.n) for r in rs], [vec(r.mu) for r in rs], [vec(r.s2) for r in rs])
    open(joinpath(out, "truth.tsv"), "w") do io
        println(io, "target\tmean\tsd\tn\tse")
        for j in eachindex(ids)
            println(io, ids[j], '\t', f17(m[j]), '\t', f17(sqrt(max(v[j], 0))), '\t', N[j], '\t',
                    f17(sqrt(max(v[j], 0) / N[j])))
        end
    end
    open(joinpath(out, "status.tsv"), "w") do io
        println(io, "arm\tintervened\tstatus\tcount")
        for arm in ARMS, (iv, key) in ((false, :statuses), (true, :statuses_intervened))
            d = Dict{Symbol, Int}()
            for r in parts[arm], (s, c) in getproperty(r, key)
                d[s] = get(d, s, 0) + c
            end
            for (s, c) in _sorted(d)
                println(io, arm, '\t', iv, '\t', s, '\t', c)
            end
        end
    end
    open(joinpath(out, "collect.md"), "w") do io
        println(io, "# KIM-95 collect\n\ncommit `", COMMIT, "`, job ", JOB, SMOKE ? ", SMOKE" : "", "\n")
        for arm in ARMS
            rs = parts[arm]
            println(io, @sprintf("- %s: %d chunks, %d cells (both arms of do(ptsG)), %.1f CPU-h",
                                 arm, length(rs), sum(r.cells for r in rs),
                                 sum(r.seconds * r.threads for r in rs) / 3600))
        end
        println(io, "- ", length(seeds), " cell seeds: unique, disjoint from the dataset's ",
                length(pp.dataset_seeds), " and the reference ensemble's ", length(pp.reference_seeds))
        println(io, "- ", length(ids), " targets, forecast saves ", Int(t[1]), "–", Int(t[end]), " s")
    end
    say("collect: ", length(ids), " targets; ", length(seeds), " cell seeds unique and disjoint")
end

# ---------------------------------------------------------------- main

function main(args)
    stage, out = args[1], args[2]
    say("KIM-95 ", stage, " at ", COMMIT, ", job ", JOB, SMOKE ? ", SMOKE" : "",
        ", ", Threads.nthreads(), " threads")
    if stage == "prep"
        run_prep(out)
    elseif stage == "sizing"
        run_sizing(out)
    elseif stage == "run"
        arm = Symbol(args[3])
        arm in ARMS || error("arm $(args[3]) is not one of $ARMS")
        run_chunk(out, arm, parse(Int, args[4]))
    elseif stage == "collect"
        run_collect(out)
    else
        error("unknown stage $stage")
    end
    say("KIM-95 ", stage, " done")
end

if abspath(PROGRAM_FILE) == @__FILE__
    main(ARGS)
end
