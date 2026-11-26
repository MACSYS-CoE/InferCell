# KIM-116 item 8: D3 scores from the plug-in variant's cut draws (prereg/kim-116 §6).
# The plug-in's 2000 cut draws are pushed forward by `predict` exactly as KIM-95 pushed
# the joint: KIM-95's driver is included, so the call (`predict_arm`), the rows
# (`pair_rows`), the targets, the horizon, the cut at COREA_FIT_S and the frozen sizes
# (100 draws evenly spaced, 973 cells a draw, 20 chunks) are KIM-95's own code. The
# truth runs and the prior arm are reused from KIM-95's archive (prereg §6), copied in
# after their SHA256SUMS are checked, so that KIM-95's scorer runs on this directory
# unchanged and its scores sit beside KIM-95's on the same truth and prior.
#
# Usage: julia --project dev/scripts/corea_modular_kim116_predict.jl <stage> <out> [<k>]
# Submit through dev/scripts/corea_modular_kim116_predict.sh, from a worktree at a
# pushed commit.
# Stages:
#   prep <out>      the plug-in's cut draws from MODULAR's plugin/stage2_*.jls, as
#                   PosteriorDraws (chain = the stage-1 draw k, iteration = its inner
#                   draw), with KIM-95's θ*; to prep.jls and prep.md;
#   run <out> <k>   chunk k of the posterior arm, on rng seed PLUGIN_SEED + k, to
#                   rows_posterior_<k>.jls;
#   collect <out>   the chunks joined; every cell seed checked unique and disjoint from
#                   the dataset's, the reference ensemble's and every KIM-95 cell's;
#                   KIM-95's targets.tsv, truth.tsv and predictive_prior_*.tsv copied in
#                   (SHA256-checked); this arm's predictive_posterior_*.tsv, status.tsv
#                   and collect.md written. Then corea_predict_kim95_score.jl <out>.
# SMOKE=1: KIM-95's smoke sizes (a 120 s run cut at 60 s, 2 draws × 2 cells), on the
# KIM-116 smoke pipeline's cut draws and KIM-95's smoke directory: plumbing only.

include(joinpath(@__DIR__, "corea_predict_kim95.jl"))


const MODULAR = get(ENV, "MODULAR", SMOKE ?
    "/fred/oz022/tkimpson/infercell-kim116-smoke-all-7b6c0b1-2026-11-19" :
    "/fred/oz022/tkimpson/infercell-kim116-run-7b6c0b1-2026-11-19")
const K95 = SMOKE ? "/fred/oz022/tkimpson/infercell-kim95-smoke-2026-11-18" :
                    "/fred/oz022/tkimpson/infercell-kim95-archive-2026-11-19"
const KIM92 = "/fred/oz022/tkimpson/infercell-kim92-archive-2026-11-17"
const VARIANT = "plugin"                           # prereg §6: the plug-in only
const K2, NKEEP2 = SMOKE ? (4, 2) : (200, 10)      # prereg §3, as in corea_modular_kim116.jl
const STAGE2 = [:log_kcatF_R_ENO, :log_kcatF_R_FBA, :log_sigma_b]
# predict's rng seeds: a block of their own, beside KIM-95's 17_95x_xxx. The cell seeds
# they give are checked against every KIM-95 cell in `collect`.
const PLUGIN_SEED = 17_961_000
const REUSED = ["targets.tsv", "truth.tsv", "predictive_prior_mu.tsv", "predictive_prior_s2.tsv",
                "predictive_prior_n.tsv"]

# ---------------------------------------------------------------- prep

function run_prep116(out)
    mkpath(out)
    d = joinpath(MODULAR, VARIANT)
    gate = joinpath(d, SMOKE ? "gate_5.md" : "gate_400.md")
    SMOKE || occursin("**The gate passes**", read(gate, String)) ||
        error("$gate does not record a passing gate")
    ds = SMOKE ? corea_dataset(17_200_001; ncells = 2, horizon = 120.0, targets = COREA_TARGETS) :
                 load(KIM92, "dataset")
    names = mixed_parameter_names(mixed_setup(ds; models = MODELS, targets = COREA_TARGETS))
    names == InferCell._log_names(COREA_TARGETS) || error("the setup's coordinates are not COREA_TARGETS'")
    psi = [findfirst(==(n), names) for n in STAGE2]
    phi = setdiff(eachindex(names), psi)
    X = zeros(K2 * NKEEP2, length(names))
    chain, iteration = zeros(Int, size(X, 1)), zeros(Int, size(X, 1))
    for k in 1:K2
        r = load(d, "stage2_$k")
        r.k == k && size(r.z) == (NKEEP2, 3) && length(r.phi) == length(phi) && (SMOKE || r.extended) ||
            error("stage2_$k is not inner chain $k of the extended plug-in fit")
        rows = (k - 1) * NKEEP2 .+ (1:NKEEP2)
        X[rows, phi] .= permutedims(r.phi)
        X[rows, psi] .= r.z
        chain[rows] .= k
        iteration[rows] .= 1:NKEEP2
    end
    all(isfinite, X) || error("a cut draw is not finite")
    k95 = load(K95, "prep")
    isig = findfirst(==(:log_sigma_b), names)
    icol = [findfirst(==(Symbol(:log_, n)), names) for n in COREA_TARGETS]
    pd = PosteriorDraws(COREA_TARGETS, exp.(X[:, icol]), exp.(X[:, isig]), chain, iteration, nothing,
                        k95.posterior.truth,
                        (; kind = :cut, variant = VARIANT, route = :mixed, modular = MODULAR,
                         k2 = K2, nkeep2 = NKEEP2))
    isempty(held_parameters(MODELS, pd).held) || error("the cut draws leave a parameter held")
    save(out, "prep", (; posterior = pd, log_draws = X, names, dataset_seeds = k95.dataset_seeds,
                       reference_seeds = k95.reference_seeds); stage = "prep")
    open(joinpath(out, "prep.md"), "w") do io
        println(io, "# KIM-116 item 8 prep\n\ncommit `", COMMIT, "`, job ", JOB, SMOKE ? ", SMOKE" : "", "\n")
        println(io, "- posterior: the ", VARIANT, " variant's ", size(X, 1), " cut draws (", K2, " stage-1 draws × ",
                NKEEP2, "), from ", d, SMOKE ? "" : " (gate passed: gate_400.md)")
        println(io, "- θ*, the dataset's and the reference ensemble's cell seeds: KIM-95's prep, ", K95)
        println(io, "- truth and prior arm: reused from ", K95, " in `collect`")
    end
    say("prep: ", size(X, 1), " ", VARIANT, " cut draws")
end

# ---------------------------------------------------------------- run

function run_chunk116(out, k)
    pd = load(out, "prep").posterior
    draws = chunk_draws(:posterior, length(pd), k)
    n_cells = N_CELLS[:posterior]
    t = @elapsed p = predict_arm(pd, draws, n_cells, PLUGIN_SEED + k)
    r = pair_rows(p)
    save(out, "rows_posterior_$k", (; r..., arm = :posterior, k, seconds = t, threads = Threads.nthreads(),
                                    cells = 2 * length(draws) * n_cells); stage = "run")
    say(@sprintf("run plug-in chunk %d: %d draws × %d cells × 2 arms in %.0f s on %d threads",
                 k, length(draws), n_cells, t, Threads.nthreads()))
end

# ---------------------------------------------------------------- collect

function sums_of(dir)
    p = joinpath(dir, "SHA256SUMS")
    isfile(p) || return nothing
    return Dict(String(l[67:end]) => String(l[1:64]) for l in readlines(p))
end

function run_collect116(out)
    pp = load(out, "prep")
    K = CHUNKS[:posterior]
    rs = [load(out, "rows_posterior_$k") for k in 1:K]
    # KIM-95's files, copied in after their checksums are checked (the smoke dir has none).
    sums = sums_of(K95)
    SMOKE || sums !== nothing || error("$K95 has no SHA256SUMS")
    for f in REUSED
        src = joinpath(K95, f)
        sums === nothing || first(split(read(`sha256sum $src`, String))) == sums[f] || error("$src fails its SHA256SUMS")
        cp(src, joinpath(out, f); force = true)
    end
    ids = String.(split(readline(joinpath(out, "predictive_prior_mu.tsv")), '\t')[4:end])
    all(r -> r.ids == ids, rs) || error("the chunks' targets are not KIM-95's")
    # Every cell seed: unique here, and disjoint from the dataset's, the reference's and
    # every cell KIM-95 ran (its truth cells above all, which this arm is scored against).
    seeds = reduce(vcat, [r.seeds for r in rs])
    allunique(seeds) || error("a cell seed repeats across chunks")
    isdisjoint(seeds, pp.dataset_seeds) || error("a cell seed is one of the dataset's")
    isdisjoint(seeds, pp.reference_seeds) || error("a cell seed is one of the reference ensemble's")
    n95 = 0
    for arm in ARMS, k in 1:(SMOKE ? 2 : CHUNKS[arm])
        s95 = load(K95, "rows_$(arm)_$k").seeds
        n95 += length(s95)
        isdisjoint(seeds, s95) || error("a cell seed is one of KIM-95's $arm chunk $k")
    end
    pd = pp.posterior
    draws = reduce(vcat, [r.draws for r in rs])
    draws == unique(draws) || error("a draw repeats across chunks")
    o = sortperm(draws)
    lead = [(d, pd.chain[d], pd.iteration[d]) for d in draws[o]]
    for (name, f) in (("mu", r -> r.mu), ("s2", r -> r.s2), ("n", r -> r.n))
        Xs = reduce(vcat, [f(r) for r in rs])[o, :]
        write_wide(joinpath(out, "predictive_posterior_$name.tsv"), ["draw", "chain", "iteration"],
                   (; ids, lead), Xs)
    end
    open(joinpath(out, "status.tsv"), "w") do io
        println(io, "arm\tintervened\tstatus\tcount")
        for (iv, key) in ((false, :statuses), (true, :statuses_intervened))
            dct = Dict{Symbol, Int}()
            for r in rs, (s, c) in getproperty(r, key)
                dct[s] = get(dct, s, 0) + c
            end
            for (s, c) in _sorted(dct)
                println(io, "posterior\t", iv, '\t', s, '\t', c)
            end
        end
    end
    open(joinpath(out, "collect.md"), "w") do io
        println(io, "# KIM-116 item 8 collect\n\ncommit `", COMMIT, "`, job ", JOB, SMOKE ? ", SMOKE" : "", "\n")
        println(io, @sprintf("- posterior (the %s variant's cut draws): %d chunks, %d draws, %d cells (both arms of do(ptsG)), %.1f CPU-h",
                             VARIANT, K, length(draws), sum(r.cells for r in rs),
                             sum(r.seconds * r.threads for r in rs) / 3600))
        println(io, "- ", length(seeds), " cell seeds: unique, disjoint from the dataset's ", length(pp.dataset_seeds),
                ", the reference ensemble's ", length(pp.reference_seeds), " and KIM-95's ", n95)
        println(io, "- reused from ", K95, sums === nothing ? " (no SHA256SUMS: smoke)" : ", SHA256-checked",
                ": ", join(REUSED, ", "))
        println(io, "- ", length(ids), " targets")
    end
    say("collect: ", length(ids), " targets; ", length(seeds), " cell seeds unique and disjoint, ",
        "KIM-95's ", n95, " included")
end

# ---------------------------------------------------------------- main

function main116(args)
    stage, out = args[1], args[2]
    say("KIM-116 item 8 ", stage, " at ", COMMIT, ", job ", JOB, SMOKE ? ", SMOKE" : "",
        ", ", Threads.nthreads(), " threads; modular run ", MODULAR)
    if stage == "prep"
        run_prep116(out)
    elseif stage == "run"
        run_chunk116(out, parse(Int, args[3]))
    elseif stage == "collect"
        run_collect116(out)
    else
        error("unknown stage $stage")
    end
    say("KIM-116 item 8 ", stage, " done")
end

if abspath(PROGRAM_FILE) == @__FILE__
    main116(ARGS)
end
