# Probabilistic trajectories, for paper 1's figure. No new inference and no new scored
# claim: two stages replay published runs with `trajectory` on, and must reproduce their
# endpoints exactly, and one releases 62 fixed constants to their priors, unscored.
# Design: docs/superpowers/specs/2026-11-07-probabilistic-trajectories-design.md. Its §12
# entry fixes everything below before the run.
#
# Usage (dev/scripts/trajectories_predict.sh submits every stage):
#   julia --project dev/scripts/trajectories_predict.jl <stage> <outdir> [<index>]
# Stages:
#   a i          M0's replay on V7 replicate 1, exact posterior, single-core: i = 1…10 is
#                block i of M0's run, and i = 0 is its truth. Baseline only;
#   b            the reduced model's replay of KIM-42's width 1, baseline only;
#   c_precheck   the release, 4 draws × 5 held × (4 + 4) cells, for its acceptance;
#   c            the release, 40 draws × 50 held × (4 + 4) cells, under M4's frozen rule;
#   data         the regenerated M0 dataset 1, and its truth against the fit's;
#   replay       the exact checks on a and b, in replay.md, and stops on a mismatch;
#   report       the same checks, the TSVs and checks.md.
# Outputs, never overwritten: <stage>.jls files (whole `Prediction`s, which open only at
# this commit), and tsv/. Each is written under tmp/ and then moved into place.
#
# Seeds, from M5's range. Panel A blocks 19_540_000 + (b − 1), truth 19_541_501, M0's
# own for dataset r = 1. Panel B 19_542_000, KIM-42's. Panel C 19_543_000, its pre-check
# 19_543_001, the new-bulk band 19_543_002. The smoke pass 19_543_900 on.
#
# SMOKE=1 is a plumbing pass, never a result: B7's held-out v6 fit, 630 s on the reduced
# model and 60 s on M0, 4 draws × 3 held × 2 + 2 cells, a floor on V so that it rejects,
# and every check reported without stopping. `data` runs in full, since it costs little
# and reads only the V7 file.

include(joinpath(@__DIR__, "pilot_common.jl"))
using InferCell: AbstractSubModel, LAC_EXPORT_METER
using Random
using Printf
using Statistics: median, std

const STAGE = Symbol(ARGS[1])
const OUT = ARGS[2]
const IDX = length(ARGS) >= 3 ? parse(Int, ARGS[3]) : -1

const V7_R1 = "/fred/oz022/tkimpson/infercell-16a9b-archive-2026-10-17/v7/sbc_v7_16a9b/replicate1.jls"
const V7_R1_SHA256 = "d1a4500b59433352076d8b3d28a44df6fde992343f92fb88960ebee4f1ebf759"
const ARCHIVE = "/fred/oz022/tkimpson/infercell-predict-archive-2026-11-06/runs"
const DATASET_SEED = 16_200_001

const M_ = LAC_EXPORT_METER
const SPECIES = [M_, :M_lac__L_e, :M_gtp_c, :M_pyr_c]
const BULK_SPECIES = [:M_lac__L_e, :M_gtp_c]        # in M0_PANEL; M_pyr_c is not
const LEVELS = (0.5, 0.9)
const EVERY = 60.0

const H_M0 = SMOKE ? 60.0 : M0_HORIZON_S
const H_R = SMOKE ? COREA_CYCLE_S / 10 : COREA_CYCLE_S
# M0: draws, blocks, cells per draw, and the truth's cells.
const A_DRAWS, A_BLOCKS, A_CELLS, A_TRUTH = SMOKE ? (4, 2, 3, 4) : (100, 10, 887, 9_849)
# The reduced model's replay: draws, held, cells.
const B_DRAWS, B_HELD, B_CELLS = SMOKE ? (4, 3, 2) : (100, 20, 4)
# The release: draws, held, cells and V's own, and the pre-check's draws and held.
const C_DRAWS, C_HELD, C_CELLS, C_V = SMOKE ? (4, 3, 2, 2) : (40, 50, 4, 4)
const P_DRAWS, P_HELD = SMOKE ? (2, 2) : (4, 5)
const N_NEW_BULK = SMOKE ? 200 : 20_000
const NCELLS_BULK = 50

const SEED_A(b) = SMOKE ? 19_543_900 + b : 19_540_000 + (b - 1)
const SEED_TRUTH = SMOKE ? 19_543_910 : 19_541_501
const SEED_B = SMOKE ? 19_543_911 : 19_542_000
const SEED_C = SMOKE ? 19_543_912 : 19_543_000
const SEED_P = SMOKE ? 19_543_913 : 19_543_001
const SEED_NEW_BULK = SMOKE ? 19_543_914 : 19_543_002

const RULE = ViabilityRule(; label = SMOKE ? "smoke: a floor, so that V rejects" : "M4 freeze, KIM-36",
                           allow_clipped = true, min_viable_fraction = 0.5, starved_fraction = 0.5,
                           floors = SMOKE ? [:M_gtp_c => 3800.0] : Pair{Symbol, Float64}[])

# ---- files
function put(name, x)
    dest = joinpath(OUT, "$name.jls")
    isfile(dest) && error("$dest exists; this run never overwrites an output")
    tmp = joinpath(OUT, "tmp")
    mkpath(tmp)
    save(tmp, name, x; stage = string(STAGE))
    mv(joinpath(tmp, "$name.jls"), dest)
    say("wrote $dest")
end
function read_out(name)
    x = deserialize(joinpath(OUT, "$name.jls"))
    startswith(COMMIT_FULL, x.commit) || error("$name was written at $(x.commit), not $COMMIT")
    x.smoke == SMOKE || error("$name was written with SMOKE = $(x.smoke), not $SMOKE")
    return x.result
end
# A published run's plain record, which opens at any commit.
archived(parts...) = deserialize(joinpath(ARCHIVE, parts...)).result

function posterior()
    SMOKE && return posterior_draws(V6)
    got = first(split(read(`sha256sum $V7_R1`, String)))
    got == V7_R1_SHA256 || error("$V7_R1 has SHA-256 $got, not the one M0 pinned")
    pd = posterior_draws(V7_R1)
    length(pd) == 556 && pd.names == M0_TARGETS ||
        error("V7 replicate 1: $(length(pd)) draws of $(pd.names)")
    return pd
end
even_draws(pd, n) = round.(Int, range(1, length(pd); length = n))
block_draws(pd, b) = even_draws(pd, A_DRAWS)[((b - 1) * (A_DRAWS ÷ A_BLOCKS) + 1):(b * (A_DRAWS ÷ A_BLOCKS))]

# ---- the runs
function stage_a(i)
    0 <= i <= A_BLOCKS || error("a's index $i is outside 0…$A_BLOCKS")
    ms = collect(AbstractSubModel, m0_models())
    if i == 0
        truth = (SMOKE ? posterior_draws(V6) : posterior()).truth
        θ = reshape([truth[n] for n in M0_TARGETS], 1, :)
        pd = PosteriorDraws(M0_TARGETS, θ, nothing, [1], [1], nothing, truth,
                            (kind = :truth, dataset = 1))
        p = predict(ms, pd; horizon = H_M0, n_cells = A_TRUTH, rng = Xoshiro(SEED_TRUTH),
                    threads = false, trajectory = SPECIES, every = EVERY)
        put("a_truth", p)
    else
        pd = posterior()
        p = predict(ms, pd; horizon = H_M0, n_cells = A_CELLS, draws = block_draws(pd, i),
                    rng = Xoshiro(SEED_A(i)), threads = false, trajectory = SPECIES, every = EVERY)
        put("a_block_$i", p)
    end
end

function stage_b()
    pd = posterior()
    p = predict(d11_models(), pd; horizon = H_R, draws = even_draws(pd, B_DRAWS),
                completion = :prior, n_held = B_HELD, n_cells = B_CELLS, width = 1,
                rng = Xoshiro(SEED_B), trajectory = SPECIES, every = EVERY)
    say("b: statuses $(_sorted(p.statuses))")
    put("b", p)
end

# The released composition: every fixed kcatF, kcatR and Km of the two metabolic modules,
# with ENO's and FBA's constants, which d11_models frees already.
function released_models()
    kin(m) = [x.name for x in parameters(m)
              if x.fixed && any(s -> startswith(String(x.name), s), ("kcatF_", "kcatR_", "km_"))]
    free_now = [:kcatF_R_ENO, :kcatR_R_ENO, :kcatF_R_FBA, :kcatR_R_FBA]
    ms = d11_models()
    ms[1] = CentralGlycolysis(enzymes = :translated,
                              free = vcat(free_now, kin(CentralGlycolysis(enzymes = :translated))))
    ms[3] = NucleotideRecycling(enzymes = :translated,
                                free = kin(NucleotideRecycling(enzymes = :translated)))
    return ms
end

function check_released(ms, pd)
    h = held_parameters(metered(ms), pd)
    n = 5 + length(h.held)
    n == 82 || error("the release draws $n parameters, not the fixed 82")
    rev = filter(x -> startswith(String(x), "kcatR_"), h.held)
    isempty(rev) || error("the release would draw reverse constants: $rev")
    return h
end

function stage_c(name, draws, held, seed)
    pd = posterior()
    ms = released_models()
    check_released(ms, pd)
    p = predict(ms, pd; horizon = H_R, draws = even_draws(pd, draws),
                completion = :viability_conditioned_prior, viability = RULE, n_held = held,
                n_cells = C_CELLS, n_viability_cells = C_V, width = 1, rng = Xoshiro(seed),
                trajectory = SPECIES, every = EVERY)
    say("$name: acceptance $(p.conditioning.per_theta_I.acceptance), statuses $(_sorted(p.statuses))")
    put(name, p)
end

function stage_data()
    ds = m0_dataset(DATASET_SEED; ncells = 50)
    fit = posterior_draws(V7_R1).truth
    dev = maximum(abs(v / fit[n] - 1) for (n, v) in ds.truth.values if n in M0_TARGETS)
    put("data", (; times = ds.times, bulk = ds.bulk, panel = collect(M0_PANEL), sigma = ds.sigma,
                 truth = ds.truth.values, truth_rel_dev = dev))
end

# ---- the report
struct Checks
    io::IO
    failed::Vector{String}
end
function check!(c::Checks, name, ok::Bool, detail = "")
    println(c.io, "- ", ok ? "PASS" : "FAIL", ": ", name, isempty(detail) ? "" : " ($detail)")
    flush(c.io)
    ok || push!(c.failed, name)
    say(ok ? "PASS " : "FAIL ", name, " ", detail)
    return ok
end

first_diff(a, b) = size(a) != size(b) ? "sizes $(size(a)), $(size(b))" :
                   (i = findfirst(j -> !isequal(a[j], b[j]), eachindex(a));
                    i === nothing ? "" : "first at $i: $(a[i]) against $(b[i])")

function published_truth()
    rows = split.(readlines(joinpath(ARCHIVE, "m0-kim40-kim41", "score_datasets.tsv")), '\t')
    hdr = rows[1]
    col(n) = findfirst(==(n), hdr)
    r = only(x for x in rows[2:end]
             if x[col("target")] == "q" && x[col("dataset")] == "1" && x[col("posterior")] == "v7" &&
                x[col("arm")] == "hybrid")
    return parse(Float64, r[col("truth")])
end

# Every exact check on the replays. A mismatch stops the work.
function replay_checks!(c::Checks)
    fin = archived("m0-kim40-kim41", "finish", "000_r1_v7.jls")
    for b in 1:A_BLOCKS
        p = read_out("a_block_$b")
        rows = ((b - 1) * (A_DRAWS ÷ A_BLOCKS) + 1):(b * (A_DRAWS ÷ A_BLOCKS))
        arch = fin.hybrid.mu[rows, :]
        check!(c, "panel A block $b: μ̂ equals M0's published row", isequal(p.mu, arch), first_diff(p.mu, arch))
    end
    t = read_out("a_truth")
    pub = published_truth()
    check!(c, "panel A truth: the mean equals M0's published truth", t.mu[1, 1] == pub,
           "$(t.mu[1, 1]) against $pub")
    b = read_out("b")
    w1 = archived("budget-kim42", "w1.jls")
    check!(c, "panel B: μ̂ equals KIM-42's width 1", isequal(b.mu, w1.mu), first_diff(b.mu, w1.mu))
    check!(c, "panel B: the seeds equal KIM-42's width 1", b.seeds == w1.seeds)
    check!(c, "panel B: θ_H equals KIM-42's width 1", b.theta_H == w1.theta_H)
    return b, w1
end

function stage_replay()
    dest = joinpath(OUT, "replay.md")
    isfile(dest) && error("$dest exists; this run never overwrites an output")
    open(dest, "w") do io
        println(io, "# Replay checks at $COMMIT, job $JOB, $(now())\n")
        c = Checks(io, String[])
        replay_checks!(c)
        report_precheck(io)
        finish_checks(c)
    end
end

function finish_checks(c::Checks)
    println(c.io, isempty(c.failed) ? "\nALL PASS" : "\nFAILED: $(c.failed)")
    flush(c.io)
    !isempty(c.failed) && !SMOKE && error("checks failed: $(c.failed)")
end

function report_precheck(io)
    p = read_out("c_precheck")
    println(io, "\n## The release's pre-check (panel C goes ahead whatever this says)\n")
    println(io, "$(size(p.mu)) pairs, acceptance $(p.conditioning.per_theta_I.acceptance), ",
            "statuses $(_sorted(p.statuses)), drawn parameters $(5 + length(p.held))")
end

function band_rows(io, panel, arm, species, b, wname)
    for (layer, L) in pairs((; rows = b.rows, pairs = b.pairs, cells = b.cells)), i in eachindex(b.t)
        println(io, join((panel, arm, species, layer, wname, b.t[i], L.mean[i], L.lo[0.5][i],
                          L.hi[0.5][i], L.lo[0.9][i], L.hi[0.9][i]), '\t'))
    end
end
const BAND_HEADER = join(("panel", "arm", "species", "layer", "weights", "t", "mean", "lo50", "hi50",
                          "lo90", "hi90"), '\t')
_q_like(s) = s === M_     # the meter's baseline is its start, so that it reads Q(t)

# A new bulk measurement: the mean of 50 cells, resampled from a posterior row's, with
# that row's σ_b, as `observe_bulk` makes one.
function new_bulk_rows(io, a, s)
    rng = Xoshiro(SEED_NEW_BULK)
    t, vals = InferCell._trajectory_values(a, s, false)
    M = size(a.mu, 1)
    for i in eachindex(t)
        t[i] > 0 || continue
        x = map(1:N_NEW_BULK) do _
            m = rand(rng, 1:M)
            X = vals[m, 1]
            c = rand(rng, 1:size(X, 1), NCELLS_BULK)
            max(mean(view(X, c, i)), 1.0) * exp(a.sigma_b[m] * randn(rng))
        end
        q(v) = InferCell._weighted_quantile(x, ones(length(x)), v)
        println(io, join(("a", "posterior", s, "new_bulk", "uniform", t[i], mean(x), q(0.25), q(0.75),
                          q(0.05), q(0.95)), '\t'))
    end
end

function stage_report()
    dest = joinpath(OUT, "checks.md")
    isfile(dest) && error("$dest exists; this run never overwrites an output")
    tsv = joinpath(OUT, "tsv")
    mkpath(tsv)
    open(dest, "w") do io
        println(io, "# Trajectories: checks at $COMMIT, job $JOB, $(now())\n")
        c = Checks(io, String[])
        b, w1 = replay_checks!(c)

        a = combine_rows([read_out("a_block_$i") for i in 1:A_BLOCKS])
        t = read_out("a_truth")
        cp = read_out("c")
        data = read_out("data")

        # The dataset, and the rescaling of its lactate onto Q's axis.
        check!(c, "the regenerated dataset's truth equals replicate 1's", data.truth_rel_dev < 1e-12,
               "largest relative difference $(data.truth_rel_dev)")
        ratios = Float64[]
        for cell in t.cells[1, 1][1:min(10, end)]
            tr = cell.trajectory
            dq = tr.values[2:end, 1] .- tr.values[1, 1]
            dl = tr.values[2:end, 2] .- tr.values[1, 2]
            append!(ratios, dq ./ dl)
        end
        R = median(ratios)
        spread = maximum(abs.(ratios ./ R .- 1))
        check!(c, "export over M_lac__L_e's change is constant", spread < 1e-9,
               "R = $R, largest relative deviation $spread, over 10 cells")
        println(io, "\nR = $R\n")

        # TSVs: panel A.
        open(joinpath(tsv, "traj_a.tsv"), "w") do f
            println(f, BAND_HEADER)
            for s in SPECIES
                band_rows(f, "a", "posterior", s, trajectory_bands(a, s; baseline = _q_like(s), levels = LEVELS), "uniform")
            end
            foreach(s -> new_bulk_rows(f, a, s), BULK_SPECIES)
        end
        open(joinpath(tsv, "traj_a_truth.tsv"), "w") do f
            println(f, join(("species", "t", "mean", "se"), '\t'))
            for s in SPECIES
                tt, vals = InferCell._trajectory_values(t, s, _q_like(s))
                X = vals[1, 1]
                for i in eachindex(tt)
                    println(f, join((s, tt[i], mean(view(X, :, i)), std(view(X, :, i)) / sqrt(size(X, 1))), '\t'))
                end
            end
        end
        open(joinpath(tsv, "traj_a_data.tsv"), "w") do f
            println(f, join(("species", "t", "value"), '\t'))
            for s in BULK_SPECIES, (i, x) in enumerate(data.times)
                println(f, join((s, x, data.bulk[findfirst(==(s), data.panel), i]), '\t'))
            end
        end

        # Panel B, and the split over time.
        open(joinpath(tsv, "traj_b.tsv"), "w") do f
            println(f, BAND_HEADER)
            for s in SPECIES
                band_rows(f, "b", "held20", s, trajectory_bands(b, s; baseline = _q_like(s), levels = LEVELS), "uniform")
            end
        end
        sp = trajectory_split(b, M_; baseline = true)
        open(joinpath(tsv, "traj_split.tsv"), "w") do f
            println(f, join(("t", "held_share_i_first", "held_share_h_first", "total", "mc", "flagged",
                             "mc_undefined", "dropped"), '\t'))
            foreach(i -> println(f, join((sp.t[i], sp.held_share.i_first[i], sp.held_share.h_first[i],
                                          sp.total[i], sp.mc[i], sp.flagged[i], sp.mc_undefined[i],
                                          sp.dropped[i]), '\t')), eachindex(sp.t))
        end

        # Panel C. The conditioned arms exist only if V kept a pair.
        M, K = size(cp.mu)
        kept = count(cp.viable)
        open(joinpath(tsv, "traj_c.tsv"), "w") do f
            println(f, BAND_HEADER)
            for s in SPECIES
                bd(w) = trajectory_bands(cp, s; baseline = _q_like(s), levels = LEVELS, weights = w)
                band_rows(f, "c", "released82", s, bd(ones(M, K)), "unconditioned")
                kept > 0 || continue
                band_rows(f, "c", "released82", s, bd(:per_theta_I), "per_theta_I")
                band_rows(f, "c", "released82", s, bd(:joint), "joint")
            end
        end

        # What the result note reports.
        println(io, "\n## Layer 2's Monte Carlo share at the end of the run (panel B)\n")
        sb = b.summary.uniform
        println(io, "mc = $(sb.mc), var = $(sb.var), mc / var = $(sb.mc / sb.var)")
        println(io, "final held share, I first: $(sp.held_share.i_first[end]); KIM-42's headline at width 1 is 0.908")

        println(io, "\n## Panel C\n")
        acc = kept / (M * K)
        statuses = _sorted(cp.statuses)
        n_all = sum(last, statuses)
        starved = cp.statuses[:starved] / n_all
        med(w) = InferCell._weighted_quantile(vec(cp.mu), vec(w), 0.5)
        unc = med(ones(M, K))
        println(io, "$M × $K pairs, drawn parameters $(5 + length(cp.held)), kept pairs $kept, acceptance $acc")
        println(io, "per_theta_I: ", cp.conditioning.per_theta_I)
        println(io, "joint: ", cp.conditioning.joint)
        println(io, "statuses of μ̂'s cells: $statuses; starved fraction $starved")
        println(io, "median of pair means at the end: panel B ", median(vec(b.mu)), ", released unconditioned ", unc,
                kept > 0 ? ", released per_theta_I $(med(cp.weights.per_theta_I)), joint $(med(cp.weights.joint))" : "")
        acc < 0.1 && println(io, "\nBELOW 10% KEPT: the figure drops the filled conditioned band and shows the unconditioned one, labelled.")
        finish_checks(c)
    end
end

function main()
    mkpath(OUT)
    say("trajectories, stage $STAGE $IDX, at $COMMIT, job $JOB, $(Threads.nthreads()) thread(s), ",
        "GC $(GC_THREADS), on $(RAN_ON.host) ($(RAN_ON.cpu_name))")
    STAGE === :a ? stage_a(IDX) :
    STAGE === :b ? stage_b() :
    STAGE === :c_precheck ? stage_c("c_precheck", P_DRAWS, P_HELD, SEED_P) :
    STAGE === :c ? stage_c("c", C_DRAWS, C_HELD, SEED_C) :
    STAGE === :data ? stage_data() :
    STAGE === :replay ? stage_replay() :
    STAGE === :report ? stage_report() : error("unknown stage $STAGE")
    say("done")
end

main()
