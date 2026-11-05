# KIM-42: the reduced model's uncertainty budget (M5). UNSCORED: there is no truth
# here, so nothing below is a calibration result. The design is the M4 freeze's
# (`93abbce`), with KIM-70's K_d diagnostic (`ed179c2`); this header and KIM-42's §12
# entry declare the rest before the run.
#
# Usage (dev/scripts/reduced_budget_kim42.sh submits every stage):
#   julia --project dev/scripts/reduced_budget_kim42.jl <stage> <outdir> [<index>]
# Stages:
#   width i   i = 0, 1, 2: one `predict` call at prior width ½, 1 or 2, threaded on
#             one node (a crossed matrix cannot be split over calls);
#   report    the tables, in budget.md and budget_rows.tsv;
#   secheck   the SE method's check on made-up matrices (below). No cell runs.
# Outputs, never overwritten, in the order a width writes them: prediction_w<i>.jls (the
# whole `Prediction`, which opens only at this commit), plain_w<i>.jls (`plain`, every
# cell with its minima and accounts) and w<i>.jls (the numbers `report` reads, plain).
#
# The design (frozen):
# - `d11_models()` at 6,300 s (`COREA_CYCLE_S`), the hybrid only, Q = :lactate_out;
# - θ_I: V7 replicate 1, 100 draws evenly spaced, round.(Int, range(1, 556; length =
#   100)), as M0 draws them. Its smallest bulk ESS is 39, so they are about 39
#   effective. Its stored truth is never read;
# - θ_H: K = 20 prior draws at each width, crossed with every row;
# - 4 cells per pair for μ̂, and V's own 4;
# - V: clipping allowed, starved at 0.5, no floors, θ viable at a fraction ≥ 0.5;
# - the ptsG knockdown at dose 0.3, on the pairs V keeps;
# - the headline weights are `per_theta_I`, with `joint` beside them.
#
# Declared here:
# - **Seeds:** rng Xoshiro(19_542_000) at every width (smoke 19_542_900; secheck
#   19_542_500 + its case). So the widths run the same θ_H draws scaled, and the same
#   cells, as M0's V7 and B7 blocks do. `report` checks the seeds are equal. A width's
#   difference from another is then read without most of the cells' noise, but the
#   three do not average to one less noisy estimate.
# - **The headline:** the held share of Var(μ) in the I-first order,
#   E_I[Var(μ | θ_I)] / Var(μ), under `per_theta_I` at width 1, with H first and the
#   gap beside it. ΔQ's split is secondary.
# - **The SE,** the freeze's "a method the script declares": a delete-a-group
#   jackknife over both factors of the crossed matrix, the two variances added.
#   Rows are deleted in 10 blocks of 10 consecutive draws, because V7's chain is
#   autocorrelated. Columns are deleted one at a time. Each deletion recomputes V's
#   weights and the whole `_variance_split`. Bootstrapping is not used: a repeated
#   row or column would repeat cell noise, which the U-statistics assume is
#   independent between entries.
# - **The check** (`secheck`), before the run: on KIM-66's kind of made-up matrix at
#   100 × 20, for each shape, mask and row law, the ratio of the mean SE to the
#   empirical SD of the estimate over 2,000 replicates. **The gate:** if the held
#   share's ratio is outside [0.7, 1.5] in any case, stop and ask the user. It adds
#   a shape in θ_I's expected regime, a share near 3%, which KIM-66 did not cover,
#   and AR(1) rows at ρ = 0.44, the ρ an ESS of 39 in 100 gives.
# - **Descriptive only:** correlations of log θ_I with log θ_H under each weight set
#   (zero by construction under uniform weights), and each parameter's first-order
#   attribution, the correlation of the line means with its log value. Twenty
#   columns cannot separate the held parameters.
# - σ_b is not used: it enters only for a future noisy observation, not for Q.
#
# SMOKE=1 is a plumbing pass, never a result. It runs B7's held-out `v6` fit, at
# 630 s, 4 draws × 3 held × 2 + 2 cells, on 8 threads, with 2 row blocks, and reads
# no evidence posterior. Its V alone has a floor, M_gtp_c ≥ 3,800, so that V rejects
# pairs and those paths run; the run's V has none.

include(joinpath(@__DIR__, "pilot_common.jl"))
using Random
using Printf
using Statistics: std

const STAGE = Symbol(ARGS[1])
const OUT = ARGS[2]
const IDX = length(ARGS) >= 3 ? parse(Int, ARGS[3]) : -1

const V7_R1 = "/fred/oz022/tkimpson/infercell-16a9b-archive-2026-10-17/v7/sbc_v7_16a9b/replicate1.jls"
# As pinned by M0's prep, /fred/oz022/tkimpson/infercell-m5-kim40-2026-11-01/inputs.sha256.
const V7_R1_SHA256 = "d1a4500b59433352076d8b3d28a44df6fde992343f92fb88960ebee4f1ebf759"
const WIDTHS = (0.5, 1.0, 2.0)
const H = SMOKE ? COREA_CYCLE_S / 10 : COREA_CYCLE_S
const N_DRAWS, N_HELD, N_CELLS, N_V, N_ROW_BLOCKS = SMOKE ? (4, 3, 2, 2, 2) : (100, 20, 4, 4, 10)
const SEED = SMOKE ? 19_542_900 : 19_542_000
const KD = Knockdown(; dose = 0.3)
const RULE = ViabilityRule(; label = SMOKE ? "smoke: a floor, so that V rejects" : "M4 freeze, KIM-36",
                           allow_clipped = true, min_viable_fraction = 0.5, starved_fraction = 0.5,
                           floors = SMOKE ? [:M_gtp_c => 3800.0] : Pair{Symbol, Float64}[])
const K_D_THRESHOLDS = (2018.0, 2160.0)     # §12 KIM-70: a lower and an upper bound
const LAWS = (:per_theta_I, :joint)

# ---- files: written in full under tmp/, then moved; never replaced.
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
    startswith(COMMIT_FULL, x.commit) ||
        error("$name was written at $(x.commit), not $COMMIT")
    x.smoke == SMOKE || error("$name was written with SMOKE = $(x.smoke), not $SMOKE")
    return x.result
end

# ---- the run
function posterior_and_rows()
    if SMOKE
        pd = posterior_draws(V6)
    else
        got = first(split(read(`sha256sum $V7_R1`, String)))
        got == V7_R1_SHA256 || error("$V7_R1 has SHA-256 $got, not the one M0 pinned")
        pd = posterior_draws(V7_R1)
        length(pd) == 556 && pd.names == M0_TARGETS ||
            error("V7 replicate 1: $(length(pd)) draws of $(pd.names)")
    end
    return pd, round.(Int, range(1, length(pd); length = N_DRAWS))
end

function kd_row(cells)
    m = [get(c.minima, :M_gtp_c, NaN) for c in cells]
    n = length(m)
    below = [count(<(x), m) for x in K_D_THRESHOLDS]
    return (; n, n_no_minimum = count(isnan, m), below_2018 = below[1], below_2160 = below[2],
            fraction_2018 = n == 0 ? NaN : below[1] / n, fraction_2160 = n == 0 ? NaN : below[2] / n)
end
_all(cells) = reduce(vcat, vec(cells); init = CellRun[])

function width_stage(i)
    0 <= i < length(WIDTHS) || error("width index $i is outside 0…$(length(WIDTHS) - 1)")
    w = WIDTHS[i + 1]
    for name in ("w$i", "plain_w$i", "prediction_w$i")
        isfile(joinpath(OUT, "$name.jls")) && error("$name is already written")
    end
    pd, rows = posterior_and_rows()
    say("width $w: $(length(rows)) draws × $N_HELD held × ($N_CELLS + $N_V) cells, $H s, ",
        "knockdown $(KD.param) at $(KD.dose), rng $SEED")
    sec = @elapsed p = predict(d11_models(), pd; horizon = H, draws = rows,
                               completion = :viability_conditioned_prior, n_held = N_HELD,
                               width = w, viability = RULE, n_cells = N_CELLS,
                               n_viability_cells = N_V, quantity = :lactate_out,
                               intervention = KD, rng = Xoshiro(SEED))
    say("predict took $(round(sec / 3600; digits = 2)) h; statuses $(_sorted(p.statuses))")
    # Saved first, so that nothing below can lose the hours above.
    put("prediction_w$i", p)
    put("plain_w$i", plain(p; cells = true, full = true))
    effect_n = [InferCell._moments(InferCell._paired_delta(p.cells[j], p.cells_intervened[j]))[3]
                for j in CartesianIndices(p.cells)]
    kept = p.viable
    kd = (; baseline_all = kd_row(vcat(_all(p.cells), _all(p.viability_cells))),
          baseline_kept = kd_row(vcat(_all(p.cells[kept]), _all(p.viability_cells[kept]))),
          knockdown = kd_row(_all(p.cells_intervened)))
    reasons = Dict(s => sum(d.counts[s] for d in p.decisions) for s in CELL_STATUSES)
    put("w$i", (; width = w, draws = p.draws, seed = SEED, seeds = p.seeds,
                viability_seeds = p.viability_seeds, theta_I = p.theta_I, held = p.held,
                theta_H = p.theta_H, mu = p.mu, sigma2 = p.sigma2, n_ok = p.n_ok,
                effect = p.effect, effect_sigma2 = p.effect_sigma2, effect_n, viable = p.viable,
                conditioning = p.conditioning, summary = p.summary,
                effect_summary = p.effect_summary, statuses = _sorted(p.statuses),
                statuses_intervened = _sorted(p.statuses_intervened),
                v_cell_counts = _sorted(reasons), kd, seconds = sec))
end

# ---- the numbers
# One matrix's mean and split, with V's weights made from `viable`.
function split_stat(mu, τ, viable, law)
    w = getproperty(viability_weights(viable), law)
    e = InferCell._variance_split(mu, τ, w; held = true)
    ok = isfinite.(mu) .& (w .> 0)
    return (; mean = sum(w[ok] .* mu[ok]) / sum(w[ok]), total = e.total,
            i_between = e.i_first.between, i_within = e.i_first.within,
            h_between = e.h_first.between, h_within = e.h_first.within,
            share_i = e.held_share.i_first, share_h = e.held_share.h_first, gap = e.gap,
            flagged = count(values(e.flagged)), dropped = e.dropped)
end

# The two-way delete-a-group jackknife declared above.
function jackknife_se(mu, τ, viable, law; n_blocks)
    M, K = size(mu)
    at(r, c) = split_stat(mu[r, c], τ[r, c], viable[r, c], law)
    rows = [at(setdiff(1:M, b), 1:K) for b in Iterators.partition(1:M, M ÷ n_blocks)]
    cols = [at(1:M, setdiff(1:K, k)) for k in 1:K]
    jv(x) = (length(x) - 1) / length(x) * sum(abs2, x .- mean(x))
    fields = (:mean, :total, :i_between, :i_within, :h_between, :h_within, :share_i, :share_h, :gap)
    return NamedTuple{fields}(map(f -> sqrt(jv(getfield.(rows, f)) + jv(getfield.(cols, f))), fields))
end

estimate(mu, τ, viable, law; n_blocks) =
    (; est = split_stat(mu, τ, viable, law), se = jackknife_se(mu, τ, viable, law; n_blocks))

log_matrix(θ) = reduce(vcat, [[log(last(p)) for p in r]' for r in θ])

# The weighted correlation of x over the rows with y over the columns, under `w`.
function crossed_corr(x, y, w)
    W = w ./ sum(w)
    rx, cy = vec(sum(W; dims = 2)), vec(sum(W; dims = 1))
    mx, my = sum(rx .* x), sum(cy .* y)
    sx, sy = sqrt(sum(rx .* (x .- mx) .^ 2)), sqrt(sum(cy .* (y .- my) .^ 2))
    return (x' * W * y - mx * my) / (sx * sy)
end

# Each parameter's first-order attribution: the line means against its log value.
function attribution(mu, w, X, names)
    ok = isfinite.(mu)
    W = w .* ok
    mass = vec(sum(W; dims = 2))
    use = mass .> 0
    line = vec(sum(W .* ifelse.(ok, mu, 0.0); dims = 2))[use] ./ mass[use]
    return [begin
                x = X[use, j]
                c = cov_(x, line)
                (; name = names[j], r = c / sqrt(cov_(x, x) * cov_(line, line)), slope = c / cov_(x, x))
            end for j in eachindex(names)]
end
cov_(a, b) = mean((a .- mean(a)) .* (b .- mean(b)))

function report()
    isfile(joinpath(OUT, "budget.md")) && error("budget.md exists; this run never overwrites an output")
    recs = [read_out("w$i") for i in eachindex(WIDTHS) .- 1]
    all(r -> r.seeds == recs[1].seeds && r.viability_seeds == recs[1].viability_seeds, recs) ||
        error("the widths' cell seeds differ, against the declared design")
    io = open(joinpath(OUT, "budget.md"), "w")
    tsv = open(joinpath(OUT, "budget_rows.tsv"), "w")
    out(x...) = (println(io, x...); println(x...))
    println(tsv, join(["width", "target", "law", "quantity", "estimate", "se"], '\t'))
    out("# KIM-42: the reduced model's uncertainty budget (UNSCORED)\n")
    out("At $COMMIT, job $JOB, $(now()). No truth is compared with anything here. ",
        "$N_DRAWS × $N_HELD pairs, $N_CELLS + $N_V cells each, $H s, Q = lactate_out, ",
        "knockdown dose $(KD.dose). SEs: the two-way jackknife the script declares.\n")

    out("## Acceptance and the correlations V induces\n")
    out("| width | acceptance | failed rows | ESS, per_theta_I | ESS, joint | largest abs corr(log θ_I, log θ_H), per_theta_I | joint |")
    out("|---|---|---|---|---|---|---|")
    for r in recs
        vw = viability_weights(r.viable)
        X, Y = log_matrix(r.theta_I), log_matrix(r.theta_H)
        mx(w) = maximum(abs(crossed_corr(X[:, a], Y[:, b], w)) for a in axes(X, 2), b in axes(Y, 2))
        out(@sprintf("| %g | %.4f | %d | %.0f | %.0f | %.3f | %.3f |", r.width, vw.acceptance,
                     length(vw.failed_rows), vw.ess_per_theta_I, vw.ess_joint,
                     mx(vw.per_theta_I), mx(vw.joint)))
    end

    for (target, title) in ((:q, "Q = lactate_out"), (:dq, "ΔQ under the knockdown (secondary)"))
        out("\n## $title\n")
        out("| width | weights | mean ± SE | √Var(μ) ± SE | I first: between, within | H first: between, within | gap | held share, I first | held share, H first |")
        out("|---|---|---|---|---|---|---|---|---|")
        notes = String[]
        for r in recs, law in LAWS
            mu, τ = target === :q ? (r.mu, r.sigma2 ./ r.n_ok) : (r.effect, r.effect_sigma2 ./ r.effect_n)
            b = estimate(mu, τ, r.viable, law; n_blocks = N_ROW_BLOCKS)
            e, s = b.est, b.se
            out(@sprintf("| %g | %s | %.5g ± %.2g | %.4g ± %.2g | %.3g ± %.2g, %.3g ± %.2g | %.3g ± %.2g, %.3g ± %.2g | %.3g ± %.2g | %.3f ± %.3f | %.3f ± %.3f |",
                         r.width, law, e.mean, s.mean, sqrt(e.total), s.total / (2sqrt(e.total)),
                         e.i_between, s.i_between, e.i_within, s.i_within, e.h_between,
                         s.h_between, e.h_within, s.h_within, e.gap, s.gap, e.share_i, s.share_i,
                         e.share_h, s.share_h))
            for f in keys(s)
                println(tsv, join([r.width, target, law, f, getfield(e, f), getfield(s, f)], '\t'))
            end
            e.flagged + e.dropped > 0 &&
                push!(notes, "Width $(r.width), $law: $(e.flagged) components came out negative " *
                             "and are reported as 0; $(e.dropped) pairs dropped.")
        end
        isempty(notes) || out("\n", join(notes, "\n"))
    end

    out("\n## The law of total variance, Q (per_theta_I)\n")
    out("Parameter uncertainty against the cells' own spread: a new cell's variance is Var(μ) + E[σ²].\n")
    out("| width | Var(μ) | E[σ²] | √E[σ²] | Var(μ) / (Var(μ) + E[σ²]) | MC removed, per pair |")
    out("|---|---|---|---|---|---|")
    for r in recs
        s = r.summary.per_theta_I
        t = split_stat(r.mu, r.sigma2 ./ r.n_ok, r.viable, :per_theta_I).total
        out(@sprintf("| %g | %.4g | %.4g | %.4g | %.4f | %.4g |", r.width, t, s.intrinsic_var,
                     sqrt(s.intrinsic_var), t / (t + s.intrinsic_var), s.mc))
    end

    out("\n## Attribution, Q (per_theta_I; descriptive: 20 columns cannot separate the held parameters)\n")
    for r in recs
        w = viability_weights(r.viable).per_theta_I
        X, Y = log_matrix(r.theta_I), log_matrix(r.theta_H)
        a = attribution(r.mu, w, X, first.(r.theta_I[1]))
        c = attribution(r.mu', w', Y, first.(r.theta_H[1]))
        fmt(v) = join([@sprintf("%s r = %.2f", x.name, x.r) for x in sort(v; by = x -> -abs(x.r))[1:min(5, end)]], ", ")
        out("- width $(r.width). θ_I: ", fmt(a), ". θ_H, largest 5: ", fmt(c), ".")
    end

    out("\n## The knockdown's cells, and the cells by status\n")
    out("| width | baseline, unfiltered | knockdown, kept pairs | V's cells | n_failed_intervened, per_theta_I |")
    out("|---|---|---|---|---|")
    for r in recs
        out("| $(r.width) | $(r.statuses) | $(r.statuses_intervened) | $(r.v_cell_counts) | ",
            r.effect_summary.per_theta_I.n_failed_intervened, " |")
    end

    out("\n## K_d diagnostic: `minima[:M_gtp_c]` below 2,018 (a lower bound) and 2,160 (an upper bound)\n")
    out("For KIM-47's model-form caveat; not scored.\n")
    out("| width | cells | n | below 2,018 | fraction | below 2,160 | fraction | no minimum |")
    out("|---|---|---|---|---|---|---|---|")
    for r in recs, (lab, x) in (("baseline, all 8 of every pair", r.kd.baseline_all),
                                ("baseline, kept pairs", r.kd.baseline_kept),
                                ("knockdown, kept pairs", r.kd.knockdown))
        out(@sprintf("| %g | %s | %d | %d | %.4f | %d | %.4f | %d |", r.width, lab, x.n,
                     x.below_2018, x.fraction_2018, x.below_2160, x.fraction_2160, x.n_no_minimum))
    end
    close(io)
    close(tsv)
end

# ---- the SE method's check, on made-up matrices (KIM-66's kind)
const SHAPES = (balanced = (a = 1.0, b = 0.7, g = 0.5, q = 0.3),
                held_dominant = (a = 0.3, b = 1.0, g = 0.3, q = 0.0),
                theta_I_near_3pc = (a = 0.17, b = 1.0, g = 0.1, q = 0.0))
const MASKS = (none = (x, y) -> 1.0, dependent = (x, y) -> 1 / (1 + exp(-(1 + x + y))))
const ROW_LAWS = (iid = 0.0, ar1 = 0.44)

function secheck()
    R = SMOKE ? 20 : 2_000
    M, K, n_cells = 100, 20, 4
    tsv = open(joinpath(OUT, "secheck.tsv"), "w")
    println(tsv, "# commit $COMMIT, job $JOB, $(now()), R = $R, $M × $K, $n_cells cells of sd 1")
    println(tsv, join(["shape", "mask", "rows", "quantity", "mean_estimate", "empirical_sd",
                       "mean_se", "ratio", "within_2se"], '\t'))
    case = 0
    gate = true
    for (sname, s) in pairs(SHAPES), (mname, π) in pairs(MASKS), (rname, ρ) in pairs(ROW_LAWS)
        case += 1
        rng = Xoshiro(19_542_500 + case)
        res = map(1:R) do _
            x = randn(rng, M)
            for m in 2:M
                x[m] = ρ * x[m - 1] + sqrt(1 - ρ^2) * x[m]
            end
            y = randn(rng, K)
            cells = [s.a * x[m] + s.b * y[k] + s.g * x[m] * y[k] + s.q * (x[m]^2 - 1) .+ randn(rng, n_cells)
                     for m in 1:M, k in 1:K]
            v = [rand(rng) < π(x[m], y[k]) for m in 1:M, k in 1:K]
            estimate(mean.(cells), var.(cells) ./ n_cells, v, :per_theta_I; n_blocks = 10)
        end
        for f in (:mean, :total, :i_between, :i_within, :h_between, :h_within, :share_i, :share_h)
            e = [getfield(r.est, f) for r in res]
            se = [getfield(r.se, f) for r in res]
            sd, ratio = std(e), mean(se) / std(e)
            within = mean(abs.(e .- mean(e)) .<= 2se)
            println(tsv, join([sname, mname, rname, f, mean(e), sd, mean(se), ratio, within], '\t'))
            @printf("%-17s %-10s %-4s %-10s mean %8.4f sd %7.4f se %7.4f ratio %.2f within 2 SE %.3f\n",
                    sname, mname, rname, f, mean(e), sd, mean(se), ratio, within)
            f === :share_i && (gate &= 0.7 <= ratio <= 1.5)
        end
        flush(stdout)
    end
    close(tsv)
    println(gate ? "GATE PASSES: the held share's ratio is within [0.7, 1.5] in every case" :
                   "GATE FAILS: stop and ask the user")
end

function main()
    mkpath(OUT)
    say("KIM-42 (unscored), stage $STAGE, at $COMMIT, job $JOB, $(Threads.nthreads()) thread(s), ",
        "GC $(GC_THREADS), on $(RAN_ON.host) ($(RAN_ON.cpu_name)), BLAS threads ",
        get(ENV, "OPENBLAS_NUM_THREADS", "unset"))
    STAGE === :width ? width_stage(IDX) :
    STAGE === :report ? report() :
    STAGE === :secheck ? secheck() : error("unknown stage $STAGE")
    say("done")
end

main()
