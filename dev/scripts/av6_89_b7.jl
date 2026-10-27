# Issue #89 B7: AV6, the mixed posterior against the exact sampler's on M0
# (spec/phases/approximate-method.md §3 AV6; the rules are §12 2026-10-24,
# "B7 as planned", fixed before any run).
#
# Fits, each a directory under <out>:
#   r1 … r50   V7's datasets, m0_dataset(16_200_000 + r): SBC, and against V7;
#   v6         m0_v6, V6's replicate 16085: against V6, beside F10's matched nulls;
#   ks_v6, ks_low1, ks_low2
#              KS: the base fit (v6, and the two V7 replicates with the lowest
#              kcatF_R_ENO truths) refitted with S re-estimated at its posterior
#              mean, δ held;
#   tx_v6, bulk_v6
#              the stream arms: m0_v6 with transcripts only and metabolites only,
#              at v6's reference (θ_ref, δ, S).
#   null_v6, null_low1
#              KS's matched null (§12 2026-10-27, added after the result): v6 and
#              ks_low1's base refitted at the same (θ_ref, δ, S) on the base arm's
#              chain seeds 3 and 4. No prep of their own; chains only.
#
# Stages:
#   prep <out> <fit>       a base fit (r*, v6): the dataset, θ_ref (the bounded MAP at
#                          δ = S = 0), the 800-cell reference ensemble, δ and S;
#   ksprep <out> <fit>     a KS fit: S from 800 cells at the base's posterior mean;
#   chain <out> <fit> <c>  one NUTS chain, NADAPT + NDRAW, from a jittered start;
#   merge <out>            every rule of "B7 as planned", and the TSVs the figures read;
#   nullcmp <out>          W1(base, null) beside KS's W1, to <out>/av6_ks_null.md.
#
# Usage: julia --project dev/scripts/av6_89_b7.jl <stage> <out> [<fit> [<c>]]
# Regenerate: dev/scripts/av6_89_b7.sh <out outside the worktrees>.
# The matched null is not in the launcher. It was submitted by hand after the
# merge, from a worktree at 09be2f5 (jobs 19334375-19334380), as below; <the four>
# is the chain jobs' ids, colon-separated, which --parsable prints
#   s="sbatch --parsable --partition=milan --exclude=dave[148-149] --export=ALL,SMOKE=0
#      --output=<out>/logs/%x_%j.log"
#   for f in null_v6 null_low1; do for c in 1 2; do
#     $s --time=3-00:00:00 --job-name=av6_chain_$f dev/scripts/av6_89_b7.slurm chain <out> $f $c
#   done; done
#   $s --time=01:00:00 --dependency=afterok:<the four> --job-name=av6_nullcmp \
#      dev/scripts/av6_89_b7.slurm nullcmp <out>
# SMOKE=1: 3 cells, a 120 s horizon, 8 reference cells, 10 + 10 NUTS steps, and
# only the replicates in AV6_REPS (default 1:50). AV6_COMPARE=1: see COMPARE.

using InferCell
using InferCell: _mixed_solve
using Random
using Serialization
using Statistics
using LinearAlgebra
using Printf
using Dates
using Distributions: Binomial, cdf, quantile
using Turing: sample, NUTS
using MCMCChains
using MCMCChains: MCMCDiagnosticTools
const MDT = MCMCDiagnosticTools

const SMOKE = get(ENV, "SMOKE", "0") == "1"
# The merge's comparisons with V7 and V6 (rules 4, 5 and KS's V6 floor) are off in a
# smoke, whose datasets are not V7's; AV6_COMPARE=1 runs them anyway, to test the code.
const COMPARE = !SMOKE || get(ENV, "AV6_COMPARE", "0") == "1"
const COMMIT = strip(read(`git rev-parse --short HEAD`, String))
const NCELLS = SMOKE ? 3 : 50
const HORIZON = SMOKE ? 120.0 : M0_HORIZON_S
const NREF = SMOKE ? 8 : 800
envint(k, default) = (v = strip(get(ENV, k, "")); isempty(v) ? default : parse(Int, v))
const NADAPT = envint("NADAPT", SMOKE ? 10 : 200)
const NDRAW = envint("NDRAW", SMOKE ? 10 : 200)
const NCHAIN = 2
const TARGET_ACCEPT = 0.8
const JITTER_SD = 0.25                       # prior SDs per coordinate
const MAP_ITERS = SMOKE ? 5 : 200
# The rules (§12 2026-10-24, "B7 as planned").
const L = 100                                # SBC: draws per replicate after thinning
const RHAT_MAX = 1.02
const ESS_MIN = 100
const KC_BOUND_CPUH = 50.75                  # a tenth of F10's 507.5 CPU-h
const DATA_SEED0 = 16_200_000                # V7's datasets
const V6_SEED = 16_085
const REPS = let v = strip(get(ENV, "AV6_REPS", ""))
    isempty(v) ? (1:50) : eval(Meta.parse(v))
end
const V7_DIR = get(ENV, "V7_DIR",
                   "/fred/oz022/tkimpson/infercell-16a9b-archive-2026-10-17/v7/sbc_v7_16a9b")
const V7_BURN, V7_SWEEPS = 100, 656
const F10_DRAWS = joinpath(@__DIR__, "f10_16a10_draws.tsv")
slurm_job() = get(ENV, "SLURM_ARRAY_JOB_ID", get(ENV, "SLURM_JOB_ID", "none"))

# Each fit's arm, which picks its chain seeds (mixed_chain_seeds' `arm`).
const ARMS = Dict(:base => 0, :ks => 1, :tx => 2, :bulk => 3)

fitdir(out, fit) = joinpath(out, String(fit))
prep_file(out, fit) = joinpath(fitdir(out, fit), "prep.jls")
chain_file(out, fit, c) = joinpath(fitdir(out, fit), "chain_$c.jls")

base_fits() = vcat([Symbol(:r, r) for r in REPS], :v6)
replicate_seed(fit::Symbol) = fit === :v6 ? V6_SEED : DATA_SEED0 + parse(Int, String(fit)[2:end])

# The two V7 replicates with the lowest kcatF_R_ENO truths. The truth is
# m0_dataset's first draw, so it is read without simulating the cells.
const LOWEST_ENO = Ref{Vector{Symbol}}()
function lowest_eno()
    isassigned(LOWEST_ENO) && return LOWEST_ENO[]
    eno = map(REPS) do r
        t = draw_truth(Xoshiro(DATA_SEED0 + r), m0_models(); names = M0_TARGETS,
                       purpose = :calibration)
        Dict(t.values)[:kcatF_R_ENO]
    end
    return LOWEST_ENO[] = [Symbol(:r, REPS[i]) for i in sortperm(eno)[1:2]]
end

# What a fit runs: its kind, the base fit it reads its dataset and reference from.
function fit_kind(fit::Symbol)
    f = String(fit)
    fit === :ks_v6 && return (; kind = :ks, base = :v6)
    fit === :ks_low1 && return (; kind = :ks, base = lowest_eno()[1])
    fit === :ks_low2 && return (; kind = :ks, base = lowest_eno()[2])
    fit === :tx_v6 && return (; kind = :tx, base = :v6)
    fit === :bulk_v6 && return (; kind = :bulk, base = :v6)
    fit === :null_v6 && return (; kind = :null, base = :v6)
    fit === :null_low1 && return (; kind = :null, base = lowest_eno()[1])
    (fit === :v6 || occursin(r"^r\d+$", f)) && return (; kind = :base, base = fit)
    error("unknown fit $fit")
end

# ---------------------------------------------------------------- prep

function run_prep(out, fit)
    fit_kind(fit).kind === :base || error("prep is for base fits; $fit is not one")
    mkpath(fitdir(out, fit))
    t_ds = @elapsed ds = m0_dataset(replicate_seed(fit); ncells = NCELLS, horizon = HORIZON)
    s0 = mixed_setup(ds)
    t_map = @elapsed ref = reference_theta(s0; maxiters = MAP_ITERS)
    values = mixed_theta_values(s0, ref.x)
    m_ref = mean_field_log_bulk(s0, _mixed_solve(s0, exp.(ref.x[1:length(s0.names)])).sol)
    seeds = reference_seeds(ds, NREF)
    t_ens = @elapsed cells = reference_ensemble(values, seeds; horizon = HORIZON)
    d = delta_and_S(cells, m_ref; C = NCELLS)
    serialize(prep_file(out, fit),
              (; fit, ds, ref, values, m_ref, seeds, cells, delta = d.delta, S = d.S,
               commit = COMMIT, smoke = SMOKE, job = slurm_job(),
               seconds = (dataset = t_ds, map = t_map, ensemble = t_ens)))
    @printf("%s prep %s: θ_ref by %s (log density %.3f) in %.0f s; %d reference cells in %.0f s\n",
            now(), fit, ref.method, ref.logdensity, t_map, length(seeds), t_ens)
end

# S re-estimated at the base fit's posterior mean, δ held (KS). The ensemble reuses
# the base's reference seeds, so the two S differ by θ, not by the cells' noise.
function run_ksprep(out, fit)
    k = fit_kind(fit)
    k.kind === :ks || error("ksprep is for KS fits; $fit is not one")
    mkpath(fitdir(out, fit))
    p = deserialize(prep_file(out, k.base))
    X = posterior_draws(out, k.base)
    x̄ = vec(mean(pooled(X); dims = 1))
    s0 = mixed_setup(p.ds)
    values = mixed_theta_values(s0, x̄)
    m_at = mean_field_log_bulk(s0, _mixed_solve(s0, exp.(x̄[1:length(s0.names)])).sol)
    t_ens = @elapsed cells = reference_ensemble(values, p.seeds; horizon = HORIZON)
    S = delta_and_S(cells, m_at; C = NCELLS).S
    serialize(prep_file(out, fit),
              (; fit, base = k.base, xbar = x̄, values, S, delta = p.delta, commit = COMMIT,
               smoke = SMOKE, job = slurm_job(), seconds = (ensemble = t_ens,)))
    @printf("%s ksprep %s (base %s): %d cells at the posterior mean in %.0f s; tr S %.4g against %.4g\n",
            now(), fit, k.base, length(p.seeds), t_ens, tr(S), tr(p.S))
end

# ---------------------------------------------------------------- chains

# The setup and the start's centre a fit's chains run from.
function fit_setup(out, fit)
    k = fit_kind(fit)
    p = deserialize(prep_file(out, k.base))
    S = k.kind === :ks ? deserialize(prep_file(out, fit)).S : p.S
    streams = k.kind === :tx ? (:tx,) : k.kind === :bulk ? (:bulk,) : InferCell.MIXED_STREAMS
    return (; s = with_reference(mixed_setup(p.ds; streams), p.delta, S), p, k)
end

function run_chain(out, fit, c)
    f = fit_setup(out, fit)
    # A null fit reruns its base's setup on the base arm's next NCHAIN seeds.
    seed = f.k.kind === :null ? mixed_chain_seeds(f.p.ds, 2NCHAIN; arm = ARMS[:base])[NCHAIN + c] :
           mixed_chain_seeds(f.p.ds, NCHAIN; arm = ARMS[f.k.kind])[c]
    rng = Xoshiro(seed)
    st = mixed_jittered_start(rng, f.s, f.p.ref.x; scale = JITTER_SD)
    t = @elapsed chain = sample(rng, mixed_model(f.s), NUTS(NADAPT, TARGET_ACCEPT), NDRAW;
                                initial_params = st.x, progress = false)
    mkpath(fitdir(out, fit))                         # the stream arms have no prep of their own
    serialize(chain_file(out, fit, c),
              (; fit, c, seed, start = st.x, start_tries = st.tries, chain, seconds = t,
               nadapt = NADAPT, ndraw = NDRAW, commit = COMMIT, job = slurm_job()))
    @printf("%s chain %s/%d: seed %d, %d + %d steps in %.0f s\n", now(), fit, c, seed, NADAPT,
            NDRAW, t)
end

# ---------------------------------------------------------------- merge

# draws × parameters × chains, on the log scale, in mixed_parameter_names' order.
function posterior_draws(out, fit)
    runs = chains_of(out, fit)
    P = count(nm -> startswith(String(nm), "x["), names(runs[1].chain))
    return cat((Array(r.chain[:, [Symbol("x[$i]") for i in 1:P], 1]) for r in runs)...; dims = 3)
end

chains_of(out, fit) = [deserialize(chain_file(out, fit, c)) for c in 1:NCHAIN]

# A fit's chains pooled: (draws · chains) × parameters.
pooled(X) = reshape(permutedims(X, (1, 3, 2)), :, size(X, 2))

# Rank-normalised split R̂ and bulk ESS per parameter, as MCMCChains reports them.
function diagnostics(X)
    rh = MDT.rhat(permutedims(X, (1, 3, 2)))
    es = MDT.ess(permutedims(X, (1, 3, 2)); kind = :bulk)
    return rh, es
end

# W1 between two samples: ∫ |F_a − F_b|.
function w1(a::AbstractVector, b::AbstractVector)
    xs = sort(vcat(a, b))
    sa, sb = sort(a), sort(b)
    s = 0.0
    for i in 1:(length(xs) - 1)
        Fa = searchsortedlast(sa, xs[i]) / length(sa)
        Fb = searchsortedlast(sb, xs[i]) / length(sb)
        s += abs(Fa - Fb) * (xs[i + 1] - xs[i])
    end
    return s
end

# count(<(i), ranks) for i in 1:L, in one pass over ranks in 0:L.
function ecdf_counts(ranks, L)
    h = zeros(Int, L + 1)
    for r in ranks
        h[r + 1] += 1
    end
    return cumsum(h)[1:L]
end

in_band(counts, lo, hi) = all(lo .<= counts .<= hi)

# V7's simultaneous bands, as sbc_v7_16a9b.jl draws them (same generator), so the
# two methods' rank-ECDFs are read against the same bands.
function ecdf_bands(R, L; α = 0.05, nsim = 50_000, rng = Xoshiro(16_299_999))
    z = (1:L) ./ (L + 1)
    tail(k, p) = (d = Binomial(R, p); 2 * min(cdf(d, k), 1 - cdf(d, k - 1)))
    mins = map(1:nsim) do _
        counts = ecdf_counts(rand(rng, 0:L, R), L)
        minimum(tail(counts[i], z[i]) for i in 1:L)
    end
    γ = quantile(mins, α)
    lo = [quantile(Binomial(R, p), γ / 2) for p in z]
    hi = [quantile(Binomial(R, p), 1 - γ / 2) for p in z]
    return z, lo, hi, γ
end

function truth_x(ds, s)
    tv = Dict(ds.truth.values)
    return vcat([log(tv[n]) for n in s.names], log(ds.sigma))
end

# V7's kept draws for replicate r, columns in mixed_parameter_names' order.
function v7_draws(r, pnames)
    x = deserialize(joinpath(V7_DIR, "replicate$(r).jls"))
    kept = x.state.trace[(V7_BURN + 1):V7_SWEEPS]
    X = reduce(hcat, [vcat(log.(t.cme), t.u, log(t.σ)) for t in kept])'
    vnames = [n === :sigma_b ? :log_sigma_b : Symbol(:log_, n) for n in x.meta.names]
    return X[:, [findfirst(==(n), vnames) for n in pnames]]
end

# V6 (chains 1–4, sweeps 101–1,650) and F10's four matched nulls, from F10's TSV.
function v6_reference(pnames)
    lines = readlines(F10_DRAWS)
    head = Symbol.(split(lines[1], '\t'))
    rows = reduce(hcat, [parse.(Float64, split(l, '\t')) for l in lines[2:end]])'
    col(n) = findfirst(==(n === :log_sigma_b ? :sigma_b : Symbol(String(n)[5:end])), head)
    cols = [col(n) for n in pnames]
    pick(chs, lo, hi) = rows[[r for r in axes(rows, 1)
                              if rows[r, 1] in chs && lo <= rows[r, 2] <= hi], cols]
    last_ = Int(maximum(rows[rows[:, 1] .== 1, 2]))
    ref = pick(1:4, 101, last_)
    sd = vec(std(ref; dims = 1))
    null = [[w1(pick([k], 101, 656)[:, p], pick(setdiff(1:4, k), 101, last_)[:, p]) / sd[p]
             for p in eachindex(pnames)] for k in 1:4]
    return (; ref, sd, null = reduce(hcat, null)')            # null: 4 × parameters
end

# KS's unit and floor for marginal q (rule 6): on m0_v6, V6's SD and F10's largest
# matched null; elsewhere the base's SD and the larger chain-vs-chain W1.
function ks_scale_floor(Xb, Xk, q, v6)
    v6 === nothing || return v6.sd[q], maximum(v6.null[:, q])
    sd = std(pooled(Xb)[:, q])
    return sd, max(w1(Xb[:, q, 1], Xb[:, q, 2]), w1(Xk[:, q, 1], Xk[:, q, 2])) / sd
end

# CPU-h by sacct for the jobs a fit ran (jobs.tsv, written by the launcher). Falls
# back to the stages' own seconds, labelled, if sacct is unavailable.
function sacct_cpuh(out, fit)
    file = joinpath(out, "jobs.tsv")
    isfile(file) || return nothing
    ids = [split(l, '\t')[2] for l in readlines(file) if split(l, '\t')[1] == String(fit)]
    isempty(ids) && return nothing
    total = 0.0
    for id in ids
        txt = try
            read(`sacct -j $id -X -n -P --format=CPUTimeRAW`, String)
        catch
            return nothing
        end
        total += sum(parse(Float64, strip(l)) for l in split(strip(txt), '\n') if !isempty(strip(l)); init = 0.0)
    end
    return total / 3600
end

stage_cpuh(out, fit) = (p = deserialize(prep_file(out, fit));
                        (p.seconds.map + p.seconds.ensemble +
                         sum(r.seconds for r in chains_of(out, fit))) / 3600)

function run_merge(out)
    io = IOBuffer()
    say(args...) = (println(io, args...); println(args...))
    fits = base_fits()
    missing_ = [f for f in fits if !all(c -> isfile(chain_file(out, f, c)), 1:NCHAIN)]
    isempty(missing_) || error("AV6 merges only a complete set (§12 2026-10-24, " *
                               "\"B7 as planned\"); missing chains for: $(join(missing_, ", "))")
    s_any = fit_setup(out, :v6).s
    pnames = mixed_parameter_names(s_any)
    P = length(pnames)
    say("# AV6 (#89 B7): the mixed posterior against the exact sampler's on M0\n")
    say("Merge at `$(COMMIT)` (job $(slurm_job())). SMOKE = $(SMOKE). $(NCELLS) cells, ",
        "$(HORIZON) s, $(NREF) reference cells, $(NCHAIN) chains of $(NADAPT) + $(NDRAW) per fit, ",
        "starts jittered by $(JITTER_SD) prior SD. Rules: §12 2026-10-24, \"B7 as planned\".\n")

    # Per replicate: diagnostics, ranks, against V7, cost.
    reps = [Symbol(:r, r) for r in REPS]
    R = length(reps)
    ranks = zeros(Int, R, P)
    rh = zeros(R, P); es = zeros(R, P)
    d = fill(NaN, R, P); mcse = fill(NaN, R, P); sdratio = fill(NaN, R, P)
    cost = zeros(R); cost_src = String[]
    rows_tsv = String[]
    for (j, f) in enumerate(reps)
        p = deserialize(prep_file(out, f))
        X = posterior_draws(out, f)
        rh[j, :], es[j, :] = diagnostics(X)
        flat = pooled(X)
        idx = round.(Int, range(1, size(flat, 1); length = L))
        tx = truth_x(p.ds, s_any)
        ranks[j, :] = [count(<(tx[q]), flat[idx, q]) for q in 1:P]
        r = parse(Int, String(f)[2:end])
        if isfile(joinpath(V7_DIR, "replicate$(r).jls")) && COMPARE
            V = v7_draws(r, pnames)
            ev = MDT.ess(reshape(V, size(V, 1), 1, P); kind = :bulk)
            sdv = vec(std(V; dims = 1))
            d[j, :] = (vec(mean(flat; dims = 1)) .- vec(mean(V; dims = 1))) ./ sdv
            mcse[j, :] = 1 ./ sqrt.(ev)
            sdratio[j, :] = vec(std(flat; dims = 1)) ./ sdv
        end
        c = sacct_cpuh(out, f)
        cost[j] = c === nothing ? stage_cpuh(out, f) : c
        push!(cost_src, c === nothing ? "stages" : "sacct")
        push!(rows_tsv, join(vcat(String(f), string(p.ds.replicate_seed), string.(ranks[j, :]),
                                  [@sprintf("%.4f", v) for v in rh[j, :]],
                                  [@sprintf("%.1f", v) for v in es[j, :]],
                                  [@sprintf("%.6g", v) for v in d[j, :]],
                                  [@sprintf("%.6g", v) for v in mcse[j, :]],
                                  [@sprintf("%.6g", v) for v in sdratio[j, :]],
                                  @sprintf("%.3f", cost[j]), cost_src[end]), '\t'))
    end
    open(joinpath(out, "av6_replicates.tsv"), "w") do o
        hdr = vcat("fit", "seed", ["rank_$n" for n in pnames], ["rhat_$n" for n in pnames],
                   ["ess_$n" for n in pnames], ["d_$n" for n in pnames], ["mcse_$n" for n in pnames],
                   ["sdratio_$n" for n in pnames], "cpuh", "cpuh_source")
        println(o, join(hdr, '\t'))
        foreach(l -> println(o, l), rows_tsv)
    end

    # Rule 2: convergence, reported, not used to drop.
    bad = [reps[j] for j in 1:R if !(all(rh[j, :] .<= RHAT_MAX) && all(es[j, :] .>= ESS_MIN))]
    say("## Convergence (rule 2: reported, not used to drop)\n")
    say(@sprintf("Worst R̂ %.4f; smallest bulk ESS %.0f; median smallest ESS %.0f.",
                 maximum(rh), minimum(es), median(vec(minimum(es; dims = 2)))))
    say("Below R̂ ≤ $(RHAT_MAX) and ESS ≥ $(ESS_MIN): ", isempty(bad) ? "none" : join(bad, ", "),
        " ($(length(bad)) of $R).\n")

    # Rule 3: SBC, the pass line.
    α = 0.05 / P
    z, lo, hi, γ = ecdf_bands(R, L; α)
    within = falses(P)
    say("## SBC (rule 3, AV6's pass line)\n")
    say("Rank-ECDF of $R replicates, L = $L, against simultaneous bands at 1 − 0.05/$P ",
        "(pointwise γ = $(round(γ; sigdigits = 3))), V7's construction.\n")
    say("| Parameter (log) | Outside the band at | Largest |ECDF − z| | Verdict |")
    say("|---|---|---|---|")
    for q in 1:P
        counts = ecdf_counts(ranks[:, q], L)
        outside = [i for i in 1:L if !(lo[i] <= counts[i] <= hi[i])]
        within[q] = isempty(outside)
        say(@sprintf("| %s | %s | %.3f | %s |", pnames[q],
                     isempty(outside) ? "nowhere" : join(round.(z[outside]; digits = 2), ", "),
                     maximum(abs.(counts ./ R .- z)), within[q] ? "within" : "OUTSIDE"))
    end
    say("\n**AV6's SBC ", all(within) ? "passes" : "FAILS", "**: ", count(within), " of $P within.")
    if !isempty(bad) && length(bad) < R
        keep = [j for j in 1:R if !(reps[j] in bad)]
        _, lo2, hi2, _ = ecdf_bands(length(keep), L; α)
        n2 = count(q -> in_band(ecdf_counts(ranks[keep, q], L), lo2, hi2), 1:P)
        say("Without the $(length(bad)) below the convergence target (sensitivity, not the verdict): ",
            "$n2 of $P within, on $(length(keep)).")
    end
    serialize(joinpath(out, "av6_ranks.jls"), (; pnames, reps, ranks, z, lo, hi, γ, α, within, bad))
    open(joinpath(out, "av6_bands.tsv"), "w") do o
        println(o, "z\tlo\thi")
        for i in 1:L
            println(o, z[i], '\t', lo[i], '\t', hi[i])
        end
    end

    # Rule 5: against V7, information.
    if COMPARE
        say("\n## Against V7 over the $R (rule 5, information)\n")
        say("d_r = (mean_mixed − mean_V7) / sd_V7 per replicate; V7's own MCSE is 1/√ESS_r.\n")
        say("| Parameter (log) | mean d ± SE | median |d| | median V7 MCSE | replicates with |d| > 3 MCSE | median SD ratio |")
        say("|---|---|---|---|---|---|")
        for q in 1:P
            dq = d[:, q]
            say(@sprintf("| %s | %+.3f ± %.3f | %.3f | %.3f | %d | %.2f |", pnames[q], mean(dq),
                         std(dq) / sqrt(R), median(abs.(dq)), median(mcse[:, q]),
                         count(abs.(dq) .> 3 .* mcse[:, q]), median(sdratio[:, q])))
        end
    end

    # Rule 4: m0_v6 against V6, beside F10's matched nulls.
    v6 = COMPARE ? v6_reference(pnames) : nothing
    if COMPARE
        Xv = pooled(posterior_draws(out, :v6))
        say("\n## m0_v6 against V6 (rule 4, information; F10's matched nulls)\n")
        say("| Parameter (log) | mixed mean ± sd | V6 mean ± sd | W1 / SD_V6 | matched nulls | reading |")
        say("|---|---|---|---|---|---|")
        for q in 1:P
            w = w1(Xv[:, q], v6.ref[:, q]) / v6.sd[q]
            mx = maximum(v6.null[:, q])
            say(@sprintf("| %s | %.4f ± %.4f | %.4f ± %.4f | %.3f | %s | %s |", pnames[q],
                         mean(Xv[:, q]), std(Xv[:, q]), mean(v6.ref[:, q]), v6.sd[q], w,
                         join((@sprintf("%.3f", v) for v in v6.null[:, q]), ", "),
                         w > mx ? @sprintf("resolved: above all four, by %.3f", w - mx) : "not resolved"))
        end
    end

    # Rule 6: KS.
    say("\n## KS (rule 6): S re-estimated at the posterior mean, δ held\n")
    say("| Fit | base | Parameter (log) | W1(fixed, re-estimated) | floor | over the floor |")
    say("|---|---|---|---|---|---|")
    fires = false
    for kf in (:ks_v6, :ks_low1, :ks_low2)
        all(c -> isfile(chain_file(out, kf, c)), 1:NCHAIN) || (say("| $kf | — | not run | | | |"); continue)
        base = fit_kind(kf).base
        Xb = posterior_draws(out, base)
        Xk = posterior_draws(out, kf)
        fb, fk = pooled(Xb), pooled(Xk)
        v6k = base === :v6 ? v6 : nothing
        unit = v6k === nothing ? "base SD; floor = larger chain-vs-chain W1" :
               "V6 SD; floor = F10's largest matched null"
        for q in 1:P
            sd, floor_ = ks_scale_floor(Xb, Xk, q, v6k)
            w = w1(fb[:, q], fk[:, q]) / sd
            over = w > floor_
            fires |= over
            say(@sprintf("| %s | %s | %s | %.3f | %.3f (%s) | %s |", kf, base, pnames[q], w, floor_,
                         unit, over ? "yes" : "no"))
        end
    end
    say("\n**KS ", fires ? "fires" : "does not fire", "**.")

    # Rule 7: KC.
    say("\n## KC (rule 7)\n")
    say(@sprintf("CPU-h per posterior (MAP, ensemble, chains; source: %s): median %.2f, min %.2f, max %.2f, against %.2f. ",
                 join(unique(cost_src), "/"), median(cost), minimum(cost), maximum(cost), KC_BOUND_CPUH),
        "**KC ", median(cost) > KC_BOUND_CPUH ? "fires" : "does not fire", "**.")

    # Rule 8: the stream arms, the mixed side (16b.4's exact side is B7.8).
    say("\n## Stream arms on m0_v6 (rule 8; the mixed side)\n")
    say("| Fit | Parameter (log) | mean ± sd | shrinkage 1 − var_post / var_prior | R̂ | ESS |")
    say("|---|---|---|---|---|---|")
    psd = mixed_prior_sd(s_any)
    for sf in (:v6, :tx_v6, :bulk_v6)
        all(c -> isfile(chain_file(out, sf, c)), 1:NCHAIN) || (say("| $sf | not run | | | | |"); continue)
        X = posterior_draws(out, sf)
        r_, e_ = diagnostics(X)
        flat = pooled(X)
        open(joinpath(out, "av6_$(sf)_draws.tsv"), "w") do o
            println(o, join(String.(pnames), '\t'))
            foreach(i -> println(o, join((@sprintf("%.17g", v) for v in flat[i, :]), '\t')), axes(flat, 1))
        end
        for q in 1:P
            say(@sprintf("| %s | %s | %.4f ± %.4f | %.3f | %.4f | %.0f |", sf, pnames[q], mean(flat[:, q]),
                         std(flat[:, q]), 1 - var(flat[:, q]) / psd[q]^2, r_[q], e_[q]))
        end
    end

    write(joinpath(out, "av6_summary.md"), String(take!(io)))
end

# KS's matched null (§12 2026-10-27, added after the result): two fits with the same
# θ_ref, δ and S and different chain seeds, base against null, in KS's units and
# beside KS's own W1 and floor. Information only; KS's verdict stands.
function run_nullcmp(out)
    pnames = mixed_parameter_names(fit_setup(out, :v6).s)
    P = length(pnames)
    v6ref = v6_reference(pnames)
    io = IOBuffer()
    say(xs...) = (println(io, xs...); println(xs...))
    say("# AV6 (#89 B7): KS's matched null\n")
    say("Commit `", COMMIT, "` (job ", slurm_job(), "). Two fits per base, same θ_ref, δ and S, ",
        "different chain seeds. W1 in KS's units.\n")
    say("| Base | Parameter (log) | W1(base, null) | W1(fixed, re-estimated) | KS floor | KS over the null |")
    say("|---|---|---|---|---|---|")
    for (nf, kf) in ((:null_v6, :ks_v6), (:null_low1, :ks_low1))
        base = fit_kind(nf).base
        Xb, Xn, Xk = posterior_draws(out, base), posterior_draws(out, nf), posterior_draws(out, kf)
        fb, fn, fk = pooled(Xb), pooled(Xn), pooled(Xk)
        for q in 1:P
            sd, floor_ = ks_scale_floor(Xb, Xk, q, base === :v6 ? v6ref : nothing)
            wn, wk = w1(fb[:, q], fn[:, q]) / sd, w1(fb[:, q], fk[:, q]) / sd
            say(@sprintf("| %s | %s | %.3f | %.3f | %.3f | %s |", base, pnames[q], wn, wk, floor_,
                         wk > wn ? "yes" : "no"))
        end
    end
    write(joinpath(out, "av6_ks_null.md"), String(take!(io)))
end

function main(args)
    stage, out = args[1], args[2]
    if stage == "prep"
        run_prep(out, Symbol(args[3]))
    elseif stage == "ksprep"
        run_ksprep(out, Symbol(args[3]))
    elseif stage == "chain"
        run_chain(out, Symbol(args[3]), parse(Int, args[4]))
    elseif stage == "nullcmp"
        run_nullcmp(out)
    elseif stage == "merge"
        run_merge(out)
    elseif stage == "lowest_eno"
        println(join(lowest_eno(), " "))
    else
        error("unknown stage $stage")
    end
end

main(ARGS)
