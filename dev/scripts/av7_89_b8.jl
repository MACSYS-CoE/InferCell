# Issue #89 B8: AV7, the mixed posterior's calibration and efficiency on Core A′
# (spec/phases/approximate-method.md §3 AV7; the rules are §12 2026-10-29,
# "B8 as planned", fixed before any code).
#
# Fits, each a directory under <out>:
#   r1 … r200     the injections, corea_dataset(17_200_000 + r): SBC;
#   tx_r1 … tx_r20, bulk_r1 … bulk_r20
#                 the stream arms on the subset r = 1–20, at their base's
#                 (θ_ref, δ, S);
#   ks_ptsg_lo, ks_eno_lo, ks_ptsg_mid
#                 KS: the base fit chosen by rule from the 200 drawn truths, refitted
#                 with S re-estimated at its posterior mean, δ held;
#   null_ptsg_lo, null_eno_lo, null_ptsg_mid
#                 KS's matched null: the same base refitted at the same (θ_ref, δ, S)
#                 on the base arm's chain seeds 3 and 4.
#
# Stages:
#   prep <out> <fit>       a base fit: the dataset and θ_ref (the bounded MAP at
#                          δ = S = 0);
#   ens <out> <fit> <k>    chunk k of NCHUNK of the 800-cell reference ensemble, at
#                          θ_ref for a base fit and at the base's posterior mean
#                          (log scale) for a KS fit, on the base's reference seeds;
#   ref <out> <fit>        δ and S from the chunks (a KS fit keeps its base's δ);
#   chain <out> <fit> <c>  one NUTS chain, NADAPT + NDRAW, from a jittered start;
#   merge <out> [<bound>]  every rule of "B8 as planned", and the TSVs the figures
#                          read; <bound> is bound_89_b8.jl's output directory;
#   ks_bases               print the three KS bases.
#
# Usage: julia --project dev/scripts/av7_89_b8.jl <stage> <out> [<fit> [<k|c>]]
# Regenerate: dev/scripts/av7_89_b8.sh <out outside the worktrees>.
# SMOKE=1: 3 cells, a 120 s horizon, 8 reference cells in 2 chunks, 10 + 10 NUTS
# steps, and only the injections in AV7_REPS (default 1:200).

using InferCell
using InferCell: _mixed_solve
using Random
using Serialization
using Statistics
using LinearAlgebra
using Printf
using Dates
using Distributions: Binomial, cdf, quantile, params
using Turing: sample, NUTS
using MCMCChains
using MCMCChains: MCMCDiagnosticTools
const MDT = MCMCDiagnosticTools

const SMOKE = get(ENV, "SMOKE", "0") == "1"
const COMMIT = strip(read(`git rev-parse --short HEAD`, String))
const NCELLS = SMOKE ? 3 : 200
const HORIZON = SMOKE ? 120.0 : COREA_CYCLE_S
const NREF = SMOKE ? 8 : 800
const NCHUNK = SMOKE ? 2 : 16
envint(k, default) = (v = strip(get(ENV, k, "")); isempty(v) ? default : parse(Int, v))
const NADAPT = envint("NADAPT", SMOKE ? 10 : 200)
const NDRAW = envint("NDRAW", SMOKE ? 10 : 200)
const NCHAIN = 2
const TARGET_ACCEPT = 0.8
const JITTER_SD = 0.25                       # prior SDs per coordinate
const MAP_ITERS = SMOKE ? 5 : 200
# The rules (§12 2026-10-29, "B8 as planned").
const L = 100                                # SBC: draws per replicate after thinning
const RHAT_MAX = 1.02
const ESS_MIN = 100
const RHO_REACHES = 1.25                     # median posterior SD / bound SD
const DATA_SEED0 = 17_200_000
const REPS = let v = strip(get(ENV, "AV7_REPS", ""))
    isempty(v) ? (1:200) : eval(Meta.parse(v))
end
const SUBSET = [r for r in 1:20 if r in REPS]   # the stream arms and the bound
const KS_FITS = (:ks_ptsg_lo, :ks_eno_lo, :ks_ptsg_mid)
const NULL_FITS = (:null_ptsg_lo, :null_eno_lo, :null_ptsg_mid)
slurm_job() = get(ENV, "SLURM_ARRAY_JOB_ID", get(ENV, "SLURM_JOB_ID", "none"))

# Each fit's arm, which picks its chain seeds (mixed_chain_seeds' `arm`).
const ARMS = Dict(:base => 0, :ks => 1, :tx => 2, :bulk => 3)

fitdir(out, fit) = joinpath(out, String(fit))
prep_file(out, fit) = joinpath(fitdir(out, fit), "prep.jls")
ens_file(out, fit, k) = joinpath(fitdir(out, fit), "ens_$k.jls")
ref_file(out, fit) = joinpath(fitdir(out, fit), "ref.jls")
chain_file(out, fit, c) = joinpath(fitdir(out, fit), "chain_$c.jls")

base_fits() = [Symbol(:r, r) for r in REPS]
replicate_seed(fit::Symbol) = DATA_SEED0 + parse(Int, String(fit)[2:end])

corea_setup(ds; kw...) = mixed_setup(ds; models = d11_models(), targets = D11_TARGETS, kw...)

# The three KS bases, by rule on the drawn truths, read without simulating cells:
# the lowest ptsG truth, the lowest ENO truth (the next lowest if that is the same
# injection), and the ptsG truth nearest its prior median on the log scale (the
# next nearest if that is one of the other two).
const KS_BASES = Ref{Vector{Symbol}}()
function ks_bases()
    isassigned(KS_BASES) && return KS_BASES[]
    ms = d11_models()
    ptsg, eno = promoter_param(:JCVISYN3A_0779), :kcatF_R_ENO
    tv = [Dict(corea_draws(DATA_SEED0 + r; ncells = NCELLS, models = ms).truth.values) for r in REPS]
    pick(order, taken) = REPS[first(i for i in order if !(REPS[i] in taken))]
    lo_p = pick(sortperm([t[ptsg] for t in tv]), Int[])
    lo_e = pick(sortperm([t[eno] for t in tv]), [lo_p])
    pr = only(p for m in ms for p in parameters(m) if p.name === ptsg).prior
    mid = pick(sortperm([abs(log(t[ptsg]) - params(pr)[1]) for t in tv]), [lo_p, lo_e])
    return KS_BASES[] = [Symbol(:r, r) for r in (lo_p, lo_e, mid)]
end

# What a fit runs: its kind, and the base fit it reads its dataset and reference from.
function fit_kind(fit::Symbol)
    f = String(fit)
    i = findfirst(==(fit), KS_FITS)
    i === nothing || return (; kind = :ks, base = ks_bases()[i])
    i = findfirst(==(fit), NULL_FITS)
    i === nothing || return (; kind = :null, base = ks_bases()[i])
    m = match(r"^(tx|bulk)_(r\d+)$", f)
    m === nothing || return (; kind = Symbol(m[1]), base = Symbol(m[2]))
    occursin(r"^r\d+$", f) && return (; kind = :base, base = fit)
    error("unknown fit $fit")
end

# ---------------------------------------------------------------- prep

function run_prep(out, fit)
    fit_kind(fit).kind === :base || error("prep is for base fits; $fit is not one")
    mkpath(fitdir(out, fit))
    t_ds = @elapsed ds = corea_dataset(replicate_seed(fit); ncells = NCELLS, horizon = HORIZON)
    s0 = corea_setup(ds)
    t_map = @elapsed ref = reference_theta(s0; maxiters = MAP_ITERS)
    values = mixed_theta_values(s0, ref.x)
    m_ref = mean_field_log_bulk(s0, _mixed_solve(s0, exp.(ref.x[1:length(s0.names)])).sol)
    seeds = reference_seeds(ds, NREF)
    serialize(prep_file(out, fit),
              (; fit, ds, ref, values, m_ref, seeds, commit = COMMIT, smoke = SMOKE,
               job = slurm_job(), seconds = (dataset = t_ds, map = t_map)))
    @printf("%s prep %s: %d cells in %.0f s; θ_ref by %s (log density %.3f) in %.0f s\n",
            now(), fit, NCELLS, t_ds, ref.method, ref.logdensity, t_map)
end

# The θ a fit's ensemble runs at, with the bulk mean field there: θ_ref for a base
# fit, and the base's posterior mean on the log scale for a KS fit.
function ensemble_at(out, fit)
    k = fit_kind(fit)
    p = deserialize(prep_file(out, k.base))
    k.kind in (:base, :ks) || error("$fit has no ensemble of its own")
    x = k.kind === :base ? p.ref.x : vec(mean(pooled(posterior_draws(out, k.base)); dims = 1))
    s0 = corea_setup(p.ds)
    m_at = mean_field_log_bulk(s0, _mixed_solve(s0, exp.(x[1:length(s0.names)])).sol)
    return (; p, k, x, values = mixed_theta_values(s0, x), m_at)
end

chunk(seeds, k) = seeds[k:NCHUNK:end]

function run_ens(out, fit, k)
    e = ensemble_at(out, fit)
    mkpath(fitdir(out, fit))
    seeds = chunk(e.p.seeds, k)
    t = @elapsed cells = reference_ensemble(e.values, seeds; models = d11_models(),
                                            horizon = HORIZON)
    serialize(ens_file(out, fit, k), (; fit, k, seeds, x = e.x, cells, commit = COMMIT,
                                      job = slurm_job(), seconds = t))
    @printf("%s ens %s/%d: %d cells in %.0f s\n", now(), fit, k, length(seeds), t)
end

function run_ref(out, fit)
    e = ensemble_at(out, fit)
    parts = [deserialize(ens_file(out, fit, k)) for k in 1:NCHUNK]
    all(q -> q.x == e.x, parts) || error("$fit's chunks are at another θ")
    seeds = reduce(vcat, (q.seeds for q in parts))
    sort(seeds) == sort(e.p.seeds) || error("$fit's chunks do not hold its reference seeds")
    order = [findfirst(==(s), seeds) for s in e.p.seeds]
    cells = cat((q.cells for q in parts)...; dims = 1)[order, :, :]
    d = delta_and_S(cells, e.m_at; C = NCELLS)
    δ = e.k.kind === :ks ? deserialize(ref_file(out, e.k.base)).delta : d.delta
    serialize(ref_file(out, fit),
              (; fit, base = e.k.base, x = e.x, values = e.values, m_at = e.m_at, delta = δ,
               S = d.S, commit = COMMIT, job = slurm_job(),
               seconds = (ensemble = sum(q.seconds for q in parts),)))
    @printf("%s ref %s: %d cells; tr S %.4g; |δ|∞ %.3g\n", now(), fit, length(seeds), tr(d.S),
            maximum(abs, δ))
end

# ---------------------------------------------------------------- chains

# The setup and the start's centre a fit's chains run from.
function fit_setup(out, fit)
    k = fit_kind(fit)
    p = deserialize(prep_file(out, k.base))
    rf = deserialize(ref_file(out, k.kind === :ks ? fit : k.base))
    streams = k.kind === :tx ? (:tx,) : k.kind === :bulk ? (:bulk,) : InferCell.MIXED_STREAMS
    return (; s = with_reference(corea_setup(p.ds; streams), rf.delta, rf.S), p, k)
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
    mkpath(fitdir(out, fit))
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
has_chains(out, fit) = all(c -> isfile(chain_file(out, fit, c)), 1:NCHAIN)

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

# V7's simultaneous bands, as sbc_v7_16a9b.jl and B7 draw them (same generator).
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

# Ranks of the truth among L thinned pooled draws.
function sbc_ranks(X, tx)
    flat = pooled(X)
    idx = round.(Int, range(1, size(flat, 1); length = L))
    return [count(<(tx[q]), flat[idx, q]) for q in eachindex(tx)]
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

function stage_cpuh(out, fit)
    p, rf = deserialize(prep_file(out, fit)), deserialize(ref_file(out, fit))
    return (p.seconds.dataset + p.seconds.map + rf.seconds.ensemble +
            sum(r.seconds for r in chains_of(out, fit))) / 3600
end

function write_draws(out, fit, pnames)
    flat = pooled(posterior_draws(out, fit))
    open(joinpath(out, "av7_$(fit)_draws.tsv"), "w") do o
        println(o, join(String.(pnames), '\t'))
        foreach(i -> println(o, join((@sprintf("%.17g", v) for v in flat[i, :]), '\t')), axes(flat, 1))
    end
end

function run_merge(out, bound_dir)
    io = IOBuffer()
    say(args...) = (println(io, args...); println(args...))
    fits = base_fits()
    missing_ = [f for f in fits if !has_chains(out, f)]
    isempty(missing_) || error("AV7 merges only a complete set (§12 2026-10-29, " *
                               "\"B8 as planned\"); missing chains for: $(join(missing_, ", "))")
    s_any = fit_setup(out, first(fits)).s
    pnames = mixed_parameter_names(s_any)
    P = length(pnames)
    psd = mixed_prior_sd(s_any)
    say("# AV7 (#89 B8): the mixed posterior's calibration and efficiency on Core A′\n")
    say("Merge at `$(COMMIT)` (job $(slurm_job())). SMOKE = $(SMOKE). $(NCELLS) cells, ",
        "$(HORIZON) s, $(NREF) reference cells, $(NCHAIN) chains of $(NADAPT) + $(NDRAW) per fit, ",
        "starts jittered by $(JITTER_SD) prior SD. Rules: §12 2026-10-29, \"B8 as planned\".\n")

    # Per injection: diagnostics, ranks, cost.
    R = length(fits)
    ranks = zeros(Int, R, P)
    rh = zeros(R, P); es = zeros(R, P)
    cost = zeros(R); cost_src = String[]
    rows_tsv = String[]
    for (j, f) in enumerate(fits)
        p = deserialize(prep_file(out, f))
        X = posterior_draws(out, f)
        rh[j, :], es[j, :] = diagnostics(X)
        ranks[j, :] = sbc_ranks(X, truth_x(p.ds, s_any))
        c = sacct_cpuh(out, f)
        cost[j] = c === nothing ? stage_cpuh(out, f) : c
        push!(cost_src, c === nothing ? "stages" : "sacct")
        push!(rows_tsv, join(vcat(String(f), string(p.ds.replicate_seed), string(p.ref.method),
                                  string.(ranks[j, :]),
                                  [@sprintf("%.4f", v) for v in rh[j, :]],
                                  [@sprintf("%.1f", v) for v in es[j, :]],
                                  @sprintf("%.3f", cost[j]), cost_src[end]), '\t'))
    end
    open(joinpath(out, "av7_replicates.tsv"), "w") do o
        hdr = vcat("fit", "seed", "theta_ref", ["rank_$n" for n in pnames],
                   ["rhat_$n" for n in pnames], ["ess_$n" for n in pnames], "cpuh", "cpuh_source")
        println(o, join(hdr, '\t'))
        foreach(l -> println(o, l), rows_tsv)
    end

    # Convergence, reported, not used to drop.
    bad = [fits[j] for j in 1:R if !(all(rh[j, :] .<= RHAT_MAX) && all(es[j, :] .>= ESS_MIN))]
    say("## Convergence (reported, not used to drop)\n")
    say(@sprintf("Worst R̂ %.4f; smallest bulk ESS %.0f; median smallest ESS %.0f.",
                 maximum(rh), minimum(es), median(vec(minimum(es; dims = 2)))))
    say("Below R̂ ≤ $(RHAT_MAX) and ESS ≥ $(ESS_MIN): ", isempty(bad) ? "none" : join(bad, ", "),
        " ($(length(bad)) of $R).\n")

    # SBC over all seven, the pass line.
    α = 0.05 / P
    z, lo, hi, γ = ecdf_bands(R, L; α)
    within = falses(P)
    say("## SBC (AV7's pass line)\n")
    say("Rank-ECDF of $R injections, L = $L, against simultaneous bands at 1 − 0.05/$P ",
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
    say("\n**AV7's SBC ", all(within) ? "passes" : "FAILS", "**: ", count(within), " of $P within.")
    if !isempty(bad) && length(bad) < R
        keep = [j for j in 1:R if !(fits[j] in bad)]
        _, lo2, hi2, _ = ecdf_bands(length(keep), L; α)
        n2 = count(q -> in_band(ecdf_counts(ranks[keep, q], L), lo2, hi2), 1:P)
        say("Without the $(length(bad)) below the convergence target (sensitivity, not the verdict): ",
            "$n2 of $P within, on $(length(keep)).")
    end
    serialize(joinpath(out, "av7_ranks.jls"), (; pnames, fits, ranks, z, lo, hi, γ, α, within, bad))
    open(joinpath(out, "av7_bands.tsv"), "w") do o
        println(o, "z\tlo\thi")
        for i in 1:L
            println(o, z[i], '\t', lo[i], '\t', hi[i])
        end
    end

    # The stream arms on the subset: shrinkage, and their rank-ECDFs (drawn, not read).
    say("\n## Stream arms on r = 1–20 (F6; rank-ECDFs drawn, not read as a verdict)\n")
    arm_rows = String[]
    arm_ranks = Dict{Symbol, Matrix{Int}}()
    say("| Arm | Parameter (log) | median shrinkage | range | median R̂ | median ESS |")
    say("|---|---|---|---|---|---|")
    for arm in (:joint, :tx, :bulk)
        fs = [arm === :joint ? Symbol(:r, r) : Symbol(arm, :_r, r) for r in SUBSET]
        ok = [f for f in fs if has_chains(out, f)]
        isempty(ok) && (say("| $arm | not run | | | | |"); continue)
        shr = zeros(length(ok), P); rr = zeros(length(ok), P); ee = zeros(length(ok), P)
        rk = zeros(Int, length(ok), P)
        for (j, f) in enumerate(ok)
            X = posterior_draws(out, f)
            rr[j, :], ee[j, :] = diagnostics(X)
            flat = pooled(X)
            shr[j, :] = 1 .- vec(var(flat; dims = 1)) ./ psd .^ 2
            rk[j, :] = sbc_ranks(X, truth_x(deserialize(prep_file(out, fit_kind(f).base)).ds, s_any))
            push!(arm_rows, join(vcat(String(arm), String(fit_kind(f).base),
                                      [@sprintf("%.6g", v) for v in vec(mean(flat; dims = 1))],
                                      [@sprintf("%.6g", v) for v in vec(std(flat; dims = 1))],
                                      [@sprintf("%.6g", v) for v in shr[j, :]], string.(rk[j, :])), '\t'))
        end
        arm_ranks[arm] = rk
        for q in 1:P
            arm === :tx && pnames[q] === :log_sigma_b && continue    # its prior, not reported
            say(@sprintf("| %s | %s | %.3f | %.3f – %.3f | %.4f | %.0f |", arm, pnames[q],
                         median(shr[:, q]), minimum(shr[:, q]), maximum(shr[:, q]),
                         median(rr[:, q]), median(ee[:, q])))
        end
        length(ok) < length(fs) && say("| $arm | $(length(fs) - length(ok)) of $(length(fs)) not run | | | | |")
    end
    open(joinpath(out, "av7_arms.tsv"), "w") do o
        println(o, join(vcat("arm", "base", ["mean_$n" for n in pnames], ["sd_$n" for n in pnames],
                             ["shrink_$n" for n in pnames], ["rank_$n" for n in pnames]), '\t'))
        foreach(l -> println(o, l), arm_rows)
    end

    # The efficiency reading against the per-truth bound.
    say("\n## Efficiency against the per-truth bound\n")
    if bound_dir === nothing
        say("No bound directory given; not read.")
    else
        eff_rows = String[]
        D = length(D11_TARGETS)
        for (arm, design) in ((:bulk, "metabolites only"), (:joint, "joint"), (:tx, "transcripts only"))
            ρ = fill(NaN, length(SUBSET), D); inf = falses(length(SUBSET), D)
            for (j, r) in enumerate(SUBSET)
                bfile = joinpath(bound_dir, "r$r", "bound.jls")
                f = arm === :joint ? Symbol(:r, r) : Symbol(arm, :_r, r)
                (isfile(bfile) && has_chains(out, f)) || continue
                b = deserialize(bfile)
                b.cols == D11_TARGETS || error("r$r's bound is over $(b.cols)")
                x = b.bound[design]
                sd = vec(std(pooled(posterior_draws(out, f)); dims = 1))[1:D] ./ psd[1:D]
                ρ[j, :] = sd ./ x.ratio
                inf[j, :] = x.ratio .< x.bracket
                push!(eff_rows, join(vcat(String(arm), "r$r", [@sprintf("%.6g", v) for v in sd],
                                          [@sprintf("%.6g", v) for v in x.ratio],
                                          [@sprintf("%.6g", v) for v in x.bracket]), '\t'))
            end
            n_ok = count(j -> !isnan(ρ[j, 1]), axes(ρ, 1))
            say("### $design ($(arm) arm; $(n_ok) of $(length(SUBSET)) with a bound)\n")
            say(arm === :bulk ? "The reading: for each parameter informed at every truth, " *
                                "\"reaches the bound\" if the median ρ ≤ $(RHO_REACHES).\n" :
                "Reported, not read: the transcript rows understate per-cell counts, so ρ < 1 is possible.\n")
            say("| Parameter (log) | informed at | median ρ | range | reading |")
            say("|---|---|---|---|---|")
            for q in 1:D
                v = filter(!isnan, ρ[:, q])
                isempty(v) && continue
                ni = count(inf[:, q])
                reading = arm !== :bulk ? "reported" :
                          ni < n_ok ? "not read: not informed at every truth" :
                          median(v) <= RHO_REACHES ? "reaches the bound" : "does not reach the bound"
                say(@sprintf("| %s | %d of %d | %.3f | %.3f – %.3f | %s |", pnames[q], ni, n_ok,
                             median(v), minimum(v), maximum(v), reading))
            end
            say()
        end
        open(joinpath(out, "av7_efficiency.tsv"), "w") do o
            cols = String.(pnames[1:D])
            println(o, join(vcat("arm", "fit", ["post_sd_ratio_$c" for c in cols],
                                 ["bound_$c" for c in cols], ["bracket_$c" for c in cols]), '\t'))
            foreach(l -> println(o, l), eff_rows)
        end
    end

    # KS, against its matched null.
    say("\n## KS: S re-estimated at the posterior mean, δ held, against the matched null\n")
    say("| Fit | base | Parameter (log) | W1(fixed, re-estimated) / SD | W1(fixed, null) / SD | over the null |")
    say("|---|---|---|---|---|---|")
    fires = false
    for (kf, nf) in zip(KS_FITS, NULL_FITS)
        base = fit_kind(kf).base
        (base in fits && has_chains(out, kf) && has_chains(out, nf)) ||
            (say("| $kf | $base | not run | | | |"); continue)
        fb, fk, fn = (pooled(posterior_draws(out, f)) for f in (base, kf, nf))
        for q in 1:P
            sd = std(fb[:, q])
            wk, wn = w1(fb[:, q], fk[:, q]) / sd, w1(fb[:, q], fn[:, q]) / sd
            fires |= wk > wn
            say(@sprintf("| %s | %s | %s | %.3f | %.3f | %s |", kf, base, pnames[q], wk, wn,
                         wk > wn ? "yes" : "no"))
        end
    end
    say("\n**KS ", fires ? "fires" : "does not fire", "** (reported, not a stop; §8).")

    # Cost.
    say("\n## Cost per posterior (Table 2)\n")
    say(@sprintf("CPU-h per posterior (dataset, MAP, ensemble, chains; source: %s): median %.1f, min %.1f, max %.1f. ",
                 join(unique(cost_src), "/"), median(cost), minimum(cost), maximum(cost)),
        "The exact sampler is projected at 60k–150k CPU-h on Core A′ (framing note §1).")

    for f in vcat(fits[1:min(end, 1)], [f for f in KS_FITS if has_chains(out, f)])
        write_draws(out, f, pnames)
    end
    write(joinpath(out, "av7_summary.md"), String(take!(io)))
end

function main(args)
    stage = args[1]
    if stage == "ks_bases"
        println(join(ks_bases(), " "))
        return
    end
    out = args[2]
    if stage == "prep"
        run_prep(out, Symbol(args[3]))
    elseif stage == "ens"
        run_ens(out, Symbol(args[3]), parse(Int, args[4]))
    elseif stage == "ref"
        run_ref(out, Symbol(args[3]))
    elseif stage == "chain"
        run_chain(out, Symbol(args[3]), parse(Int, args[4]))
    elseif stage == "merge"
        run_merge(out, length(args) >= 3 ? args[3] : nothing)
    else
        error("unknown stage $stage")
    end
end

main(ARGS)
