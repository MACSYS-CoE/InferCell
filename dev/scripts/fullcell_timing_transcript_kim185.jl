# KIM-185 (Linear M11, full-cell timing 3): the transcript factor on 455 genes, by the
# data plug-in. It times one density and one dual pass of the gene-expression block over
# the full syn3A model's 455 mRNA genes, serial and threaded over genes, and projects the
# block's cost per posterior against the threshold of the Decision log entry of
# 2026-11-26 (closing KIM-176): at most 5,000 CPU-h.
#
# What "the transcript factor" is here. src/'s exact factor, `TranscriptFactor` and
# `transcript_loglik` (src/transcript_factor.jl), called as is. Its rates are the data
# plug-in's (KIM-116, `plugin_loglik` in src/mixed_likelihood.jl): each gene's
# transcription constant is the reduced model's rate law, `rate_constants` of a
# `CoreATranscription` (src/organisms/coreA/transcription.jl), at the gene's promoter
# strength S_g and at nucleotide pools fixed from the observed panel. `plugin_loglik` is
# wired to the reduced model's 17-gene `MixedSetup`, so the 455-gene case is built here,
# not in src: one single-gene `CoreATranscription` per gene, and one single-gene
# `TranscriptFactor` per gene and cell. Nothing in src changes.
#
# The block, per gene g: log S_g, plus one log krnadeg shared by every gene, with
# μ_g = krnadeg / len_g, as in the reduced model. 455 genes give 456 free coordinates.
# Given krnadeg, the genes decouple: gene g's term depends on log S_g and log krnadeg only.
#
# The AD, fixed before any timing was read (Linear comment on KIM-185):
#   per-gene   the headline. ForwardDiff, one dual pass per gene with Chunk{2} (∂/∂log S_g
#              and ∂/∂log krnadeg), summed. Exact, and the value comes out of the same pass.
#   naive      for comparison: ForwardDiff over all 456 coordinates at its default chunk
#              (pickchunksize(456) = 12), on the whole factor. One chunk pass is timed and
#              multiplied by the number of passes (ceil(456 / 12) = 38); it is labelled so.
# The reduced model's sampler (KIM-116) used neither: it called ForwardDiff.gradient over
# 18 coordinates at the default chunk (two passes of 9).
#
# Inputs.
#   genes      /fred/oz022/tkimpson/syn3a-ensemble/data/model/cme/{reactions,species}.csv:
#              each of the 455 transcription reactions that makes an mRNA gives the
#              gene's length (the ATP_trsc increment), its four base counts (ATP_mRNA,
#              CTP_mRNA, GTP_mRNA, UTP_mRNA), its first-minute rate constant, and its
#              mRNA's initial count (procid 1). First two bases: AU for every gene (the
#              ATG start; all 17 reduced genes are AU). The reference promoter strength
#              S_ref is the first-minute rate over the same law at S = 1, published base
#              mapping and the first-minute concentrations (MinCell_CMEODE.py's
#              minute 1, which the export's `rate` column is).
#   pools      the reduced model's plug-in states on KIM-92's dataset
#              (`panel_plugin_states` at θ_ref, archive infercell-kim92-archive-2026-11-17):
#              ATP and GTP from the observed panel, CTP and UTP the chemostats, at t = 0
#              and each of the 52 switches.
#   design     the reduced model's: 200 cells, saves every 60 s to 3120 s (52 windows),
#              KIM-92's dataset (corea_pilot_kim92.jl: ncells = 200, fit_to = COREA_FIT_S).
#   truth      log S_g drawn from its prior, Normal(log S_ref, log 2) (the reduced model's
#              promoter prior), seed SEED; krnadeg at KIM-92's truth.
#   counts     exact: each window is 59 s at k[:, m] then 1 s at k[:, m + 1], and a
#              linear birth-death step of length Δ from n is
#              Binomial(n, e^{−μΔ}) + Poisson(k (1 − e^{−μΔ}) / μ). Every cell starts at
#              the same y0, procid 1's counts, as the reduced model's cells share theirs.
#   support    each gene's factor is sized as `mixed_setup` sizes the reduced model's:
#              kmax = TURNOVER_CEILING / (len − 1), μmin = krnadeg prior's 1e-6 quantile / len.
#
# Usage: julia --project dev/scripts/fullcell_timing_transcript_kim185.jl <stage> <out>
# Submit through dev/scripts/fullcell_timing_transcript_kim185.sh.
# Stages:
#   data <out>      the 455-gene inputs, rates and counts, to data.jls and data.md, with
#                   the checks on the 17 reduced genes;
#   serial <out>    on one thread: the 17-gene checks against src's plug-in (density and
#                   gradient) and its timing on this node; then the 455-gene density and
#                   per-gene dual pass, median of NREP, and the naive chunk pass. serial.jls/md;
#   threaded <out>  the 455-gene density and per-gene dual pass with W workers over genes
#                   for each W in WORKERS, median of NREP. threaded.jls/md;
#   project <out>   gradients per posterior counted from KIM-116's plug-in chains, and the
#                   projection against 5,000 CPU-h. projection.md.
# SMOKE=1: 3 cells, 12 genes, 3 windows, 2 repeats. A plumbing pass, never a result.

include(joinpath(@__DIR__, "pilot_common.jl"))
using InferCell: _parameter_by_name, TranscriptionGene, rate_constants
using ForwardDiff
using Distributions
using Random
using Printf
using Statistics: median, quantile
using SHA: sha256
using LinearAlgebra: BLAS
using MCMCChains

const DR = ForwardDiff.DiffResults
const CME = "/fred/oz022/tkimpson/syn3a-ensemble/data/model/cme"
const SYN3A = "/fred/oz022/tkimpson/syn3a-ensemble"
const ARCHIVE92 = "/fred/oz022/tkimpson/infercell-kim92-archive-2026-11-17"
const ARCHIVE116 = "/fred/oz022/tkimpson/infercell-kim116-archive-2026-11-26"
const SEED = 18_500_001
const NGENES = 455
const NREP = SMOKE ? 2 : 5
const NREP_NAIVE = SMOKE ? 1 : 3
const SMOKE_SIZE = (cells = 3, genes = 12, windows = 3)
const THRESHOLD_CPUH = 5_000.0          # Decision log 2026-11-26, closing KIM-176
const STAGE2 = [:log_kcatF_R_ENO, :log_kcatF_R_FBA, :log_sigma_b]
const KIM116_PLUGIN_STAGE1_CPUH = 165.4 # KIM-116 result note, cost table
const KIM116_PLUGIN_GRADIENT_S = 11.26  # KIM-116 profile, job 21196497, at the truth

sha(path) = bytes2hex(open(sha256, path))

# ---------------------------------------------------------------- the 455 genes

function read_csv(path)
    lines = readlines(path)
    header = split(lines[1], ',')
    rows = [split(l, ',') for l in lines[2:end] if !isempty(l)]
    all(r -> length(r) == length(header), rows) || error("$path: a row does not have $(length(header)) fields")
    return header, rows
end

# Each mRNA gene's numbers from the full model's CME tables.
function full_model_genes()
    _, sp = read_csv(joinpath(CME, "species.csv"))
    name = Dict(parse(Int, r[1]) => String(r[2]) for r in sp)
    count0 = Dict(String(r[2]) => parse(Int, r[3]) for r in sp)
    idx = Dict(v => k for (k, v) in name)
    _, rx = read_csv(joinpath(CME, "reactions.csv"))
    genes = []
    for r in rx
        st = Dict{String, Int}()
        for pr in split(r[5], ';')
            a, b = split(pr, ':')
            st[name[parse(Int, a)]] = parse(Int, b)
        end
        haskey(st, "ATP_trsc") || continue
        mrna = [s for (s, v) in st if v > 0 && startswith(s, "M_RNA_")]
        isempty(mrna) && continue
        m = only(mrna)
        locus = Symbol(m[7:(end - 2)])            # M_RNA_<locus>_c
        len = st["ATP_trsc"]
        counts = (A = st["ATP_mRNA"], C = st["CTP_mRNA"], G = st["GTP_mRNA"], U = st["UTP_mRNA"])
        sum(counts) == len || error("$locus: base counts $(counts) do not sum to its length $len")
        push!(genes, (; locus, len, counts, rate1 = parse(Float64, r[3]), y0 = count0[m],
                      reactant = name[parse(Int, r[4])]))
    end
    length(genes) == NGENES || error("found $(length(genes)) mRNA genes, not $NGENES")
    allunique(g.locus for g in genes) || error("a locus is transcribed twice")
    return genes
end

# The gene as the reduced model's rate law reads it. S_ref is the first-minute rate over
# the law at S = 1 (ptn_count = PROMOTER_DIVISOR), published mapping, first-minute
# concentrations; valid while RNAPOL_KCAT · S_ref < TURNOVER_CEILING, which is checked.
function as_transcription_gene(g)
    unit = TranscriptionGene(g.locus, g.locus, g.len, g.counts, (:A, :U), PROMOTER_DIVISOR, 1.0)
    k_unit = transcription_rate_constant(unit, InferCell.SETUP_CONSTANTS; base_mapping = :published)
    s_ref = g.rate1 / k_unit
    RNAPOL_KCAT * s_ref < InferCell.TURNOVER_CEILING ||
        error("$(g.locus): S_ref = $s_ref puts kcat_mod at the ceiling; the inversion does not hold")
    tg = TranscriptionGene(g.locus, g.locus, g.len, g.counts, (:A, :U),
                           max(1, round(Int, s_ref * PROMOTER_DIVISOR)), Float64(g.y0))
    return tg, s_ref
end

# The reduced model's setup and its plug-in states on KIM-92's dataset.
function reduced()
    ds = load(ARCHIVE92, "dataset")
    p = load(ARCHIVE92, "prep")
    txo = with_reference(mixed_setup(ds; models = d11_models(), targets = COREA_TARGETS,
                                     streams = (:tx,)), p.delta, p.S)
    states = panel_plugin_states(txo, p.ref.x)
    r = only(rb for rb in txo.route.rebuilds if rb.model isa CoreATranscription)
    r.species == collect(TRANSCRIPTION_POOLS) || error("the rebuild's pools are $(r.species)")
    pools = [[k == 0 ? r.pool_held[j] : u[k] for (j, k) in enumerate(r.pool_idxs)] for u in states]
    tv = Dict(ds.truth.values)
    x_truth = vcat([log(tv[n]) for n in txo.names], log(ds.sigma))
    krnadeg_prior = _parameter_by_name(txo.models, :krnadeg).prior
    return (; ds, p, txo, states, rebuild = r, pools, x_truth, krnadeg_prior,
            krnadeg = tv[:krnadeg], times = txo.switches)
end

# One linear birth-death step of length Δ from n, exactly.
bd_step(rng, n, k, μ, Δ) = (q = exp(-μ * Δ); rand(rng, Binomial(n, q)) + rand(rng, Poisson(k * (1 - q) / μ)))

function run_data(out)
    mkpath(out)
    R = reduced()
    genes = full_model_genes()
    conv = map(as_transcription_gene, genes)
    tgenes, sref = first.(conv), last.(conv)

    # Check 1: the 17 reduced genes, as the full model's tables give them, are the reduced
    # model's extract: the same lengths and base counts, and S_ref = ptn_count / 180.
    red = read_transcription_genes()
    chk = map(red) do g
        i = findfirst(t -> t.locus == g.locus, tgenes)
        i === nothing && error("reduced gene $(g.locus) is not among the 455")
        t = tgenes[i]
        (t.length == g.length && t.counts == g.counts) ||
            error("$(g.locus): full-model length/counts $(t.length)/$(t.counts), extract $(g.length)/$(g.counts)")
        (; locus = g.locus, sref = sref[i], extract = g.ptn_count / PROMOTER_DIVISOR,
         rel = abs(sref[i] / (g.ptn_count / PROMOTER_DIVISOR) - 1))
    end
    say(@sprintf("check 1: the 17 reduced genes match the extract's lengths and base counts; S_ref against ptn_count/180, max rel. diff %.2e",
                 maximum(c.rel for c in chk)))

    G = SMOKE ? SMOKE_SIZE.genes : NGENES
    C = SMOKE ? SMOKE_SIZE.cells : size(R.ds.transcripts, 1)
    M = SMOKE ? SMOKE_SIZE.windows : length(R.ds.times)
    keep = 1:G
    times, pools = R.times[1:(M + 1)], R.pools[1:(M + 1)]
    mods = [CoreATranscription(genes = [tgenes[g]]) for g in keep]
    lens = [tgenes[g].length for g in keep]

    # Check 2: the single-gene module's rate is the 17-gene rebuild's, at KIM-92's truth.
    θ = exp.(R.x_truth[1:length(R.txo.names)])
    pr = mean_field_parameters(R.txo.route, InferCell._mixed_writes(R.txo, θ))
    k17 = transcription_rates_from_states(R.txo, pr, R.states)
    S17 = Dict(l => θ[findfirst(==(promoter_param(l)), R.txo.names)] for l in R.txo.loci)
    rel2 = maximum(enumerate(R.txo.loci)) do (i, l)
        m1 = CoreATranscription(genes = [red[findfirst(g -> g.locus == l, red)]])
        k1 = [rate_constants([S17[l]], t, m1, P)[1] for (t, P) in zip(R.times, R.pools)]
        maximum(abs.(k1 ./ k17[i, :] .- 1))
    end
    say(@sprintf("check 2: single-gene rate_constants against the 17-gene plug-in rebuild at the truth, max rel. diff %.2e", rel2))

    # The truth and the rates.
    rng = Xoshiro(SEED)
    logS_ref = log.(sref[keep])
    logS = logS_ref .+ log(2.0) .* randn(rng, G)
    kd = R.krnadeg
    μ = kd ./ lens
    k = [rate_constants([exp(logS[g])], t, mods[g], P)[1] for g in 1:G, (t, P) in zip(times, pools)]
    kmax = InferCell.TURNOVER_CEILING ./ (lens .- 1)
    μmin = quantile(R.krnadeg_prior, 1e-6) ./ lens
    all(k .<= kmax) && all(μ .>= μmin) || error("a truth rate leaves the factor's support")

    # The counts: cells × genes × windows, exactly.
    y0 = [genes[g].y0 for g in keep]
    y = zeros(Int, C, G, M)
    for c in 1:C, g in 1:G
        n = y0[g]
        for m in 1:M
            n = bd_step(rng, n, k[g, m], μ[g], 59.0)
            n = bd_step(rng, n, k[g, m + 1], μ[g], 1.0)
            y[c, g, m] = n
        end
    end

    # The caps the factor sizes for each gene and cell.
    caps = [only(TranscriptFactor([y0[g]], y[c, g:g, :]; kmax = kmax[g], μmin = μmin[g]).caps)
            for c in 1:C, g in 1:G]
    tables = Dict(f => sha(joinpath(CME, f)) for f in ("reactions.csv", "species.csv"))
    syn3a_commit = strip(read(`git -C $SYN3A rev-parse HEAD`, String))
    full_deg = Dict{Symbol, Float64}()
    _, rx = read_csv(joinpath(CME, "reactions.csv"))
    _, sp = read_csv(joinpath(CME, "species.csv"))
    nm = Dict(parse(Int, r[1]) => String(r[2]) for r in sp)
    for r in rx
        occursin("ATP_mRNAdeg", join(nm[parse(Int, split(p, ':')[1])] for p in split(r[5], ';'))) || continue
        full_deg[Symbol(nm[parse(Int, r[4])][7:(end - 2)])] = parse(Float64, r[3])
    end
    deg_ratio = [full_deg[tgenes[g].locus] / μ[g] for g in keep]
    d = (; loci = [tgenes[g].locus for g in keep], genes = tgenes[keep], sref = sref[keep], lens, y0,
         logS_ref, logS, krnadeg = kd, logkd_prior = (μ = params(R.krnadeg_prior)[1], σ = params(R.krnadeg_prior)[2]),
         μ, k, kmax, μmin, y, caps, times, pools, check1 = chk, check2 = rel2, tables, syn3a_commit,
         deg_ratio, saves = R.ds.times[1:M], seed = SEED)
    save(out, "data", d; stage = "data")
    write_data_md(out, d)
end

function write_data_md(out, d)
    io = IOBuffer()
    p(xs...) = println(io, xs...)
    C, G, M = size(d.y)
    p("# KIM-185 item 1: synthetic transcripts for the full model's genes\n")
    p("Commit `", COMMIT, "`, job ", JOB, ", on ", RAN_ON.host, " (", RAN_ON.cpu, "), Julia ", RAN_ON.julia,
      ". SMOKE = ", SMOKE, ".\n")
    p("- **Counts `y`:** ", C, " cells × ", G, " genes × ", M, " windows (`", summary(d.y), "`), saves every 60 s from ",
      d.saves[1], " to ", d.saves[end], " s. Initial counts `y0`: ", length(d.y0), " genes, shared by every cell.")
    p("- **Rates `k`:** ", size(d.k, 1), " × ", size(d.k, 2), " (genes × windows + 1); death rates `μ`: ", length(d.μ), ".")
    p("- **Free coordinates:** ", G, " log promoter strengths + log krnadeg = ", G + 1, ".")
    p("- **Truth:** log S drawn from Normal(log S_ref, log 2), seed ", d.seed, "; krnadeg = ", @sprintf("%.6g", d.krnadeg),
      " (KIM-92's truth). S / S_ref: median ", @sprintf("%.3g", median(exp.(d.logS .- d.logS_ref))), ".")
    p("- **Gene lengths:** median ", median(d.lens), ", range ", minimum(d.lens), "–", maximum(d.lens), " nt.")
    p("- **S_ref:** median ", @sprintf("%.3g", median(d.sref)), ", range ", @sprintf("%.3g–%.3g", extrema(d.sref)...), ".")
    p("- **Birth rate k (truth):** median ", @sprintf("%.3g", median(d.k)), " /s, range ", @sprintf("%.3g–%.3g", extrema(d.k)...),
      ". Counts: mean ", @sprintf("%.3g", sum(d.y) / length(d.y)), ", max ", maximum(d.y), "; y0 mean ",
      @sprintf("%.3g", sum(d.y0) / length(d.y0)), ", max ", maximum(d.y0), ".")
    cq = quantile(vec(d.caps), [0.0, 0.5, 0.9, 0.99, 1.0])
    p("- **Caps** over genes × cells: min ", cq[1], ", median ", cq[2], ", 90% ", cq[3], ", 99% ", cq[4], ", max ", cq[5], ".")
    p("- **Death rate:** the reduced model's law, krnadeg / len. The full model's own first-order mRNA degradation constant over it: median ",
      @sprintf("%.3g", median(d.deg_ratio)), ", range ", @sprintf("%.3g–%.3g", extrema(d.deg_ratio)...), " (reported, not used).")
    p("- **Check 1:** the 17 reduced genes in the full model's tables have the extract's lengths and base counts, and S_ref matches ptn_count/180 to a max relative difference of ",
      @sprintf("%.2e", maximum(c.rel for c in d.check1)), ".")
    p("- **Check 2:** a single-gene `CoreATranscription`'s `rate_constants` reproduces the 17-gene plug-in rebuild at KIM-92's truth to a max relative difference of ",
      @sprintf("%.2e", d.check2), ".")
    p("- **Tables:** syn3a-ensemble commit `", d.syn3a_commit, "`; sha256 ", join(["`$k` $(v)" for (k, v) in sort(collect(d.tables))], ", "), ".")
    write(joinpath(out, "data.md"), String(take!(io)))
    print(read(joinpath(out, "data.md"), String))
end

# ---------------------------------------------------------------- the block

# Gene g's term over every cell: its single-gene factors, module, length, and the pools.
struct GeneBlock
    fs::Vector{TranscriptFactor}
    mod::CoreATranscription
    len::Int
    times::Vector{Float64}
    pools::Vector{Vector{Float64}}
end

# Gene g's log-likelihood at z = (log S_g, log krnadeg), by src's factor.
function gene_ll(b::GeneBlock, z::AbstractVector)
    S, kd = exp(z[1]), exp(z[2])
    k = reshape([rate_constants([S], t, b.mod, P)[1] for (t, P) in zip(b.times, b.pools)], 1, :)
    return transcript_loglik(b.fs, k, [kd / b.len])
end

function blocks(d; cells = axes(d.y, 1))
    G = size(d.y, 2)
    return map(1:G) do g
        fs = [TranscriptFactor([d.y0[g]], d.y[c, g:g, :]; kmax = d.kmax[g], μmin = d.μmin[g]) for c in cells]
        GeneBlock(fs, CoreATranscription(genes = [d.genes[g]]), d.lens[g], d.times, d.pools)
    end
end

# The prior on x = (log S_1..G, log krnadeg) and its gradient.
log_prior(d, x) = sum(logpdf(Normal(m, log(2.0)), xi) for (m, xi) in zip(d.logS_ref, x[1:(end - 1)])) +
                  logpdf(Normal(d.logkd_prior.μ, d.logkd_prior.σ), x[end])
grad_prior(d, x) = vcat(-(x[1:(end - 1)] .- d.logS_ref) ./ log(2.0)^2, -(x[end] - d.logkd_prior.μ) / d.logkd_prior.σ^2)

# Each gene's term on `workers` tasks pulling genes off a shared counter, heaviest first.
# `f(g)` returns a tuple; results land in gene order, so the sum is the same at any W.
function over_genes(f, G, workers, order)
    out = Vector{Any}(undef, G)
    if workers == 1
        for g in order
            out[g] = f(g)
        end
        return out
    end
    next = Threads.Atomic{Int}(1)
    @sync for _ in 1:workers
        Threads.@spawn while true
            i = Threads.atomic_add!(next, 1)
            i > length(order) && break
            g = order[i]
            out[g] = f(g)
        end
    end
    return out
end

function density(bs, d, x; workers = 1, order = eachindex(bs))
    v = over_genes(g -> gene_ll(bs[g], [x[g], x[end]]), length(bs), workers, order)
    return sum(Float64.(v)) + log_prior(d, x)
end

# One dual pass per gene (Chunk{2}): the value and the gradient of the whole block.
function value_gradient(bs, d, x; workers = 1, order = eachindex(bs), pergene = nothing)
    G = length(bs)
    r = over_genes(G, workers, order) do g
        z = [x[g], x[end]]
        f = w -> gene_ll(bs[g], w)
        res = DR.GradientResult(z)
        t = @elapsed ForwardDiff.gradient!(res, f, z, ForwardDiff.GradientConfig(f, z, ForwardDiff.Chunk{2}()))
        (DR.value(res), DR.gradient(res)[1], DR.gradient(res)[2], t)
    end
    pergene === nothing || (pergene .= [q[4] for q in r])
    val = sum(q[1] for q in r) + log_prior(d, x)
    grad = vcat([q[2] for q in r], sum(q[3] for q in r)) .+ grad_prior(d, x)
    return val, grad
end

# The whole factor at once: one TranscriptFactor per cell over every gene, as src's
# 17-gene plug-in holds them, at x of any number type.
function whole_factors(d; cells = axes(d.y, 1))
    [TranscriptFactor(d.y0, d.y[c, :, :]; kmax = d.kmax, μmin = d.μmin) for c in cells]
end
function whole_ll(fs, bs, d, x)
    G = length(bs)
    k = reduce(vcat, [permutedims([rate_constants([exp(x[g])], t, bs[g].mod, P)[1]
                                   for (t, P) in zip(d.times, d.pools)]) for g in 1:G])
    return transcript_loglik(fs, k, exp(x[end]) ./ d.lens)
end

# One naive chunk pass: every coordinate a Dual of width N, the first N seeded.
function naive_pass(fs, bs, d, x, ::Val{N}) where {N}
    T = ForwardDiff.Tag(naive_pass, Float64)
    xd = [ForwardDiff.Dual{typeof(T)}(x[i], ForwardDiff.Partials(ntuple(j -> Float64(i == j), N))) for i in eachindex(x)]
    return whole_ll(fs, bs, d, xd)
end

function reps(f, n)
    map(1:n) do _
        GC.gc()
        st = @timed f()
        (; seconds = st.time, gib = st.bytes / 2^30, gc = st.gctime, value = st.value)
    end
end
med(v, s = :seconds) = median([getproperty(r, s) for r in v])

function environment()
    (; threads = Threads.nthreads(), blas = BLAS.get_num_threads(), gc = GC_THREADS,
     forwarddiff = string(pkgversion(ForwardDiff)), host = RAN_ON.host, cpu = RAN_ON.cpu, julia = RAN_ON.julia)
end

# ---------------------------------------------------------------- serial

function check_reduced()
    R = reduced()
    txo, names = R.txo, mixed_parameter_names(R.txo)
    x = R.x_truth
    cells = SMOKE ? (1:SMOKE_SIZE.cells) : eachindex(txo.factors)
    # A smoke keeps the first few cells, in a copy of the setup with only their factors.
    s1 = !SMOKE ? txo :
         InferCell.MixedSetup(txo.models, txo.route, txo.factors[cells], txo.loci, txo.lengths,
                              txo.tx_slots, txo.rows, txo.saves, txo.switches, txo.logz, txo.delta,
                              txo.S, txo.floor, txo.names, txo.forwards, txo.prior_mu, txo.prior_sd,
                              txo.sigma_prior, txo.abstol, txo.reltol, txo.maxiters, txo.box_sd,
                              txo.streams)
    plug = modular_stage(s1, setdiff(names, STAGE2), x; loglik = (s, y) -> plugin_loglik(s, R.states, y))
    f = z -> stage_logdensity(plug, z)
    z0 = stage_x0(plug)
    # The 17 genes as this script's per-gene blocks.
    red = read_transcription_genes()
    θ = exp.(x[1:length(txo.names)])
    pos(n) = findfirst(==(n), names)
    bs = map(enumerate(txo.loci)) do (g, l)
        fs = [TranscriptFactor([txo.factors[c].y0[g]], txo.factors[c].y[g:g, :];
                               kmax = txo.factors[c].kmax[g], μmin = txo.factors[c].μmin[g]) for c in cells]
        all(only(fs[i].caps) == txo.factors[c].caps[g] for (i, c) in enumerate(cells)) || error("a cap differs")
        GeneBlock(fs, CoreATranscription(genes = [red[findfirst(q -> q.locus == l, red)]]), txo.lengths[g],
                  txo.switches, R.pools)
    end
    iS = [pos(Symbol(:log_, promoter_param(l))) for l in txo.loci]
    ik = pos(:log_krnadeg)
    prior(z) = sum(logpdf(Normal(m, s), zi) for (m, s, zi) in zip(mixed_prior_median(txo)[plug.free],
                                                               mixed_prior_sd(txo)[plug.free], z))
    full(zf) = (xf = stage_x(plug, zf); xf)
    function mine(zf)
        xf = full(zf)
        r = [begin
                 w = [xf[iS[g]], xf[ik]]
                 fg = q -> gene_ll(bs[g], q)
                 res = DR.GradientResult(w)
                 ForwardDiff.gradient!(res, fg, w, ForwardDiff.GradientConfig(fg, w, ForwardDiff.Chunk{2}()))
                 (DR.value(res), DR.gradient(res))
             end for g in eachindex(bs)]
        gx = zeros(length(xf))
        for (g, q) in enumerate(r)
            gx[iS[g]] += q[2][1]
            gx[ik] += q[2][2]
        end
        z = zf
        gz = gx[plug.free] .+ ForwardDiff.gradient(prior, z)
        return sum(q[1] for q in r) + prior(z), gz
    end
    mine_density(zf) = (xf = full(zf); sum(gene_ll(bs[g], [xf[iS[g]], xf[ik]]) for g in eachindex(bs)) + prior(zf))
    # compile
    f(z0); ForwardDiff.gradient(f, z0); mine(z0); mine_density(z0)
    v_src, g_src = f(z0), ForwardDiff.gradient(f, z0)
    v_me, g_me = mine(z0)
    dv = abs(v_me - v_src) / abs(v_src)
    dg = maximum(abs.(g_me .- g_src)) / maximum(abs.(g_src))
    say(@sprintf("17-gene check at KIM-92's truth: density %.10f (src) vs %.10f (per gene), rel. diff %.2e; gradient max |Δ| / max |g| %.2e",
                 v_src, v_me, dv, dg))
    dv < 1e-10 && dg < 1e-8 || error("the per-gene block does not reproduce src's plug-in")
    t = (; src_density = reps(() -> f(z0), NREP), src_gradient = reps(() -> ForwardDiff.gradient(f, z0), NREP),
         me_density = reps(() -> mine_density(z0), NREP), me_gradient = reps(() -> mine(z0), NREP))
    say(@sprintf("17-gene timing: src plug-in density %.2f s, gradient %.2f s (ForwardDiff default chunk %d); per gene density %.2f s, dual pass %.2f s",
                 med(t.src_density), med(t.src_gradient), ForwardDiff.pickchunksize(length(z0)),
                 med(t.me_density), med(t.me_gradient)))
    return (; ncoord = length(z0), dv, dg, chunk_src = ForwardDiff.pickchunksize(length(z0)),
            ncells = length(cells), nwindows = size(txo.factors[1].y, 2),
            timing = map(v -> [Base.structdiff(r, (; value = 0)) for r in v], t))
end

function setup455(out)
    d = load(out, "data")
    bs = blocks(d)
    x = vcat(d.logS, log(d.krnadeg))
    order = sortperm([sum(f.caps[1] for f in b.fs) for b in bs]; rev = true)
    # compile on a one-cell copy, with the same types
    tiny = blocks(d; cells = 1:1)
    density(tiny, d, x); value_gradient(tiny, d, x)
    return d, bs, x, order, tiny
end

function run_serial(out)
    env = environment()
    say("serial: ", env)
    red = check_reduced()
    d, bs, x, order, tiny = setup455(out)
    G = length(bs)
    # The per-gene sum is the whole factor's density.
    fs = whole_factors(d)
    v_whole = whole_ll(fs, bs, d, x) + log_prior(d, x)
    v_gene = density(bs, d, x)
    say(@sprintf("455-gene check: whole factor %.10f vs per-gene sum %.10f, rel. diff %.2e", v_whole, v_gene,
                 abs(v_whole - v_gene) / abs(v_whole)))
    abs(v_whole - v_gene) / abs(v_whole) < 1e-10 || error("the per-gene sum is not the whole factor")
    # A subset gradient against ForwardDiff on the whole factor (first 6 genes + krnadeg).
    sub = 1:min(6, G)
    dsub = merge(d, (; logS_ref = d.logS_ref[sub], lens = d.lens[sub], y0 = d.y0[sub], kmax = d.kmax[sub],
                      μmin = d.μmin[sub], y = d.y[:, sub, :], genes = d.genes[sub]))
    bsub, fsub = bs[sub], whole_factors(dsub)
    xs = vcat(x[sub], x[end])
    g_whole = ForwardDiff.gradient(w -> whole_ll(fsub, bsub, dsub, w) + log_prior(dsub, w), xs)
    _, g_gene = value_gradient(bsub, dsub, xs)
    dgs = maximum(abs.(g_whole .- g_gene)) / maximum(abs.(g_whole))
    say(@sprintf("455-gene gradient check on genes %s + krnadeg: max |Δ| / max |g| %.2e", string(sub), dgs))
    dgs < 1e-8 || error("the per-gene gradient is not the whole factor's")
    # Timing.
    pergene = zeros(G)
    td = reps(() -> density(bs, d, x), NREP)
    tg = reps(() -> value_gradient(bs, d, x; pergene), NREP)
    say(@sprintf("455 genes, serial: density %s s (median %.2f); dual pass %s s (median %.2f)",
                 join((@sprintf("%.2f", r.seconds) for r in td), ", "), med(td),
                 join((@sprintf("%.2f", r.seconds) for r in tg), ", "), med(tg)))
    # The naive chunk pass, on the whole factor.
    N = ForwardDiff.pickchunksize(G + 1)
    npass = cld(G + 1, N)
    naive_pass(whole_factors(d; cells = 1:1), bs, d, x, Val(N))
    tn = reps(() -> naive_pass(fs, bs, d, x, Val(N)), NREP_NAIVE)
    say(@sprintf("naive: one chunk-%d pass %s s (median %.2f); × %d passes = %.0f s per gradient",
                 N, join((@sprintf("%.2f", r.seconds) for r in tn), ", "), med(tn), npass, npass * med(tn)))
    strip_(v) = [Base.structdiff(r, (; value = 0)) for r in v]
    s = (; env, red, G, ncells = size(d.y, 1), nwindows = size(d.y, 3), check_density = abs(v_whole - v_gene) / abs(v_whole),
         check_gradient = dgs, density = strip_(td), gradient = strip_(tg), value = td[1].value,
         gradnorm = sqrt(sum(abs2, tg[1].value[2])), pergene, naive = (; chunk = N, passes = npass, pass = strip_(tn)))
    save(out, "serial", s; stage = "serial")
    write_serial_md(out, s, d)
end

fmt(v) = join((@sprintf("%.2f", r.seconds) for r in v), ", ")

function write_serial_md(out, s, d)
    io = IOBuffer()
    p(xs...) = println(io, xs...)
    e = s.env
    p("# KIM-185 item 2: the transcript factor on ", s.G, " genes, serial\n")
    p("Commit `", COMMIT, "`, job ", JOB, ", on ", e.host, " (", e.cpu, "), Julia ", e.julia, ", ForwardDiff ",
      e.forwarddiff, ", ", e.threads, " thread(s), BLAS ", e.blas, " thread(s). SMOKE = ", SMOKE, ". ",
      s.ncells, " cells × ", s.G, " genes × ", s.nwindows, " windows. Each call ran ", NREP,
      " times after compiling on a one-cell copy; a full GC ran before each.\n")
    p("## The 17-gene plug-in, on this node\n")
    r = s.red
    p("At KIM-92's truth, ", r.ncells, " cells × ", r.nwindows, " windows, ", r.ncoord,
      " coordinates. The per-gene block reproduces src's `plugin_loglik` stage density to a relative ",
      @sprintf("%.2e", r.dv), " and its gradient to ", @sprintf("%.2e", r.dg), " (max |Δ| over max |g|).\n")
    p("| call | AD | each, s | median, s | GiB allocated (median) |")
    p("|---|---|---|---|---|")
    t = r.timing
    p("| src plug-in density | — | ", fmt(t.src_density), " | ", @sprintf("%.2f", med(t.src_density)), " | ", @sprintf("%.2f", med(t.src_density, :gib)), " |")
    p("| src plug-in gradient | ForwardDiff.gradient, default chunk ", r.chunk_src, " (", cld(r.ncoord, r.chunk_src), " passes) | ",
      fmt(t.src_gradient), " | ", @sprintf("%.2f", med(t.src_gradient)), " | ", @sprintf("%.2f", med(t.src_gradient, :gib)), " |")
    p("| per-gene density | — | ", fmt(t.me_density), " | ", @sprintf("%.2f", med(t.me_density)), " | ", @sprintf("%.2f", med(t.me_density, :gib)), " |")
    p("| per-gene dual pass | ForwardDiff, Chunk{2} per gene | ", fmt(t.me_gradient), " | ", @sprintf("%.2f", med(t.me_gradient)), " | ",
      @sprintf("%.2f", med(t.me_gradient, :gib)), " |")
    p("\nKIM-116's profile (job 21196497, dave14) measured 1.18 s and ", KIM116_PLUGIN_GRADIENT_S, " s for the first two rows.\n")
    p("## ", s.G, " genes\n")
    p("- The per-gene sum equals the whole factor (one `TranscriptFactor` per cell over every gene) to a relative ",
      @sprintf("%.2e", s.check_density), ", and the per-gene gradient equals ForwardDiff's on the whole factor, on 6 genes + krnadeg, to ",
      @sprintf("%.2e", s.check_gradient), ".")
    p("- Log density at the truth ", @sprintf("%.6f", s.value), "; gradient norm ", @sprintf("%.6g", s.gradnorm), ".\n")
    p("| call | AD | each, s | median, s | GiB allocated (median) | GC s (median) |")
    p("|---|---|---|---|---|---|")
    p("| density | — | ", fmt(s.density), " | ", @sprintf("%.2f", med(s.density)), " | ", @sprintf("%.2f", med(s.density, :gib)), " | ",
      @sprintf("%.2f", med(s.density, :gc)), " |")
    p("| dual pass (value + gradient, ", s.G + 1, " coordinates) | ForwardDiff, Chunk{2} per gene | ", fmt(s.gradient), " | ",
      @sprintf("%.2f", med(s.gradient)), " | ", @sprintf("%.2f", med(s.gradient, :gib)), " | ", @sprintf("%.2f", med(s.gradient, :gc)), " |")
    n = s.naive
    p("| naive: one chunk pass | ForwardDiff Dual width ", n.chunk, " over all ", s.G + 1, " coordinates, whole factor | ", fmt(n.pass), " | ",
      @sprintf("%.2f", med(n.pass)), " | ", @sprintf("%.2f", med(n.pass, :gib)), " | ", @sprintf("%.2f", med(n.pass, :gc)), " |")
    p("\nThe naive gradient is the chunk pass × ", n.passes, " passes = ", @sprintf("%.0f", n.passes * med(n.pass)),
      " s. It is a multiplication, not a timed gradient.\n")
    pg = s.pergene
    srt = sort(pg; rev = true)
    p("Per-gene dual-pass time (last repeat): median ", @sprintf("%.3g", median(pg)), " s, max ", @sprintf("%.3g", maximum(pg)),
      " s; the 10% most expensive genes take ", @sprintf("%.0f%%", 100 * sum(srt[1:cld(length(srt), 10)]) / sum(srt)), " of the total.")
    write(joinpath(out, "serial.md"), String(take!(io)))
    print(read(joinpath(out, "serial.md"), String))
end

# ---------------------------------------------------------------- threaded

function run_threaded(out)
    env = environment()
    say("threaded: ", env)
    d, bs, x, order, _ = setup455(out)
    ws = [parse(Int, w) for w in split(get(ENV, "WORKERS", SMOKE ? "1,2" : "1,4,16,32"), ',')]
    all(<=(Threads.nthreads()), ws) || error("WORKERS $ws exceed $(Threads.nthreads()) threads")
    v1, g1 = value_gradient(bs, d, x)
    rows = map(ws) do w
        density(bs, d, x; workers = w, order); value_gradient(bs, d, x; workers = w, order)
        vw, gw = value_gradient(bs, d, x; workers = w, order)
        (vw == v1 && gw == g1) || error("$w workers do not give the serial value and gradient bit for bit")
        td = reps(() -> density(bs, d, x; workers = w, order), NREP)
        tg = reps(() -> value_gradient(bs, d, x; workers = w, order), NREP)
        say(@sprintf("%d workers: density %s s (median %.2f); dual pass %s s (median %.2f)", w, fmt(td), med(td), fmt(tg), med(tg)))
        (; workers = w, density = [Base.structdiff(r, (; value = 0)) for r in td],
         gradient = [Base.structdiff(r, (; value = 0)) for r in tg])
    end
    t = (; env, rows, G = length(bs))
    save(out, "threaded", t; stage = "threaded")
    write_threaded_md(out, t)
end

function write_threaded_md(out, t)
    io = IOBuffer()
    p(xs...) = println(io, xs...)
    e = t.env
    p("# KIM-185 item 3: the transcript factor on ", t.G, " genes, threaded over genes\n")
    p("Commit `", COMMIT, "`, job ", JOB, ", on ", e.host, " (", e.cpu, "), Julia ", e.julia, ", ForwardDiff ", e.forwarddiff,
      ", ", e.threads, " threads, GC threads ", e.gc.total, " (sweep ", e.gc.sweep, "), BLAS ", e.blas, " thread(s). SMOKE = ", SMOKE,
      ". W tasks pull genes off a shared counter, heaviest (by summed cap) first; each gene's dual pass is ForwardDiff with Chunk{2}. ",
      "At every W the value and gradient equal the one-worker result bit for bit. Each call ran ", NREP, " times; a full GC ran before each.\n")
    b = t.rows[1]
    p("| workers | density, each s | density median s | speed-up | dual pass, each s | dual pass median s | speed-up | efficiency | GC s (median) |")
    p("|---|---|---|---|---|---|---|---|---|")
    for r in t.rows
        sd, sg = med(b.density) / med(r.density), med(b.gradient) / med(r.gradient)
        p(@sprintf("| %d | %s | %.2f | %.2f | %s | %.2f | %.2f | %.2f | %.2f |", r.workers, fmt(r.density), med(r.density), sd,
                   fmt(r.gradient), med(r.gradient), sg, sg / r.workers, med(r.gradient, :gc)))
    end
    p("\nSpeed-up is against the first row, on the same node in the same job.")
    write(joinpath(out, "threaded.md"), String(take!(io)))
    print(read(joinpath(out, "threaded.md"), String))
end

# ---------------------------------------------------------------- the projection

function kim116_steps()
    d = joinpath(ARCHIVE116, "run", "plugin")
    rows = map(1:4) do c
        ch = deserialize(joinpath(d, "chain_$c.jls"))
        ex = deserialize(joinpath(d, "extend_$c.jls"))
        nc, ne = Int.(vec(ch.result.chain[:n_steps])), Int.(vec(ex.result.chain[:n_steps]))
        (; c, chain_iters = length(nc), chain_steps = sum(nc), ext_iters = length(ne), ext_steps = sum(ne),
         seconds = ch.result.seconds + ex.result.seconds, depth_max = maximum(vcat(Int.(vec(ch.result.chain[:tree_depth])),
                                                                                   Int.(vec(ex.result.chain[:tree_depth])))))
    end
    m = deserialize(joinpath(d, "metric.jls")).result
    return (; rows, metric_gradients = 2 * size(m.H, 1), ncoord = size(m.H, 1), metric_seconds = m.seconds)
end

function run_project(out)
    K = kim116_steps()
    s, t = load(out, "serial"), load(out, "threaded")
    steps = sum(r.chain_steps + r.ext_steps for r in K.rows)
    secs = sum(r.seconds for r in K.rows)
    D = s.G + 1
    ngrad = steps + 2 * D                       # the chains' leapfrog steps, plus the metric's 2·D gradients
    tg = med(s.gradient)
    cpuh = ngrad * tg / 3600
    dimfac = (D / K.ncoord)^(1 / 4)
    cpuh_dim = (steps * dimfac + 2 * D) * tg / 3600
    naive = ngrad * s.naive.passes * med(s.naive.pass) / 3600
    ratio = KIM116_PLUGIN_STAGE1_CPUH * tg / med(s.red.timing.src_gradient)
    verdict = cpuh <= THRESHOLD_CPUH ? "PASS" : "FAIL"
    io = IOBuffer()
    p(xs...) = println(io, xs...)
    p("# KIM-185: the gene-expression block's cost by the data plug-in\n")
    p("Commit `", COMMIT, "`, job ", JOB, ". Timings from `serial.md` (job ", load_job(out, "serial"), ") and `threaded.md` (job ",
      load_job(out, "threaded"), "). Gradients per posterior from KIM-116's plug-in stage 1 (", ARCHIVE116, "/run/plugin).\n")
    p("## Gradients per posterior, counted\n")
    p("| chain | iterations (chain + extension) | leapfrog steps (chain + extension) | max tree depth | hours |")
    p("|---|---|---|---|---|")
    for r in K.rows
        p(@sprintf("| %d | %d + %d | %d + %d | %d | %.1f |", r.c, r.chain_iters, r.ext_iters, r.chain_steps, r.ext_steps, r.depth_max, r.seconds / 3600))
    end
    p(@sprintf("\nTotal %d leapfrog steps in %.1f h of sampling (%.2f s a step, against the profile's %.2f s a gradient). The warm metric took %d gradients at %d coordinates (%.0f s). At %d coordinates it takes 2·%d = %d. So one posterior is %d gradients, before any scaling with dimension.\n",
               steps, secs / 3600, secs / steps, KIM116_PLUGIN_GRADIENT_S, K.metric_gradients, K.ncoord, K.metric_seconds, D, D, 2D, ngrad))
    p("## The projection\n")
    p("| projection | CPU-h |")
    p("|---|---|")
    p(@sprintf("| **headline:** %d gradients × %.2f s (serial per-gene dual pass) | **%.0f** |", ngrad, tg, cpuh))
    p(@sprintf("| the same, with the leapfrog steps × (%d / %d)^¼ = %.2f for NUTS's step size at higher dimension | %.0f |", D, K.ncoord, dimfac, cpuh_dim))
    p(@sprintf("| KIM-116's stage-1 %.1f CPU-h × (%.2f s ÷ %.2f s, src's 17-gene plug-in gradient on this node) | %.0f |",
               KIM116_PLUGIN_STAGE1_CPUH, tg, med(s.red.timing.src_gradient), ratio))
    p(@sprintf("| naive AD: %d gradients × %d chunk passes × %.2f s | %.0f |", ngrad, s.naive.passes, med(s.naive.pass), naive))
    for r in t.rows
        w = r.workers
        tw = med(r.gradient)
        p(@sprintf("| threaded, %d workers: %d gradients × %.2f s × %d CPUs (wall %.1f h per chain of %d) | %.0f |",
                   w, ngrad, tw, w, (ngrad / 4) * tw / 3600, ngrad ÷ 4, ngrad * tw * w / 3600))
    end
    p("\n**Verdict: ", verdict, " against the 2026-11-26 threshold** (at most ", Int(THRESHOLD_CPUH),
      " CPU-h for the gene-expression block over 455 genes by the data plug-in; Decision log 2026-11-26, closing KIM-176). ",
      "The verdict is on the headline row, fixed before the timing was read.\n")
    write(joinpath(out, "projection.md"), String(take!(io)))
    save(out, "projection", (; K, steps, ngrad, tg, cpuh, cpuh_dim, dimfac, naive, ratio, verdict); stage = "project")
    print(read(joinpath(out, "projection.md"), String))
end

load_job(out, name) = deserialize(joinpath(out, "$name.jls")).job

function main()
    length(ARGS) == 2 || error("usage: julia --project $(PROGRAM_FILE) <data|serial|threaded|project> <out>")
    stage, out = ARGS[1], abspath(ARGS[2])
    stage == "data" ? run_data(out) :
    stage == "serial" ? run_serial(out) :
    stage == "threaded" ? run_threaded(out) :
    stage == "project" ? run_project(out) : error("unknown stage $stage")
end

main()
