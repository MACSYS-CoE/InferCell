# spec/phases/16-recovery.md task 16a.9b: V7, simulation-based calibration of the
# M0 Gibbs chain at production settings (PG, N = 20, C = 50 cells, bulk
# metabolites), sized by §12 2026-10-08.
#
# Replicate r's dataset is `m0_dataset(16_200_000 + r)`. That seed draws the
# truth of M0's five targets from their priors (`draw_truth`, purpose
# :calibration), then σ_b from `SIGMA_MET_PRIOR`, then the cells' and noise seeds.
# The chain starts from θ and σ_b drawn from their priors with
# Xoshiro(16_300_000 + r), so the start is independent of the truth. It runs
# SWEEPS sweeps and checkpoints after each one.
#
# The merge drops BURN sweeps, thins the rest to L draws, and ranks each
# parameter's truth among them. The ranks then get a rank-ECDF with 95%
# simultaneous bands (Säilynoja, Bürkner and Vehtari 2022, by simulation).
#
# Stages:
#   run <r>     replicate r, resumable
#   merge       ranks, bands and per-replicate diagnostics over every replicate
#               that has a checkpoint; missing or short ones are named
# Regenerate: sbc_v7_16a9b.sh, from a worktree at a pushed commit.

using InferCell
using Random
using Serialization
using Statistics
using Printf
using Distributions: Binomial, cdf, quantile
using LinearAlgebra: BLAS
using MCMCChains: MCMCDiagnosticTools
const MDT = MCMCDiagnosticTools

BLAS.set_num_threads(1)

const DIR = get(ENV, "V7_DIR", joinpath(@__DIR__, "sbc_v7_16a9b"))   # elsewhere for a smoke run
mkpath(DIR)
const NCELLS = parse(Int, get(ENV, "V7_NCELLS", "50"))       # 3 for a smoke run
const HORIZON = parse(Float64, get(ENV, "V7_HORIZON", string(M0_HORIZON_S)))
const SWEEPS = parse(Int, get(ENV, "V7_SWEEPS", "656"))      # §12 2026-10-08
const BURN = parse(Int, get(ENV, "V7_BURN", "100"))
const L = parse(Int, get(ENV, "V7_L", "100"))                # draws per replicate after thinning
const N = 20
const DATA_SEED0 = 16_200_000
const START_SEED0 = 16_300_000

chainfile(r) = joinpath(DIR, "replicate$(r).jls")

function run_replicate(r)
    ds = m0_dataset(DATA_SEED0 + r; ncells = NCELLS, horizon = HORIZON)
    g = gibbs_setup(ds; N, horizon = HORIZON)
    @printf("V7 replicate %d: seed %d, %d cells, N = %d, %d threads, %d sweeps\n",
            r, DATA_SEED0 + r, NCELLS, N, Threads.nthreads(), SWEEPS)
    if isfile(chainfile(r))
        saved = deserialize(chainfile(r))
        s = gibbs_resume(saved.state, g)
        meta = saved.meta
        println("  resumed at sweep $(s.sweep)")
    else
        rng = Xoshiro(START_SEED0 + r)
        t0 = draw_truth(rng, g.models; names = M0_TARGETS, purpose = :calibration)
        v = Dict(t0.values)
        σ0 = rand(rng, SIGMA_MET_PRIOR)
        cme0 = [v[n] for n in g.cme]
        u0 = [log(v[f]) for f in g.forwards]
        tinit = @elapsed s = gibbs_init(rng, g; cme = cme0, u = u0, σ = σ0)
        truth = Dict(ds.truth.values)
        meta = (; replicate = r, seed = DATA_SEED0 + r, names = vcat(g.cme, g.forwards, :sigma_b),
                truth = vcat([log(truth[n]) for n in g.cme], [log(truth[f]) for f in g.forwards],
                             log(ds.sigma)),
                start = (cme = cme0, u = u0, σ = σ0), init_seconds = tinit, ncells = NCELLS, N,
                threads = Threads.nthreads())
        @printf("  truth (log) = %s; start: cme = %s, u = %s, σ_b = %.3g; init took %.0f s\n",
                round.(meta.truth; digits = 3), round.(cme0; sigdigits = 3),
                round.(u0; digits = 3), σ0, tinit)
    end
    tmp = chainfile(r) * ".tmp"
    save() = (serialize(tmp, (; state = gibbs_checkpoint(s), meta)); mv(tmp, chainfile(r); force = true))
    save()
    while s.sweep < SWEEPS
        gibbs_sweep!(s, g)
        save()
        t = s.trace[end]
        @printf("  sweep %d: %.0f s (B3 %.0f, B1 %.0f); ln θ = %s; σ_b = %.4f; windows changed %.2f\n",
                s.sweep, sum(t.seconds), t.seconds.b3, t.seconds.b1,
                round.(vcat(log.(t.cme), t.u); digits = 3), t.σ, mean(t.changed))
        flush(stdout)
    end
end

draws(t) = vcat(log.(t.cme), t.u, log(t.σ))

# Simultaneous 1 − α bands for the ECDF of R ranks uniform on 0:L, on the grid
# z_i = i / (L + 1): the pointwise level γ is the α quantile, over simulated
# uniform rank sets, of the smallest two-sided binomial tail probability.
function ecdf_bands(R, L; α = 0.05, nsim = 10_000, rng = Xoshiro(16_299_999))
    z = (1:L) ./ (L + 1)
    tail(k, p) = (d = Binomial(R, p); 2 * min(cdf(d, k), 1 - cdf(d, k - 1)))
    mins = map(1:nsim) do _
        ranks = rand(rng, 0:L, R)
        minimum(tail(count(<(i), ranks), z[i]) for i in 1:L)
    end
    γ = quantile(mins, α)
    lo = [quantile(Binomial(R, p), γ / 2) for p in z]
    hi = [quantile(Binomial(R, p), 1 - γ / 2) for p in z]
    return z, lo, hi, γ
end

function merge_replicates()
    files = filter(f -> occursin(r"^replicate\d+\.jls$", f), readdir(DIR))
    saved = sort([deserialize(joinpath(DIR, f)) for f in files]; by = x -> x.meta.replicate)
    done = filter(x -> length(x.state.trace) >= SWEEPS, saved)
    short = [(x.meta.replicate, length(x.state.trace)) for x in saved if length(x.state.trace) < SWEEPS]
    isempty(done) && error("no replicate has reached $SWEEPS sweeps")
    names = string.(first(done).meta.names)
    P = length(names)
    R = length(done)
    println("V7: $R replicates at $SWEEPS sweeps, $BURN dropped, thinned to L = $L; ",
            "$(first(done).meta.ncells) cells, N = $N.")
    isempty(short) || println("Not yet complete (replicate, sweeps): ", short)

    ranks = zeros(Int, R, P)
    ess_min = zeros(R)
    println("\n| Replicate | Seed | Min bulk ESS (kept) | Ranks (of $L) |")
    println("|---|---|---|---|")
    for (j, x) in enumerate(done)
        kept = x.state.trace[(BURN + 1):SWEEPS]
        X = reduce(hcat, draws.(kept))'                       # sweeps × parameters
        ess_min[j] = minimum(MDT.ess(reshape(X, size(X, 1), 1, P); kind = :bulk))
        idx = round.(Int, range(1, size(X, 1); length = L))
        for p in 1:P
            ranks[j, p] = count(<(x.meta.truth[p]), X[idx, p])
        end
        @printf("| %d | %d | %.0f | %s |\n", x.meta.replicate, x.meta.seed, ess_min[j],
                join(ranks[j, :], ", "))
    end
    nlow = count(e -> !(e >= L), ess_min)                    # NaN (too few draws) counts as low
    println("\nReplicates whose smallest bulk ESS over the kept sweeps is below L = $L: $nlow of $R ",
            "(their thinned draws are autocorrelated, which widens the rank distribution).")

    z, lo, hi, γ = ecdf_bands(R, L)
    println("\nRank-ECDF against 95% simultaneous bands (pointwise level γ = $(round(γ; sigdigits = 3))):\n")
    println("| Parameter (log) | Outside the band at | Largest |ECDF − z| | Verdict |")
    println("|---|---|---|---|")
    for p in 1:P
        counts = [count(<(i), ranks[:, p]) for i in 1:L]
        out = [i for i in 1:L if !(lo[i] <= counts[i] <= hi[i])]
        dev = maximum(abs.(counts ./ R .- z))
        @printf("| %s | %s | %.3f | %s |\n", names[p],
                isempty(out) ? "nowhere" : join(round.(z[out]; digits = 2), ", "), dev,
                isempty(out) ? "within" : "OUTSIDE")
    end
    serialize(joinpath(DIR, "ranks.jls"), (; names, ranks, ess_min, R, L, z, lo, hi, γ,
                                          replicates = [x.meta.replicate for x in done]))
end

if ARGS[1] == "run"
    run_replicate(parse(Int, ARGS[2]))
elseif ARGS[1] == "merge"
    merge_replicates()
else
    error("usage: run <r> | merge")
end
