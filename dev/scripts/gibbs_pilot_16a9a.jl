# spec/phases/16-recovery.md task 16a.9a: a pilot of the three-block Gibbs chain
# on M0, to size 16a.9b (V6 and V7).
#
# Four chains on M0 replicate 16085 (16a.7c's), 50 cells observed in bulk, particle
# Gibbs at N = 20. Chain k starts from θ and σ_b drawn from their priors with
# Xoshiro(16_090_000 + k), and from paths drawn by transcript-weighted SMC at
# that θ (`gibbs_init`), so the starts are overdispersed.
#
# Stages:
#   run <chain> <sweeps>   resumable: it checkpoints after every sweep
#   merge [burn]           diagnostics on the sweeps after `burn` (default: half)
# Regenerate: gibbs_pilot_16a9a.sh, from a worktree at a pushed commit.

using InferCell
using Random
using Serialization
using Statistics
using Printf
using LinearAlgebra: BLAS
using MCMCChains: MCMCDiagnosticTools
const MDT = MCMCDiagnosticTools

BLAS.set_num_threads(1)

const DIR = get(ENV, "PILOT_DIR", joinpath(@__DIR__, "gibbs_pilot_16a9a"))   # elsewhere for a smoke run
mkpath(DIR)
const REPLICATE = 16085
const NCELLS = parse(Int, get(ENV, "PILOT_NCELLS", "50"))     # 3 for a smoke run
const HORIZON = parse(Float64, get(ENV, "PILOT_HORIZON", string(M0_HORIZON_S)))
const N = 20
const THRESHOLD = 0.10

setup() = (ds = m0_dataset(REPLICATE; ncells = NCELLS, horizon = HORIZON);
           (ds, gibbs_setup(ds; N, horizon = HORIZON)))

chainfile(k) = joinpath(DIR, "chain$(k).jls")

# A plain cycle of one cell, timed warm: the cost unit.
function cycle_seconds(g)
    run1() = (d = restore(g.pristine[1]); Random.seed!(1);
              for _ in 1:round(Int, HORIZON); handshake_step!(d); end)
    run1()
    return @elapsed run1()
end

function save(k, s, meta)
    tmp = chainfile(k) * ".tmp"
    serialize(tmp, (; state = gibbs_checkpoint(s), meta))
    mv(tmp, chainfile(k); force = true)
end

function run_chain(k, sweeps)
    ds, g = setup()
    cyc = cycle_seconds(g)
    @printf("pilot chain %d: M0 replicate %d, %d cells, N = %d, %d threads; a plain cycle takes %.2f s\n",
            k, REPLICATE, NCELLS, N, Threads.nthreads(), cyc)
    if isfile(chainfile(k))
        saved = deserialize(chainfile(k))
        s = gibbs_resume(saved.state, g)
        meta = saved.meta
        println("  resumed at sweep $(s.sweep)")
    else
        rng = Xoshiro(16_090_000 + k)
        t0 = draw_truth(rng, g.models; names = M0_TARGETS, purpose = :calibration)
        v = Dict(t0.values)
        σ0 = rand(rng, SIGMA_MET_PRIOR)
        cme0 = [v[n] for n in g.cme]
        u0 = [log(v[f]) for f in g.forwards]
        tinit = @elapsed s = gibbs_init(rng, g; cme = cme0, u = u0, σ = σ0)
        meta = (; start = (cme = cme0, u = u0, σ = σ0), init_seconds = tinit,
                names = vcat(g.cme, g.forwards, :sigma_b), truth = ds.truth.values,
                sigma_true = ds.sigma, ncells = NCELLS, N, replicate = REPLICATE)
        @printf("  start: cme = %s, u = %s, σ_b = %.3g; init took %.0f s\n",
                round.(cme0; sigdigits = 3), round.(u0; digits = 3), σ0, tinit)
    end
    meta = merge(meta, (; cycle = cyc, threads = Threads.nthreads()))
    save(k, s, meta)
    while s.sweep < sweeps
        gibbs_sweep!(s, g)
        save(k, s, meta)
        t = s.trace[end]
        @printf("  sweep %d: %.0f s (B3 %.0f, B1 %.0f, B2 %.2f, rebuild %.1f); ln θ = %s; σ_b = %.4f; windows changed %.2f\n",
                s.sweep, sum(t.seconds), t.seconds.b3, t.seconds.b1, t.seconds.b2,
                t.seconds.rebuild, round.(vcat(log.(t.cme), t.u); digits = 3), t.σ, mean(t.changed))
        flush(stdout)
    end
end

function merge_chains(burn_arg)
    files = sort(filter(f -> occursin(r"^chain\d+\.jls$", f), readdir(DIR)))
    saved = [deserialize(joinpath(DIR, f)) for f in files]
    meta = first(saved).meta
    S = minimum(length(r.state.trace) for r in saved)
    burn = burn_arg === nothing ? S ÷ 2 : burn_arg
    kept = S - burn
    kept >= 4 || error("only $kept sweeps after burn-in")
    names = string.(meta.names)
    P = length(names)
    # draws × chains × parameters, all on the log scale.
    x = zeros(kept, length(saved), P)
    for (c, r) in enumerate(saved), (i, t) in enumerate(r.state.trace[(burn + 1):S])
        x[i, c, :] = vcat(log.(t.cme), t.u, log(t.σ))
    end
    truth = Dict(meta.truth)
    tv = vcat([log(truth[Symbol(n)]) for n in names[1:end - 1]], log(meta.sigma_true))
    println("Pilot of 16a.9a: $(length(saved)) chains, $S sweeps each, the first $burn dropped; ",
            "$(meta.ncells) cells, N = $(meta.N), $(meta.threads) threads.\n")

    rhat = MDT.rhat(x; kind = :rank)
    essb = MDT.ess(x; kind = :bulk)
    esst = MDT.ess(x; kind = :tail)
    total = kept * length(saved)
    println("| Parameter (log) | R̂ (rank) | Bulk ESS | Tail ESS | Bulk ESS per sweep | Tail ESS per sweep | ",
            "MCSE(q05)/SD | MCSE(q95)/SD | Posterior mean ± SD | Truth | Truth's quantile |")
    println("|---|---|---|---|---|---|---|---|---|---|---|")
    mq = zeros(P, 2)
    for p in 1:P
        v = vec(x[:, :, p])
        sd = std(v)
        for (j, q) in enumerate((0.05, 0.95))
            mq[p, j] = MDT.mcse(x[:, :, p:p]; kind = Base.Fix2(quantile, q))[1] / sd
        end
        @printf("| %s | %.3f | %.0f | %.0f | %.3f | %.3f | %.3f | %.3f | %.3f ± %.3f | %.3f | %.2f |\n",
                names[p], rhat[p], essb[p], esst[p], essb[p] / total, esst[p] / total,
                mq[p, 1], mq[p, 2], mean(v), sd, tv[p], mean(v .<= tv[p]))
    end

    # Update rates per cell and window after burn-in. A window whose observed
    # counts are 0 → 0 for every gene is reported apart: its posterior sits near
    # the empty path, so a low rate is expected there (16a.7c).
    ds, g = setup()
    T = size(g.bulk, 2)
    rate = zeros(T, meta.ncells)
    for r in saved, t in r.state.trace[(burn + 1):S]
        rate .+= t.changed
    end
    rate ./= total
    y0 = [g.pristine[c].jump.u[s] for s in g.tmap.states, c in 1:meta.ncells]
    live = falses(T, meta.ncells)
    for c in 1:meta.ncells, m in 1:T
        prev = m == 1 ? y0[:, c] : g.transcripts[c][:, m - 1]
        live[m, c] = any(prev .> 0) || any(g.transcripts[c][:, m] .> 0)
    end
    lv = rate[live]
    @printf("\nUpdate rates: over %d live cell-windows the lowest is %.2f, and %d are below %.2f; ",
            length(lv), minimum(lv), count(<(THRESHOLD), lv), THRESHOLD)
    @printf("the %d empty cell-windows average %.2f.\n", count(.!live),
            count(.!live) > 0 ? mean(rate[.!live]) : NaN)

    # Cost, from every sweep after each chain's first.
    secs = [t.seconds for r in saved for t in r.state.trace[2:S]]
    b3, b1, b2, rb = (mean(getfield.(secs, f)) for f in (:b3, :b1, :b2, :rebuild))
    sweep_s = b3 + b1 + b2 + rb
    threads = meta.threads
    @printf("\nCost per sweep, wall-clock at %d threads: %.0f s (B3 %.0f, B1 %.0f, B2 %.2f, rebuild %.1f), ",
            threads, sweep_s, b3, b1, b2, rb)
    @printf("which is %.1f plain cycles of %.2f s; %.2f CPU-h per sweep at %d cores. Init: %.0f s.\n",
            sweep_s / meta.cycle, meta.cycle, sweep_s * threads / 3600, threads, meta.init_seconds)

    # Projections. V6 needs bulk and tail ESS ≥ 1,000 for every parameter and
    # every 5% and 95% quantile's MCSE ≤ 0.05 SD, on four chains. ESS grows with
    # the kept sweeps, and an MCSE falls as their square root.
    emin = minimum(min(essb[p], esst[p]) / total for p in 1:P)
    need_ess = 1000 / emin
    need_mcse = maximum(mq)^2 / 0.05^2 * total
    v6_kept = max(need_ess, need_mcse)
    v6_sweeps = v6_kept + 4burn            # four chains' burn-in, at the pilot's
    per_chain = v6_sweeps / 4
    # Reference at N_ref = 50: B3 is linear in N (16a.7c: about 0.84N cycles per
    # cell-sweep); ESS per sweep is taken as the pilot's, which is conservative.
    sweep50 = b3 * 50 / meta.N + b1 + b2 + rb
    println("\nProjection for V6 (4 chains):")
    @printf("- the slowest ESS per sweep is %.3f, so ESS ≥ 1,000 needs %.0f kept sweeps; the quantile MCSE needs %.0f;\n",
            emin, need_ess, need_mcse)
    for (lab, sw) in (("N_ref = 20", sweep_s), ("N_ref = 50", sweep50))
        @printf("- %s: %.0f sweeps in all, %.0f per chain: %.0f CPU-h, and %.1f days of wall-clock with the chains in parallel at %d threads each.\n",
                lab, v6_sweeps, per_chain, v6_sweeps * sw * threads / 3600, per_chain * sw / 86400, threads)
    end
    # V7: each replicate a chain with the pilot's burn-in and enough sweeps for
    # 100 effective draws of its slowest parameter, at production settings.
    rep_sweeps = burn + 100 / emin
    println("\nProjection for V7 (one chain per replicate, burn-in $(burn) plus 100 effective draws):")
    for R in (50, 100)
        @printf("- R = %d: %.0f sweeps per replicate, %.0f CPU-h in all, %.1f days of wall-clock with every replicate in parallel.\n",
                R, rep_sweeps, R * rep_sweeps * sweep_s * threads / 3600, rep_sweeps * sweep_s / 86400)
    end
end

if ARGS[1] == "run"
    run_chain(parse(Int, ARGS[2]), parse(Int, ARGS[3]))
elseif ARGS[1] == "merge"
    merge_chains(length(ARGS) >= 2 ? parse(Int, ARGS[2]) : nothing)
else
    error("stage must be run or merge")
end
