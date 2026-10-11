# spec/phases/16-recovery.md §12 2026-10-07, option 2: where block 3's time goes.
#
# The pilot's B3 takes about 130 s of a 205 s sweep at 20 threads, with core use
# near 62%. This splits one cell's CSMC sweep into its serial and parallel parts,
# window by window, on chain 1's state at sweep 700 of V6 stage 1:
#   1. threaded against serial, and GC time, for whole cell sweeps (`csmc_sweep`);
#   2. a copy of `csmc_sweep`'s loop, timing the resampling copies (serial), the
#      parallel region, and each particle's own time, so the barrier's idle share
#      is the gap between the slowest particle and the mean;
#   3. a serial profile of one cell sweep, by self time and by inclusive time.
# It changes nothing in the library. Run: profile_b3_16a9b.slurm <checkpoint>.

using InferCell
using InferCell: _gibbs_spec, advance_window!, _normalise, _categorical, _panel_latent,
                 restore, CellData, bulk_window_observations
using Random
using Serialization
using Statistics
using Printf
using Profile
using LinearAlgebra: BLAS

BLAS.set_num_threads(1)

const REPLICATE = 16085
const NCELLS = 50
const CELLS = (1, 17, 34)        # three cells, spread over the population

# Block 3's inputs for cell c, as `gibbs_sweep!` forms them.
function cell_inputs(s, g, c)
    spec = _gibbs_spec(g, s.σ)
    total = sum(s.latents)
    others = total .- s.latents[c]
    bulk = bulk_window_observations(g.bulk, others, NCELLS)
    return spec, CellData(g.transcripts[c], g.bulk), bulk
end

# `csmc_sweep`'s particle-Gibbs loop (lag 0, no schedule, bulk weights), timed.
function timed_sweep(rng, base, spec, data, ref, bulk; N, threads)
    T = size(data.transcripts, 2)
    y0 = [base.jump.u[s] for s in spec.tmap.states]
    tcopy = zeros(T)
    tpar = zeros(T)
    tpart = zeros(N, T)
    particles = [restore(base) for _ in 1:N]
    lw = zeros(N)
    for m in 1:T
        yprev = m == 1 ? y0 : data.transcripts[:, m - 1]
        if m > 1
            w = _normalise(lw)
            a = [_categorical(rng, w) for _ in 1:N]
            a[1] = 1
            tcopy[m] = @elapsed particles .= [restore(particles[a[i]]) for i in 1:N]
        end
        seeds = rand(rng, UInt64, N)
        body(i) = (tpart[i, m] = @elapsed begin
            prng = Xoshiro(seeds[i])
            fixed = i == 1 ? ref[m] : nothing
            lw[i], _ = advance_window!(prng, particles[i], spec, yprev,
                                       data.transcripts[:, m], bulk[m]; fixed)
            _panel_latent(particles[i], spec)
        end)
        tpar[m] = @elapsed if threads
            Threads.@threads for i in 1:N
                body(i)
            end
        else
            foreach(body, 1:N)
        end
    end
    return (; tcopy, tpar, tpart)
end

function main(ckpt)
    ds = m0_dataset(REPLICATE; ncells = NCELLS, horizon = M0_HORIZON_S)
    g = gibbs_setup(ds; N = 20, horizon = M0_HORIZON_S)
    saved = deserialize(ckpt)
    s = gibbs_resume(saved.state, g)
    N = g.N
    nt = Threads.nthreads()
    plain = let d = restore(g.pristine[1])
        f() = (e = restore(d); Random.seed!(1); for _ in 1:round(Int, M0_HORIZON_S); handshake_step!(e); end)
        f(); @elapsed f()
    end
    tr = @elapsed restore(s.bases[1])
    println("B3 profile on chain state at sweep $(s.sweep), $nt threads, N = $N, ",
            "$(size(g.bulk, 2)) windows. A plain cycle takes $(round(plain; digits = 3)) s, ",
            "one restore $(round(1e3 * tr; digits = 2)) ms.\n")

    rng = Xoshiro(16_097_000)
    # Warm up both paths.
    spec, data, bulk = cell_inputs(s, g, 1)
    csmc_sweep(rng, s.bases[1], spec, data, s.paths[1]; N, bulk, threads = true)
    csmc_sweep(rng, s.bases[1], spec, data, s.paths[1]; N, bulk, threads = false)
    timed_sweep(rng, s.bases[1], spec, data, s.paths[1], bulk; N, threads = true)

    println("## 1. Whole cell sweeps (`csmc_sweep`), 3 repeats each\n")
    println("| Cell | Threaded wall (s) | GC share | Serial wall (s) | GC share | Speed-up | Parallel efficiency |")
    println("|---|---|---|---|---|---|---|")
    for c in CELLS
        spec, data, bulk = cell_inputs(s, g, c)
        th = [@timed csmc_sweep(rng, s.bases[c], spec, data, s.paths[c]; N, bulk, threads = true)
              for _ in 1:3]
        se = [@timed csmc_sweep(rng, s.bases[c], spec, data, s.paths[c]; N, bulk, threads = false)
              for _ in 1:3]
        wt, ws = median(t.time for t in th), median(t.time for t in se)
        gt, gs = median(t.gctime / t.time for t in th), median(t.gctime / t.time for t in se)
        @printf("| %d | %.2f | %.1f%% | %.2f | %.1f%% | %.1f | %.0f%% |\n",
                c, wt, 100gt, ws, 100gs, ws / wt, 100 * ws / wt / min(nt, N))
    end

    println("\n## 2. Inside a threaded cell sweep, per window (copy of the loop), 3 repeats per cell\n")
    println("| Cell | Wall (s) | Resampling copies, serial | Parallel region | of which: slowest particle | mean particle | reference particle (mean) | Idle share of the parallel region |")
    println("|---|---|---|---|---|---|---|---|")
    for c in CELLS
        spec, data, bulk = cell_inputs(s, g, c)
        for _ in 1:3
            local r
            wall = @elapsed r = timed_sweep(rng, s.bases[c], spec, data, s.paths[c], bulk; N,
                                            threads = true)
            slow = sum(maximum(r.tpart; dims = 1))
            mean_ = sum(mean(r.tpart; dims = 1))
            ref = mean(r.tpart[1, :])
            idle = 1 - sum(r.tpart) / (sum(r.tpar) * min(nt, N))
            @printf("| %d | %.2f | %.2f | %.2f | %.2f | %.2f | %.3f | %.0f%% |\n",
                    c, wall, sum(r.tcopy), sum(r.tpar), slow, mean_, ref, 100idle)
        end
    end
    spec, data, bulk = cell_inputs(s, g, CELLS[1])
    r = timed_sweep(rng, s.bases[CELLS[1]], spec, data, s.paths[CELLS[1]], bulk; N, threads = false)
    println("\nSerial, cell $(CELLS[1]): particle seconds per window, summed over particles 2 to $N, ",
            "and the reference particle's:")
    for m in axes(r.tpart, 2)
        @printf("- window %d: proposals %.2f s (max %.3f, min %.3f), reference %.3f s\n",
                m, sum(r.tpart[2:end, m]), maximum(r.tpart[2:end, m]), minimum(r.tpart[2:end, m]),
                r.tpart[1, m])
    end

    println("\n## 3. Serial profile of one cell sweep (cell $(CELLS[1]), two sweeps)\n")
    Profile.init(n = 10^7, delay = 0.001)
    Profile.clear()
    @profile for _ in 1:2
        csmc_sweep(rng, s.bases[CELLS[1]], spec, data, s.paths[CELLS[1]]; N, bulk, threads = false)
    end
    println("### By self time\n```")
    Profile.print(IOContext(stdout, :displaysize => (200, 220)); format = :flat,
                  sortedby = :overhead, C = false, mincount = 20, noisefloor = 0)
    println("```\n### Tree, inclusive time\n```")
    Profile.print(IOContext(stdout, :displaysize => (400, 220)); format = :tree, C = false,
                  mincount = 100, maxdepth = 40, noisefloor = 2)
    println("```")
end

main(ARGS[1])
