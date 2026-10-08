# spec/phases/16-recovery.md §12 2026-10-07, option 2: what a Gibbs sweep costs
# at each thread count. Block 3's cell sweeps run at 40 to 51% parallel efficiency
# at 20 threads (profile_b3_16a9b.jl), and threads do not change the chain (the
# 16a.9a tests), so fewer threads per chain may cost fewer CPU-h per sweep.
#
# Resumes chain 1 at sweep 700 of V6 stage 1, runs one warm-up sweep and
# `nsweeps` timed ones, and saves nothing. The final θ is printed to full
# precision, so runs at different thread counts can be compared bitwise.
# Run: thread_scaling_16a9b.slurm <checkpoint> <nsweeps>, one job per count.

using InferCell
using Serialization
using Printf
using LinearAlgebra: BLAS

BLAS.set_num_threads(1)

function main(ckpt, nsweeps)
    ds = m0_dataset(16085; ncells = 50, horizon = M0_HORIZON_S)
    g = gibbs_setup(ds; N = 20, horizon = M0_HORIZON_S)
    s = gibbs_resume(deserialize(ckpt).state, g)
    nt = Threads.nthreads()
    println("Thread scaling at $nt threads from sweep $(s.sweep); one warm-up sweep, then $nsweeps.\n")
    println("| Sweep | Wall (s) | B3 (s) | B1 (s) | CPU-h at $nt cores | GC share |")
    println("|---|---|---|---|---|---|")
    for k in 0:nsweeps
        t = @timed gibbs_sweep!(s, g)
        sec = s.trace[end].seconds
        @printf("| %d%s | %.1f | %.1f | %.1f | %.3f | %.1f%% |\n", s.sweep, k == 0 ? " (warm-up)" : "",
                t.time, sec.b3, sec.b1, t.time * nt / 3600, 100 * t.gctime / t.time)
        flush(stdout)
    end
    println("\nFinal state: cme = ", repr(s.cme), ", u = ", repr(s.u), ", σ = ", repr(s.σ))
end

main(ARGS[1], parse(Int, ARGS[2]))
