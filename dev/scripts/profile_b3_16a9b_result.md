# 16a.9b: where block 3's time goes

§12 2026-10-07, option 2. Job 18255035 at `3a43f7d`, on chain 1's state at sweep
700 of V6 stage 1, at 20 threads. N = 20, and M0 has 10 windows. A plain cycle
takes 1.835 s, and one `restore` (a deep copy of a driver) takes 2.03 ms.

## 1. Whole cell sweeps: parallel efficiency is 40 to 51%

These are medians of 3 runs each.

| Cell | Threaded wall (s) | GC share | Serial wall (s) | GC share | Speed-up | Parallel efficiency |
|---|---|---|---|---|---|---|
| 1 | 3.51 | 13.8% | 36.06 | 1.3% | 10.3 | 51% |
| 17 | 2.58 | 12.5% | 23.41 | 1.3% | 9.1 | 45% |
| 34 | 2.15 | 9.9% | 17.25 | 1.4% | 8.0 | 40% |

**GC goes from about 1% of the serial time to 10 to 14% of the threaded time.**
Julia's collector stops every thread, so the allocation that is cheap serially
becomes a shared pause.

## 2. Inside a threaded cell sweep

The harness reruns a copy of `csmc_sweep`'s loop, with timers added. Times are
seconds per cell sweep, summed over its 10 windows.

| Cell | Wall | Resampling copies (serial) | Parallel region | Slowest particle | Mean particle | Idle share of the parallel region |
|---|---|---|---|---|---|---|
| 1 | 3.17 to 3.26 | 0.33 to 0.42 | 2.76 to 2.86 | = parallel region | 2.39 to 2.42 | 14 to 15% |
| 17 | 2.48 to 2.73 | 0.32 to 0.42 | 2.12 to 2.27 | = parallel region | 1.58 to 1.66 | 22 to 30% |
| 34 | 1.98 to 2.48 | 0.32 to 0.38 | 1.60 to 2.11 | = parallel region | 1.17 to 1.24 | 27 to 41% |

The wall time splits into three losses:
- **Resampling copies, about 12%.** After each window, all 20 particles are
  deep-copied from their ancestors on the main task: 9 resamplings × 20 copies ×
  about 2 ms.
- **The barrier, 14 to 41% of the parallel region.** Each window waits for its
  slowest particle. A particle's window takes from 0.03 to 0.26 s, depending on
  how many events it proposes. The spread suggests some proposals end early, at
  weight zero, but that was not measured directly.
- **Contention, about 1.5×.** A particle takes about 1.5 times longer when it
  runs alongside the others than when it runs alone. Cell 1's mean is 2.4 s
  threaded against about 1.6 s serially. GC pauses (§1) account for part of
  this.

## 3. Serial profile: the stiff ODE dominates

The profile covers two serial sweeps of cell 1. Nearly all of the time (97%) is
in the integrator's `step!`: Rodas5P on M0's 32-state ODE, out of place with
`SVector`s, and with ForwardDiff Jacobians. Within `step!`:
- **The Jacobian is about 39%.** It is ForwardDiff over the right-hand side,
  every step, because Rosenbrock methods need a fresh one.
- **The LU of `W` is at least 13%.** A 32 × 32 `SMatrix` is beyond StaticArrays'
  unrolled sizes, so `lu` falls back to LAPACK. That fallback allocates heap
  arrays, which is a GC source.
- **Right-hand-side evaluations and the stage algebra** take the rest. That
  includes `_fold_contributions` and `_accumulate` in `orchestrator.jl`.

Everything outside the ODE step comes to about 3% together: jump replay, bridges,
proposals, weights and copies.

## What this suggests

| Lever | Expected gain | Changes results? | Effort |
|---|---|---|---|
| Fewer threads per chain | **Measured** (`thread_scaling_16a9b_result.md`): at 5 threads a sweep costs 0.78 CPU-h against 1.22 at 20, 36% less, at 2.6× the wall-clock | No: the final state is bitwise equal at 1, 5, 10 and 20 threads | None |
| Do the resampling copies inside the parallel region | About 12% of B3's wall | No: each copy is independent | Small |
| Cut allocation in the ODE step, for example an in-place 32-state solve with a cached LU | Part of the 10–14% GC share and the 1.5× contention | Roundoff, so bitwise replays against recorded paths need checking | Moderate |
| A cheaper Jacobian (analytic or sparse), or a solver that reuses it | Up to about 39% of a particle's time | Yes: a different integrator changes the paths' latents at tolerance level, so V2 and V4 to V5 would need rerunning | Large |
