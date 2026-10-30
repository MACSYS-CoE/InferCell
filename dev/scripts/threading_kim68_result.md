# KIM-68: why `predict`'s threaded cells scale poorly (M5 prep)

**Run:** `dev/scripts/threading_kim68.sh` at `6e2a45c` (pushed; tag
`kim68-run-6e2a45c`), into `/fred/oz022/tkimpson/infercell-m5-kim68-2026-10-30/threading/`
(`jobs.tsv`, `logs/`, `summary.md`). The archive's `SHA256SUMS` is made at the end of
KIM-68 and its digest added here.
- **Jobs and nodes:**
  - precompile 19583235 (dave29); serial 19583236 (dave3);
  - 16 threads 19583237 (dave105) and 64 threads 19583238 (dave108);
  - at 64 threads, `--gcthreads=32,1` 19583239 (dave13), `--gcthreads=1` 19583240
    (dave50) and `OPENBLAS_NUM_THREADS=1` 19583241 (dave67);
  - 4 processes × 16 threads on one node 19583242 (dave114);
  - allocations 19583243 (dave8); report 19583244 (dave8).
- **The report stage was rerun at `c25b659`** (tag `kim68-run-c25b659`, job 19586524).
  It replaces the first summary's "tail after the median stop" with the idle share and
  GC's pauses. That column was `wall − median(stop)`, which is about half the wall for
  any balanced run, so it measured nothing. The first summary is `summary_6e2a45c.md`;
  its other columns are unchanged.
- **A SMOKE pass** ran end to end first, at `d9261a2` (jobs 19581856–19581866,
  `threading-smoke/`).
- Julia 1.10.5, milan, `--exclude=dave[148-149]`.

**The matrix** is KIM-38's: d11, B7's held-out `v6`, `completion = :prior`, M 4 × K 8 ×
4 cells = 128 cells at 6,300 s, `rng = Xoshiro(18_438_301)` (M4's block, reused so the
before and after compare like for like). Nothing reads an evidence posterior.
- **"chunked"** is the frozen `_each`, `Threads.@threads` over the cells. **"greedy"**
  is the new `_each`. **"predict"** is one full `predict` call on the matrix, KIM-38's
  measurement.
- At 16 and at 64 threads (default GC), both passes' Q equals `predict`'s. No cell
  failed in any run.
- **Efficiency** = 128 × the node's serial median / (threads × wall). The node's serial
  median is that of the matrix's first 2 cells, run one at a time on the same node
  first.
- **The idle share** = 1 − Σ cell time / (threads × wall). **GC's share of the wall**
  is the process's own GC clock over its wall.

## Results

| variant | threads | GC threads | loop | wall (s) | efficiency | serial median (s) | per-cell slowdown | GC share of cell time | GC share of wall | GC pauses | mean pause (s) | idle share | MaxRSS (GiB) |
|---|---|---|---|---|---|---|---|---|---|---|---|---|---|
| default | 16 | 8 | chunked | 564 | 0.36 | 25.16 | 2.7 | 0.46 | 0.45 | 1,915 | 0.13 | 0.06 | 5.2 |
| default | 16 | 8 | greedy | 551 | 0.37 | 25.16 | 2.7 | 0.46 | 0.46 | 1,873 | 0.14 | 0.03 | 5.2 |
| default | 16 | 8 | predict | 553 | 0.36 | 25.16 | | | | | | | 5.2 |
| default | 64 | 32 | chunked | 453 | 0.12 | 26.09 | 8.1 | 0.78 | 0.77 | 293 | 1.19 | 0.05 | 29.9 |
| default | 64 | 32 | greedy | 464 | 0.11 | 26.09 | 8.6 | 0.80 | 0.79 | 277 | 1.33 | 0.04 | 31.4 |
| default | 64 | 32 | predict | 477 | 0.11 | 26.09 | | | | | | | 31.7 |
| `--gcthreads=32,1` | 64 | 32 + 1 sweep | greedy | **166** | **0.31** | 26.17 | 3.0 | 0.30 | 0.30 | 284 | 0.17 | 0.05 | 29.9 |
| `--gcthreads=1` | 64 | 1 | greedy | 441 | 0.11 | 24.86 | 8.5 | 0.77 | 0.77 | 286 | 1.18 | 0.04 | 31.4 |
| `OPENBLAS_NUM_THREADS=1` | 64 | 32 | greedy | 446 | 0.13 | 28.11 | 7.6 | 0.77 | 0.77 | 291 | 1.18 | 0.04 | 31.0 |
| 4 processes × 16 | 64 | 8 each | greedy | **135** | **0.37** | 24.83 | 2.6 | 0.45 | 0.44 | 2,238 | 0.11 | 0.06 | 4.3 each |

**Serial** (dave3, 8 cells, one at a time): 25.90–31.15 s per cell, median 28.88 s; GC
0.027 of the time; 24.76 GiB allocated per cell (median).

## Reading

1. **The time is lost inside the cells, to GC. Idle threads and the scheduler are not
   the cause.**
   - Every pass keeps its threads busy: the idle share is 0.03–0.06.
   - A cell takes 2.7× its serial time at 16 threads and 8.1–8.6× at 64. GC is 0.46 and
     0.78–0.80 of the cells' time, against 0.027 serial.
   - The process spends the same share of its wall in GC (0.45–0.46; 0.77–0.79). A GC
     pause stops every thread, so each pause costs all of them.
   - The chunked loop and the greedy `_each` are within 3% of each other, one way at 16
     threads and the other at 64. The new scheduler is about errors (the lowest failing
     index's, thrown as itself; the rest of the work stopped), not speed. `predict`'s
     own wall matches the passes'.
2. **BLAS is ruled out**, as the PR #123 review read it: one BLAS thread gives 446 s
   against 464 s.
3. **Parallel marking is not the cost.** One GC thread (441 s) does as well as 32
   (464 s).
4. **The sweep is.** `--gcthreads=32,1` adds one concurrent sweep thread to the default
   32 mark threads, and nothing else.
   - At 64 threads the wall falls from 464 s to 166 s, and GC's share of it from 0.79
     to 0.30.
   - The pause count stays the same (284 against 277). The mean pause falls from 1.33 s
     to 0.17 s.
   - This is one run, on one node. The node baselines span 24.9–28.1 s, about 13%,
     against a 2.8× change.
   - Not measured: `,1` at 16 threads (where the default pause is already 0.14 s, but
     GC is still 0.46 of the wall); fewer mark threads with `,1` (64 + 33 threads
     oversubscribe 64 cores); and `predict`'s own call under the flag.
5. **Separate processes avoid the 64-thread heap, not GC.** 4 processes × 16 threads
   on one node finish in 135 s, efficiency 0.37. Each process's GC is 0.44 of its wall,
   as in one 16-thread process, and 551 s / 4 = 138 s. Each heap peaks at 4.3 GiB,
   against about 30 GiB for one 64-thread process.
6. **Allocation is the source.** A 6,300-s cell allocates 24.76 GiB, and a 128-cell
   pass 3,251.5 GiB in every variant. One 60-s cell (job 19583243, `Profile.Allocs` at
   sample rate 0.01) allocates 358.4 MiB, nearly all at three sites:
   - `_exchange!` at `handshake.jl:1833`, `step!(d.ode, …)`: 172.0 MiB, mostly
     `Matrix{Float64}` (168.3 MiB in about 23,400 allocations, about 7.5 KiB each);
   - `handshake_step!` at `handshake.jl:1785`, `reset_aggregated_jumps!`: 93.1 MiB;
   - `handshake_step!` at `handshake.jl:1787`, `step!(d.jump, …)`: 88.6 MiB.

   The last two are mostly closures over a `Reaction` and its `_JumpSlot`, the jump
   closures that `_global_jump` builds (`orchestrator.jl`), about 105 MiB. These sites
   are reported only (the user's call): KIM-71, related to KIM-68 and KIM-42.
   `handshake.jl` is not edited here.

**Against KIM-38.**
- The walls for the same 128 cells agree: 453–477 s at 64 threads against KIM-38's
  454 s, and 551–564 s at 16 against 619 s.
- The efficiencies differ because the serial baselines do. KIM-38 timed 10 cells at one
  parameter set, `θ_ref`, at `a742224` on dave7: a median of 36.8 s. This run times the
  matrix's own cells on each node: 24.9–28.1 s.
- With 36.8 s, these walls give 0.52–0.53 at 16 threads and 0.15–0.16 at 64, KIM-38's
  figures.
- Whether the gap is the parameter set, the commit or the node is not tested here.

## What it means for M5 (for the gate)

The freeze's "about 24 h for each width on one node" (spec §12, 2026-10-30) is KIM-38's
64-thread wall, scaled. This run reproduces that wall at the default (464 s, greedy).
Scaled by the measured walls, on the assumption that KIM-42's cells cost like these:

| setting | wall, 128 cells | per width, one node | core-h, three widths (frozen: 3,000–4,500) | needs |
|---|---|---|---|---|
| default, 64 threads | 464 s | about 24 h | 3,000–4,500 | nothing |
| `--gcthreads=32,1`, 64 threads | 166 s | about 8.6 h | about 1,100–1,600 | a launcher flag |
| 4 processes × 16 threads | 135 s | about 7 h | about 870–1,300 | the matrix split over processes and combined (KIM-40's combiner), and `predict`'s rules on splitting rows under completion |

- M0 (KIM-40, KIM-41) runs single-core tasks, at the serial GC share of 0.027. This
  diagnosis does not change M0's estimate.
- The user's call at the gate, and any confirmation run, are recorded in spec §12.
