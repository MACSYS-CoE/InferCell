# KIM-38: the cost of the evidence package, and the gate (M4 pilot)

**Run:** `dev/scripts/cost_pilot_kim38.sh` at `a742224` (pushed), into
`/fred/oz022/tkimpson/infercell-m4-kim38-2026-10-30/` (`SHA256SUMS`, `jobs.tsv`).
- Jobs: precompile 19560132; cells 19560133; threads 19560134 (16) and
  19560135 (64); sizing 19560136; projection 19560250.
- The projection reads KIM-35's arms. It is `projection.md` in that directory.

**Held out:** B7's `v6` fit, and seeds from M4's block.

## Per run (warm, one thread, median)

| run | median | p90 |
|---|---|---|
| M0 hybrid, 600 s | 1.74 s | 2.25 s |
| M0 hybrid, 600 s, 25% knockdown | 1.67 s | 2.23 s |
| M0 mean field, 600 s | 0.39 s | 0.41 s |
| reduced-model hybrid, 6,300 s | 36.8 s | 37.8 s |
| reduced-model hybrid, 600 s | 3.99 s | 4.14 s |
| reduced-model mean field, 6,300 s | 82.2 s | 83.1 s |

**Threading** (KIM-62), on a 128-cell reduced-model matrix at 6,300 s:
- 16 threads: 619 s, an efficiency of 0.48, with MaxRSS 4.9 GiB;
- 64 threads: 454 s, an efficiency of 0.16, with MaxRSS 22.1 GiB.

The results are the serial ones, but threading scales poorly. A likely cause is
allocation and GC contention, or BLAS threads on top of Julia's. It is not
diagnosed here.

## Sizing (declared rules; M0 on v6, 20 draws × 10 cells)

| Q | Var_θ[μ] | σ² (cell to cell) | n_cells (MC ≤ 10%) | n_ref = n_truth (SE ≤ 10% of SD) |
|---|---|---|---|---|
| `lactate_out` | 1.87e8 | 2.18e9 | 106 | 1,167 |
| `exported` | 8.9e−7 | 1.7e−5 | 174 | 1,927 |

On M0 the posterior predictive is tight next to the cells' spread. The SD over
draws is 13.7k on a mean of 480k, and the SD between cells is 46.7k. So the rules
ask for many cells. On the reduced model at 6,300 s (KIM-35) the spread over
parameters dominates, and n_cells = 4 suffices.

## Projection (Q = `lactate_out`)

| part | CPU-h | wall |
|---|---|---|
| M0 crossed comparison (KIM-40) | 598 | 6.0 h at 100 concurrent tasks |
| intervention panel (KIM-41) | 575 | 5.7 h at 100 tasks |
| reduced-model budget, 6,300 s, n_cells 4, 3 widths (KIM-42) | 245 | 7.9 h per width, one 64-thread node each |
| … with V reading its own cells | 491 | 15.8 h per width |
| **total** | **1,420 to 1,660** | **under a day, if the widths run on 3 nodes at once** |

Options not costed above:
- an intervened arm on the reduced model, which adds up to the budget's baseline
  again;
- `exported`'s larger sizes, which add about 750 CPU-h.

## The gate

- **The scale is confirmed. No cut is needed.**
- Started by about 13 November, M5's package finishes in about two days of
  wall-clock, with a margin for the queue and a smoke pass. That is well before
  the gate's 20 November. Neither the intervention panel nor a prior width is
  cut.
- **Caveats:**
  - The sizes come from one held-out posterior (v6). Posterior widths vary across
    the 50 datasets, and the M0 MC rule may be loose or tight on some.
  - The concurrency of 100 tasks is an assumption. B7 ran about 150 at once.
