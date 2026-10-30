# KIM-38: the cost of the evidence package, and the gate (M4 pilot)

**Run:** `dev/scripts/cost_pilot_kim38.sh` at `a742224` (pushed), into
`/fred/oz022/tkimpson/infercell-m4-kim38-2026-10-30/` (`SHA256SUMS`, `jobs.tsv`).
- Jobs: precompile 19560132; cells 19560133; threads 19560134 (16) and
  19560135 (64); sizing 19560136; projection 19560250.
- The projection reads KIM-35's arms. It is `projection.md` in that directory.
- **The projection was submitted by hand**, after its inputs, with
  `sbatch … --dependency=afterok:19560133:19560134:19560135:19560136:19560138:19560139:19560140:19560141 … --wrap "… julia --project dev/scripts/cost_pilot_kim38.jl project /fred/oz022/tkimpson/infercell-m4-kim38-2026-10-30 /fred/oz022/tkimpson/infercell-m4-kim35-2026-10-30"`.
  The launcher's `project` mode submits the same thing without the dependency.
- **The nodes:** cells on dave7, 16 threads on dave23, 64 threads on dave309,
  sizing on dave16.
- **The archive's `SHA256SUMS`** has SHA256 `892e3c4b3ffcf995…`. It was
  remade after `sizing_resolved` added its files; the first run's was
  `1c543c4f87054dd1…`.
- **`sizing_resolved`** (job 19568624, at `82c14d6`): added after the PR #123
  review. Its output, `sizing_resolved.jls`, holds plain numbers and opens at any
  later commit.
- **`sizing.jls` opens only at `a742224`.** It stores whole `Prediction`s, and
  KIM-36 added fields to that type afterwards. So the projection stage reruns
  from a worktree at `a742224`. `projection.md` carries its numbers.

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

These runs only time the work. That a threaded result equals the serial one
rests on `test_predict.jl`'s threaded-equals-serial testset, at the suite's 4
threads (job 19560130).
- Threading scales poorly.
- The PR #123 review read the accounting: user CPU is flat from 16 to 64
  threads, while system CPU grows. That points to allocation and GC, with an idle
  tail, rather than BLAS.
- It is not diagnosed here. KIM-68 does that, as M5's first job.

## Sizing (declared rules; M0 on v6)

**The first sizing, 20 draws × 10 cells (job 19560136): superseded.**

| Q | Var_θ[μ] | σ² (cell to cell) | n_cells (MC ≤ 10%) | n_ref = n_truth (SE ≤ 10% of SD) |
|---|---|---|---|---|
| `lactate_out` | 1.87e8 | 2.18e9 | 106 | 1,167 |
| `exported` | 8.9e−7 | 1.7e−5 | 174 | 1,927 |

The script reported no SE for Var_θ. The PR #123 review bootstrapped it: 1.87e8
is 1.2 SE from zero. So the n it gave was not determined by its own rule.

**The resolved sizing, 200 draws × 40 cells (job 19568624, `sizing_resolved`),
with the mean field on the same draws:**

| Q | Var_θ[μ] ± SE | mean field's var | σ² | n_cells | n_ref = n_truth | M0 CPU-h (comparison + panel) |
|---|---|---|---|---|---|---|
| `lactate_out` | 2.81e7 ± 0.91e7 (3.1 SE) | 1.65e7 | 2.77e9 | 887 | 9,849 | 4,999 + 4,808 |
| `exported` | 6.6e−8 ± 4.7e−8 (1.4 SE) | 5.6e−8 | 1.6e−5 | 2,171 | 24,121 | 12,235 + 11,767 |

- **On M0 the posterior predictive of μ is very tight next to the cells'
  spread.** The SD over draws is about 5.3k on a mean of about 480k, and the SD
  between cells is about 52.6k. So the declared rules ask for many cells.
- **The hybrid's spread over θ is close to the mean field's** on the same draws.
- On the reduced model at 6,300 s (KIM-35) the spread over parameters dominates,
  and n_cells = 4 suffices.

## Projection (Q = `lactate_out`)

| part | CPU-h | wall |
|---|---|---|
| M0 comparison (KIM-40), at the resolved sizes | 4,999 | about 1 day at 200 single-core tasks (rows split over tasks) |
| M0 knockdown panel (KIM-41), run in the same calls | 4,808 | about 1 day at 200 tasks |
| reduced-model budget, 6,300 s, n_cells 4 + V's 4, 3 widths (KIM-42) | 491 | about 16 h per width, one 64-thread node each |
| reduced-model knockdown arm, on the kept pairs | ≤ 245 | about 8 h more per width |
| **total** | **about 10,500** | **about 2 days at 200 cores, 4 at 100** |

At the first, superseded sizing it was about 1,900.

## The gate

- **Not decided: M0's sizes are a todo (KIM-69).** On 2026-10-30 the user chose
  to draft the paper first and revisit, rather than wait on the run.
- **The options and their cost:**
  - keep the declared rules (about 10,500 CPU-h);
  - cut the panel (about 5,700);
  - halve the draws (about 6,400);
  - relax the rules (about 2,600, as a recorded post-pilot change).
- **Against the gate's 20 November:** at the declared rules, the package
  finishes before it if M5 starts by about 12–14 November.
- **Caveats:**
  - The sizes come from one held-out posterior (v6). Posterior widths vary across
    the 50 datasets, and the rules may be loose or tight on some.
  - The concurrency is an assumption. B7 ran about 150 tasks at once.
