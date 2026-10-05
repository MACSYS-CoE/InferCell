# 16a.7c: the path update under bulk metabolites

Code at `50e5af4`, with runs at `22d1f3e`.

**The tests pass,** 211 of 211 (job 18049605). They include:
- the bulk record recomputing from the latent and its seed;
- block 1 scoring the cells' mean;
- the CSMC bulk weight equalling the bulk likelihood with the particle in place;
- a sweep's returned latent equalling its path's replay.

## Update rates: KP no longer fires (jobs 18052532 to 18052536; per-cell merge 18098390 at `516b3b5`)

**The setup:**
- **The population:** M0 replicate 16085's 50 cells, observed in bulk at the
  replicate's σ = 0.0727 as σ_b.
- **The update:** cells 1 to 5 are updated in turn by PG, each conditioned on the
  other 49. The rest are held at their true paths.
- **The run:** 16 scans, so 16 sweeps per cell and 80 cell-sweeps per N, from the
  true paths.
- **Cost** is per cell-sweep, in plain M0 cycles timed in the same job, with the
  first scan dropped.

**Corrected after the #81 review.** The first table pooled the rate over the five
cells. The per-cell rates show windows that never change: cell 2's windows 7 to
10, and cell 5's windows 4 to 8.
- **Each of those has an empty true path,** transcripts 0 → 0, whose posterior sits
  on the empty path, so a near-zero rate is expected there.
- **D16.3's criterion is applied per cell, over the other 41 cell-windows.**

| N | Cycles per cell-sweep | Pooled: lowest window rate | Per cell: lowest rate over live windows | Live cell-windows below 10% |
|---|---|---|---|---|
| 5 | 4.2 | 0.05 | 0.00 | 12 of 41 |
| 10 | 8.4 | 0.28 | 0.12 | 0 of 41 |
| 20 | 16.5 | 0.55 | 0.50 | 0 of 41 |
| 50 | 40.3 | 0.59 | 0.62 | 0 of 41 |

**Chosen: N = 20.** N = 10's lowest live rate, 0.12 from 16 sweeps (SE about 0.08),
does not robustly clear 10%. N = 20's, 0.50, does. Cost is about 0.84N cycles per
cell-sweep. The per-cell, per-window rates are in the merge log.

## V5's bulk "on" half: passes, with a watch item (prep, kernel and merge jobs 18052537 to 18052542)

**The setup:** three cells share the data's transcripts. The reference is 200,000
random triples of kept runs, weighted by the bulk likelihood of their mean. The
kernel draws are 2,000 chains of five scans at N = 5. The functionals are averaged
over the three cells.

| Case | σ_b | Reference ESS | z: births w1, w2, w3 | z: translation w1, w2, w3 |
|---|---|---|---|---|
| Toy | 0.1 (the data's) | 23,353 | 2.51, 0.02, 0.33 | −0.19, −0.19, 1.60 |
| M0 | 0.8 (inflated 8×) | 11,588 | 0.89, **2.85**, **−2.81** | −0.81, 2.29, −0.79 |

**Both pass at 3 SE.** These z-scores used an ESS-based reference SE, which the
fresh-reference check below found too small. `bulkmerge` now bootstraps over runs,
from `516b3b5`.

**Watch item: births in window 2 on M0.** In every M0 V5 comparison so far, the
kernel puts more births in window 2 than the reference does.

| Run | Kernel | Reference | z |
|---|---|---|---|
| Per-cell, metabolites off | 0.018 | 0.013 | about 1.6 |
| Per-cell, metabolites off, annealed | 0.018 | 0.013 | about 1.6 |
| Per-cell, metabolites on | 0.015 | 0.011 | 1.38 |
| Per-cell, metabolites on, annealed | 0.013 | 0.009 | 1.24 |
| Bulk | 0.013 | 0.009 | 2.85 |

- **The comparisons are not independent.** Every one uses the same rejection
  reference set, so one low draw of that set would make them all agree.
- **The bulk run's −2.81 for births in window 3** is the mirror image: windows 2
  and 3 trade a birth.
- **What would decide it:** a fresh rejection reference with new seeds, about 30
  CPU-h.

## Still to do for 16a.7c

- **V4 on block 1's bulk form.** Block 1's bulk term is implemented and unit-tested,
  but the grid check has not been rerun.
- **The full-scale half of parent 15.7b.** 15.7's bulk observations need the latent
  record regenerated, about 1.8 CPU-h.

## The watch item resolved: a fresh, independent V5 reference (2026-10-05)

**The fresh reference** (`V5_FRESH=1`, jobs 18060205 and 18060210 at `cfd41dd`):
a second rejection set from seeds offset by 50,000,000. M0 kept 8,884 runs of
100,000, and the toy 11,232 of 500,000.

**The checks** (`refcheck`, jobs 18060206 to 18060213) tested the earlier kernel
draws, which do not depend on the reference. The pooled checks (`V5_POOL`, jobs
18061130 to 18061135 at `1190cdc`) then ran against the fresh and original sets
combined: 17,594 M0 runs and 22,156 toy runs.

**A flaw in the first bulk merge, fixed at `1190cdc`.**
- **The flaw:** its reference SE was ESS-based over 200,000 triples drawn from
  about 4,400 runs. The triples share runs, so that SE was far too small. The first
  reference read window-2 births at 0.009 with an SE under 0.0005; the fresh one
  read 0.013.
- **So the first bulk z-scores were overstated,** including 2.85 and −2.81.
- **The fix:** the reference SE is now a bootstrap over runs.

**Against the pooled reference:**

| Sampler | M0 | Toy |
|---|---|---|
| PG per cell, metabolites off | χ² p ≥ 0.049 | χ² p ≥ 0.058 |
| PG per cell, metabolites on (σ inflated as before) | \|z\| ≤ 1.22 | \|z\| ≤ 1.89 |
| **PG bulk, the production kernel** (bootstrap SE) | \|z\| ≤ 1.79 | \|z\| ≤ 1.16 |
| Annealed, K = 5 (off by default) | \|z\| ≤ 1.48, except translation w2 on at **−2.90** | \|z\| ≤ 1.05 |

- **Window-2 births on M0 agree with the pooled reference:** z = 0.09 per cell and
  0.81 bulk. The watch item was the first reference set coming out low.
- **The production kernels pass V5 cleanly** on both cases.
- **The pooled reference contains the original kept set,** whose first half supplied
  the kernel's starts. So the pooled checks are not fully independent; z is
  understated by about 10% at most.
- **The fresh-only checks of the production kernels (PG per cell and bulk) are fully independent** (jobs 18060206, 18060209, 18060211 and 18060213):
  - M0: off χ² p ≥ 0.032, per-cell on |z| ≤ 1.64, bulk |z| ≤ 2.07 (ESS-based SE);
  - toy: off χ² p ≥ 0.047, per-cell on |z| ≤ 1.28, bulk |z| ≤ 1.78.
- **To regenerate the pooled checks:**
  `POOL=<original run>/csmc_v5_16a7/m0 dev/scripts/csmc_v5_fresh_16a7c.sh m0 <original>/m0 <annealed>/m0 <bulk>/m0`.
  The same for toy, from `516b3b5`. The three directories are the kernel runs of
  `csmc_v5_16a7.sh`, `csmc_v5_anneal_16a7.sh` and `csmc_v5_bulk_16a7c.sh`.
- **New watch item, annealed kernel only:** M0's metabolite-on translation in
  window 2 reads low against every reference: −1.97, −3.14 fresh and −2.90 pooled.
  The annealed window is off by default and not used in production. It must be
  investigated before it is ever switched on.

## V4 on block 1's bulk form: passes (jobs 18060214 to 18060217, 18061281, 18061355 to 18061363)

**The setup:** M0 replicate 16085, three cells. Their panel is observed in bulk, one
measurement per pool per save of their mean, at σ = 0.0727, held. The grid's
likelihood averages the three replays and scores the mean by its own code (`V4_BULK=1`).

**What the grids needed:**
- **Grids 1 and 2** ran at `cfd41dd`.
- **Merge 2 found mass at both axes' upper edges** (ENO 7.1e-6, FBA 1.5e-5). FBA's
  posterior has a long right tail under bulk data. The script refused a truncated
  FBA axis, so it was generalised at `5c40835`: grid 3 is recomputed whole, 81 × 81,
  on ln ENO 3.80 to 4.63 and ln FBA 1.08 to 6.92.
- **Merge 3** (`4e4a9b6`): ln ENO 4.2147 ± 0.0264 and ln FBA 4.0007 ± 0.1562, with
  edge mass 6e-47 and 5e-17.
- **The widened grid has about two points per ENO SD.** So ENO's exact quantiles come
  from a new ENO-fine grid (161 × 41, edge mass at most 1.9e-7), as FBA's come from
  the FBA-fine grid.

**The target check:** block 1's bulk target and the grid's log posterior differ by a
constant to 1.7e-12 over 15 points.

**The final test:** 20 chains × 45 sweeps after burn-in, σ held. Chain fraction below
each exact quantile:

| Constant | 5% | 25% | 50% | 75% | 95% |
|---|---|---|---|---|---|
| ln ENO | 0.042 ± 0.007 | 0.262 ± 0.014 | 0.491 ± 0.018 | 0.742 ± 0.014 | 0.941 ± 0.010 |
| ln FBA | 0.062 ± 0.009 | 0.281 ± 0.012 | 0.490 ± 0.016 | 0.742 ± 0.013 | 0.934 ± 0.008 |

**All ten are within 3 SE. V4 passes.**
- **The largest deviation is FBA at its 25% point, +2.57 SE.** It was misquoted at
  first as the 95% point's −1.95.
- **FBA's chain fractions sit above their levels at 5% and 25%,** as in the per-cell
  V4. That watch item carries to V6 and V7.
- About 45 CPU-h.

## 16a.7c: met

| Item | Result |
|---|---|
| M0 bulk datasets | `m0_dataset` |
| Full-scale bulk dataset | parent 15.7b |
| B1's bulk form | V4 passes |
| B3's bulk weight | unit-tested |
| V5, both cases | passes against the pooled reference |
| Update rates | every window above 10% from N = 10 |

**Chosen: PG, lag 0, N = 20** (corrected from N = 10 after the #81 review).
