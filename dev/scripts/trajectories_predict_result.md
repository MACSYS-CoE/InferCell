# KIM-87: probabilistic trajectories for paper 1's figure. UNSCORED

**No new scored claim.** Panel A replays M0's scored prediction with the time series recorded.
Panels B and C are on the reduced model, with no truth compared. The runs add no inference.
The design is §12's 2026-11-08 entry, "Probabilistic trajectories: runs and figure,
pre-registered".

**Run:** `dev/scripts/trajectories_predict.sh` at `554ee2b` (main, #158's merge), from a
detached worktree, into `/fred/oz022/tkimpson/infercell-trajectories-2026-11-08/`. That
folder is the archive: `replay.md`, `checks.md`, `bulk_check.md`, `tsv/` (the figure's six
TSVs), the `.jls` outputs, `logs/`, `jobs.tsv`, `sacct.tsv` and `SHA256SUMS`. The `.jls`
files open only at `554ee2b`. Julia came from the cluster's `julia/1.10.5` module.

- **Jobs:** first half precompile 20271496, A 20271497 (array 0–10: the 10 blocks and the truth),
  B 20271498, C's pre-check 20271499, data 20271500, replay 20271501. Second half precompile
  20293796, C 20293797, report 20293798. Bulk check precompile 20293799, bulk 20293800, run
  from a worktree at `dfe570a`. All completed with exit 0, on `znver3` nodes.
- **Order:** `rest` was launched only after `replay.md` ended in `ALL PASS`.

## The checks

All passed. `replay.md` (job 20271501) and `checks.md` (job 20293798) agree line for line.

- **Panel A:** each of the 10 blocks' μ̂ equals its rows in M0's `finish/000_r1_v7.jls`
  exactly, and the truth's mean equals M0's published truth, 597951.252683243 in both.
- **Panel B:** μ̂, the seeds and θ_H equal KIM-42's width 1 `w1.jls`, exactly.
- **Dataset:** the regenerated `m0_dataset(16_200_001; ncells = 50)` has the truth and σ_b of
  replicate 1, with a largest relative difference of 0.0 for both.
- **Ratio:** exported lactate over `M_lac__L_e`'s change is R = 99,921, with a largest
  relative deviation of 0.245% over 100 cell-times (10 truth cells × 10 times), against the
  1% bound. That is about 0.14 of Q's predictive SD, using §12's 1% = 0.57 SD. It is larger
  than the smoke pass's 0.09% because that pass ran 240 s, not 600 s. The overlay is
  rescaled by the single R, so its points can be off by up to this much.
- **Bulk data at V7's commit:** `bulk_check.md` reports the regenerated bulk data equal
  `dfe570a`'s exactly, largest relative difference 0.0. The caption needs no caveat on this.
- **Bands:** every row of `traj_a.tsv`, `traj_b.tsv` and `traj_c.tsv` is finite and ordered
  lo90 ≤ lo50 ≤ hi50 ≤ hi90. The mean falls outside the 90% band in 25 rows, all at t = 0
  for `M_gtp_c` and `M_pyr_c`, where the band has zero width and the mean differs from it
  only in rounding.

## Panel B: layer 2 and the held share

- **Layer 2's Monte Carlo share at 6,300 s:** mc = 1.418e11 against var = 2.213e12, so
  mc / var = 0.064.
- **The final held share, I first:** 0.9080, which reproduces KIM-42's headline of 0.908 at
  width 1. With H first it is 0.898.
- **What the figure's rule shows:** of the 105 times after t = 0, 99 pass (`flagged` and
  `mc_undefined` false, `total > mc`). The hidden ones are 60–360 s, where the parameter part
  is smaller than the cells' Monte Carlo noise. Nothing is flagged and nothing is dropped. The
  shown share runs from 0.838 to 0.908.

## Panel C: the released recipe

82 drawn parameters (5 inferred, 15 held, 62 released) on 40 × 50 pairs.

- **Pre-check** (4 × 5, seed `19_543_001`): acceptance 0.6, 48 of 80 pairs. It took 16 min
  49 s on 16 threads.
- **Acceptance:** 0.623, 1,246 of 2,000 pairs kept. Every θ_I row keeps between 0.60 and
  0.66 of its pairs, and no row failed. That is far above the 10% cut, so the figure shows
  the filled conditioned band.
- **`per_theta_I` beside `joint`:** the same acceptance, with ESS 1,245.6 against 1,246. With
  row acceptance this even, the two weightings are nearly the same.
- **Starved:** 3,030 of μ̂'s 8,000 cells are starved, a fraction of 0.379. The rest are 3,175
  clean and 1,795 clipped, with none failed or below the floor.
- **Median move,** of pair means at 6,300 s:

| | median | against panel B |
|---|---|---|
| Panel B, KIM-42's width 1, `:prior` | 9.545e6 | |
| Released, unconditioned | 8.454e6 | −11.4% |
| Released, `per_theta_I` | 9.816e6 | +2.8% |
| Released, `joint` | 9.814e6 | +2.8% |

Releasing 62 more parameters with no rule moves the median down by about a tenth. The
viability rule puts it back to within 3% of panel B. Panel B's median is a reference, not a
truth.

## Time and cost

| Stage | wall | core-h allocated | CPU-h used | §12's estimate, core-h |
|---|---|---|---|---|
| A, 11 single-core jobs | 5.07–5.96 h each | 60.6 | 60.1 | 75 |
| B, 64 threads | 3.65 h | 233.9 | 154.0 | 200 |
| C's pre-check, 16 threads | 0.28 h | 4.5 | 3.6 | |
| C, 64 threads | 5.75 h | 367.7 | 237.0 | 400 |
| The rest | | 0.3 | 0.2 | |
| **Total** | | **666.9** | **454.9** | **about 680** |

The allocation is within the estimate, and CPU time used is below the design's 500 CPU-h.
§12's "75 for A" had no stated source. M0's numbers give about 60, and 60.6 is what A took.
The largest MaxRSS was 46.5 GiB, panel B's, against 192G. Core-h are AllocCPUS × elapsed
and CPU-h are `TotalCPU`, both from `sacct.tsv`.
