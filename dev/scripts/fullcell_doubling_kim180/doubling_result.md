# KIM-180: the doubling-time definition on the default run

Linear: [KIM-180](https://linear.app/kimpson/issue/KIM-180/freeze-the-block-selection-rule-on-doubling-time).
The definition is `spec/prereg/kim-180.md` §1. The check (§5) is that at the published
parameter values it gives a doubling time inside Thornburg et al. 2022's range, 88 to 112
min, with a mean of 97 (Cell 185:345, Fig. 3F).

**The check holds.** Both default ensembles give a median of 99 min.

| ensemble | cells | healthy | of them reach t_SA2 | doubling time |
|---|---|---|---|---|
| syn3a-ensemble's default, 120 min (`results/cycle/cells.csv`, `2075809`) | 10,000 | 8,469 | 8,434 | 99 min |
| KIM-180's default, 180 min (job 21839731) | 100 | 84 | 84 | 99 min |
| the same, cut at 120 min | 100 | 84 | 84 | 99 min |

**The run.**
- Job 21839731, after the precompile job 21839730, on dave99 (AMD EPYC 7543).
- InferCell commit `d376bdf` and syn3a-ensemble `2075809`, whose code and data were
  clean.
- KIM-183's driver, `dev/scripts/fullcell_timing_cycle/timing_cycle.jl`, unchanged: procids
  1 to 100, seed 1, `MINUTES=181`, 32 threads. It took 8 min 38 s of wall time.
- The outputs are in `/fred/oz022/tkimpson/infercell-kim180-2026-11-28/ensemble/`. The
  per-cell table and the run's metadata are copied to `results/ensemble/` here; the
  per-minute record (670 KB) is not.

**What the run shows.**
- No cell threw or ended in an ODE failure.
- 25 cells collapsed. 16 of those collapsed before their t_SA2, so they are unhealthy.
- 5 cells reached t_SA2 after minute 120. All 5 had collapsed first.
- **These are not new cells.** Up to minute 120 their t_SA2 and t_short match KIM-183's
  ensemble (job 21664422) and syn3a-ensemble's rows for procids 1 to 100, seed 1, exactly.
  So the second row checks determinism and minutes 120 to 180, not a fresh sample.
- At the published values, the 180-min window therefore censors no healthy cell that the
  120-min window keeps. The window matters only under the prior, where syn3a-ensemble's
  demonstration found 58 of 245 draw medians past minute 120.

**Reproduce:**
- `dev/scripts/fullcell_doubling_kim180/doubling.sh <outdir>` from the worktree root at a
  pushed commit;
- then `python3 -I dev/scripts/fullcell_doubling_kim180/doubling_reads.py <outdir>`, or with
  `dev/scripts/fullcell_doubling_kim180/results` in place of `<outdir>`. It refuses a
  syn3a-ensemble `cells.csv` other than the one at `2075809`.
