# KIM-37: the ptsG knockdown dose screen (M4 pilot)

**Run:** `dev/scripts/dose_screen_kim37.sh` at `a742224` (pushed), into
`/fred/oz022/tkimpson/infercell-m4-kim37-2026-10-30/` (`SHA256SUMS`, `jobs.tsv`).
- Jobs: precompile 19560145; M0 19560146; reduced model 19560147; merge 19560148.
- The merge's tables are in `dose_screen.md` in that directory.

**Held out:** θ comes from fresh prior draws, the population the SBC truths came
from. The seeds are in the script's header. The screen computes no ΔQ. The
archive's `SHA256SUMS` has SHA256 `81fa1126e07ebb08…`.

**Design and criterion** (declared in the script's header before the run):
- **Draws:** M0 at 600 s, 50 draws × 4 cells; the reduced model at 6,300 s,
  20 draws × 2 cells.
- **Doses** 0 to 0.5, on the same seeds at every dose.
- **Watched:** ATP and GTP against K_d = 0.1 mM, and charged tRNA against its
  per-amino-acid K_d (0.02 mM lumped). Each is read at every 60-s rebuild instant.
- **A dose is admissible** when, for each watched pool:
  1. its e95, the 95th percentile of the pool's largest elasticity, is at most
     max(0.1, 1.1 × the dose-0 e95);
  2. its fraction of instants below K_d is at most the dose-0 fraction + 0.01.
- **The pick** is the largest admissible dose in {0.1, 0.2, 0.3}, judged on M0.
- **Added in the script, not in the approved plan:** judging the reduced model
  with and without its baseline-starving draws, for information. It changes no
  pick.

## Results

**M0** (1,200 cells; none failed, none starved). Every dose up to 0.5 is
admissible. **The rule picks 0.3.**

| dose | ATP e95 | GTP e95 | charged tRNA e95 | GTP below K_d | charged tRNA below K_d |
|---|---|---|---|---|---|
| 0 | 0.0230 | 0.0543 | 0.1001 | 0.0005 | 0 |
| 0.1 | 0.0231 | 0.0551 | 0.1007 | 0.0005 | 0 |
| 0.2 | 0.0230 | 0.0547 | 0.0986 | 0.0005 | 0.0005 |
| 0.3 | 0.0231 | 0.0541 | 0.0983 | 0 | 0 |

The limits are 0.1, 0.1 and 0.1101 for the e95s, and 0.0105 and 0.01 for the
below-K_d fractions.

- **The knockdown hardly moves the watched elasticities on M0.**
  - Up to 0.3, they change by under 2% within 600 s. At 0.4 and 0.5 they fall by
    up to 3.1% (GTP) and 4.6% (charged tRNA).
  - Against the measured 0.04–0.09 band: ATP's (0.023) lies below it, GTP's
    within it, and charged tRNA's (0.096–0.101) just above its top.

**The reduced model** (240 cells, 40 per dose; none failed; 13 starved, from 4
draws). Every dose up to 0.5 is admissible, with or without the starving draws.
**The rule would also pick 0.3.**
- **Which draws starved:**
  - draw 19 at baseline and at every dose but 0.2 (9 cells);
  - draw 14 at baseline (one cell) and at 0.2 (the other);
  - draws 16 and 20 at 0.4 only.

  Draws 14 and 19 are the two excluded in the second table.
- **No trend with dose.** On the same seeds, the starved counts per dose are 3,
  2, 1, 2, 3 and 2. At 0.3 only draw 19 starves, and it starves at baseline. So
  the screen shows cells near the starvation threshold, with no dose trend
  detectable at 2 cells per draw per dose. *(Corrected 2026-10-31, spec §12
  KIM-70.)*
  Under V at baseline, intervened failures are part of the prediction, and are
  counted.
- **The reduced model's baseline is already near K_d, at dose 0.** GTP is below
  K_d at 19% of instants (17.5% without the starving draws), with an e95
  elasticity of 0.27. Charged tRNA is below at 6.6% (3.5%), with an e95 of 0.31.
  - These draws take θ_I from its prior too, not from the M0 posterior M5
    carries.
  - The fraction varies strongly with θ: GTP is below K_d at 0% of instants in
    6 of 20 draws, 1% in 2, and 21–55% in 12.
  - So it is not measured for M5's population.
- The criterion passes because it is relative to the baseline. The absolute
  basis of the one-way reduction, elasticities of 0.04–0.09, does not hold
  across the reduced model's full cycle under prior draws, even before any
  knockdown.

## Recommendation (the freeze decides)

- **Freeze the dose at 0.3**, the rule's pick on M0. The criterion's numbers are
  above.
- **State it in the paper as a model-form caveat** for the reduced-model
  results (KIM-47's table), whatever the dose. Under prior draws of every
  parameter, the one-way reduction's basis is weaker at baseline: GTP is below
  K_d at about a fifth of the rebuild instants.
- **What M5 can check.** `predict` keeps each cell's smallest GTP, not the
  fraction of instants, so KIM-42 reports the fraction of its cells whose GTP
  ever fell below about K_d (the freeze, §12; its thresholds and cells are
  declared in §12 2026-10-31, KIM-70).
