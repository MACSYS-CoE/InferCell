# T3 at M0: the kill-criteria scoreboard (16b.8)

Published whether or not everything passed (parent `spec/spec.md` §6, T3). Every
verdict is M0's unless labelled otherwise. Verdicts forced by the cap are labelled
budget-bound (16-recovery §8). The two scoring choices for K2 and the forward
direction are the user's, made 2026-10-28 after the results (16-recovery §12
"T3 at M0").

| | Criterion | Threshold | Measured | Verdict | Source |
|---|---|---|---|---|---|
| **K1** | Calibration through feasibility | One production posterior above 1,000 CPU-h. **Set 2026-10-22, after the measurement** | M0: 507.5 to 678 CPU-h per posterior. Full scale: about 21k, projected | **Does not fire at M0. Fires at full scale: budget-bound** | 16-recovery §12 "16a.11"; `f10_16a10_result.md`; `v7_16a9b_result.md` |
| **K2** | The bidirectional claim | Under transcripts only, at 1,000 cells, KL(posterior ‖ prior) below 0.05 nat for every metabolic parameter | Not run at 1,000 cells. **M0's read-off at 50 cells:** ENO 0.0041 and FBA 0.0081 nat, each below its floor | **Budget-bound, not scored.** The firing condition is met in M0's read-off, which is not K2's verdict | `m0_16b_result.md` (job 19395555) |
| **K3** | A calibration failure the architecture causes | Both conditional checks pass and the joint fails its simultaneous band at 95% | The conditional checks are not run in phase 16 | **Not scored here.** Phase 17's | 16-recovery §7, §8 |
| **K4** | The controls | Tight-prior control's coverage at 95% outside the binomial band. **Restated at R = 50: 44 to 50 of 50** (from 90 to 99 of 100). SD / prior SD above 0.9. Equal widths for ENO and FBA | ENO covers 46 of 50 (V7). ENO's SD / prior SD is 0.267 (V6) and 0.258 (16b.1). Widths: ENO 0.042 against FBA 0.184 (V6) | **Does not fire,** on any trigger | `m0_16b_result.md` (jobs 19378382, 19399322) |
| **K5** | The non-smooth drain is in the operating regime | Any deficit at published parameters, or more than 5% of prior draws clipping, at the published 1 s drain | 14b: deficits on 10 of 10 seeds at published parameters (0.33% to 4.0% of drains); 146 of 200 prior draws clip (73.0%). 14c.2, rescored over D11's targets: 163 of 200 (81.5%) | **Fired** (14b, 2026-09-26) | parent spec §11 phases 14b and 14c |
| **K6** | Structural non-identifiability | Rank below six, or condition number above 1e6 | Rank 6 of 6. The smallest singular value is 119, against a split-half floor of 72.3. The condition number is 7.7, capped by the floor at 12.8 | **Does not fire** | parent spec §11 15.8 (job 17704779) |
| **K7** | The reduction is doing the statistical work | PYK's posterior width changes more than twofold with the dropped reactions stubbed in | Not measured. 16b.6 was not funded | **Budget-bound** | 16-recovery §11 16b.6 |

## Notes on the rows

- **K1** was set after the measurement, and the table says so, as the parent
  requires. Calibration ran on the fallback ladder's last rung, M0, and V7
  calibrated across 50 truths, so the headline is not restated.
- **K2.** The M0 read-off is 50 cells against K2's 1,000 and covers ENO and FBA,
  M0's two metabolic targets. At ENO this is the measured null that the
  information bound anticipated (§12 2026-10-05). FBA's evidence was called mixed
  then. At M0 it is not resolved from the prior either.
- **K4** is read at R = 50, not 100, and carries 16a.11 call 1's caveats:
  - 23 of V7's 50 replicates have bulk ESS below 100;
  - L = 100 makes the 99% level coarse;
  - it covers M0's five parameters plus σ_b, not D11's six.
- **K5 was answered twice.**
  - 14c built and validated the smoothed model, whose coupled half failed on GTP
    (3.75% against 1%) because of its continuous pools.
  - D16.1 then made phase 16's block 1 gradient-free on the clipped model. So no
    phase 16 posterior depends on the smoothing, and K5's "gradient-based sampling
    is invalid as posed" does not apply to it.
- **K7.** With no factor measured, the recovery claims carry no conditioning factor
  and none can be assumed to be 1 (`phase16_m0_license.md`).

## Recorded beside T3: the forward direction

F6's load-bearing cell, ptsG under metabolites only, is not a kill criterion. It is
recorded here because K2's response ("report the forward direction as the
headline") assumes the forward direction is measured.
- **It is not measured on the exact sampler at M0.** The transcript parameters did
  not mix under metabolites only (ESS 6 and 7, split-R̂ 1.13 and 1.16), and the
  user chose to report them as unmixed (§12 2026-10-28).
- **Until it is measured,** the forward direction rests on track A's linearised
  estimate (ptsG at 0.47 of prior SD; `approximate-method.md` A1) and on the
  approximate route.
- **A redesigned move and a rerun** (about 550 CPU-h) is kept as an option, to
  revisit when the paper's framing is set.
