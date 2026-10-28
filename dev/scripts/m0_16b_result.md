# 16b.1 and 16b.4 at M0: the perturbed truth recovers; F6 is drawn, with ptsG under metabolites only unmixed

The rules were fixed before any production draw existed (`spec/phases/16-recovery.md`
§12 2026-10-22, approved 2026-10-23). How 16b.4 is held is the user's decision of
2026-10-28, made after the results (§12 2026-10-28).

**The runs.** Three production chains of the three-block Gibbs sampler, side by
side, at F10's settings: 50 cells, 600 s, PG at N = 20, 656 sweeps, 5 threads,
and 100 sweeps of burn-in. Each started from the prior with `Xoshiro(16_400_000 + i)`.

| Mode | Task | Job | At | Wall | CPU-h |
|---|---|---|---|---|---|
| `perturbed` | 16b.1 | 19034901 | `03b4c4f` | 4 d 10.6 h | 532.8 |
| `tx_only` | 16b.4 | 19034902 | `03b4c4f` | 4 d 21.3 h | 586.5 |
| `met_only` | 16b.4 | 19034903 | `03b4c4f` | 4 d 14.1 h | 550.5 |

- All three were COMPLETED with exit code 0, from the detached worktree
  `.worktrees/16b-03b4c4f`. The CPU-h are `sacct`'s CPUTimeRAW.
- The checkpoints, logs and `sacct` record are archived, with checksums, at
  `/fred/oz022/tkimpson/infercell-16b-archive-2026-10-28/`.
- **The reference** is V6's chains 1 to 4 (`v6_16a9b_result.md`), sweeps 101 to
  1,650 of each. It is the joint fit.

**Reproduce** (Slurm job 19378382 at `626a1cc`; its log is in the archive):
```
RUN_MODE=perturbed PILOT_DIR=<any dir> julia --project dev/scripts/m0_16b.jl summary \
    /fred/oz022/tkimpson/infercell-16a9b-archive-2026-10-17/v6/gibbs_pilot_16a9a \
    <archive>/m0_16b_perturbed <archive>/m0_16b_tx_only <archive>/m0_16b_met_only
RUN_MODE=met_only PILOT_DIR=<archive>/m0_16b_met_only julia --project dev/scripts/m0_16b.jl rates
```
The `rates` line is from job 19376253 at `03b4c4f`, run once on met_only after the
chains ended. The check §12 set at about sweep 20 on each chain was not run.

## 16b.1: recovers, so K4 does not fire

ENO's truth moved by +0.1568 (one prior SD), and the posterior mean moved by
+0.1491 (MC SE 0.0033).

| Rule | Bound | Measured | Holds |
|---|---|---|---|
| Shrinkage (SD / prior SD) | below 0.9 | 0.258 | yes |
| The moved truth (ln ENO = 4.3386) in the central 95% | quantile in 0.025 to 0.975 | 0.099 | yes |
| Mean shift ÷ truth shift | 0.5 to 1.5 | 0.951 ± 0.021 | yes |

- The ratio is near the conjugate reading's 0.93 at V6's shrinkage.
- V6's own truth sits at quantile 0.000 of the perturbed posterior, so the
  posterior left it.
- S_0607's bulk ESS is 97 in this chain, just below 100. The rule reads ENO only,
  whose ESS is 164.

## 16b.4: transcripts only is measured; metabolites only is held open

| Fit | Parameter (log) | Mean ± SD | SD / prior SD | Bulk ESS | Split-R̂ | Truth | Truth's quantile |
|---|---|---|---|---|---|---|---|
| joint (V6, 4 chains) | S_JCVISYN3A_0607 | 1.545 ± 0.116 | 0.167 | 1081 | 1.007 | 1.489 | 0.31 |
| joint (V6, 4 chains) | S_JCVISYN3A_0779 | 0.700 ± 0.199 | 0.288 | 3441 | 1.001 | 1.012 | 0.95 |
| joint (V6, 4 chains) | krnadeg | 2.225 ± 0.076 | 0.110 | 1088 | 1.006 | 2.243 | 0.60 |
| joint (V6, 4 chains) | kcatF_R_ENO | 4.242 ± 0.042 | 0.267 | 1531 | 1.002 | 4.182 | 0.07 |
| joint (V6, 4 chains) | kcatF_R_FBA | 4.408 ± 0.184 | 0.331 | 5750 | 1.001 | 4.207 | 0.13 |
| joint (V6, 4 chains) | sigma_b | -2.625 ± 0.060 | 0.060 | 4248 | 1.000 | -2.622 | 0.52 |
| perturbed | S_JCVISYN3A_0607 | 1.522 ± 0.127 | 0.183 | 97 | 1.002 | 1.489 | 0.38 |
| perturbed | S_JCVISYN3A_0779 | 0.711 ± 0.210 | 0.303 | 252 | 0.998 | 1.012 | 0.94 |
| perturbed | krnadeg | 2.215 ± 0.077 | 0.111 | 135 | 1.001 | 2.243 | 0.64 |
| perturbed | kcatF_R_ENO | 4.391 ± 0.040 | 0.258 | 164 | 1.000 | 4.339 | 0.10 |
| perturbed | kcatF_R_FBA | 4.403 ± 0.202 | 0.362 | 485 | 0.999 | 4.207 | 0.16 |
| perturbed | sigma_b | -2.621 ± 0.059 | 0.059 | 288 | 1.009 | -2.622 | 0.46 |
| tx_only | S_JCVISYN3A_0607 | 1.502 ± 0.117 | 0.169 | 142 | 1.000 | 1.489 | 0.43 |
| tx_only | S_JCVISYN3A_0779 | 0.684 ± 0.200 | 0.288 | 305 | 1.001 | 1.012 | 0.96 |
| tx_only | krnadeg | 2.202 ± 0.073 | 0.105 | 160 | 1.009 | 2.243 | 0.72 |
| tx_only | kcatF_R_ENO | 4.131 ± 0.147 | 0.934 | 405 | 0.999 | 4.182 | 0.63 |
| tx_only | kcatF_R_FBA | 4.121 ± 0.562 | 1.007 | 496 | 1.003 | 4.207 | 0.58 |
| tx_only | sigma_b | held at its start, not inferred | — | — | — | — | — |
| met_only | S_JCVISYN3A_0607 | 1.480 ± 0.314 | 0.453 | 7 | 1.129 | 1.489 | 0.46 |
| met_only | S_JCVISYN3A_0779 | 0.980 ± 0.376 | 0.542 | 6 | 1.164 | 1.012 | 0.51 |
| met_only | krnadeg | 2.343 ± 0.107 | 0.154 | 39 | 1.003 | 2.243 | 0.18 |
| met_only | kcatF_R_ENO | 4.251 ± 0.042 | 0.267 | 140 | 1.000 | 4.182 | 0.06 |
| met_only | kcatF_R_FBA | 4.417 ± 0.194 | 0.347 | 568 | 1.000 | 4.207 | 0.13 |
| met_only | sigma_b | -2.622 ± 0.053 | 0.053 | 385 | 1.001 | -2.622 | 0.50 |

**Transcripts only.**
- **ENO's measured null is not flagged, by 0.00014.** SD / prior SD is 0.93445 at
  bulk ESS 404.7, so SE = 0.03285 and the flag line 1 − 2·SE is 0.93431. The point
  estimate is 2 SE below 1. The bound it is read against
  (`eno_path_information_16a7_result.md`) was measured at one truth, on Core A′ at
  200 cells, not on M0.
  - **The verdict is not robust either way.** The margin is 0.4% of one SE, far
    smaller than the uncertainty in the ESS estimate that the SE rests on. A small
    change in that estimate, or in the burn-in, would move the line across the
    ratio.
- FBA's SD / prior SD is 1.007.
- The three transcript parameters shrink as far as in the joint fit (0.169, 0.288
  and 0.105, against V6's 0.167, 0.288 and 0.110). The two agree within Monte
  Carlo error (about 0.01 at these ESS), so at M0 the bulk data add no shrinkage
  to them that this run can resolve.
- For 16b.5, M0's prior-to-posterior divergence for ENO and FBA is to be read off
  this chain. It is not computed here.

**Metabolites only.**
- **The metabolic parameters are measured.** ENO (ESS 140), FBA (568) and σ_b
  (385) have split-R̂ ≤ 1.001. Their posteriors match the joint fit's to within
  0.01 in the mean and 0.02 in the SD ratio.
- **The transcript parameters are not.** S_0607 and ptsG (S_JCVISYN3A_0779) have
  bulk ESS 7 and 6, and split-R̂ 1.129 and 1.164. krnadeg has ESS 39.
  - **ptsG's shrinkage under metabolites only, the load-bearing cell, is therefore
    not measured.** The table's 0.542 is one wandering chain's spread, and is not
    read against A1's 0.47 or against anything else.
- **The symptom named in advance did not occur; the cause is not tested.** §12
  named met_only's bootstrap proposal as the mixing risk, checked by a live
  cell-window updating at a rate below 0.10. Over 414 live cell-windows the lowest
  rate is 0.37 and the median 0.54, so none is below 0.10. That shows the paths are
  resampled, not that they decorrelate from S_g: a window counts as changed if any
  event in it differs, metabolic ones included. One reading is that each S_g is
  updated given the path's transcription counts and the path given S_g, so the two
  move together slowly when no transcript data separate them. That coupling comes
  from updating S_g and the path in separate blocks. Neither that reading nor any part
  the proposal plays is tested here.
- **What reaching ESS 100 would cost.** ptsG's ESS of 6 comes from 556 kept draws,
  so ESS 100 needs about 9,270 kept draws, or about 8,710 more sweeps. At 0.839
  CPU-h a sweep (550.47 over 656), that is about 7,300 CPU-h more. It assumes ESS
  grows in proportion to the draws, which an ESS of 6 at split-R̂ 1.16 does not
  establish. About 4.6k remain under phase 16's cap.

**Held open, then decided (the user's decisions, 2026-10-28).** 16b.1 closed
first, with 16b.4 held open. The user then chose to report the metabolites-only
transcript parameters as unmixed (§12 "16b.4 at M0: metabolites only reported as
unmixed"). Redesigning the sampler's move under metabolites only and rerunning it
is kept as an option to revisit.

**It followed that** F6 waited on that choice, which is now made (below). B7.8's
ptsG comparison has no exact counterpart. B7.8 also needs a floor:
each 16b.4 arm is one chain, so "16b.4's chain-vs-chain W1" does not exist as
written. A split-half floor at ESS 6 would make the ptsG comparison vacuous.

## F6 at M0

The rules were fixed before any divergence was computed (§12 2026-10-28 "F6 at M0:
the divergence, its floor and how each cell reads"). The shrinkages are the ones
in job 19378382.

**Reproduce** (Slurm job 19395555 at `f0dd5b2`, 0.15 CPU-h; its log is in the archive):
```
RUN_MODE=tx_only PILOT_DIR=<archive>/m0_16b_tx_only julia --project dev/scripts/m0_16b.jl f6 \
    /fred/oz022/tkimpson/infercell-16a9b-archive-2026-10-17/v6/gibbs_pilot_16a9a \
    <archive>/m0_16b_tx_only <archive>/m0_16b_met_only
```
The divergence is KL(posterior ‖ prior) in nats, per marginal, on the log scale.
"Floor" is the 95th percentile of 20 prior samples, each the size of the cell's
bulk ESS, scored the same way.

| Parameter | Data | SD / prior SD | KL (nat) | KL at 8,000 | Gaussian KL | Floor | Bulk ESS | Split-R̂ | Reads |
| S_JCVISYN3A_0607 | transcripts only | 0.169 | 1.5515 | 1.5515 | 1.5683 | 0.0322 | 142 | 1.000 | resolved |
| S_JCVISYN3A_0607 | metabolites only | — | — | — | — | — | 7 | 1.129 | unmixed, not read |
| S_JCVISYN3A_0607 | joint (V6) | 0.167 | 1.5329 | 1.5329 | 1.5355 | 0.0097 | 1081 | 1.007 | resolved |
| S_JCVISYN3A_0779 | transcripts only | 0.288 | 1.4957 | 1.4957 | 1.5305 | 0.0192 | 305 | 1.001 | resolved |
| S_JCVISYN3A_0779 | metabolites only | — | — | — | — | — | 6 | 1.164 | unmixed, not read |
| S_JCVISYN3A_0779 | joint (V6) | 0.288 | 1.4867 | 1.4867 | 1.5035 | 0.0033 | 3441 | 1.001 | resolved |
| krnadeg | transcripts only | 0.105 | 2.6590 | 2.6590 | 2.6944 | 0.0521 | 160 | 1.009 | resolved |
| krnadeg | metabolites only | 0.154 | 2.5850 | 2.5850 | 2.6176 | 0.1768 | 39 | 1.003 | resolved, ESS below 100 |
| krnadeg | joint (V6) | 0.110 | 2.6700 | 2.6700 | 2.6933 | 0.0072 | 1088 | 1.006 | resolved |
| kcatF_R_ENO | transcripts only | 0.934 | 0.0041 | 0.0041 | 0.0044 | 0.0194 | 405 | 0.999 | not resolved from the prior |
| kcatF_R_ENO | metabolites only | 0.267 | 1.1273 | 1.1273 | 1.1513 | 0.0280 | 140 | 1.000 | resolved |
| kcatF_R_ENO | joint (V6) | 0.267 | 1.0991 | 1.0991 | 1.1109 | 0.0048 | 1531 | 1.002 | resolved |
| kcatF_R_FBA | transcripts only | 1.007 | 0.0081 | 0.0081 | 0.0017 | 0.0138 | 496 | 1.003 | not resolved from the prior, KDE and Gaussian differ by > 20% |
| kcatF_R_FBA | metabolites only | 0.347 | 0.7903 | 0.7903 | 0.7910 | 0.0140 | 568 | 1.000 | resolved |
| kcatF_R_FBA | joint (V6) | 0.331 | 0.8399 | 0.8399 | 0.8244 | 0.0023 | 5750 | 1.001 | resolved |
| sigma_b | transcripts only | held at its start | — | — | — | — | — | — | held |
| sigma_b | metabolites only | 0.053 | 2.9196 | 2.9196 | 2.9483 | 0.0161 | 385 | 1.001 | resolved |
| sigma_b | joint (V6) | 0.060 | 2.8183 | 2.8183 | 2.8297 | 0.0027 | 4248 | 1.000 | resolved |

**What it shows.**
- **Transcripts only:**
  - the transcript parameters are resolved (1.50 to 2.66 nat);
  - ENO and FBA are not resolved from the prior: 0.0041 and 0.0081 nat, under
    floors of 0.0194 and 0.0138.
  - ENO's cell is the measured null, read in 16b.4 above.
  - FBA's KDE and Gaussian values differ (0.0081 and 0.0017). Both are below the
    floor, so the reading does not depend on which one is used.
- **Metabolites only:**
  - ENO, FBA and σ_b are resolved: 1.13, 0.79 and 2.92 nat, against the joint
    fit's 1.10, 0.84 and 2.82;
  - krnadeg is resolved at 2.59 nat, with ESS 39;
  - S_0607 and ptsG are unmixed and not read. **The load-bearing cell is therefore
    not measured on the exact sampler at M0.**
- **Joint (V6):** every parameter is resolved.
- **The grid check holds in every cell.** The value at 8,000 points matches 4,000
  to four decimals. The KDE and Gaussian values agree within 20% in every read
  cell except FBA under transcripts only.

**16b.5's M0 read-off.** Under transcripts only, ENO's and FBA's divergences
(0.0041 and 0.0081 nat) are below 0.05 nat. That is K2's firing condition, met at
M0 with 50 cells. It is M0's read-off, not K2's verdict, which is defined at
1,000 cells and is budget-bound. 16b.5 is not ticked.

## 16b.2 and 16b.3 at M0: read off V6 and V7

Neither needs new compute (§12 "16a.11", call 1). Both read results that already
exist, against the clauses in §11.

### 16b.2: the full parameter set recovers (V6)

M0 has five of D11's six targets; `S_PGI` is not in M0. The numbers are V6's joint
rows from job 19378382 (also in `v6_16a9b_result.md`).

| Parameter | SD / prior SD | Truth's quantile |
|---|---|---|
| S_JCVISYN3A_0607 | 0.167 | 0.31 |
| S_JCVISYN3A_0779 (ptsG) | 0.288 | 0.95 |
| krnadeg | 0.110 | 0.60 |
| kcatF_R_ENO | 0.267 | 0.07 |
| kcatF_R_FBA | 0.331 | 0.13 |
| sigma_b | 0.060 | 0.52 |

- **Every marginal contains its truth.** All six quantiles lie inside the central
  95% (0.025 to 0.975). ptsG's 0.95 is the closest to an edge.
- **At least one target shrinks below 0.5 (K4):** all six do, from 0.060 to 0.331.
- **Cost per posterior, against 16a.11's measured 507.5 to 678 CPU-h:** the three
  16b production posteriors cost 532.8, 586.5 and 550.5 CPU-h (`sacct`), inside
  that range. V6 itself is a reference posterior (4 × 1,650 sweeps), at about
  6,809.7 CPU-h (`v6_16a9b_result.md`). It is not a production posterior.

### 16b.3: coverage over repeated datasets (V7)

Read off V7's ranks by `dev/scripts/v7_coverage_16a11.py` (Slurm job 19399322 at
`e84a44e`). R = 50 replicates and L = 100 thinned draws. Bold marks a count
outside the binomial 95% band.

| Level | S_JCVISYN3A_0607 | S_JCVISYN3A_0779 | krnadeg | kcatF_R_ENO | kcatF_R_FBA | sigma_b | Band |
| 50% | 23 | 27 | 22 | 31 | 24 | 26 | 18–32 |
| 80% | 36 | 41 | 37 | 42 | 40 | 42 | 34–45 |
| 90% | 42 | 47 | 45 | 45 | 42 | 48 | 41–49 |
| 95% | 46 | 47 | 46 | 46 | 44 | 49 | 44–50 |
| 99% | 49 | 48 | **47** | 50 | 49 | 49 | 48–50 |

- **29 of 30 cells are in band.** One or two misses are expected by chance over 30.
  The miss is krnadeg at 99%, where L = 100 is coarsest.
- **The tight-prior control (ENO) at nominal 95%** covers 46 of 50. That is inside
  K4's band restated at R = 50 (44 to 50; §12 "16a.11", call 1), so K4's first
  trigger does not fire.
- **K4's third trigger (equal widths for ENO and FBA) does not fire.** In V6, ENO's
  posterior SD is 0.042 and FBA's 0.184, which is 0.267 and 0.331 of their priors.
  The same holds under 16b.1's perturbation: 0.040 and 0.202.
- **Carried into T3, as call 1 required:**
  - 23 of the 50 replicates have bulk ESS below 100;
  - L = 100 makes the 99% level coarse;
  - it covers M0's five parameters plus σ_b, not D11's six;
  - K4 asks for 100 replications, and V7 has 50, so the band is restated at R = 50.
