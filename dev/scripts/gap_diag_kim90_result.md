# KIM-90: diagnosing the M0 gap

Paper §4 (text and Limits); Spec risk "The M0 gap's cause may bear on the reduced model."

On m0_v6 the mixed posterior sits away from V6, the long exact reference, by
0.468 SD on S_0779, 0.564 on krnadeg and 0.419 on ENO
(`av6_89_b7_result.md`). Re-estimating S did not close the gap. This note
tests the candidate causes one at a time.

**Verdict.** No single cause is named.
- **Shared cells** account for most of krnadeg: predicted +0.37 of its +0.56.
- **Ruled out:**
  - target identity;
  - Monte Carlo error on either side;
  - the Jensen correction.
- **Partial:** varying the bulk reference (more cells, or δ and S
  recomputed at the posterior mean) moves S_0779 by about 0.17 SD and krnadeg
  by up to 0.16.
- **Unexplained:** about 0.3 SD of S_0779 and 0.4 SD of ENO.
- **Untested:**
  - the mean field on this dataset;
  - the transcript factor's population-mean constants;
  - a direct test of shared cells (a disjoint-cell exact run, 508–678 CPU-h).
  This last one was not run (Tom's call, 2026-11-23).

All distances are in V6 posterior SDs: shift = mean(mixed) − mean(V6). The
archive is `/fred/oz022/tkimpson/infercell-kim90-archive-2026-11-23` (`run/`
the refits and cross-covariance check, `reads/` items 1 and 2, SHA256SUMS).

## What was run

| Check | Script | Slurm jobs | Commit |
|---|---|---|---|
| 1. Target identity | `gap_diag_identity.jl` | 21386523 | 289a1d7 |
| 2. Monte Carlo error | `gap_diag_mixing.jl` | 21386523 | 289a1d7 |
| 3. Ensemble size (n1600) | `gap_diag_refit.jl` | prep 21098790, chains 21098791/92 | 66431db |
| 4a. Jensen correction off (nojensen) | `gap_diag_refit.jl` | prep 21098793, chains 21098794/95 | 66431db |
| control (archived δ and S) | `gap_diag_refit.jl` | prep 21098787, chains 21098788/89 | 66431db |
| 4b. Shared cells (cross-covariance) | `gap_diag_crosscov.jl` | ensemble 21098797, analysis 21098798 | 66431db |
| 4c. δ at the posterior mean (postmean) | `gap_diag_refit.jl` | prep 21186223, chains 21186224/25 | 66431db |
| merge of all refit arms | `gap_diag_refit.jl merge` | 21186226 | 66431db |

The design of every refit is AV6's: 2 NUTS chains of 200 + 200 at target acceptance 0.8.
- **Chains:** v6's chain seeds and jittered starts.
- **Input:** B7's archived prep (`infercell-89-b7-archive-2026-10-27/run/v6`).
- **What varies:** arms differ only in the bulk reference's δ and S.
- **The 4c arm** was added after 4b came out indeterminate (Decision log 2026-11-19).

## 1. Target identity: ruled out

Job 21386523 prints "RESULT: target identity holds."
- **Fingerprints across commits:** the dataset (replicate 16085, 50 cells, every field), the priors of the
  five targets and σ_b, the panel and the horizon hash the same at V6's commit
  (281d548), B7's (7a102fc) and 289a1d7.
- **B7's archived bulk:** matches the regenerated one exactly in this job (max |Δ log| = 0;
  on the login node it differs by 1.8e-10, roundoff between machines).
- **Both setups built from the dataset have the same:**
  - free parameters;
  - priors on the log;
  - bulk record;
  - per-cell transcripts;
  - floor (1.0);
  - window (60 s, 10 saves);
  - panel.
- **Truncations:**
  - No draw of either posterior is near the mixed route's μmin on krnadeg. V6's lowest draw is 5.7 prior SDs above it and m0_v6's lowest is 6.0.
  - All 400 m0_v6 draws have a finite log density, so the kmax truncation never bound.

## 2. Monte Carlo error: ruled out

| Parameter | shift | V6 split-half W1 | V6 MCSE | mixed chain 1 vs 2 W1 | mixed MCSE | z |
|---|---|---|---|---|---|---|
| S_0779 | −0.468 | 0.079 | 0.017 | 0.142 | 0.046 | −9.5 |
| krnadeg | +0.564 | 0.087 | 0.030 | 0.104 | 0.040 | +11.2 |
| ENO | +0.419 | 0.071 | 0.026 | 0.068 | 0.054 | +7.1 |

V6 is 4 chains × sweeps 101–1650 (jobs 18276511–14, 18420986–89). m0_v6 is 2 ×
200 draws (jobs 19177810/11). z is the shift over the two sides' MCSEs added in
quadrature. The gap is a shift in where the posterior sits, beyond both sides'
Monte Carlo error. The full table is in `reads/mixing.md`.

## 3–4. The refits and the cross-covariance check

| Arm | S_0779 | krnadeg | ENO | max R̂ |
|---|---|---|---|---|
| m0_v6 (the gap) | −0.468 | +0.564 | +0.419 | 1.010 |
| control | −0.468 | +0.564 | +0.419 | 1.010 |
| n1600 | −0.304 | +0.523 | +0.314 | 1.008 |
| nojensen | −0.002 | +2.187 | +1.546 | 1.008 |
| postmean | −0.292 | +0.408 | +0.396 | 1.010 |
| cross-covariance prediction | −0.143 [−0.327, +0.046] | +0.369 [+0.158, +0.579] | +0.085 [−0.044, +0.253] | — |

The prediction's brackets are 90% bootstrap intervals. All six parameters,
with W1 and the shift against control, are in `run/refit.md` and
`run/crosscov/crosscov.md`.

**control reproduces m0_v6 exactly.** It uses the same seeds, starts and inputs, so it is a check
that this branch's code reproduces B7, not a noise floor.
- **The noise floor for an arm against control:** a changed target makes the
  same-seed chains diverge, so each arm is in effect an independent pair of
  chains.
- **Its size:** the difference of two such posterior means has an SD of
  about √2 × the mixed MCSE, so 0.065 on S_0779, 0.057 on krnadeg and 0.076
  on ENO.

**3. Ensemble size: a small part.** Doubling the reference to 1600 cells moves δ by 0.189 (norm; prep 21098790) and leaves tr S at 1.715 against 1.718.
- **S_0779:** the gap narrows by 0.164, 2.5× the noise floor.
- **ENO:** narrows by 0.105, 1.4× the floor.
- **krnadeg:** narrows by 0.041, within noise.
- So the 800-cell reference's Monte Carlo error accounts for at most a third of S_0779 and a quarter of ENO.

**4a. Jensen correction: ruled out.** With δ = 0 (S held), the gap grows
three- to fourfold:
- krnadeg +2.19 and ENO +1.55;
- S_0607 and σ_b move by 1.29 and 0.94 where they were close before.

Only S_0779 lands on V6 (−0.002), so for that one parameter the correction
may overshoot. It is needed, and it is not the cause.

**4b. Shared cells: named for krnadeg, indeterminate for S_0779 and ENO.** The method:
- **Regression:** across the 800 reference cells (re-simulated on the archived seeds;
  the panel matches the archive's exactly), each cell's panel is regressed on
  its transcripts (20 features: 2 genes × 10 saves).
- **Δm:** this gives the shift in the bulk's expected log that the observed
  transcripts imply, which the mixed route ignores. ‖Δm‖ = 0.514, or 0.135 of
  the bulk noise per row.
- **Displacement:** carried through V6's posterior covariance by linear response, Δm predicts
  a displacement of the mixed posterior.
- **Reading rule:** fixed in the script's header before the run.
- **Result:**
  - krnadeg's predicted +0.369 has an interval that contains the observed +0.564 and excludes 0: named.
  - For S_0779 and ENO the predictions have the gap's sign but are a third
    and a fifth of it, and their intervals include 0: indeterminate.
- **Caveat:** the linear response is first order and uses V6's covariance for the mixed
  route's information. It sizes the effect but is not a refit.

**4c. δ at the posterior mean: a small part, like 3.** θ_ref sits about 1.3 SD from the posterior in krnadeg
and S_0607. Recomputing δ and S on the same 800 seeds at the posterior mean
moves δ by 0.344 and tr S to 1.704 (prep 21186223); nojensen's δ = 0 is 0.813 away.
- **S_0779:** narrows by 0.176 (2.7× the floor).
- **krnadeg:** narrows by 0.156 (2.7×).
- **ENO:** narrows by 0.023 (within noise).

S_0779 moves about as far as under n1600, whose δ moved about half as far
(0.189). So δ's dependence on θ and its sampling noise cannot be separated
here. Either way, the reference's
placement does not close the gap.

## What this leaves

- **krnadeg:** mostly the shared cells. The mixed route treats the transcript and
  bulk streams as independent, but in the data the bulk is the mean of the same 50
  cells whose transcripts are observed.
- **S_0779 and ENO:** about 0.3 and 0.4 SD remain after every check.
  Candidates not tested here:
  - **the mean field on this dataset:** m(θ) = log x̄_MF against the ensemble mean.
    δ corrects this at θ_ref only, to first order.
  - **the transcript factor's population-mean constants.**
  - **shared cells beyond first order.** A disjoint-cell exact run would test this directly, at
    508–678 CPU-h. Not run: Tom's call, 2026-11-23.

**For the reduced model.** The reduced model uses the same mixed route:
- δ and S at θ_ref;
- transcripts and bulk from the same cells, there 200 cells and all 17 genes.

So the shared-cell mechanism that moves krnadeg on M0 applies there too.
Its size on the reduced model is not measured. The pilot's promoters sit
about 0.8 SD low (`corea_pilot_kim92_result.md`). This note does not show
that the two are linked.
