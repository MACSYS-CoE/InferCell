# KIM-90: diagnosing the M0 gap

Paper §4 (text and Limits); Spec risk "The M0 gap's cause may bear on the reduced model."

On m0_v6 the mixed posterior sits away from V6, the long exact reference, by
0.468 SD on S_0779, 0.564 on krnadeg and 0.419 on ENO
(`av6_89_b7_result.md`). Re-estimating S did not close the gap. This note
tests the candidate causes one at a time.

**Verdict.** No cause is named. The checks narrow it down without closing it.
- **Ruled out:**
  - target identity;
  - MCMC error on either side.
- **The bulk reference (δ and S): partial.**
  - Switching the Jensen correction off makes krnadeg's and ENO's gaps about four
    times larger, so the correction is needed. This does not test whether δ's value
    is accurate.
  - Varying the reference (more cells, or δ and S recomputed at m0_v6's posterior
    mean) closes about 0.17 SD of S_0779 and up to 0.16 of krnadeg, and leaves most
    of the gap.
- **Shared cells: indeterminate.** The cross-covariance check predicted +0.37 of
  krnadeg's +0.56. That prediction omitted first-order terms that are not small. With
  them, krnadeg's first-order shift is +0.12 with the bulk covariance held and −0.14
  with it conditioned on the transcripts, and S_0779's changes sign (+0.15, −0.24).
  The one-step response these rest on is not valid here either (§4b), so the check
  names nothing.
- **A new lead: the mixed route counts the transcripts' information twice.** The bulk
  is the mean of the same 50 cells whose transcripts are observed. Given the
  transcripts, the bulk carries, to first order, 0.19–0.35 of the information the
  mixed route gives it on S_0779, 0.32 on krnadeg and 0.36–0.66 on S_0607 (§4b).
  This bears on width: the mixed posterior should be too narrow there. It is
  narrower on krnadeg (0.85) and S_0779 (0.93), and on S_0607 barely (0.98, about the
  SD's sampling error). Its effect on where the posterior sits is not measured.
- **Against shared cells, and double counting, for krnadeg:** under metabolites
  only, with no transcripts, the mixed fit already sits 0.032 above the exact on log
  krnadeg (AV6 B7.8), about 0.4 V6 SD. Neither mechanism can act there. That exact
  fit has an ESS of 39, and whether this is the same gap was not tested.
- **Untested:**
  - the mean field on this dataset beyond δ's first-order correction;
  - the transcript factor's population-mean constants;
  - a direct test of shared cells (a disjoint-cell exact run, 508–678 CPU-h). This
    was not run (Tom's call, 2026-11-23).

All distances are in V6 posterior SDs: shift = mean(mixed) − mean(V6). The
archive is `/fred/oz022/tkimpson/infercell-kim90-archive-2026-11-23`:
- `run/` holds the refits and the first cross-covariance run;
- `reads/` holds items 1 and 2;
- `review/` holds the round run after the PR #189 review;
- `SHA256SUMS` covers all of it.

## What was run

| Check | Script | Slurm jobs | Commit |
|---|---|---|---|
| 1. Target identity | `gap_diag_identity.jl` | 21397711 (first run 21386523 at 289a1d7) | 64f9da0 |
| 2. MCMC error | `gap_diag_mixing.jl` | 21397711 (first run 21386523) | 64f9da0 |
| 3. Ensemble size (n1600) | `gap_diag_refit.jl` | prep 21098790, chains 21098791/92 | 66431db |
| 4a. Jensen correction off (nojensen) | `gap_diag_refit.jl` | prep 21098793, chains 21098794/95 | 66431db |
| control (archived δ and S) | `gap_diag_refit.jl` | prep 21098787, chains 21098788/89 | 66431db |
| 4b. Shared cells, first run (`analyse`) | `gap_diag_crosscov.jl` | ensemble 21098797, analysis 21098798; rerun 21397710 at 64f9da0 | 66431db |
| 4b. Shared cells, full first-order terms (`full`) | `gap_diag_crosscov.jl` | 21403381 | 89c0134 |
| 4c. δ at m0_v6's posterior mean (postmean) | `gap_diag_refit.jl` | prep 21186223, chains 21186224/25 | 66431db |
| merge of all refit arms | `gap_diag_refit.jl merge` | 21186226 | 66431db |

The refit script has changed since 66431db only in its comments.
- **The 64f9da0 runs:** items 1 and 2 ran again there because the identity script's
  truncation checks changed. The `analyse` rerun at 64f9da0 reproduces the 66431db
  table exactly.
- **Submission:** `gap_diag_refit.sh` documents how each round was submitted. The
  review round was `sbatch` calls from worktrees at the commits, into
  `/fred/oz022/tkimpson/infercell-kim90-review-2026-11-23`:
  - at 64f9da0: `full` (21397709, `full_64f9da0.md`) and `analyse` (21397710) through
    `gap_diag_refit.slurm`, and `gap_diag_reads.slurm <out>/reads` (21397711);
  - at 0a2a79b: `full` again (21397829, `full_0a2a79b.md`), adding the θ_ref line;
  - at a02eb6e: `full` again (21398305, `full_a02eb6e.md`), adding the conditioned
    bulk covariance to the displacement and the Hessian (the fix gate's finding);
  - at 89c0134: `full` again (21403381, `full.md`), taking both displacement columns
    from one formula. 21398305's "sum, W held" column, summed term by term, dropped
    GᵀWΔm(x̄) and is withdrawn. 21403381 is the one cited; its other outputs match
    the earlier runs.
- **Wording:** the mixing script and its table say "Monte Carlo error" for sampler
  error; this note says MCMC error. "Monte Carlo error" here means the reference
  ensemble's.

The design of every refit is AV6's: 2 NUTS chains of 200 + 200 at target acceptance 0.8.
- **Chains:** v6's chain seeds and jittered starts.
- **Input:** B7's archived prep (`infercell-89-b7-archive-2026-10-27/run/v6`).
- **What varies:** arms differ only in the bulk reference's δ and S.
- **The 4c arm** was added after 4b's first run came out indeterminate for S_0779
  and ENO (Decision log 2026-11-19).

## 1. Target identity: ruled out

Job 21397711 prints "RESULT: target identity holds."
- **Fingerprints across commits:** these hash the same at V6's commit (281d548), at
  B7's (7a102fc) and at 64f9da0:
  - the dataset (replicate 16085, 50 cells, every field);
  - the priors of the five targets and σ_b;
  - the panel and the horizon.
- **B7's archived bulk:** matches the regenerated one exactly (max |Δ log| = 0).
- **Both setups built from the dataset have the same:**
  - free parameters;
  - priors on the log (the CME and σ_b priors are checked; ENO's and FBA's are read
    from the model in both routes, so they agree by construction);
  - bulk record;
  - per-cell transcripts;
  - floor (1.0);
  - window (60 s, 10 saves);
  - panel.
- **Truncations:**
  - No draw of either posterior is near the mixed route's μmin on krnadeg. V6's
    lowest draw is 5.7 prior SDs above it and m0_v6's lowest is 6.0.
  - m0_v6's 400 stored draws have no divergences. A −Inf density on any leapfrog
    step (past kmax, μmin or the box) would show as one.
  - Every m0_v6 draw is inside the ±6 prior-SD box; the closest is 3.88 prior SDs
    from an edge.
  - The exact route's block-2 bound is not checked.

## 2. MCMC error: ruled out

| Parameter | shift | V6 split-half W1 | V6 MCSE | mixed chain 1 vs 2 W1 | mixed MCSE | z |
|---|---|---|---|---|---|---|
| S_0779 | −0.468 | 0.079 | 0.017 | 0.142 | 0.046 | −9.5 |
| krnadeg | +0.564 | 0.087 | 0.030 | 0.104 | 0.040 | +11.2 |
| ENO | +0.419 | 0.071 | 0.026 | 0.068 | 0.054 | +7.1 |

V6 is 4 chains × sweeps 101–1650 (jobs 18276511–14, 18420986–89). m0_v6 is 2 ×
200 draws (jobs 19177810/11). z is the shift over the two sides' MCSEs added in
quadrature. The gap is a shift in where the posterior sits, beyond both samplers'
error. The full table is in `review/reads/mixing.md` (job 21397711; the same as
`reads/mixing.md` from 21386523).

## 3–4. The refits and the cross-covariance check

| Arm | S_0779 | krnadeg | ENO | max R̂ |
|---|---|---|---|---|
| m0_v6 (the gap) | −0.468 | +0.564 | +0.419 | 1.010 |
| control | −0.468 | +0.564 | +0.419 | 1.010 |
| n1600 | −0.304 | +0.523 | +0.314 | 1.008 |
| nojensen | −0.002 | +2.187 | +1.546 | 1.008 |
| postmean | −0.292 | +0.408 | +0.396 | 1.010 |

Every entry is the arm's shift against V6. All six parameters, with W1 and the
shift against control, are in `run/refit.md` (job 21186226).

**control reproduces m0_v6 exactly.** It uses the same seeds, starts and inputs, so
it checks that this branch's code reproduces B7. It is not a noise floor.
- **The noise floor for an arm against control:** a changed target should make the
  same-seed chains diverge, so each arm is in effect an independent pair of chains.
  The difference of two such posterior means then has an SD of about √2 × the mixed
  MCSE: 0.065 on S_0779, 0.057 on krnadeg and 0.076 on ENO.
- **This floor is assumed, not measured.** Nothing here tests that the same-seed
  arms decorrelate from control. AV6's matched null measured two identical fits on
  different chain seeds (job 19334380) at W1 0.061 / 0.090 / 0.098 in that note's
  units, so it is a distance, not a difference of means.

**3. Ensemble size: a part.** Doubling the reference to 1600 cells moves δ by 0.189
(norm; prep 21098790) and leaves tr S at 1.715 against 1.718.
- **S_0779:** the gap narrows by 0.164, 2.5× the floor and 0.35 of the gap.
- **ENO:** narrows by 0.105, 1.4× the floor and 0.25 of the gap.
- **krnadeg:** narrows by 0.041, within noise.
- **σ_b** moves by +0.153 against control (`run/refit.md`). A changed δ changes the
  bulk residual, so σ_b can move for a real reason, but no floor is computed for
  it here.
- **Not an upper bound.** The 1600 cells contain the archived 800, so the doubling
  removes only part of the 800-cell reference's Monte Carlo error. More of the gap
  may come from it than these moves show.

**4a. Jensen correction off: makes the gap worse.** With δ = 0
(S held), krnadeg's and ENO's gaps grow about fourfold:
- krnadeg sits at +2.19 and ENO at +1.55;
- S_0607 and σ_b sit at +1.29 and +0.94 against V6, where they were close before.

Only S_0779 lands on V6 (−0.002), so for that one parameter the correction may
overshoot. Switching the correction off is wrong. This does not test whether δ's
value is accurate: δ corrects the mean field at θ_ref, to first order. 4c and §3
both move the gap by moving δ.

**4b. Shared cells: indeterminate.**

*The first run (`analyse`, job 21098798, reproduced by 21397710).*
- **Regression:** across the 800 reference cells (re-simulated on the archived
  seeds; the panel matches the archive's exactly), each cell's panel is regressed on
  its transcripts (20 features: 2 genes × 10 saves).
- **Δm:** this gives the shift in the bulk's expected log that the observed
  transcripts imply, which the mixed route ignores. ‖Δm‖ = 0.514, or 0.135 of the
  bulk noise per row.
- **Response:** carried through V6's posterior covariance as Σ Jᵀ W Δm, it predicts
  −0.143 on S_0779, +0.369 on krnadeg and +0.085 on ENO. By the reading rule fixed
  before the run, that named shared cells for krnadeg.

*Why that reading does not stand (job 21403381 at 89c0134, PR #189 review and its fix gate).*
- **What `analyse` omitted.** Two first-order terms:
  - Δm depends on θ through the transcripts' expected value, with
    G = ∂Δm/∂log θ ≈ −B ∂t̄/∂log θ / x̄. This enters through the bulk residual r,
    whose RMS in units of the bulk noise is 1.01, against 0.135 for Δm.
  - Given the transcripts, the exact route's bulk covariance is smaller, σ²I + S_c,
    where S_c is the panel's covariance given the transcripts. tr S_c / tr S = 0.62,
    and S dominates the noise (tr S / tr σ²I = 1.91).
- **The full displacement.** To first order, mixed − exact is V6's posterior
  covariance times the difference of the two routes' bulk gradients at m0_v6's
  posterior mean x̄: Σ (JᵀW r − (J+G)ᵀW* (r − Δm(x̄))). Δm(x̄) = Δm + G (x̄ − x_ref) is
  Δm carried from θ_ref to x̄. W* is the mixed route's W = (σ²I + S)⁻¹ (W held), or
  W_c = (σ²I + S_c)⁻¹ (S conditioned). Both columns come from this one formula.

  | Parameter | gap | Σ JᵀWΔm, Δm at θ_ref (`analyse`) | full, W held | full, S conditioned |
  |---|---|---|---|---|
  | S_0607 | −0.129 | −0.006 | −0.152 | −0.756 |
  | S_0779 | −0.468 | −0.143 | +0.147 | −0.243 |
  | krnadeg | +0.564 | +0.369 | +0.117 | −0.142 |
  | ENO | +0.419 | +0.085 | +0.093 | +0.181 |

  ∂t̄ is taken from the mean field, whose transcripts differ from the ensemble mean
  by up to 13%; G carries that error. An earlier version of the W-held column
  (job 21398305) was summed term by term, dropped GᵀWΔm(x̄) and gave +0.635 on
  krnadeg; it is withdrawn.
- **Unstable to what is kept.** krnadeg's first-order shift is +0.369 from `analyse`,
  +0.117 in full with W held and −0.142 with S conditioned. S_0779's changes sign
  between the two full columns (+0.147, −0.243), and S_0607's moves from −0.15 to
  −0.76.
- **The response is not valid here.** H is the exact route's bulk precision minus the
  mixed route's. Its largest-magnitude eigenvalue, relative to V6's posterior
  precision, is
  −0.734 with W held and −0.702 with S conditioned. So in one direction the change
  is most of the posterior's own curvature, and a one-step linear response cannot
  size the shift.
- **No size is given.** None of +0.369, +0.117 or −0.142 is a reliable size for
  krnadeg, and no interval is given for the full terms.

*What the check does measure: the bulk's information.* Each entry is the diagonal of
the bulk's Fisher information as a share of V6's posterior precision (job 21403381):

| Parameter | ‖J+G‖/‖J‖ | mixed route | given the transcripts, W held | given the transcripts, S conditioned |
|---|---|---|---|---|
| S_0607 | 0.321 | 0.261 | 0.093 | 0.171 |
| S_0779 | 0.185 | 0.325 | 0.062 | 0.115 |
| krnadeg | 0.139 | 0.758 | 0.240 | 0.242 |
| ENO | 1.001 | 0.853 | 0.853 | 0.863 |
| FBA | 1.000 | 0.989 | 0.989 | 0.989 |

The two right-hand columns are first order. The last also replaces S by the
panel's covariance given the transcripts.
- **What is double-counted.** On S_0607, S_0779 and krnadeg, much of the bulk's
  sensitivity to θ comes through the transcripts. Once they are known, the bulk adds
  0.19–0.35 of what the mixed route credits it with on S_0779, 0.32 on krnadeg and
  0.36–0.66 on S_0607 (the two columns' ratios to the mixed share).
- **What this implies.** The mixed route should be too narrow on these three. That is
  a statement about width, and it is computed, not measured. m0_v6's posterior SDs
  against V6's are narrower on log krnadeg (0.0645 against 0.0763, 0.85) and log
  S_0779 (0.1864 against 0.1995, 0.93). On log S_0607, whose information is
  double-counted too (0.36–0.66), they are barely narrower: 0.1137 against 0.1159
  (0.98), about the SD's own sampling error, roughly 1/√(2 ESS) ≈ 3% at
  ESS 547 (`av6_89_b7_result.md`; V6 SDs also in job
  21403381). The direction agrees on all three; the size on S_0607 is not borne out.
- **What is not measured.** Whether this also moves where the posterior sits, and by
  how much, would need a refit with the bulk conditioned on the transcripts.
- **ENO and FBA** are learned from the bulk independently of the transcripts, and are
  unaffected.

*Against shared cells for krnadeg (AV6 B7.8, job 19411150 at db30305).*
- **The result.** Under metabolites only, the mixed fit sits above the exact met_only
  fit (job 19034903) by 0.032 on log krnadeg, 0.30 exact SD.
- **Why it cuts against shared cells.** That fit has no transcripts, so shared cells
  cannot act there.
- **In V6 units.** It is about 0.4 V6 SD (V6's SD on log krnadeg is 0.0763; job
  21403381), against the joint gap of 0.564. So the bulk factor alone may carry much
  of krnadeg's offset.
- **Caveat.** The exact met_only fit has an ESS of 39, and whether this is the same
  gap was not tested there.

**4c. δ at m0_v6's posterior mean: a part, like 3.** θ_ref sits 1.27 V6 SDs from
m0_v6's posterior mean in S_0607 and 1.17 in krnadeg (job 21403381). In σ_b it sits
11.9 SDs away, but δ and S come from simulated cells, which σ_b does not enter.
Recomputing δ and S on the same 800 seeds at m0_v6's posterior mean moves δ by 0.344
and tr S to 1.704 (prep 21186223); nojensen's δ = 0 is 0.813 away.
- **S_0779:** narrows by 0.176 (2.7× the floor).
- **krnadeg:** narrows by 0.156 (2.7×).
- **ENO:** narrows by 0.023 (within noise).

S_0779 moves about as far as under n1600, whose δ moved about half as far (0.189).
So δ's dependence on θ and its sampling noise cannot be separated here. Either way,
the reference's placement does not close the gap.

## What this leaves

- **krnadeg:** two leads, neither measured as a shift.
  - **Double counting:** the bulk's information on krnadeg is about three times what
    it is given the transcripts. It bears on width; whether it also moves the mean is
    open.
  - **The bulk factor itself:** B7.8 shows about 0.4 V6 SD with no transcripts, where
    neither shared cells nor double counting can act.

  Changing the reference closes up to 0.16 of the 0.56.
- **S_0779:** about 0.3 SD remains after §3 or §4c. Its information is double-counted
  too, which is a width statement and does not account for the remainder.
- **ENO:** about 0.3 SD remains after §3, and §4c does not move it. ENO's
  information is not double-counted, and the first-order cross-covariance
  displacement there is +0.093 (W held) or +0.181 (S conditioned), with the same
  validity caveat.
  It is the least explained of the three.
- **Candidates not tested here:**
  - **a mixed refit with the bulk conditioned on the transcripts:** this would test
    the double counting directly, but it needs a new likelihood in `src/`;
  - **the mean field on this dataset:** m(θ) = log x̄_MF against the ensemble mean,
    which δ corrects at θ_ref only, to first order;
  - **the transcript factor's population-mean constants;**
  - **a disjoint-cell exact run**, which would test shared cells directly, at
    508–678 CPU-h. Not run: Tom's call, 2026-11-23.

**For the reduced model.** The reduced model uses the same mixed route:
- δ and S at θ_ref;
- transcripts and bulk from the same cells, there 200 cells and all 17 genes.

So the double counting measured here can arise there too, on any parameter the bulk
learns mainly through the transcripts. Its size on the reduced model is not
measured. The pilot's promoters sit about 0.8 SD low (`corea_pilot_kim92_result.md`).
This note does not show that the two are linked.
