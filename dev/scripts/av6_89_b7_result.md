# Issue #89 B7: AV6, the mixed posterior against the exact sampler's on M0

Paper section 5.2. Rules fixed before any run in §12 2026-10-24, "B7 as planned";
the MAP's stopping rule changed during the run (§12 2026-10-25 and 2026-10-26).

**Verdicts:**
- **AV6's pass line, SBC over V7's 50 datasets: passes**, 6 of 6 parameters within
  the simultaneous bands at 1 − 0.05/6. Also 6 of 6 on the 33 replicates that meet
  the convergence target.
- **KS fires** (rule 6): on 5 of 18 marginals, each by at most 0.027 SD. A matched
  null run after the result (information only; the verdict stands) puts two of the
  five at about twice the noise between two identical fits (v6 S_0607, r12 ENO),
  one at 1.3× it (v6 FBA) and two within it (v6 σ_b, r12 S_0607).
  - On the two datasets with a null (v6 and r12; r19 had none), holding S fixed
    shows as a W1 of about 0.2 SD on one marginal each. That raw W1 includes
    about 0.1 SD of chain noise.
  - KS's floors are not matched to its 400-against-400 comparison. Against the
    matched null, the floor is the higher of the two on 5 of 6 marginals on v6
    (F10's largest null; σ_b is the exception) and on 4 of 6 on r12 (the larger
    200-against-200 chain W1; S_0607 and FBA are the exceptions). So the floors
    mostly make KS harder to fire. The exceptions include two of the five
    firings, v6 σ_b and r12 S_0607: the two that are within the null.
  - The response §8 names, the linear noise approximation, is not built in B7
    (the user's decision, 2026-10-27).
- **KC does not fire** (rule 7): median 6.57 CPU-h per posterior against 50.75.
- **On m0_v6** (rule 4, information): five of six marginals are resolved from V6
  under the reading fixed in advance (W1 above all four of F10's matched nulls).
  - Three by about 0.3 SD over the largest null: S_0779, krnadeg and ENO, at
    0.47, 0.56 and 0.42 SD.
  - Two narrowly: FBA by 0.005 and σ_b by 0.013.
  - S_0607 is not resolved.
  - Re-estimating S at the posterior mean moves the three large ones by
    0.08–0.09 SD and does not close the gap; S's form was not varied. Open
    question, §9.
- **Against V7 over the 50** (rule 5, information): the largest mean bias is
  krnadeg's, +0.095 ± 0.050 SD. Per parameter, 9–22 of the 50 replicates differ
  from V7 by more than 3× V7's own MCSE. The median SD ratio is 0.92–1.07.
- **Stream arms** (rule 8): the mixed side is below. The comparison with 16b.4's
  exact fits (B7.8, §12 2026-10-28 "B7.8 as planned") is in "B7.8" below:
  - under transcripts only, no marginal is resolved from the exact fit;
  - under metabolites only, ENO, σ_b and krnadeg are resolved, by small margins.
    FBA is not resolved. ptsG has no exact counterpart.

## Runs

- **Code:** the merge ran at `8da6045` (job 19292232). Stages ran at three
  commits, identical except for `reference_theta`'s stopping rule:
  - `7a102fc`: the launch, jobs 19177653–19177825 (39 replicates and v6, ks_v6,
    the stream arms);
  - `6526b45`: the second resume, jobs 19238017–19238056 (ten base fits);
  - `8da6045`: the third resume, jobs 19292223–19292232 (r33, ks_low1, ks_low2,
    the merge).
- **θ_ref's rule per fit:** Optim's tests alone for 40 base fits; the relative-
  change stop (`reltol = 1e-10`) for r5, r7, r18, r23, r25, r28, r29, r37, r42,
  r50; the x-step stop (`xtol = 1e-6`) for r33. The KS refits reuse their base's
  θ_ref (v6, r12, r19: Optim's tests alone). The rules
  agree to at most 2.1e-5 prior SD wherever compared (jobs 19232720–19232726,
  19290049–19290055).
- **Matched null:** at `09be2f5`, chains 19334375, 19334376, 19334378, 19334379
  and comparison 19334380 (`av6_ks_null.md`); smoke 19334381 on copied chains gave
  W1 = 0 and reproduced the merge's KS numbers. The launcher does not submit
  these; `dev/scripts/av6_89_b7.jl`'s header has the commands, run against the
  run directory after the merge.
- **Run size:** 50 cells, 600 s, an 800-cell reference ensemble, 2 chains of
  200 + 200 per fit from starts jittered by 0.25 prior SD.
- **Archive:** `/fred/oz022/tkimpson/infercell-89-b7-archive-2026-10-27/`, with
  `SHA256SUMS` (529 files):
  - the run directory, the MAP traces, the MAP checks and the null's smoke;
  - `profile/`: B7.2's snapshots and the profile logs (19159993, 19160196,
    19161734, 19162281, 19162737, 19163602);
  - `tests/`: the test-file logs (19160188, 19163839, 19163842, 19164153,
    19232727, 19232740, 19290056) and the suites' (19342349; 19344907 after
    #104's review);
  - `replay/`: #104's review's replay of the merge and nullcmp (19344905,
    19344906);
  - B7.5's two smoke directories, r1 and r2.
- **Figures:** `dev/notes/figures/av6-89b7/` (`dev/scripts/av6_89_b7.py <run>`):
  Figure 4 (`f4_mixed_vs_exact`) and Figure 5's middle panel (`f5_rank_ecdf_m0`).
- **Cost:** 391.9 CPU-h of record across the three batches' completed jobs, and
  274.1 CPU-h lost to the MAP fault (`dev/scripts/phase16_ledger.md`).

The tables below are the merge's own (`av6_summary.md`, job 19292232), copied
verbatim, then the null's (job 19334380).

## Convergence (rule 2: reported, not used to drop)

Worst R̂ 1.0463; smallest bulk ESS 81; median smallest ESS 365.
Below R̂ ≤ 1.02 and ESS ≥ 100: r3, r4, r8, r10, r13, r17, r19, r20, r22, r25, r26, r27, r32, r34, r36, r44, r50 (17 of 50).

## SBC (rule 3, AV6's pass line)

Rank-ECDF of 50 replicates, L = 100, against simultaneous bands at 1 − 0.05/6 (pointwise γ = 0.00076), V7's construction.

| Parameter (log) | Outside the band at | Largest |ECDF − z| | Verdict |
|---|---|---|---|
| log_S_JCVISYN3A_0607 | nowhere | 0.144 | within |
| log_S_JCVISYN3A_0779 | nowhere | 0.084 | within |
| log_krnadeg | nowhere | 0.183 | within |
| log_kcatF_R_ENO | nowhere | 0.136 | within |
| log_kcatF_R_FBA | nowhere | 0.113 | within |
| log_sigma_b | nowhere | 0.135 | within |

**AV6's SBC passes**: 6 of 6 within.
Without the 17 below the convergence target (sensitivity, not the verdict): 6 of 6 within, on 33.

## Against V7 over the 50 (rule 5, information)

d_r = (mean_mixed − mean_V7) / sd_V7 per replicate; V7's own MCSE is 1/√ESS_r.

| Parameter (log) | mean d ± SE | median |d| | median V7 MCSE | replicates with |d| > 3 MCSE | median SD ratio |
|---|---|---|---|---|---|
| log_S_JCVISYN3A_0607 | +0.019 ± 0.030 | 0.109 | 0.089 | 9 | 0.99 |
| log_S_JCVISYN3A_0779 | -0.064 ± 0.031 | 0.137 | 0.053 | 22 | 0.92 |
| log_krnadeg | +0.095 ± 0.050 | 0.153 | 0.083 | 18 | 0.99 |
| log_kcatF_R_ENO | +0.064 ± 0.031 | 0.087 | 0.052 | 17 | 1.07 |
| log_kcatF_R_FBA | -0.029 ± 0.018 | 0.071 | 0.045 | 12 | 0.93 |
| log_sigma_b | +0.032 ± 0.018 | 0.064 | 0.046 | 13 | 1.01 |

## m0_v6 against V6 (rule 4, information; F10's matched nulls)

| Parameter (log) | mixed mean ± sd | V6 mean ± sd | W1 / SD_V6 | matched nulls | reading |
|---|---|---|---|---|---|
| log_S_JCVISYN3A_0607 | 1.5305 ± 0.1137 | 1.5454 ± 0.1159 | 0.137 | 0.084, 0.117, 0.120, 0.186 | not resolved |
| log_S_JCVISYN3A_0779 | 0.6067 ± 0.1864 | 0.6999 ± 0.1995 | 0.468 | 0.044, 0.040, 0.052, 0.115 | resolved: above all four, by 0.353 |
| log_krnadeg | 2.2676 ± 0.0645 | 2.2246 ± 0.0763 | 0.564 | 0.110, 0.033, 0.161, 0.210 | resolved: above all four, by 0.355 |
| log_kcatF_R_ENO | 4.2599 ± 0.0454 | 4.2423 ± 0.0419 | 0.419 | 0.082, 0.048, 0.094, 0.118 | resolved: above all four, by 0.301 |
| log_kcatF_R_FBA | 4.3906 ± 0.1870 | 4.4077 ± 0.1844 | 0.105 | 0.043, 0.100, 0.056, 0.078 | resolved: above all four, by 0.005 |
| log_sigma_b | -2.6217 ± 0.0557 | -2.6250 ± 0.0601 | 0.091 | 0.042, 0.038, 0.079, 0.076 | resolved: above all four, by 0.013 |

## KS (rule 6): S re-estimated at the posterior mean, δ held

| Fit | base | Parameter (log) | W1(fixed, re-estimated) | floor | over the floor |
|---|---|---|---|---|---|
| ks_v6 | v6 | log_S_JCVISYN3A_0607 | 0.213 | 0.186 (V6 SD; floor = F10's largest matched null) | yes |
| ks_v6 | v6 | log_S_JCVISYN3A_0779 | 0.086 | 0.115 (V6 SD; floor = F10's largest matched null) | no |
| ks_v6 | v6 | log_krnadeg | 0.076 | 0.210 (V6 SD; floor = F10's largest matched null) | no |
| ks_v6 | v6 | log_kcatF_R_ENO | 0.093 | 0.118 (V6 SD; floor = F10's largest matched null) | no |
| ks_v6 | v6 | log_kcatF_R_FBA | 0.111 | 0.100 (V6 SD; floor = F10's largest matched null) | yes |
| ks_v6 | v6 | log_sigma_b | 0.084 | 0.079 (V6 SD; floor = F10's largest matched null) | yes |
| ks_low1 | r12 | log_S_JCVISYN3A_0607 | 0.115 | 0.111 (base SD; floor = larger chain-vs-chain W1) | yes |
| ks_low1 | r12 | log_S_JCVISYN3A_0779 | 0.063 | 0.154 (base SD; floor = larger chain-vs-chain W1) | no |
| ks_low1 | r12 | log_krnadeg | 0.061 | 0.104 (base SD; floor = larger chain-vs-chain W1) | no |
| ks_low1 | r12 | log_kcatF_R_ENO | 0.216 | 0.197 (base SD; floor = larger chain-vs-chain W1) | yes |
| ks_low1 | r12 | log_kcatF_R_FBA | 0.104 | 0.158 (base SD; floor = larger chain-vs-chain W1) | no |
| ks_low1 | r12 | log_sigma_b | 0.107 | 0.139 (base SD; floor = larger chain-vs-chain W1) | no |
| ks_low2 | r19 | log_S_JCVISYN3A_0607 | 0.101 | 0.195 (base SD; floor = larger chain-vs-chain W1) | no |
| ks_low2 | r19 | log_S_JCVISYN3A_0779 | 0.093 | 0.178 (base SD; floor = larger chain-vs-chain W1) | no |
| ks_low2 | r19 | log_krnadeg | 0.080 | 0.199 (base SD; floor = larger chain-vs-chain W1) | no |
| ks_low2 | r19 | log_kcatF_R_ENO | 0.095 | 0.146 (base SD; floor = larger chain-vs-chain W1) | no |
| ks_low2 | r19 | log_kcatF_R_FBA | 0.049 | 0.173 (base SD; floor = larger chain-vs-chain W1) | no |
| ks_low2 | r19 | log_sigma_b | 0.070 | 0.194 (base SD; floor = larger chain-vs-chain W1) | no |

**KS fires**.

### KS's matched null (added after the result; information only)

Two fixed-S fits per base with the same θ_ref, δ and S and different chain seeds
(the base arm's chain seeds 3 and 4). W1(base, null) is the noise between two
identical fits, in KS's units. Job 19334380 at `09be2f5`.

| Base | Parameter (log) | W1(base, null) | W1(fixed, re-estimated) | KS floor | KS over the null |
|---|---|---|---|---|---|
| v6 | log_S_JCVISYN3A_0607 | 0.113 | 0.213 | 0.186 | yes |
| v6 | log_S_JCVISYN3A_0779 | 0.061 | 0.086 | 0.115 | yes |
| v6 | log_krnadeg | 0.090 | 0.076 | 0.210 | no |
| v6 | log_kcatF_R_ENO | 0.098 | 0.093 | 0.118 | no |
| v6 | log_kcatF_R_FBA | 0.086 | 0.111 | 0.100 | yes |
| v6 | log_sigma_b | 0.128 | 0.084 | 0.079 | no |
| r12 | log_S_JCVISYN3A_0607 | 0.147 | 0.115 | 0.111 | no |
| r12 | log_S_JCVISYN3A_0779 | 0.100 | 0.063 | 0.154 | no |
| r12 | log_krnadeg | 0.085 | 0.061 | 0.104 | no |
| r12 | log_kcatF_R_ENO | 0.101 | 0.216 | 0.197 | yes |
| r12 | log_kcatF_R_FBA | 0.186 | 0.104 | 0.158 | no |
| r12 | log_sigma_b | 0.056 | 0.107 | 0.139 | yes |

Reading of the five firings:
- v6 S_0607 (0.213 against a null of 0.113) and r12 ENO (0.216 against 0.101) are
  about twice the null.
- v6 FBA (0.111 against 0.086) is 1.3× the null.
- v6 σ_b (0.084 against 0.128) and r12 S_0607 (0.115 against 0.147) are below it.
- Each null is one pair of 400-draw fits, so it is itself noisy; "about twice" is
  the precision it supports.

## KC (rule 7)

CPU-h per posterior (MAP, ensemble, chains; source: sacct): median 6.57, min 5.04, max 13.08, against 50.75. **KC does not fire**.

Per replicate, from `av6_replicates.tsv`: 5.04 (r7) to 13.08 CPU-h (r38). The
three dearest, r38, r9 and r4 (13.08, 12.90, 12.03), carry MAPs of 6.41, 7.66
and 3.50 h under Optim's tests alone. Under the x-step stop the same MAPs take
0.06–0.07 h. KC counts only completed jobs: the 274.1 CPU-h lost to timed-out and
cancelled preps are not in any replicate's cost.

## Stream arms on m0_v6 (rule 8; the mixed side)

| Fit | Parameter (log) | mean ± sd | shrinkage 1 − var_post / var_prior | R̂ | ESS |
|---|---|---|---|---|---|
| v6 | log_S_JCVISYN3A_0607 | 1.5305 ± 0.1137 | 0.973 | 1.0101 | 547 |
| v6 | log_S_JCVISYN3A_0779 | 0.6067 ± 0.1864 | 0.928 | 1.0084 | 411 |
| v6 | log_krnadeg | 2.2676 ± 0.0645 | 0.991 | 1.0029 | 438 |
| v6 | log_kcatF_R_ENO | 4.2599 ± 0.0454 | 0.916 | 1.0027 | 406 |
| v6 | log_kcatF_R_FBA | 4.3906 ± 0.1870 | 0.888 | 0.9959 | 668 |
| v6 | log_sigma_b | -2.6217 ± 0.0557 | 0.997 | 0.9984 | 414 |
| tx_v6 | log_S_JCVISYN3A_0607 | 1.5119 ± 0.1166 | 0.972 | 0.9988 | 635 |
| tx_v6 | log_S_JCVISYN3A_0779 | 0.6779 ± 0.1997 | 0.917 | 0.9979 | 658 |
| tx_v6 | log_krnadeg | 2.2083 ± 0.0740 | 0.989 | 0.9982 | 615 |
| tx_v6 | log_kcatF_R_ENO | 4.1221 ± 0.1443 | 0.153 | 0.9989 | 572 |
| tx_v6 | log_kcatF_R_FBA | 4.1031 ± 0.5444 | 0.047 | 1.0087 | 353 |
| tx_v6 | log_sigma_b | -1.5924 ± 1.0205 | -0.041 | 1.0016 | 477 |
| bulk_v6 | log_S_JCVISYN3A_0607 | 1.3291 ± 0.3780 | 0.703 | 1.0006 | 288 |
| bulk_v6 | log_S_JCVISYN3A_0779 | 0.9707 ± 0.4683 | 0.544 | 0.9974 | 295 |
| bulk_v6 | log_krnadeg | 2.3749 ± 0.1084 | 0.976 | 0.9995 | 312 |
| bulk_v6 | log_kcatF_R_ENO | 4.2565 ± 0.0454 | 0.916 | 1.0176 | 420 |
| bulk_v6 | log_kcatF_R_FBA | 4.4119 ± 0.1942 | 0.879 | 1.0096 | 455 |
| bulk_v6 | log_sigma_b | -2.6126 ± 0.0597 | 0.996 | 1.0034 | 421 |

The exact side, 16b.4's fits, had not finished on 2026-10-27. B7.8 compares them
below.

## B7.8: the stream arms against 16b.4's exact fits

The rules were fixed before any distance was computed (§12 2026-10-28, "B7.8 as
planned").
- **W1** is in the exact posterior's SD.
- **The floor** is the larger of the exact chain's split-half W1 and the mixed
  fit's chain-vs-chain W1.
- **The exact fits** are 16b.4's single chains: tx_only (job 19034902) and met_only
  (job 19034903), sweeps 101 to 656.

**Reproduce** (Slurm job 19411150 at `db30305`, 0.05 CPU-h; identical to a
login-node run):
```
julia --project dev/scripts/av6_89_b78.jl /fred/oz022/tkimpson/infercell-89-b7-archive-2026-10-27/run \
    /fred/oz022/tkimpson/infercell-16b-archive-2026-10-28/m0_16b_tx_only \
    /fred/oz022/tkimpson/infercell-16b-archive-2026-10-28/m0_16b_met_only
```

**Transcripts only** (`tx_v6` against the exact tx_only fit):

| Parameter | Exact mean ± SD | Mixed mean ± SD | W1 / exact SD | Exact split-half | Mixed chain-vs-chain | Floor | Exact ESS | Exact split-R̂ | Reads |
|---|---|---|---|---|---|---|---|---|---|
| S_JCVISYN3A_0607 | 1.502 ± 0.117 | 1.512 ± 0.117 | 0.096 | 0.103 | 0.112 | 0.112 | 142 | 1.000 | not resolved |
| S_JCVISYN3A_0779 | 0.684 ± 0.200 | 0.678 ± 0.200 | 0.055 | 0.093 | 0.064 | 0.093 | 305 | 1.001 | not resolved |
| krnadeg | 2.202 ± 0.073 | 2.208 ± 0.074 | 0.097 | 0.146 | 0.126 | 0.146 | 160 | 1.009 | not resolved |
| kcatF_R_ENO | 4.131 ± 0.147 | 4.122 ± 0.144 | 0.066 | 0.067 | 0.126 | 0.126 | 405 | 0.999 | not resolved |
| kcatF_R_FBA | 4.121 ± 0.562 | 4.103 ± 0.544 | 0.051 | 0.145 | 0.144 | 0.145 | 496 | 1.003 | not resolved |
| sigma_b | held at its start | -1.592 ± 1.021 | — | — | — | — | — | — | not compared |

**Metabolites only** (`bulk_v6` against the exact met_only fit):

| Parameter | Exact mean ± SD | Mixed mean ± SD | W1 / exact SD | Exact split-half | Mixed chain-vs-chain | Floor | Exact ESS | Exact split-R̂ | Reads |
|---|---|---|---|---|---|---|---|---|---|
| S_JCVISYN3A_0607 | unmixed | 1.329 ± 0.378 | — | — | — | — | 7 | 1.129 | not read |
| S_JCVISYN3A_0779 | unmixed | 0.971 ± 0.468 | — | — | — | — | 6 | 1.164 | not read |
| krnadeg | 2.343 ± 0.107 | 2.375 ± 0.108 | 0.302 | 0.147 | 0.141 | 0.147 | 39 | 1.003 | resolved, by 0.155, exact ESS below 100 |
| kcatF_R_ENO | 4.251 ± 0.042 | 4.256 ± 0.045 | 0.153 | 0.098 | 0.116 | 0.116 | 140 | 1.000 | resolved, by 0.037 |
| kcatF_R_FBA | 4.417 ± 0.194 | 4.412 ± 0.194 | 0.069 | 0.097 | 0.227 | 0.227 | 568 | 1.000 | not resolved |
| sigma_b | -2.622 ± 0.053 | -2.613 ± 0.060 | 0.179 | 0.112 | 0.149 | 0.149 | 385 | 1.001 | resolved, by 0.030 |

**What it shows.**
- **Under transcripts only, the mixed route is not resolved from the exact
  sampler on any of the five marginals.** W1 runs from 0.051 to 0.097 exact SD,
  each under its floor (0.093 to 0.146).
  - This includes §3's block-2 check (B1): the transcript factor's S_g against the
    exact sampler. S_0607 is at 0.096 against 0.112, and ptsG at 0.055 against
    0.093.
  - ENO's posterior SD is 0.144 in the mixed fit and 0.147 in the exact. Neither
    fit learns ENO from transcripts, which is the measured null.
- **Under metabolites only, three of the four read marginals are resolved, by
  small margins.**
  - ENO is at 0.153 against a floor of 0.116, a margin of 0.037.
  - σ_b is at 0.179 against 0.149, a margin of 0.030.
  - krnadeg is at 0.302 against 0.147, a margin of 0.155. But its exact fit has
    ESS 39, and a split-half of so few effective draws can understate the exact
    side's noise.
  - FBA is not resolved: 0.069 against 0.227.
  - In each resolved cell the mixed mean sits above the exact one: krnadeg by
    0.032, ENO by 0.005 and σ_b by 0.009 on the log scale. That is 0.30, 0.12 and
    0.17 exact SD.
  - The joint fit on m0_v6 also sits apart from V6 on krnadeg and ENO (rule 4,
    "the m0_v6 gap", §9). Whether this is the same gap is not tested here.
- **The floor is one draw.** Each half and each mixed chain is a single sample of
  the noise. Like F10, a margin of a few hundredths of an SD is not large against
  that.
- **ptsG under metabolites only has no exact counterpart.** The exact fit did not
  mix (ESS 6), so it is not read.
  - The mixed fit gives SD / prior SD 0.676, which is 1 − var_post / var_prior of
    0.544.
  - Track A's linearised estimate (A1) was 0.47, against a split-half noise
    bracket of 0.68 (`sensitivity_88_m0_result.md`). A1 treats save points as
    independent, which overstates the information.
  - So the approximate route shows ptsG shrinking under metabolites only, by less
    than A1 projected. **Figure 4's "ptsG still shrinks under metabolites only"
    has the approximate route's support and no exact comparison at M0** (T3, "the
    forward direction").

## Open

- **The m0_v6 gap.** S_0779, krnadeg and ENO sit 0.42–0.56 SD from V6, about
  0.3 SD above every matched null, while over the 50 the largest mean bias
  against V7 is 0.095 ± 0.050 SD and SBC passes. Re-estimating S at the posterior
  mean does not close it. Untested candidates: δ's sampling noise
  (one 800-cell ensemble), the mean field itself at this dataset, and the
  transcript factor's population-mean constants (§9).
- **Convergence:** 17 of 50 replicates fall short of R̂ ≤ 1.02 and ESS ≥ 100 at
  2 × (200 + 200); the worst R̂ is 1.046. SBC is unchanged without them.
