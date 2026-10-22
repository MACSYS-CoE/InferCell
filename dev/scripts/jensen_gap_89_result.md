# Issue #89 B3: the mean field fails KG at every truth

**Run:** ensembles 18886224–18886226 (3 × 10 tasks), mean fields 18886227–18886231,
merge 18886232, all at `9686014`, submitted from a detached worktree at that pushed
commit. Smoke 18885264–18885272 at the same commit passed first.
**Rounding check:** mean fields 18888681–18888685 and merge 18888686 at `2d15e83`,
with `MF_ROUNDING=fractional_carry`, reusing the M0 ensembles above.
**Diagnosis:** the `diagnose` stage at `c8c9ceb`, job 18892098 over the run and
18892099 over the rounding check. It reads the archived outputs and runs nothing new.
**Script:** `dev/scripts/jensen_gap_89.jl`. **Gate:** `spec/phases/approximate-method.md`
§8 KG (AV5), agreed 2026-10-20 and amended for its noise before this run (§12).
**Cost:** 2.20 CPU-h for the run, 0.86 for the smoke, 0.46 for the rounding check.
**Raw outputs:** `/fred/oz022/tkimpson/infercell-89-b3-archive-2026-10-20/` and
`/fred/oz022/tkimpson/infercell-89-b3-rounding-archive-2026-10-20/` (mean fields,
`merged.md`, `merged.jls`, `diagnosis.md`, Slurm logs; ensembles in the first).
**Reference for the reduced model:** #88 A3's `live60` ensembles (200 cells, jobs
18825531 and 18825539), made by `dev/scripts/eno_feedback_88.jl` at `8813dca` on
PR #92's branch, not yet on `main`. Archive:
`/fred/oz022/tkimpson/infercell-88-a3-archive-2026-10-19/`. `src` at `8813dca`
differs from `main` only by `override_rebuild_pools!`, unset in `live60`. The merge
checks that each reference was drawn at the same truth as its mean field. Merge
#92 and this PR with merge commits, so `8813dca` and `9686014` stay reachable.
**Regenerate:** `dev/scripts/jensen_gap_89.sh <outdir outside the worktrees>`.

## Setup

- **The mean-field cell.** This is the full composition with transcription, decay and
  translation silenced: the states, parameters, edges and rebuilds are kept, and the
  reactions are removed. B2's `MeanField` mean is advanced over each 1 s jump step at
  the driver's current rate constants. So the drains (clip included), the catalytic
  edges, the volume and the 60 s rebuild run as published. The rebuild reads this
  cell's own pools, which are the mean-field pools under θ. The ODE pools are written
  back continuously.
- **Writing the mean into the integer jump state:**
  - counters with fractional carry, exact in the drained total to one particle;
  - proteins rounded to the nearest count. At 100 copies or more that is within 0.27%.
    The worst case is cytosolic ptsG at about one copy, at 0.5, but the driver never
    reads it; the mean field carries it.
- **The reference.** The 200-cell mean of the full two-way model at the same truth,
  under fractional carry, with cells seeded 40000 + k.
- **Truths:**
  - M0 (600 s, 10 saves): V6's replicate 16085, and V7's replicates 13 and 2 (ptsG
    2.75, 1.13 and 33.5).
  - The reduced model (6,300 s, 105 saves): A3's truths 1 and 2 (ptsG 4.36 and 0.80).

## Verdict: the mean field fails KG at all five truths. Stop before B4

| Truth | Total (nats, σ_b = 10%) | Noise | Total ≤ noise | Pools above 3× their noise |
|---|---|---|---|---|
| M0, V6's 16085 | 70.8 | 29.2 | no | 8 |
| M0, V7's 13 (low ptsG) | 12.4 | 14.9 | yes | 7 |
| M0, V7's 2 (high ptsG) | 25.8 | 17.4 | no | 9 |
| Reduced model, truth 1 | 99.9 | 115 | yes | 5 |
| Reduced model, truth 2 (low ptsG) | 202 | 55.7 | no | 9 |

- **The per-pool clause fails everywhere.** The total fails at three of five truths.
- **The worst gap is never large in σ_b units.** The largest single-save gap is
  3.6 σ_b (GTP, M0 16085). Most pools stay under 1.5 σ_b. But GTP's gap is systematic
  across the whole trajectory, so it adds up to many nats.
- **It is not where the sub-spec expected.** RB1 expected the gap at low ptsG. M0's
  worst total is at ptsG 2.75, and its low-ptsG truth has the smallest.
- **Per the sub-spec, this stops B4.** The mean field, as built, is not yet the second
  description, and σ_b is not to be tuned to make it pass.
- **What AV5 measures.** The mean-field cell differs from the population in two
  declared reductions at once: the Jensen gap, and the mean-field feedback (the
  rebuild reads mean pools, not each cell's). This run does not separate them.

## The rounding check: no effect

The mean-field cell writes its ODE pools back continuously, and the reference cells
use fractional carry. Rerun under the reference's policy, every number barely moves,
and each truth flags the same pools:

| Truth | Total, continuous | Total, fractional carry | GTP shift, continuous → fractional carry |
|---|---|---|---|
| M0, 16085 | 70.8 | 70.7 | 22.8 → 22.8 |
| M0, 13 | 12.4 | 12.3 | 7.46 → 7.39 |
| M0, 2 | 25.8 | 25.9 | 6.98 → 7.02 |
| Reduced model, truth 1 | 99.9 | 100 | 31.0 → 31.2 |
| Reduced model, truth 2 | 202 | 203 | 116 → 117 |

So the rounding policy does not explain the gap.

## What the gap is (diagnosis on the same outputs)

Mean of log(x̄_MF / x̄_SSA) over the saves, floored at one particle as KG is, with
its range (from `diagnosis.md`, which also gives every panel pool's CV and its
fraction of cell-saves below one particle):

| Pool | M0 16085 | M0 13 | M0 2 | Reduced 1 | Reduced 2 |
|---|---|---|---|---|---|
| GTP | −0.17 (−0.36 to 0) | −0.11 | −0.11 | −0.075 | −0.14 (−0.24 to −0.02) |
| ATP | +0.05 | +0.02 | +0.03 | +0.02 | −0.00 |
| GDP | +0.13 | +0.04 | +0.06 | +0.03 | +0.02 |
| GMP | +0.10 | +0.02 | +0.03 | +0.00 | +0.01 |
| AMP | −0.10 | −0.03 | −0.05 | −0.03 | −0.02 |
| G3P | −0.01 | −0.01 | −0.06 | −0.03 | −0.03 |
| G6P | −0.01 | −0.00 | −0.01 | −0.01 | +0.02 |

- **It is systematic and has one sign.** The mean field's GTP sits about 7–16% below
  the population mean at every truth (exp(−0.075) to exp(−0.172)). GDP sits above it
  at every truth, and GMP too, though barely on the reduced model. That is a guanylate imbalance, not noise.
- **It is not the one-particle floor.** No cell-save is below one GTP particle except
  2% at reduced truth 2. No mean-field drain clipped.
  - **The SSA side's clipping was not captured.** The clip fires when a 1 s drain's
    accrual exceeds the pool, not when a save is below one particle. Neither this
    script's ensembles nor A3's saved `clipping_census` per cell. It costs about
    7 CPU-h to add.
- **The adenylates point at the response, not the demand.**
  - If the mean field over-drained GTP, ATP would be drained harder too and would sit
    below the SSA mean. Instead it sits above at four of five truths (+0.05 at M0
    16085), with AMP below at all five.
  - GTP down, GDP up and ATP up looks like the nucleotide recycling at a different
    operating point.
- **The hypothesis, untested:** a Jensen gap in the nucleotide kinetics, or the
  mean-field feedback, or both. A pool's own spread does not predict its offset:
  GTP's per-cell CV is 0.24–0.45, but Pi, 3PG and 2PG reach 1.5 at M0 16085 and
  carry no consistent offset, and GDP's CV (0.43) is above GTP's (0.36) there.
- **It varies with θ by about 2×** across the five truths (GTP −0.075 to −0.17).
- **Not yet excluded:**
  - the mean-field feedback: the rate constants, translation's included, are
    rebuilt from the mean pools, not per cell. A3's `popmean` arm has this
    reduction alone, with stochastic gene expression but frozen population-mean
    pools (#92, `eno_feedback_88_result.md`). It moves GTP by 1.17 nats at reduced
    truth 1 and 13.2 at truth 2 (noise 2.32 and 7.88), against the mean field's 31
    and 116. At truth 2 that is above its noise, but far short of the mean field's
    shift: the feedback, as A3 froze it, does not account for most of GTP's gap. A
    live version of that arm on M0 would test this at all five truths;
  - the SSA side's clipping (above).

## Options for the rethink (none taken)

1. **Correct the mean with an offset measured at a reference θ,** δ(θ_ref) = log x̄_SSA − log x̄_MF.
   - It comes from the same ensemble that already gives S (D-S), and is held fixed
     like S.
   - **Tested on these outputs (table below): it does not transfer across prior-wide
     truths.** A δ from one truth applied to another passes KG at no pair.
     - Between M0's truths 13 and 2, the farthest apart in ptsG (1.13 and 33.5),
       it brings GTP within its noise one way only: 0.153 against 0.194 from 13 to
       2, but 0.153 against 0.136 from 2 to 13.
     - Elsewhere it often makes the total worse than no correction (16085 → 13:
       44.2 against 12.4 uncorrected).
   - **Still open:** whether δ transfers across a posterior's width, which is far
     narrower than these truths' spread. That needs truths a posterior SD apart.
   - It would be a third declared reduction beside a fixed S.
2. **A second-order description:**
   - moment closure or the linear noise approximation through the metabolism, so
     the mean carries its own Jensen term (D-S's alternative for S);
   - more framework work, but θ-dependent and differentiable.
3. **Leave the guanylates out of the Gaussian factor** (GTP, GDP and GMP; possibly
   AMP).
   - It changes 5.2's data, so the exact fits would drop them too.
   - It would cut out the pools that carry most of R.
4. **Route (i) instead of (ii):** a stochastic ensemble per θ, with no gradients. The
   sub-spec rules it out on cost.

## Cross-applied offsets (option 1's test)

δ is taken per pool and save from the first truth, added to the second truth's mean
field, and scored as KG (σ_b = 10%). From the `diagnose` stage (job 18892098), over
the archived outputs, with no new simulation.

| Offset from → applied to | Total | Noise | Pools above 3× noise | GTP shift (noise) |
|---|---|---|---|---|
| M0 16085 → 13 | 44.2 | 14.9 | 7 | 7.63 (0.136) |
| M0 16085 → 2 | 63.4 | 17.4 | 11 | 6.75 (0.194) |
| M0 13 → 16085 | 44.2 | 29.2 | 6 | 7.63 (1.05) |
| M0 13 → 2 | 19.6 | 17.4 | 4 (G3P, DHAP, FDP, ATP) | 0.153 (0.194) |
| M0 2 → 16085 | 63.4 | 29.2 | 9 | 6.75 (1.05) |
| M0 2 → 13 | 19.6 | 14.9 | 6 | 0.153 (0.136) |
| Reduced 1 → 2 | 161 | 55.7 | 10 | 32.5 (7.88) |
| Reduced 2 → 1 | 161 | 115 | 5 | 32.5 (2.32) |

The uncorrected totals are 70.8, 12.4 and 25.8 for M0 (16085, 13, 2), and 99.9 and
202 for the reduced model (1, 2).

## The full merge output

Merged at 9686014. σ_b = 10%. Scores in nats; they scale as 1/σ_b².

## m0_v6 (M0, truth seed 16085)

S_JCVISYN3A_0607 4.433, S_JCVISYN3A_0779 2.751, krnadeg 9.424, kcatF_R_ENO 65.48, kcatF_R_FBA 67.13, kcatR_R_FBA 0.6297, kcatR_R_ENO 24.7.

- Reference: 200 cells at 9686014 (job 18886224). Mean field at 9686014 (job 18886227).
- 10 saves. Protein rounding at most 0.5 relative, in `Pcyto_JCVISYN3A_0779` (0.0006 at ≥ 100 copies). 0 of 600 mean-field drains clipped.
- **Verdict (KG): FAIL**. Total 70.8 against noise 29.2; above 3× their noise: `M_g3p_c`, `M_dhap_c`, `M_gtp_c`, `M_atp_c`, `M_amp_c`, `M_gmp_c`, `M_gdp_c`, `M_adp_c`.

| Pool | noise | shift | max \|Δ log x̄\| / σ_b |
|---|---|---|---|
| `M_g6p_c` | 0.304 | 0.263 | 0.41 |
| `M_f6p_c` | 0.25 | 0.292 | 0.45 |
| `M_lac__L_c` | 3.97 | 3.69 | 1.72 |
| `M_g3p_c` | 0.0337 | **0.126** | 0.25 |
| `M_dhap_c` | 0.00711 | **0.0329** | 0.13 |
| `M_trna_chg_c` | 0.0101 | 0.00737 | 0.07 |
| `M_pi_c` | 7.52 | 5.63 | 1.99 |
| `M_3pg_c` | 7.08 | 6.23 | 1.85 |
| `M_2pg_c` | 6.52 | 6.68 | 2.18 |
| `M_gtp_c` | 1.05 | **22.8** | 3.57 |
| `M_fdp_c` | 0.00088 | 0.00241 | 0.04 |
| `M_atp_c` | 0.0373 | **1.73** | 0.76 |
| `M_lac__L_e` | 0.102 | 0.0225 | 0.12 |
| `M_amp_c` | 0.122 | **5.75** | 1.37 |
| `M_gmp_c` | 1.02 | **5.88** | 1.38 |
| `M_gdp_c` | 1.22 | **11.5** | 2.10 |
| `M_adp_c` | 0.00437 | **0.139** | 0.24 |
| **total** | **29.2** | **70.8** | |

## m0_low (M0, truth seed 16200013)

S_JCVISYN3A_0607 12.35, S_JCVISYN3A_0779 1.131, krnadeg 2.356, kcatF_R_ENO 74.58, kcatF_R_FBA 77.55, kcatR_R_FBA 0.7274, kcatR_R_ENO 28.13.

- Reference: 200 cells at 9686014 (job 18886225). Mean field at 9686014 (job 18886228).
- 10 saves. Protein rounding at most 0.5 relative, in `Pcyto_JCVISYN3A_0779` (0.0006 at ≥ 100 copies). 0 of 600 mean-field drains clipped.
- **Verdict (KG): FAIL**. Total 12.4 against noise 14.9; above 3× their noise: `M_gtp_c`, `M_fdp_c`, `M_atp_c`, `M_amp_c`, `M_gmp_c`, `M_gdp_c`, `M_adp_c`.

| Pool | noise | shift | max \|Δ log x̄\| / σ_b |
|---|---|---|---|
| `M_g6p_c` | 0.311 | 0.0802 | 0.29 |
| `M_f6p_c` | 0.275 | 0.0684 | 0.27 |
| `M_lac__L_c` | 1.34 | 0.388 | 0.47 |
| `M_g3p_c` | 0.0751 | 0.116 | 0.20 |
| `M_dhap_c` | 0.00987 | 0.0195 | 0.08 |
| `M_trna_chg_c` | 0.0366 | 0.0117 | 0.10 |
| `M_pi_c` | 3.46 | 0.512 | 0.65 |
| `M_3pg_c` | 4.83 | 0.773 | 0.79 |
| `M_2pg_c` | 4.28 | 0.648 | 0.55 |
| `M_gtp_c` | 0.136 | **7.46** | 1.80 |
| `M_fdp_c` | 9.55e-05 | **0.000882** | 0.02 |
| `M_atp_c` | 0.00358 | **0.14** | 0.27 |
| `M_lac__L_e` | 0.0257 | 0.00622 | 0.05 |
| `M_amp_c` | 0.0183 | **0.663** | 0.59 |
| `M_gmp_c` | 0.0486 | **0.394** | 0.41 |
| `M_gdp_c` | 0.0681 | **1.09** | 0.70 |
| `M_adp_c` | 0.00224 | **0.024** | 0.12 |
| **total** | **14.9** | **12.4** | |

## m0_high (M0, truth seed 16200002)

S_JCVISYN3A_0607 7.613, S_JCVISYN3A_0779 33.52, krnadeg 14.66, kcatF_R_ENO 52.29, kcatF_R_FBA 21.15, kcatR_R_FBA 0.1984, kcatR_R_ENO 19.72.

- Reference: 200 cells at 9686014 (job 18886226). Mean field at 9686014 (job 18886229).
- 10 saves. Protein rounding at most 0.49 relative, in `Pcyto_JCVISYN3A_0779` (0.00059 at ≥ 100 copies). 0 of 600 mean-field drains clipped.
- **Verdict (KG): FAIL**. Total 25.8 against noise 17.4; above 3× their noise: `M_g3p_c`, `M_dhap_c`, `M_gtp_c`, `M_fdp_c`, `M_atp_c`, `M_amp_c`, `M_gmp_c`, `M_gdp_c`, `M_adp_c`.

| Pool | noise | shift | max \|Δ log x̄\| / σ_b |
|---|---|---|---|
| `M_g6p_c` | 0.608 | 0.103 | 0.29 |
| `M_f6p_c` | 0.548 | 0.0953 | 0.27 |
| `M_lac__L_c` | 2.04 | 0.252 | 0.39 |
| `M_g3p_c` | 0.363 | **2.13** | 0.84 |
| `M_dhap_c` | 0.142 | **0.964** | 0.55 |
| `M_trna_chg_c` | 0.0813 | 0.0202 | 0.14 |
| `M_pi_c` | 4.95 | 2 | 0.86 |
| `M_3pg_c` | 4.4 | 6.04 | 1.73 |
| `M_2pg_c` | 3.62 | 2.62 | 1.21 |
| `M_gtp_c` | 0.194 | **6.98** | 1.46 |
| `M_fdp_c` | 6.48e-05 | **0.00204** | 0.02 |
| `M_atp_c` | 0.0187 | **0.354** | 0.34 |
| `M_lac__L_e` | 0.0255 | 0.0114 | 0.07 |
| `M_amp_c` | 0.118 | **1.49** | 0.74 |
| `M_gmp_c` | 0.111 | **0.585** | 0.52 |
| `M_gdp_c` | 0.201 | **2.13** | 0.92 |
| `M_adp_c` | 0.0084 | **0.0338** | 0.13 |
| **total** | **17.4** | **25.8** | |

## corea1 (reduced syn3A, truth seed 1507)

S_JCVISYN3A_0607 2.973, S_JCVISYN3A_0445 0.5672, S_JCVISYN3A_0779 4.364, krnadeg 8.775, kcatF_R_ENO 57.46, kcatF_R_FBA 32.11, kcatR_R_FBA 0.3012, kcatR_R_ENO 21.67.

- Reference: 200 cells at 8813dca (job 18825531), #88 A3's `live60`. Mean field at 9686014 (job 18886230).
- 105 saves. Protein rounding at most 0.5 relative, in `Pcyto_JCVISYN3A_0779` (0.0027 at ≥ 100 copies). 0 of 6300 mean-field drains clipped.
- **Verdict (KG): FAIL**. Total 99.9 against noise 115; above 3× their noise: `M_dhap_c`, `M_gtp_c`, `M_atp_c`, `M_amp_c`, `M_gdp_c`.

| Pool | noise | shift | max \|Δ log x̄\| / σ_b |
|---|---|---|---|
| `M_g6p_c` | 5.03 | 1.86 | 0.57 |
| `M_f6p_c` | 4.63 | 1.74 | 0.55 |
| `M_lac__L_c` | 14.5 | 5.48 | 0.90 |
| `M_g3p_c` | 2.45 | 6.76 | 0.62 |
| `M_dhap_c` | 0.788 | **2.45** | 0.36 |
| `M_trna_chg_c` | 0.62 | 0.273 | 0.19 |
| `M_pi_c` | 25.4 | 10.4 | 1.39 |
| `M_3pg_c` | 30.4 | 16 | 1.78 |
| `M_2pg_c` | 26.8 | 11.1 | 1.43 |
| `M_gtp_c` | 2.32 | **31** | 1.15 |
| `M_fdp_c` | 0.0073 | 0.0033 | 0.01 |
| `M_atp_c` | 0.233 | **1.92** | 0.27 |
| `M_lac__L_e` | 0.604 | 0.162 | 0.10 |
| `M_amp_c` | 0.833 | **6.35** | 0.52 |
| `M_gmp_c` | 0.355 | 0.187 | 0.25 |
| `M_gdp_c` | 0.225 | **4.2** | 0.49 |
| `M_adp_c` | 0.0556 | 0.0397 | 0.07 |
| **total** | **115** | **99.9** | |

## corea2 (reduced syn3A, truth seed 8802)

S_JCVISYN3A_0607 5.455, S_JCVISYN3A_0445 1.931, S_JCVISYN3A_0779 0.7992, krnadeg 2.698, kcatF_R_ENO 85.53, kcatF_R_FBA 73.14, kcatR_R_FBA 0.6861, kcatR_R_ENO 32.26.

- Reference: 200 cells at 8813dca (job 18825539), #88 A3's `live60`. Mean field at 9686014 (job 18886231).
- 105 saves. Protein rounding at most 0.5 relative, in `Pcyto_JCVISYN3A_0779` (0.0027 at ≥ 100 copies). 0 of 6300 mean-field drains clipped.
- **Verdict (KG): FAIL**. Total 202 against noise 55.7; above 3× their noise: `M_g6p_c`, `M_f6p_c`, `M_lac__L_c`, `M_g3p_c`, `M_dhap_c`, `M_gtp_c`, `M_fdp_c`, `M_lac__L_e`, `M_gdp_c`.

| Pool | noise | shift | max \|Δ log x̄\| / σ_b |
|---|---|---|---|
| `M_g6p_c` | 1.56 | **7.39** | 0.84 |
| `M_f6p_c` | 1.59 | **7.8** | 0.85 |
| `M_lac__L_c` | 4.95 | **16** | 1.25 |
| `M_g3p_c` | 1.22 | **4.83** | 0.59 |
| `M_dhap_c` | 0.119 | **0.637** | 0.22 |
| `M_trna_chg_c` | 3.1 | 7.39 | 0.89 |
| `M_pi_c` | 9.2 | 14.5 | 1.35 |
| `M_3pg_c` | 11.9 | 9.86 | 1.16 |
| `M_2pg_c` | 10.2 | 9.72 | 1.28 |
| `M_gtp_c` | 7.88 | **116** | 2.41 |
| `M_fdp_c` | 0.00309 | **0.0479** | 0.07 |
| `M_atp_c` | 0.649 | 0.152 | 0.12 |
| `M_lac__L_e` | 0.117 | **1.79** | 0.37 |
| `M_amp_c` | 1.94 | 2.17 | 0.39 |
| `M_gmp_c` | 0.927 | 0.644 | 0.26 |
| `M_gdp_c` | 0.19 | **1.99** | 0.34 |
| `M_adp_c` | 0.0907 | 0.207 | 0.12 |
| **total** | **55.7** | **202** | |

