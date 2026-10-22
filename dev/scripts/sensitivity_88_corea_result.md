# Issue #88 A2: 15.8 rerun on Core A′, with the raw outputs kept

**Run:** array 18821442 (294 tasks, 5 seeds each) and merge 18821443, both at `8813dca`,
submitted from a detached worktree at that pushed commit. Driver
`dev/scripts/sensitivity_88.jl corea`, merge `dev/scripts/sensitivity_88_merge.jl`.
**Cost:** 165.83 CPU-h for the array (the longest task 37 min), 0.02 for the merge.
**Re-merged** after the #92 review's simplifications, by job 18836375 at `2b039a0`: every
number is identical. The original merge's output is kept as `merged_8813dca.md`.
**Raw outputs:** `/fred/oz022/tkimpson/infercell-88-a2-archive-2026-10-19/`.
- 294 task files with `MD5SUMS`, 421 MiB on disk.
- `merged.md`, `merged.jls` and the Slurm logs.
- Every run's full species × times matrix, so 15.8's raw outputs, deleted in the
  2026-10-13 cleanup, exist again outside every worktree.

**Regenerate:**
`dev/scripts/sensitivity_88.sh corea 1470 294 <outdir outside the worktrees>`.

## Setup

- **15.8's design:** the sampler model (`d11_models` with `COREA_SAMPLER`), one 6300 s
  cycle, saves every 60 s, D11's six freed.
- **Seeds:** 1470 paired seeds (20001–21470, as 15.8), each run nominal and
  ln θ ± 0.2 under common random numbers.
- **One change from 15.8:** the polymerase column is dropped.
- **The bound** is `bulk_identifiability_16a7.jl`'s calculation at C = 200.

## 15.8 reproduces

The six-target identifiability check, done 15.8's way (every observable and save,
1% floor, 200 cells), gives singular values **923, 400, 254, 215, 184, 119 against a
floor of 72.3, rank 6 of 6**. That is identical, to the three figures quoted, to 15.8
(job 17704779 at `7835c20`). The 32 `src` files changed since then leave this
measurement unchanged.

## The information bound for D11's six (paper section 5.3's yardstick)

Posterior over prior SD, linearised, with the split-half bracket in parentheses. A
target is informed only if its ratio is below its bracket. All at C = 200, σ_b = 10%:

| Design | Rank (floor) | S_GAPD | S_PGI | S_ptsG | krnadeg | ENO | FBA |
|---|---|---|---|---|---|---|---|
| metabolites only | 6 of 6 (9.51) | 0.06 (0.18) | 0.14 (0.16) | 0.10 (0.19) | 0.03 (0.16) | 0.13 (0.75) | 0.04 (0.51) |
| transcripts only | 4 of 6 (54.6) | 0.01 (0.03) | 0.01 (0.03) | 0.01 (0.03) | 0.00 (0.03) | 0.17 (0.18) | 0.10 (0.11) |
| joint | 5 of 6 (55.4) | 0.01 (0.03) | 0.01 (0.03) | 0.01 (0.03) | 0.00 (0.03) | 0.10 (0.18) | 0.04 (0.10) |

- **Metabolites only:** every one of the six is informed.
  - ptsG at 0.10 (0.19) and ENO at 0.13 (0.75) are D9's numbers exactly.
  - PGI's margin is thin: 0.14 against 0.16.
- **Transcripts only:**
  - The promoters and krnadeg are pinned, at 0.01 or below.
  - ENO and FBA sit at their brackets (0.17 against 0.18, 0.10 against 0.11), so the
    transcripts' ensemble means carry no ENO or FBA information beyond the noise.
    This is the linearised counterpart of the ENO path-information null (A3).
- **Joint:** ENO 0.10 (0.18) and FBA 0.04 (0.10). The metabolites supply both.
- **σ_b changes the metabolite rows only.** At 5% every ratio roughly halves, and at
  20% roughly doubles. The full tables are below.
- **Caveats:**
  - Save points are treated as independent, which overstates the information.
  - Transcript rows are 200-cell ensemble means at a 1% floor, which understates what
    per-cell counts carry.
  - Both are as in D9.

## R per pool

R ≈ Σ Var(log x̄_C)/σ_b², over the 17 panel pools × 105 saves (n = 1785):

| C | σ_b = 5% | 10% | 20% |
|---|---|---|---|
| 200 | 279 | **70** | 17.5 |
| 50 | 1,120 | 280 | 70 |

- **These are upper bounds.** The cells are not conditioned on the transcripts, and
  the note's R is.
- **They are at the nominal parameter values.** R depends on θ, so this is not a
  prior-wide number. A3 found the reverse feedback's effect on the bulk metabolites
  differs strongly between two truths.
- **Every setting is in the middle regime** (1 ≲ R ≲ n): the metabolites constrain
  the paths collectively. So one path draw per θ is not a usable estimator anywhere
  here, and the Gaussian convolution of note §7, with the C-cell covariance S, is
  the method.
- **R is about 8× the note's guess.** Discussion §7 guessed about 9 at C = 200 and 10%,
  assuming a per-cell CV of 0.1. The measured CVs run to 0.52, with the largest in
  3PG, 2PG, Pi and GTP.
- **Four pools carry 49 of the 70 at C = 200, 10%:** 3PG, 2PG, Pi and GTP.
  - GTP is below one particle in 3% of cell-saves, the one pool here near the floor.
  - 3PG, 2PG and Pi never are on Core A′, unlike on M0 (A1).
- **On M0, A1 measured R = 64 at C = 50, 10%,** over n = 170.

## The full merge output

### Information bound (200 cells)

#### σ_b = 1%

1035 of 1785 metabolite rows have their resolution set by the 200-cell SE, not σ_b.

| Design | Rows | Rank (floor) | S_0607 (GAPD) | S_0445 (PGI) | S_0779 (ptsG) | krnadeg | ENO | FBA |
|---|---|---|---|---|---|---|---|---|
| metabolites only | 1785 | 5 of 6 (45) | 0.02 (0.04) | 0.03 (0.04) | 0.02 (0.05) | 0.01 (0.04) | 0.03 (0.26) | 0.00 (0.14) |
| transcripts only | 1785 | 4 of 6 (54.6) | 0.01 (0.03) | 0.01 (0.03) | 0.01 (0.03) | 0.00 (0.03) | 0.17 (0.18) | 0.10 (0.11) |
| joint | 3570 | 6 of 6 (69.9) | 0.01 (0.02) | 0.01 (0.02) | 0.01 (0.03) | 0.00 (0.02) | 0.03 (0.15) | 0.00 (0.08) |

#### σ_b = 5%

0 of 1785 metabolite rows have their resolution set by the 200-cell SE, not σ_b.

| Design | Rows | Rank (floor) | S_0607 (GAPD) | S_0445 (PGI) | S_0779 (ptsG) | krnadeg | ENO | FBA |
|---|---|---|---|---|---|---|---|---|
| metabolites only | 1785 | 6 of 6 (19) | 0.03 (0.09) | 0.07 (0.08) | 0.05 (0.10) | 0.01 (0.08) | 0.07 (0.49) | 0.02 (0.29) |
| transcripts only | 1785 | 4 of 6 (54.6) | 0.01 (0.03) | 0.01 (0.03) | 0.01 (0.03) | 0.00 (0.03) | 0.17 (0.18) | 0.10 (0.11) |
| joint | 3570 | 6 of 6 (57.7) | 0.01 (0.03) | 0.01 (0.03) | 0.01 (0.03) | 0.00 (0.03) | 0.06 (0.17) | 0.02 (0.10) |

#### σ_b = 10%

0 of 1785 metabolite rows have their resolution set by the 200-cell SE, not σ_b.

| Design | Rows | Rank (floor) | S_0607 (GAPD) | S_0445 (PGI) | S_0779 (ptsG) | krnadeg | ENO | FBA |
|---|---|---|---|---|---|---|---|---|
| metabolites only | 1785 | 6 of 6 (9.51) | 0.06 (0.18) | 0.14 (0.16) | 0.10 (0.19) | 0.03 (0.16) | 0.13 (0.75) | 0.04 (0.51) |
| transcripts only | 1785 | 4 of 6 (54.6) | 0.01 (0.03) | 0.01 (0.03) | 0.01 (0.03) | 0.00 (0.03) | 0.17 (0.18) | 0.10 (0.11) |
| joint | 3570 | 5 of 6 (55.4) | 0.01 (0.03) | 0.01 (0.03) | 0.01 (0.03) | 0.00 (0.03) | 0.10 (0.18) | 0.04 (0.10) |

#### σ_b = 20%

0 of 1785 metabolite rows have their resolution set by the 200-cell SE, not σ_b.

| Design | Rows | Rank (floor) | S_0607 (GAPD) | S_0445 (PGI) | S_0779 (ptsG) | krnadeg | ENO | FBA |
|---|---|---|---|---|---|---|---|---|
| metabolites only | 1785 | 6 of 6 (4.76) | 0.13 (0.34) | 0.28 (0.30) | 0.19 (0.36) | 0.05 (0.30) | 0.26 (0.91) | 0.08 (0.76) |
| transcripts only | 1785 | 4 of 6 (54.6) | 0.01 (0.03) | 0.01 (0.03) | 0.01 (0.03) | 0.00 (0.03) | 0.17 (0.18) | 0.10 (0.11) |
| joint | 3570 | 4 of 6 (54.8) | 0.01 (0.03) | 0.01 (0.03) | 0.01 (0.03) | 0.00 (0.03) | 0.14 (0.18) | 0.06 (0.11) |

#### σ_b = 10%, metabolite floor 1 particle

0 of 1785 metabolite rows have their resolution set by the 200-cell SE, not σ_b.

| Design | Rows | Rank (floor) | S_0607 (GAPD) | S_0445 (PGI) | S_0779 (ptsG) | krnadeg | ENO | FBA |
|---|---|---|---|---|---|---|---|---|
| metabolites only | 1785 | 6 of 6 (10.5) | 0.06 (0.16) | 0.13 (0.14) | 0.09 (0.17) | 0.02 (0.14) | 0.07 (0.70) | 0.04 (0.47) |
| transcripts only | 1785 | 4 of 6 (54.6) | 0.01 (0.03) | 0.01 (0.03) | 0.01 (0.03) | 0.00 (0.03) | 0.17 (0.18) | 0.10 (0.11) |
| joint | 3570 | 5 of 6 (55.6) | 0.01 (0.03) | 0.01 (0.03) | 0.01 (0.03) | 0.00 (0.03) | 0.06 (0.18) | 0.04 (0.10) |

Each cell is posterior over prior SD, linearised, with the split-half noise matrix's ratio in brackets. A target is informed only if its ratio is below its bracket. Save points are treated as independent, which overstates the information.

### R per pool

Σ over 105 save points of Var(log x̄_C) / σ_b², x̄_C the mean of C cells resampled (4000 draws) from 1470 nominal cells, floored at 1 particle. Per-cell CV is the median over save points. "< 1" is the fraction of cell-saves below one particle.

| Pool | per-cell CV | < 1 | R, C = 200, σ_b = 5% | R, C = 200, σ_b = 10% | R, C = 200, σ_b = 20% | R, C = 50, σ_b = 5% | R, C = 50, σ_b = 10% | R, C = 50, σ_b = 20% |
|---|---|---|---|---|---|---|---|---|
| `M_3pg_c` | 0.522 | 0.00 | 58.4 | 14.6 | 3.65 | 233 | 58.2 | 14.6 |
| `M_2pg_c` | 0.468 | 0.00 | 47.1 | 11.8 | 2.94 | 190 | 47.4 | 11.8 |
| `M_pi_c` | 0.467 | 0.00 | 46.5 | 11.6 | 2.91 | 187 | 46.8 | 11.7 |
| `M_gtp_c` | 0.45 | 0.03 | 44.4 | 11.1 | 2.77 | 178 | 44.6 | 11.1 |
| `M_lac__L_c` | 0.314 | 0.00 | 21.6 | 5.41 | 1.35 | 86.6 | 21.6 | 5.41 |
| `M_trna_chg_c` | 0.245 | 0.00 | 12.4 | 3.1 | 0.775 | 49.7 | 12.4 | 3.11 |
| `M_g6p_c` | 0.24 | 0.00 | 11.6 | 2.89 | 0.722 | 46.2 | 11.5 | 2.89 |
| `M_f6p_c` | 0.24 | 0.00 | 11.5 | 2.87 | 0.717 | 45.9 | 11.5 | 2.87 |
| `M_amp_c` | 0.206 | 0.00 | 8.25 | 2.06 | 0.515 | 33 | 8.25 | 2.06 |
| `M_g3p_c` | 0.157 | 0.00 | 5.22 | 1.31 | 0.327 | 20.9 | 5.23 | 1.31 |
| `M_gmp_c` | 0.144 | 0.00 | 4.1 | 1.02 | 0.256 | 16.5 | 4.12 | 1.03 |
| `M_atp_c` | 0.142 | 0.00 | 4.03 | 1.01 | 0.252 | 16.1 | 4.03 | 1.01 |
| `M_lac__L_e` | 0.1 | 0.00 | 2.22 | 0.555 | 0.139 | 8.91 | 2.23 | 0.557 |
| `M_gdp_c` | 0.0665 | 0.00 | 0.994 | 0.249 | 0.0621 | 3.97 | 0.992 | 0.248 |
| `M_dhap_c` | 0.0622 | 0.00 | 0.825 | 0.206 | 0.0516 | 3.31 | 0.828 | 0.207 |
| `M_adp_c` | 0.0395 | 0.00 | 0.322 | 0.0805 | 0.0201 | 1.29 | 0.321 | 0.0803 |
| `M_fdp_c` | 0.0133 | 0.00 | 0.0351 | 0.00877 | 0.00219 | 0.141 | 0.0351 | 0.00879 |
| **total** | | | **279** | **69.9** | **17.5** | **1.12e+03** | **280** | **70** |

n = 17 pools × 105 saves = 1785. Regimes (note §7): R ≪ 1, one path draw per θ suffices; 1 ≲ R ≲ n, Gaussian convolution; pools near the one-particle floor are non-Gaussian.

### 15.8's identifiability check, redone

Every observable and save, 1% floor, 200 cells: singular values 923, 400, 254, 215, 184, 119, floor 72.3, rank 6 of 6. 15.8 (job 17704779 at 7835c20): 923, 400, 254, 215, 184, 119, floor 72.3, rank 6 of 6.

