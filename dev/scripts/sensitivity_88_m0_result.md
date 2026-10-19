# Issue #88 A1: is there a ptsG signal under bulk metabolites at M0's size?

**Run:** array 18821440 (49 tasks) and merge 18821441, both at `8813dca`, submitted from
a detached worktree at that pushed commit. Driver `dev/scripts/sensitivity_88.jl m0`,
merge `dev/scripts/sensitivity_88_merge.jl`.
**Cost:** 11.08 CPU-h for the array, 0.02 for the merge.
**Raw outputs:** `/fred/oz022/tkimpson/infercell-88-a1-archive-2026-10-19/`: 49 task
files (`MD5SUMS`), `merged.md`, `merged.jls` and the Slurm logs.
**Regenerate:**
`dev/scripts/sensitivity_88.sh m0 1470 49 <outdir outside the worktrees>`.

## Setup

- **Model:** M0 as the V6/V7 exact sampler builds it and `m0_dataset` generates it:
  `build_m0()`, 600 s, saves every 60 s, t = 0 dropped, so 10 save points.
- **Targets:** M0's five. Each of 1470 paired seeds (20001–21470) runs nominal and
  ln θ ± 0.2 on each target, under common random numbers.
- **The bound:** `bulk_identifiability_16a7.jl`'s calculation (D9) at C = 50 cells.
  A metabolite row's resolution is max(50-cell SE, σ_b · max(|mean|, 500 particles)).
  A transcript row keeps 15.8's resolution, the SE floored at 1%, which understates
  what per-cell counts carry. The ratio is posterior over prior SD, linearised, with
  the split-half noise matrix's ratio in brackets.
- **Save points are treated as independent,** which overstates the information.

## The gate

ptsG (`S_JCVISYN3A_0779`) is informed iff its ratio is below its bracket, metabolites
only, at σ_b = 10%.

| σ_b | ptsG, metabolites only | Core A′, D9 (200 cells, 6300 s) |
|---|---|---|
| 5% | 0.34 (0.55) | 0.05 (0.10) |
| 10% | **0.47 (0.68)** | 0.10 (0.19) |
| 20% | 0.66 (0.85) | 0.19 (0.36) |
| 10%, floor 1 particle | 0.42 (0.61) | – |

- **The gate passes.** At M0's size the metabolites carry a ptsG signal: 0.47
  against a bracket of 0.68 at 10%. It is below its bracket at every σ_b and under
  either floor.
- **It is much weaker than on Core A′.** At 10% it shrinks to 0.47 of the prior SD
  against D9's 0.10, so 50 cells over 600 s keep about a fifth of the full-scale
  information in this cell. The margin over the noise also narrows as σ_b grows.
- **No larger M0 dataset is needed** for paper section 5.2 on this measurement. That
  holds only if the exact sampler realises the bound. It is an upper bound; F10
  measures the realised value.

## R at M0 (C = 50)

R ≈ Σ Var(log x̄₅₀)/σ_b², over the 17 panel pools × 10 saves:

| σ_b | 5% | 10% | 20% |
|---|---|---|---|
| R | 256 | 64 | 16 |

- **These are upper bounds.** The cells are not conditioned on the transcripts, and
  the note's R is.
- **M0 sits in the middle regime** (1 ≲ R ≲ n = 170) at 10% and 20%. There one path
  draw per θ is not a usable estimator, and the Gaussian convolution of note §7 is
  the method. At 5%, R exceeds n.
- **Three pools carry most of it:** 3PG, Pi and 2PG give 49 of the 64 at 10%.
  - Their per-cell CV is 0.87–0.92.
  - 2–6% of their cell-saves are below one particle, so they are the non-Gaussian
    pools the note warned of.
  - Without them, R is about 15 at 10%.
- **V7 draws σ_b from LogNormal(log 0.2, 1),** so its datasets span R ≈ 4 to 250 or
  more at M0. Discussion §7's table guessed 3.4 to 340, and the measurement falls
  inside it.

## The full merge output

### Information bound (50 cells)

#### σ_b = 5%

37 of 170 metabolite rows have their resolution set by the 50-cell SE, not σ_b.

| Design | Rows | Rank (floor) | S_0607 (GAPD) | S_0779 (ptsG) | krnadeg | ENO | FBA |
|---|---|---|---|---|---|---|---|
| metabolites only | 170 | 4 of 5 (6.04) | 0.34 (0.38) | 0.34 (0.55) | 0.13 (0.47) | 0.44 (1.00) | 0.12 (0.99) |
| transcripts only | 20 | 2 of 5 (2.52) | 0.44 (0.66) | 0.43 (0.81) | 0.48 (0.70) | 1.00 (1.00) | 1.00 (1.00) |
| joint | 190 | 5 of 5 (6.54) | 0.10 (0.36) | 0.10 (0.51) | 0.10 (0.42) | 0.43 (1.00) | 0.12 (0.99) |

#### σ_b = 10%

10 of 170 metabolite rows have their resolution set by the 50-cell SE, not σ_b.

| Design | Rows | Rank (floor) | S_0607 (GAPD) | S_0779 (ptsG) | krnadeg | ENO | FBA |
|---|---|---|---|---|---|---|---|
| metabolites only | 170 | 4 of 5 (3.84) | 0.44 (0.49) | 0.47 (0.68) | 0.20 (0.61) | 0.64 (1.00) | 0.23 (1.00) |
| transcripts only | 20 | 2 of 5 (2.52) | 0.44 (0.66) | 0.43 (0.81) | 0.48 (0.70) | 1.00 (1.00) | 1.00 (1.00) |
| joint | 190 | 5 of 5 (4.59) | 0.15 (0.44) | 0.15 (0.62) | 0.16 (0.53) | 0.63 (1.00) | 0.23 (0.99) |

#### σ_b = 20%

0 of 170 metabolite rows have their resolution set by the 50-cell SE, not σ_b.

| Design | Rows | Rank (floor) | S_0607 (GAPD) | S_0779 (ptsG) | krnadeg | ENO | FBA |
|---|---|---|---|---|---|---|---|
| metabolites only | 170 | 4 of 5 (2.07) | 0.60 (0.69) | 0.66 (0.85) | 0.33 (0.81) | 0.85 (1.00) | 0.43 (1.00) |
| transcripts only | 20 | 2 of 5 (2.52) | 0.44 (0.66) | 0.43 (0.81) | 0.48 (0.70) | 1.00 (1.00) | 1.00 (1.00) |
| joint | 190 | 4 of 5 (3.26) | 0.24 (0.57) | 0.23 (0.74) | 0.26 (0.63) | 0.85 (1.00) | 0.42 (1.00) |

#### σ_b = 10%, metabolite floor 1 particle

30 of 170 metabolite rows have their resolution set by the 50-cell SE, not σ_b.

| Design | Rows | Rank (floor) | S_0607 (GAPD) | S_0779 (ptsG) | krnadeg | ENO | FBA |
|---|---|---|---|---|---|---|---|
| metabolites only | 170 | 4 of 5 (4.79) | 0.36 (0.39) | 0.42 (0.61) | 0.17 (0.54) | 0.30 (1.00) | 0.23 (0.99) |
| transcripts only | 20 | 2 of 5 (2.52) | 0.44 (0.66) | 0.43 (0.81) | 0.48 (0.70) | 1.00 (1.00) | 1.00 (1.00) |
| joint | 190 | 5 of 5 (5.41) | 0.14 (0.37) | 0.14 (0.57) | 0.15 (0.47) | 0.29 (1.00) | 0.23 (0.99) |

Each cell is posterior over prior SD, linearised, with the split-half noise matrix's ratio in brackets. A target is informed only if its ratio is below its bracket. Save points are treated as independent, which overstates the information.

**A1 gate (ptsG, metabolites only, σ_b = 10%):** ratio 0.47, bracket 0.68 — informed: the gate passes.

### R per pool

Σ over 10 save points of Var(log x̄_C) / σ_b², x̄_C the mean of C cells resampled (4000 draws) from 1470 nominal cells, floored at 1 particle. Per-cell CV is the median over save points. "< 1" is the fraction of cell-saves below one particle.

| Pool | per-cell CV | < 1 | R, C = 50, σ_b = 5% | R, C = 50, σ_b = 10% | R, C = 50, σ_b = 20% |
|---|---|---|---|---|---|
| `M_3pg_c` | 0.923 | 0.04 | 71.4 | 17.9 | 4.46 |
| `M_pi_c` | 0.879 | 0.02 | 62.7 | 15.7 | 3.92 |
| `M_2pg_c` | 0.87 | 0.06 | 62.1 | 15.5 | 3.88 |
| `M_lac__L_c` | 0.607 | 0.00 | 26.7 | 6.68 | 1.67 |
| `M_gtp_c` | 0.336 | 0.00 | 8.08 | 2.02 | 0.505 |
| `M_g6p_c` | 0.317 | 0.00 | 6.99 | 1.75 | 0.437 |
| `M_f6p_c` | 0.3 | 0.00 | 6.17 | 1.54 | 0.386 |
| `M_gdp_c` | 0.153 | 0.00 | 4.19 | 1.05 | 0.262 |
| `M_gmp_c` | 0.129 | 0.00 | 2.92 | 0.73 | 0.182 |
| `M_g3p_c` | 0.135 | 0.00 | 1.31 | 0.329 | 0.0821 |
| `M_amp_c` | 0.113 | 0.00 | 1.17 | 0.293 | 0.0733 |
| `M_lac__L_e` | 0.0939 | 0.00 | 0.786 | 0.197 | 0.0491 |
| `M_trna_chg_c` | 0.0912 | 0.00 | 0.634 | 0.158 | 0.0396 |
| `M_dhap_c` | 0.0632 | 0.00 | 0.285 | 0.0713 | 0.0178 |
| `M_atp_c` | 0.0501 | 0.00 | 0.258 | 0.0645 | 0.0161 |
| `M_adp_c` | 0.0319 | 0.00 | 0.086 | 0.0215 | 0.00537 |
| `M_fdp_c` | 0.00739 | 0.00 | 0.00389 | 0.000973 | 0.000243 |
| **total** | | | **256** | **64** | **16** |

n = 17 pools × 10 saves = 170. Regimes (note §7): R ≪ 1, one path draw per θ suffices; 1 ≲ R ≲ n, Gaussian convolution; pools near the one-particle floor are non-Gaussian.

