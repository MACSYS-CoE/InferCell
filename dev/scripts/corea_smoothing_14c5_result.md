# 14c.5: the sampler model against the published one, over 5000 seeds

Merged by `dev/scripts/corea_smoothing_14c5_merge.jl` from 250 task files. Regenerate with `sbatch dev/scripts/corea_smoothing_14c5.slurm`, then this script. The gate is §12 2026-09-29's (the later entry). The tasks' headers:

- commit 87b88f1, job 17671709, Julia 1.10.5

Each seed is a full cycle under the published model (clamped drain, fractional carry) and the sampler model (drain smoothed at one particle, pools carried continuously). 50 observables: every ODE state in particles with its carried remainder, every transcript, and the volume. Floors: 500 particles for an ODE state, one copy for a transcript, none for the volume.

## Coupled half

1950 of 5000 pairs decouple within the cycle, at a median handshake of 3672 (earliest 78). Before decoupling there are 21988700 paired comparisons. The largest is 3.75% (`M_gtp_c`, seed 18227, t = 4500 s). Gate: every one within 1% — **fails**.

## Decoupled half, at the end of the cycle

Per observable: d is the difference of the ensemble means over max(|published mean|, floor), SE is from the paired differences, and R is the 200-cell dataset's standard error on the mean. The observable passes when |d| + 2·SE ≤ max(R, 1%). **50 pass, 0 fail, 0 unresolved.**

| observable | d | SE | R | (\|d\| + 2·SE) / tolerance | verdict |
|---|---|---|---|---|---|
| `mRNA_JCVISYN3A_0203` | +0.84% | 0.55% | 3.7% | 0.53 | passes |
| `mRNA_JCVISYN3A_0694` | +1% | 0.67% | 4.5% | 0.52 | passes |
| `mRNA_JCVISYN3A_0606` | -1.3% | 0.77% | 5.6% | 0.51 | passes |
| `mRNA_JCVISYN3A_0475` | -1.32% | 0.8% | 5.8% | 0.50 | passes |
| `mRNA_JCVISYN3A_0727` | +1.14% | 0.74% | 5.4% | 0.48 | passes |
| `mRNA_JCVISYN3A_0213` | +1.1% | 0.81% | 6% | 0.45 | passes |
| `mRNA_JCVISYN3A_0131` | +1.14% | 0.92% | 7% | 0.42 | passes |
| `mRNA_JCVISYN3A_0233` | +0.68% | 0.73% | 5.1% | 0.42 | passes |
| `M_3pg_c` | -0.494% | 0.52% | 3.8% | 0.40 | passes |
| `mRNA_JCVISYN3A_0729` | +0.34% | 0.75% | 4.8% | 0.39 | passes |
| `mRNA_JCVISYN3A_0607` | +0.701% | 0.67% | 5.3% | 0.38 | passes |
| `mRNA_JCVISYN3A_0651` | -0.34% | 0.58% | 3.9% | 0.38 | passes |
| `mRNA_JCVISYN3A_0234` | +0.34% | 0.67% | 4.8% | 0.35 | passes |
| `M_g3p_c` | -0.0828% | 0.15% | 1.1% | 0.34 | passes |
| `mRNA_JCVISYN3A_0779` | +0.347% | 0.95% | 6.6% | 0.34 | passes |
| `mRNA_JCVISYN3A_0221` | +0.28% | 0.89% | 6.3% | 0.33 | passes |
| `mRNA_JCVISYN3A_0344` | +0.1% | 0.55% | 3.7% | 0.32 | passes |
| `mRNA_JCVISYN3A_0445` | -0.28% | 0.55% | 4.4% | 0.31 | passes |
| `mRNA_JCVISYN3A_0220` | -0.3% | 0.72% | 5.8% | 0.30 | passes |
| `M_trna_chg_c` | -0.0375% | 0.26% | 1.9% | 0.30 | passes |
| `M_trna_c` | +0.0673% | 0.47% | 3.4% | 0.30 | passes |
| `M_pi_c` | -0.109% | 0.51% | 3.8% | 0.30 | passes |
| `M_ppi_c` | +0.0207% | 0.43% | 3.1% | 0.28 | passes |
| `M_lac__L_c` | +0.048% | 0.34% | 2.7% | 0.27 | passes |
| `M_amp_c` | -0.0993% | 0.17% | 1.7% | 0.26 | passes |
| `M_gmp_c` | -0.0729% | 0.12% | 1.2% | 0.26 | passes |
| `M_g6p_c` | -0.0141% | 0.21% | 1.7% | 0.25 | passes |
| `M_gtp_c` | +0.092% | 0.5% | 4.5% | 0.25 | passes |
| `M_f6p_c` | -0.0153% | 0.2% | 1.7% | 0.24 | passes |
| `M_2pg_c` | -0.0539% | 0.086% | 0.63% | 0.23 | passes |
| `M_atp_c` | +0.0327% | 0.12% | 1.3% | 0.21 | passes |
| `M_ptsi_c` | -0.0248% | 0.11% | 1.3% | 0.20 | passes |
| `M_ptsh_c` | +0.0993% | 0.044% | 0.52% | 0.19 | passes |
| `M_dhap_c` | -0.05% | 0.066% | 0.5% | 0.18 | passes |
| `M_ptsg_c` | +0.0143% | 0.076% | 0.91% | 0.17 | passes |
| `M_crr_c` | +0.0161% | 0.069% | 0.81% | 0.15 | passes |
| `M_gdp_c` | +0.0138% | 0.067% | 0.54% | 0.15 | passes |
| `M_pep_c` | -0.0213% | 0.035% | 0.27% | 0.09 | passes |
| `M_adp_c` | -0.0175% | 0.027% | 0.25% | 0.07 | passes |
| `M_ptsg_P_c` | -0.0124% | 0.02% | 0.15% | 0.05 | passes |
| `M_fdp_c` | +0.02% | 0.0088% | 0.12% | 0.04 | passes |
| `M_lac__L_e` | +0.0136% | 0.0097% | 0.14% | 0.03 | passes |
| `volume_litres` | +0.00104% | 0.0078% | 0.093% | 0.02 | passes |
| `M_pyr_c` | +0.00313% | 0.0035% | 0.026% | 0.01 | passes |
| `M_crr_P_c` | -0.00202% | 0.0038% | 0.029% | 0.01 | passes |
| `M_ptsh_P_c` | -0.00107% | 0.003% | 0.024% | 0.01 | passes |
| `M_nadh_c` | -0.00088% | 0.0019% | 0.014% | 0.00 | passes |
| `M_13dpg_c` | -0.000134% | 0.00053% | 0.0041% | 0.00 | passes |
| `M_ptsi_P_c` | -0.000197% | 0.00034% | 0.0025% | 0.00 | passes |
| `M_nad_c` | +9.86e-06% | 2.1e-05% | 0.00016% | 0.00 | passes |

## The null: the published model against itself on independent seeds

The published runs at the first 2500 seeds against those at the other 2500, unpaired. 50 observables vary across seeds. For the new statistic the largest |z| is 1.75, against 3.29 (5%, Bonferroni over 50) — **passes**. The old per-seed ratio, on the same pairs:

| observable | old: mean per-seed ratio | its z | new: d | its z |
|---|---|---|---|---|
| `M_gtp_c` | +88.3% ± 5.3% | 16.78 | -0.411% ± 1.7% | -0.23 |
| `mRNA_JCVISYN3A_0607` | +40.6% ± 2.7% | 15.06 | -2.4% ± 2.1% | -1.15 |
| `M_pi_c` | +46.8% ± 3.3% | 14.29 | -0.195% ± 1.5% | -0.13 |
| `mRNA_JCVISYN3A_0475` | +35.5% ± 2.5% | 13.95 | -0.527% ± 2.3% | -0.23 |
| `M_trna_c` | +28.5% ± 2.1% | 13.51 | -0.382% ± 1.3% | -0.29 |
| `mRNA_JCVISYN3A_0213` | +33.5% ± 2.5% | 13.19 | -0.275% ± 2.4% | -0.11 |
| `mRNA_JCVISYN3A_0779` | +29.8% ± 2.4% | 12.30 | +2.1% ± 2.7% | 0.79 |
| `mRNA_JCVISYN3A_0131` | +28.6% ± 2.4% | 11.97 | +3.82% ± 2.9% | 1.34 |
| `M_trna_chg_c` | +11.6% ± 1.3% | 9.00 | +0.213% ± 0.74% | 0.29 |
| `M_amp_c` | +6.42% ± 0.72% | 8.89 | +0.676% ± 0.67% | 1.00 |
| `M_f6p_c` | +7.25% ± 0.82% | 8.82 | +0.388% ± 0.69% | 0.56 |
| `M_g6p_c` | +7.32% ± 0.84% | 8.75 | +0.324% ± 0.7% | 0.46 |

## Clipping and wall-clock

| | median drains clipped | seeds carrying a deficit | median wall-clock per cycle (s) |
|---|---|---|---|
| `published` | 50.0 | 4843 of 5000 | 33.3 |
| `sampler` | 51.0 | 4849 of 5000 | 28.4 |

## Verdict

Coupled half: **fails**. Decoupled half: 50 of 50 pass, 0 fail, 0 unresolved. Null: passes. 14c.5 is **not closed**.

*Note added at merge.* The verdict above is the gate's, and it stands. §12
2026-09-29 (the entry on the coupled half) records the failure as the
smoothing's measured cost, from the seed-18227 rerun in
`corea_smoothing_14c5_diag_result.md`, and closes 14c.5 on that record.
