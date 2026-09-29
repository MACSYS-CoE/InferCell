# Tasks 15.6 and 15.8: the ensemble sensitivity

Merged by `dev/scripts/corea_sensitivity_15_merge.jl` from 1470 task files (runs at commit 7835c20, job 17704779; merged at 67792ad). Regenerate with `sbatch --array=0-1469 --export=ALL,NSEEDS=1470 dev/scripts/corea_sensitivity_15.slurm`, then this script. 1470 paired seeds, each through 15 configurations on the sampler model with D11's six freed: nominal, and ln θ ± 0.2 on each column. Rows are 42 candidate observables × 106 save points, each in units of its resolution: the 200-cell standard error of its nominal mean, floored at 1% of max(|mean|, floor) as in 14c.5's gate. The noise matrix is (J_A − J_B)/2 over two halves of 735 seeds.

## 15.8: identifiability (K6)

**D11's six** (`S_JCVISYN3A_0607`, `S_JCVISYN3A_0445`, `S_JCVISYN3A_0779`, `krnadeg`, `kcatF_R_ENO`, `kcatF_R_FBA`): singular values 923, 400, 254, 215, 184, 119. Noise floor 72.3. Rank 6 of 6 — **identifiable, K6 does not fire**. Condition number 7.74, which the floor caps at about 12.8, so it is not held to 1e6.

**With the polymerase direction** (all 17 promoters scaled together): singular values 1.4e+03, 406, 254, 221, 204, 177, 116. Noise floor 72.5. Rank 7 of 7. The seventh singular value is 116, 0.969 of the six-set's smallest. No drop is seen. The planning entry expected at most a partial ridge, since the 14 fixed promoters anchor the scale; that explanation is not tested here.

## 15.6: influence audit

Largest |∂ mean / ∂ ln θ| over the save points, in resolution units, per observable and target. An entry is marked with * when it exceeds three times the noise at the same entry. Observables are ordered by their largest entry.

| observable | S_0607 | S_0445 | S_0779 | krnadeg | kcatF_R_ENO | kcatF_R_FBA | polymerase |
|---|---|---|---|---|---|---|---|
| `M_g6p_c` | 11.26* | 3.43* | 10.86* | 36.12* | 10.07* | 1.74* | 49.22* |
| `M_f6p_c` | 11.20* | 4.21* | 10.67* | 35.65* | 10.61* | 1.73* | 49.21* |
| `M_lac__L_c` | 8.80* | 4.51 | 7.73* | 36.45* | 2.74* | 2.54* | 43.24* |
| `M_g3p_c` | 8.04* | 3.93 | 6.10* | 30.43* | 2.99* | 34.80* | 34.69* |
| `M_dhap_c` | 2.95* | 1.53 | 2.51* | 12.50* | 1.32* | 34.13* | 12.35* |
| `M_trna_chg_c` | 7.32* | 3.17* | 8.41* | 32.22* | 2.69 | 1.06* | 32.73* |
| `M_trna_c` | 7.32* | 3.17* | 8.41* | 32.22* | 2.69 | 1.06* | 32.73* |
| `M_pi_c` | 8.75* | 3.55* | 7.52* | 24.32* | 3.03* | 1.28 | 31.03* |
| `M_3pg_c` | 6.92 | 3.64 | 6.79* | 24.38* | 15.48* | 1.48 | 30.40* |
| `M_2pg_c` | 6.50* | 3.05* | 5.65* | 21.33* | 29.54* | 0.97* | 27.44* |
| `M_gtp_c` | 4.23* | 2.46* | 6.26* | 24.62* | 3.19* | 1.89* | 16.42* |
| `mRNA_JCVISYN3A_0607` | 22.72* | 2.72* | 3.33 | 24.51* | 2.55* | 1.02* | 23.11* |
| `M_fdp_c` | 3.48* | 0.76* | 1.85* | 2.90* | 0.22* | 0.87* | 21.62* |
| `M_atp_c` | 4.94* | 1.52 | 6.65* | 20.10* | 4.15* | 1.13* | 21.13* |
| `mRNA_JCVISYN3A_0475` | 4.20* | 4.02 | 3.57* | 19.83* | 2.46* | 0.88 | 21.13* |
| `mRNA_JCVISYN3A_0213` | 2.45* | 4.09* | 3.09* | 18.28* | 2.53 | 1.17* | 19.62* |
| `mRNA_JCVISYN3A_0779` | 2.47* | 2.77* | 18.90* | 16.00* | 2.91* | 0.93 | 17.44* |
| `mRNA_JCVISYN3A_0131` | 3.00 | 2.56* | 2.74 | 17.24* | 1.60* | 1.16* | 17.85* |
| `mRNA_JCVISYN3A_0221` | 3.38* | 2.74* | 3.75* | 15.17* | 2.85* | 1.14 | 16.72* |
| `M_lac__L_e` | 2.37* | 0.62* | 2.46* | 11.97* | 0.10 | 0.05* | 15.50* |
| `mRNA_JCVISYN3A_0694` | 3.94* | 4.06* | 3.07 | 14.17* | 3.35* | 1.77 | 11.34* |
| `M_amp_c` | 3.29* | 2.46 | 4.75* | 14.06* | 8.69* | 2.38* | 11.21* |
| `mRNA_JCVISYN3A_0233` | 3.04 | 4.64* | 2.96* | 14.05* | 3.53 | 1.48* | 11.00* |
| `mRNA_JCVISYN3A_0220` | 2.80 | 2.71 | 3.41* | 14.05* | 2.75 | 0.92 | 13.74* |
| `M_pep_c` | 3.20* | 1.44* | 2.38* | 10.54* | 0.93* | 0.34* | 14.04* |
| `mRNA_JCVISYN3A_0729` | 4.03 | 3.18 | 4.25* | 13.30* | 1.51 | 1.02 | 12.29* |
| `M_gmp_c` | 2.76* | 2.28* | 4.28* | 13.26* | 4.41* | 3.24* | 10.32* |
| `mRNA_JCVISYN3A_0606` | 2.86* | 3.07 | 2.30 | 12.29* | 2.82* | 1.32* | 13.23* |
| `mRNA_JCVISYN3A_0234` | 2.88* | 3.29* | 4.34* | 11.53* | 2.77* | 1.19 | 12.90* |
| `mRNA_JCVISYN3A_0727` | 3.83* | 3.03* | 3.37 | 12.68* | 2.48* | 1.39 | 12.81* |
| `mRNA_JCVISYN3A_0445` | 2.58 | 12.48* | 3.81* | 10.41* | 1.99* | 0.98 | 11.71* |
| `mRNA_JCVISYN3A_0203` | 3.35* | 3.06 | 3.81* | 8.97* | 2.84* | 1.38* | 10.97* |
| `mRNA_JCVISYN3A_0651` | 4.64* | 3.15 | 3.49* | 10.61* | 2.68 | 1.63* | 10.41* |
| `M_gdp_c` | 2.72* | 1.38* | 3.35* | 8.11* | 2.80* | 1.19* | 10.49* |
| `mRNA_JCVISYN3A_0344` | 2.96* | 3.02 | 3.68* | 10.22* | 1.82 | 1.47* | 9.11* |
| `M_adp_c` | 2.00* | 0.94 | 1.47* | 8.31* | 2.27* | 0.61* | 8.09* |
| `M_ppi_c` | 5.00* | 3.76* | 4.51* | 5.14* | 3.04* | 1.49 | 6.82 |
| `volume_litres` | 0.16* | 0.07 | 3.29* | 2.42* | 0.19* | 0.01 | 2.53* |
| `M_pyr_c` | 0.27 | 0.13* | 0.24* | 0.92* | 0.68* | 0.04 | 1.20* |
| `M_nadh_c` | 0.12 | 0.07* | 0.13* | 0.46* | 0.05* | 0.02 | 0.60* |
| `M_13dpg_c` | 0.04 | 0.02* | 0.04* | 0.14* | 0.01* | 0.01* | 0.19* |
| `M_nad_c` | 0.00 | 0.00* | 0.00* | 0.01* | 0.00* | 0.00 | 0.01* |

## The metabolite panel (15.5 into 15.7)

Kept: a pool check 1b does not exclude (`M_ppi_c`, `M_13dpg_c`, `M_nadh_c`, `M_pep_c`), other than `M_trna_c` (conserved with `M_trna_chg_c`, which is kept), where one of D11's six moves it by at least one resolution unit per unit ln θ, beyond three times its noise, at that target's largest save point. F2's rank is shown and does not select. Both tests are weak for a row near its noise, so a pool kept by a margin of less than 2× is marked.

| metabolite | F2 rank | largest entry over the six | kept |
|---|---|---|---|
| `M_g6p_c` | 11 | 36.12 | yes |
| `M_f6p_c` | 10 | 35.65 | yes |
| `M_lac__L_c` | 15 | 36.45 | yes |
| `M_g3p_c` | 17 | 34.80 | yes |
| `M_dhap_c` | 18 | 34.13 | yes |
| `M_trna_chg_c` | 16 | 32.22 | yes |
| `M_trna_c` | 13 | 32.22 | no (conserved pair) |
| `M_pi_c` | 5 | 24.32 | yes |
| `M_3pg_c` | 12 | 24.38 | yes |
| `M_2pg_c` | 9 | 29.54 | yes |
| `M_gtp_c` | 6 | 24.62 | yes |
| `M_fdp_c` | 21 | 3.48 | yes |
| `M_atp_c` | 4 | 20.10 | yes |
| `M_lac__L_e` | – | 11.97 | yes |
| `M_amp_c` | 1 | 14.06 | yes |
| `M_pep_c` | 14 | 10.54 | excluded (check 1b) |
| `M_gmp_c` | 7 | 13.26 | yes |
| `M_gdp_c` | 19 | 8.11 | yes |
| `M_adp_c` | 20 | 8.31 | yes |
| `M_ppi_c` | 2 | 5.14 | excluded (check 1b) |
| `M_pyr_c` | 22 | 0.92 | no |
| `M_nadh_c` | 8 | 0.46 | excluded (check 1b) |
| `M_13dpg_c` | 3 | 0.14 | excluded (check 1b) |
| `M_nad_c` | 23 | 0.01 | no |

Panel (17): `M_g6p_c`, `M_f6p_c`, `M_lac__L_c`, `M_g3p_c`, `M_dhap_c`, `M_trna_chg_c`, `M_pi_c`, `M_3pg_c`, `M_2pg_c`, `M_gtp_c`, `M_fdp_c`, `M_atp_c`, `M_lac__L_e`, `M_amp_c`, `M_gmp_c`, `M_gdp_c`, `M_adp_c`.
