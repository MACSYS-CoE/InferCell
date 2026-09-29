# Tasks 15.6 and 15.8: the ensemble sensitivity, 100-seed pilot

Merged by `dev/scripts/corea_sensitivity_15_merge.jl` from 100 task files (runs at commit 886a2f7, job 17703391; merged at 886a2f7). 100 paired seeds, each through 15 configurations on the sampler model with D11's six freed: nominal, and ln θ ± 0.2 on each column. Rows are 42 candidate observables × 106 save points, each in units of its resolution: the 200-cell standard error of its nominal mean, floored at 1% of max(|mean|, floor) as in 14c.5's gate. The noise matrix is (J_A − J_B)/2 over two halves of 50 seeds.

## 15.8: identifiability (K6)

**D11's six** (`S_JCVISYN3A_0607`, `S_JCVISYN3A_0445`, `S_JCVISYN3A_0779`, `krnadeg`, `kcatF_R_ENO`, `kcatF_R_FBA`): singular values 993, 412, 344, 328, 302, 266. Noise floor 316. Rank 4 of 6 — **not identifiable, K6 fires**. Condition number 3.73, which the floor caps at about 3.14, so it is not held to 1e6.

**With the polymerase direction** (all 17 promoters scaled together): singular values 1.47e+03, 422, 356, 328, 303, 299, 265. Noise floor 320. Rank 4 of 7. The seventh singular value is 265, 0.995 of the six-set's smallest. With 14 promoters fixed the scale is anchored, so a partial ridge was expected (§12).

## 15.6: influence audit

Largest |∂ mean / ∂ ln θ| over the save points, in resolution units, per observable and target. An entry is marked with * when it exceeds three times the noise at the same entry. Observables are ordered by their largest entry.

| observable | S_0607 | S_0445 | S_0779 | krnadeg | kcatF_R_ENO | kcatF_R_FBA | polymerase |
|---|---|---|---|---|---|---|---|
| `M_g6p_c` | 21.20* | 11.89* | 17.19* | 42.27* | 10.06* | 5.69 | 59.21* |
| `M_f6p_c` | 20.98* | 12.73* | 16.38* | 42.05* | 10.60* | 5.60 | 58.92* |
| `M_lac__L_c` | 16.20* | 12.21* | 16.12* | 43.89* | 9.47* | 4.96* | 49.36* |
| `M_trna_c` | 12.13* | 15.47* | 18.02* | 46.66* | 12.14* | 5.18 | 45.28* |
| `M_trna_chg_c` | 12.13* | 15.47* | 18.02* | 46.66* | 12.14* | 5.18 | 45.28* |
| `M_g3p_c` | 13.73* | 12.72 | 13.14* | 39.58* | 10.12 | 34.47* | 41.91* |
| `M_pi_c` | 14.60* | 16.57 | 18.34* | 38.53* | 14.60* | 4.79* | 37.14* |
| `M_3pg_c` | 17.33 | 13.28 | 17.21* | 37.59* | 24.24* | 5.28 | 36.46 |
| `mRNA_JCVISYN3A_0607` | 32.01* | 16.02* | 12.13* | 35.35* | 8.87* | 4.54 | 30.67* |
| `M_gtp_c` | 10.62 | 11.97* | 11.97* | 34.66* | 7.87* | 4.03* | 25.96* |
| `mRNA_JCVISYN3A_0779` | 10.46 | 9.61 | 34.17* | 20.36* | 6.61* | 6.51 | 21.74* |
| `M_dhap_c` | 5.70* | 4.70* | 4.72* | 15.97* | 4.53 | 33.86* | 13.56* |
| `mRNA_JCVISYN3A_0213` | 9.45* | 8.47 | 11.09 | 25.71* | 5.09 | 4.99 | 31.35* |
| `M_2pg_c` | 13.89 | 10.94 | 11.96* | 29.05* | 31.16* | 3.25* | 30.10* |
| `mRNA_JCVISYN3A_0475` | 10.84* | 16.79* | 9.10* | 25.58* | 8.61* | 7.51* | 29.31* |
| `mRNA_JCVISYN3A_0131` | 9.47 | 9.35* | 10.73 | 23.66 | 6.13 | 3.93 | 28.41* |
| `M_ppi_c` | 25.52* | 27.01* | 20.42* | 24.75 | 15.00* | 8.38 | 27.33 |
| `mRNA_JCVISYN3A_0234` | 9.35* | 11.66* | 10.73* | 27.32* | 9.77* | 7.51* | 17.80* |
| `mRNA_JCVISYN3A_0221` | 14.09* | 13.06 | 12.41* | 22.78* | 8.65 | 5.18* | 26.36* |
| `M_atp_c` | 11.73* | 6.97* | 9.78 | 25.43* | 4.15* | 2.76 | 25.46* |
| `mRNA_JCVISYN3A_0651` | 10.10 | 12.05 | 14.32* | 17.86* | 8.38* | 7.45* | 23.94* |
| `mRNA_JCVISYN3A_0220` | 10.88* | 8.98* | 9.50* | 22.07* | 9.66 | 3.22* | 23.22 |
| `mRNA_JCVISYN3A_0344` | 13.90* | 15.27* | 11.18* | 22.87* | 8.60* | 10.12 | 21.21* |
| `mRNA_JCVISYN3A_0606` | 10.53* | 11.52* | 13.14* | 22.83 | 14.74 | 5.01 | 16.07* |
| `mRNA_JCVISYN3A_0445` | 10.61 | 19.41* | 11.75 | 20.28* | 7.05* | 4.02 | 21.81 |
| `M_fdp_c` | 3.59* | 1.08* | 1.79 | 3.06* | 0.51 | 0.92* | 21.55* |
| `mRNA_JCVISYN3A_0233` | 10.28* | 15.28 | 9.46* | 20.92 | 13.64* | 4.20* | 20.37* |
| `mRNA_JCVISYN3A_0694` | 14.17* | 11.52 | 17.47* | 20.37 | 8.22 | 5.69* | 18.03* |
| `mRNA_JCVISYN3A_0203` | 12.40 | 9.38 | 11.56* | 20.21* | 11.32* | 7.46 | 17.68* |
| `mRNA_JCVISYN3A_0727` | 13.80* | 13.86 | 12.83* | 18.56* | 6.19 | 5.00 | 18.05* |
| `mRNA_JCVISYN3A_0729` | 10.78* | 9.98* | 13.57* | 18.09 | 7.39* | 6.32* | 14.43 |
| `M_amp_c` | 13.16* | 6.80* | 7.65 | 17.43* | 8.62* | 3.06 | 15.52* |
| `M_pep_c` | 7.42 | 4.37 | 5.23* | 15.47* | 3.42* | 1.56 | 15.96* |
| `M_gmp_c` | 9.46* | 6.21* | 7.38 | 15.58* | 4.71* | 3.31* | 12.41* |
| `M_lac__L_e` | 2.03* | 1.33* | 2.92 | 11.75* | 0.26 | 0.15 | 15.41* |
| `M_gdp_c` | 5.49 | 5.05 | 6.88* | 11.01* | 4.38 | 3.14* | 14.38* |
| `M_adp_c` | 3.59* | 2.39 | 3.63 | 9.58* | 2.23* | 0.81* | 9.04* |
| `volume_litres` | 0.40* | 0.49 | 3.97* | 1.75* | 0.18 | 0.10 | 2.59* |
| `M_pyr_c` | 0.62 | 0.41* | 0.53* | 1.29* | 0.84* | 0.17 | 1.31* |
| `M_nadh_c` | 0.30 | 0.23 | 0.27* | 0.61* | 0.16* | 0.08 | 0.63* |
| `M_13dpg_c` | 0.10 | 0.06* | 0.08* | 0.20* | 0.05* | 0.02 | 0.21 |
| `M_nad_c` | 0.00 | 0.00 | 0.00* | 0.01* | 0.00* | 0.00 | 0.01* |

## The metabolite panel (15.5 into 15.7)

Kept: a pool check 1b does not exclude (`M_ppi_c`, `M_13dpg_c`, `M_nadh_c`, `M_pep_c`) with some target moving it by at least one resolution unit per unit ln θ, beyond three times its noise.

| metabolite | F2 rank | largest entry | kept |
|---|---|---|---|
| `M_g6p_c` | 11 | 59.21 | yes |
| `M_f6p_c` | 10 | 58.92 | yes |
| `M_lac__L_c` | 15 | 49.36 | yes |
| `M_trna_c` | 13 | 46.66 | yes |
| `M_trna_chg_c` | 16 | 46.66 | yes |
| `M_g3p_c` | 17 | 41.91 | yes |
| `M_pi_c` | 5 | 38.53 | yes |
| `M_3pg_c` | 12 | 37.59 | yes |
| `M_gtp_c` | 6 | 34.66 | yes |
| `M_dhap_c` | 18 | 33.86 | yes |
| `M_2pg_c` | 9 | 31.16 | yes |
| `M_ppi_c` | 2 | 27.33 | excluded (check 1b) |
| `M_atp_c` | 4 | 25.46 | yes |
| `M_fdp_c` | 21 | 21.55 | yes |
| `M_amp_c` | 1 | 17.43 | yes |
| `M_pep_c` | 14 | 15.96 | excluded (check 1b) |
| `M_gmp_c` | 7 | 15.58 | yes |
| `M_lac__L_e` | – | 15.41 | yes |
| `M_gdp_c` | 19 | 14.38 | yes |
| `M_adp_c` | 20 | 9.58 | yes |
| `M_pyr_c` | 22 | 1.31 | yes |
| `M_nadh_c` | 8 | 0.63 | excluded (check 1b) |
| `M_13dpg_c` | 3 | 0.21 | excluded (check 1b) |
| `M_nad_c` | 23 | 0.01 | no |

Panel (19): `M_g6p_c`, `M_f6p_c`, `M_lac__L_c`, `M_trna_c`, `M_trna_chg_c`, `M_g3p_c`, `M_pi_c`, `M_3pg_c`, `M_gtp_c`, `M_dhap_c`, `M_2pg_c`, `M_atp_c`, `M_fdp_c`, `M_amp_c`, `M_gmp_c`, `M_lac__L_e`, `M_gdp_c`, `M_adp_c`, `M_pyr_c`.
