# Output F2: concentration control coefficients (spec §11 task 15.5)

Commit 129e771, job 17681807, 2026-09-29T12:12:17.126, Julia 1.10.5. Regenerate with `sbatch dev/scripts/corea_f2_15.slurm`.

The frozen-expression variant's steady state, 23 metabolite pools against 13 genes, in 57.6 s. Check 8 guards it: the worst summation deviation is 3.61e-14 (tolerance 1e-6), and every coefficient matches re-solved steady states to 2.50e-09 (concentration) and 1.51e-09 (flux), against 1e-5.

Each entry is ∂ln x/∂ln E for a gene's enzyme, summing its reactions (PGK and PYK each two). The four phosphotransferase genes set carrier totals, which multiply no rate, and are not columns. The carrier states are protein counts and are not rows (spec §4 D8).

| metabolite | steady (mM) | PGI | PFK | FBA | TPI | GAPD | PGK | PGM | ENO | PYK | LDH_L | ADK1 | GK1 | PPA |
|---|---|---|---|---|---|---|---|---|---|---|---|---|---|---|
| `M_amp_c` | 0.9524 | 0 | 0 | 0 | 0 | 0 | -0.001 | 0.001 | 0.001 | -1.381 | 0 | -0.008 | 0 | 0 |
| `M_ppi_c` | 6.474e-05 | 0 | 0 | -0.006 | -0.003 | -1.339 | 0 | 0 | 0 | 1.249 | 0 | 0.002 | 0 | -0.969 |
| `M_13dpg_c` | 3.805e-05 | 0 | 0 | 0 | 0 | 0 | -1.012 | -0.014 | -0.010 | 1.183 | 0 | 0.001 | 0 | 0 |
| `M_atp_c` | 1.356 | 0 | 0 | 0 | 0 | 0 | 0.001 | -0.001 | -0.001 | 1.127 | 0 | 0.003 | 0 | 0 |
| `M_pi_c` | 0.1736 | 0 | 0 | -0.004 | -0.002 | -1.026 | 0 | 0 | 0 | 0.620 | 0 | 0.001 | 0 | 0 |
| `M_gtp_c` | 0.4829 | 0 | 0 | 0 | 0 | 0 | -0.002 | 0.001 | 0.001 | 0.977 | 0 | 0 | 0 | 0 |
| `M_gmp_c` | 0.7147 | 0 | 0 | 0 | 0 | 0 | 0 | 0 | 0 | -0.970 | 0 | 0 | 0 | 0 |
| `M_nadh_c` | 0.0001338 | 0 | 0 | 0 | 0 | 0 | 0 | -0.001 | -0.001 | 0.605 | -0.949 | 0.001 | 0 | 0 |
| `M_2pg_c` | 0.005722 | 0 | 0 | 0 | 0 | 0 | 0 | 0.001 | -0.837 | 0.618 | 0 | 0.001 | 0 | 0 |
| `M_f6p_c` | 0.1913 | 0.002 | -0.787 | -0.009 | -0.003 | 0.004 | -0.001 | 0.001 | 0.001 | -0.814 | 0 | 0 | 0 | 0 |
| `M_g6p_c` | 1.127 | -0.059 | -0.744 | -0.008 | -0.003 | 0.004 | -0.001 | 0.001 | 0.001 | -0.733 | 0 | 0 | 0 | 0 |
| `M_3pg_c` | 0.03147 | 0 | 0 | 0 | 0 | 0 | 0 | -0.606 | -0.421 | 0.676 | 0 | 0.001 | 0 | 0 |
| `M_trna_c` | 0.1006 | 0 | 0 | 0 | 0 | 0 | 0 | 0 | 0 | -0.674 | 0 | -0.002 | 0 | 0 |
| `M_pep_c` | 0.003152 | 0 | 0 | 0 | 0 | 0 | 0 | 0.004 | 0.004 | 0.607 | 0 | 0.001 | 0 | 0 |
| `M_lac__L_c` | 0.7582 | 0 | 0 | 0 | 0 | 0 | 0 | 0 | 0 | 0.600 | 0 | 0.001 | 0 | 0 |
| `M_trna_chg_c` | 0.1494 | 0 | 0 | 0 | 0 | 0 | 0 | 0 | 0 | 0.454 | 0 | 0.001 | 0 | 0 |
| `M_g3p_c` | 1.197 | 0.001 | 0.012 | 0.305 | 0.136 | 0.002 | 0 | 0 | 0 | -0.303 | 0 | -0.001 | 0 | 0 |
| `M_dhap_c` | 0.08776 | 0.001 | 0.012 | 0.298 | -0.145 | 0.002 | 0 | 0 | 0 | -0.130 | 0 | 0 | 0 | 0 |
| `M_gdp_c` | 0.7748 | 0 | 0 | 0 | 0 | 0 | 0.001 | -0.001 | 0 | 0.286 | 0 | 0 | 0 | 0 |
| `M_adp_c` | 1.646 | 0 | 0 | 0 | 0 | 0 | 0 | 0 | 0 | -0.129 | 0 | 0.003 | 0 | 0 |
| `M_fdp_c` | 21.02 | 0.002 | 0.023 | -0.009 | -0.003 | 0.004 | 0 | 0 | 0 | -0.066 | 0 | 0 | 0 | 0 |
| `M_pyr_c` | 4.48 | 0 | 0 | 0 | 0 | 0 | 0 | 0.004 | 0.004 | -0.006 | 0 | 0 | 0 | 0 |
| `M_nad_c` | 2.21 | 0 | 0 | 0 | 0 | 0 | 0 | 0 | 0 | 0 | 0 | 0 | 0 | 0 |

## The metabolite panel

Ranked by the row's largest |coefficient| over the genes. The panel is every pool at or above 0.1, where a 10% change in some enzyme moves it by at least 1%.

| rank | metabolite | largest \|C\| | gene | in panel |
|---|---|---|---|---|
| 1 | `M_amp_c` | 1.381 | PYK | yes |
| 2 | `M_ppi_c` | 1.339 | GAPD | yes |
| 3 | `M_13dpg_c` | 1.183 | PYK | yes |
| 4 | `M_atp_c` | 1.127 | PYK | yes |
| 5 | `M_pi_c` | 1.026 | GAPD | yes |
| 6 | `M_gtp_c` | 0.977 | PYK | yes |
| 7 | `M_gmp_c` | 0.970 | PYK | yes |
| 8 | `M_nadh_c` | 0.949 | LDH_L | yes |
| 9 | `M_2pg_c` | 0.837 | ENO | yes |
| 10 | `M_f6p_c` | 0.814 | PYK | yes |
| 11 | `M_g6p_c` | 0.744 | PFK | yes |
| 12 | `M_3pg_c` | 0.676 | PYK | yes |
| 13 | `M_trna_c` | 0.674 | PYK | yes |
| 14 | `M_pep_c` | 0.607 | PYK | yes |
| 15 | `M_lac__L_c` | 0.600 | PYK | yes |
| 16 | `M_trna_chg_c` | 0.454 | PYK | yes |
| 17 | `M_g3p_c` | 0.305 | FBA | yes |
| 18 | `M_dhap_c` | 0.298 | FBA | yes |
| 19 | `M_gdp_c` | 0.286 | PYK | yes |
| 20 | `M_adp_c` | 0.129 | PYK | yes |
| 21 | `M_fdp_c` | 0.066 | PYK | no |
| 22 | `M_pyr_c` | 0.006 | PYK | no |
| 23 | `M_nad_c` | 0.000 | LDH_L | no |

Panel (20): `M_amp_c`, `M_ppi_c`, `M_13dpg_c`, `M_atp_c`, `M_pi_c`, `M_gtp_c`, `M_gmp_c`, `M_nadh_c`, `M_2pg_c`, `M_f6p_c`, `M_g6p_c`, `M_3pg_c`, `M_trna_c`, `M_pep_c`, `M_lac__L_c`, `M_trna_chg_c`, `M_g3p_c`, `M_dhap_c`, `M_gdp_c`, `M_adp_c`.
