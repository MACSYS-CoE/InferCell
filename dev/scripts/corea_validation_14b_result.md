# Phase 14b: the derivative, ensemble and census checks on the assembled Core A′

Assembled by `dev/scripts/corea_validation_14b.jl` (SECTION=report) from the section files below, each with its own provenance line. Check 7's census is in `corea_census_result.md`. Regenerate with `dev/scripts/corea_validation_14b.slurm`.

## Check 1b: the particle floor over the pinned cycle

Section `floor` at commit ebe9e9d, Slurm job 17559025, on dave14, 2026-09-26T18:42:02.231, Julia 1.10.5.

Seed 1410, one cycle of 6300 handshakes (116 s). Every ODE state, in particles at the volume it was held at, smallest minimum first. Flagged: minimum below 500.

| state | min | at (s) | median | time under 1 particle | flagged | check 1b |
|---|---|---|---|---|---|---|
| M_lac__L_e | 0 | 0 | 55.44 | 0.2% | yes | external, not assessed |
| M_gtp_c | 0 | 4914 | 5617 | 0.4% | yes | flagged |
| M_gdp_c | 0 | 20 | 1.96e+04 | 0.0% | yes | flagged |
| M_gmp_c | 0 | 3 | 1.425e+04 | 0.3% | yes | flagged |
| M_ppi_c | 0 | 47 | 1 | 41.9% | yes | **excluded**: median under 10 |
| M_ptsi_P_c | 0.0207 | 5786 | 0.3826 | 99.0% | yes | **excluded**: median under 10 |
| M_13dpg_c | 0.0309 | 5786 | 0.691 | 73.0% | yes | **excluded**: median under 10 |
| M_nadh_c | 0.176 | 5786 | 2.539 | 7.6% | yes | **excluded**: median under 10 |
| M_ptsh_P_c | 0.187 | 5786 | 3.439 | 4.0% | yes | **excluded**: median under 10 |
| M_crr_P_c | 0.217 | 5786 | 4.331 | 2.5% | yes | **excluded**: median under 10 |
| M_ptsg_P_c | 1.2 | 5786 | 20.88 | 0.0% | yes | Langevin cross-check |
| M_pep_c | 2.07 | 5786 | 54.61 | 0.0% | yes | Langevin cross-check |
| M_2pg_c | 5.34 | 5786 | 111.2 | 0.0% | yes | Langevin cross-check |
| M_trna_c | 7 | 21 | 1494 | 0.0% | yes | flagged |
| M_ptsh_c | 14.5 | 0 | 347 | 0.0% | yes | flagged |
| M_crr_c | 15.7 | 0 | 415 | 0.0% | yes | flagged |
| M_ptsi_c | 17.7 | 0 | 382 | 0.0% | yes | flagged |
| M_amp_c | 26 | 16 | 1.494e+04 | 0.0% | yes | flagged |
| M_3pg_c | 30.4 | 5786 | 631.6 | 0.0% | yes | flagged |
| M_ptsg_c | 125 | 0 | 1132 | 0.0% | yes | flagged |
| M_pi_c | 157 | 5785 | 3524 | 0.0% | yes | flagged |
| M_trna_chg_c | 162 | 2316 | 3553 | 0.0% | yes | flagged |
| M_g3p_c | 846 | 1 | 2.676e+04 | 0.0% | no |  |
| M_dhap_c | 1.11e+03 | 6 | 1986 | 0.0% | no |  |
| M_f6p_c | 1.82e+03 | 1475 | 3632 | 0.0% | no |  |
| M_lac__L_c | 2.02e+03 | 0 | 1.969e+04 | 0.0% | no |  |
| M_adp_c | 2.08e+03 | 18 | 3.25e+04 | 0.0% | no |  |
| M_g6p_c | 1.07e+04 | 1475 | 2.16e+04 | 0.0% | no |  |
| M_atp_c | 2.18e+04 | 5635 | 3.36e+04 | 0.0% | no |  |
| M_pyr_c | 3.01e+04 | 9 | 9.047e+04 | 0.0% | no |  |
| M_nad_c | 4.41e+04 | 0 | 4.461e+04 | 0.0% | no |  |
| M_fdp_c | 1.54e+05 | 0 | 4.866e+05 | 0.0% | no |  |

22 of 32 states fall below 500 particles. Excluded outright, with a median under 10 particles: `M_ppi_c`, `M_ptsi_P_c`, `M_13dpg_c`, `M_nadh_c`, `M_ptsh_P_c`, `M_crr_P_c`. Cross-checked against the Langevin ensemble: `M_ptsg_P_c`, `M_pep_c`, `M_2pg_c`.

The band mutation's trajectory, enolase's two catalytic constants halved at seed 1410, is recorded for `merge-band`.

## Check 1b: the chemical-Langevin cross-check

Section `merge-band` at commit ebe9e9d, Slurm job 17559027, on dave5, 2026-09-26T19:30:35.205, Julia 1.10.5.

The ensemble is the test-local double in `test/corea_langevin_doubles.jl`: the ODE block's reactions as a chemical Langevin equation, each reaction's forward and reverse channel separately, integrated by local linearization, with the jump path frozen at seed 1410's (spec §12, 2026-09-26 E). The band is the 5th to 95th percentile at every handshake.

**The integrator against the pinned Rodas5P over one interval**, at h = 0.01: the largest `|aₖ|` relative to the state on species no counter touches, which is this integrator's own error.

| species | h = 0.01 | h = 0.02 |
|---|---|---|
| M_13dpg_c | 0.276 | 0.536 |
| M_2pg_c | 0.263 | 0.48 |
| M_3pg_c | 0.278 | 0.53 |
| M_fdp_c | 7.65e-06 | 6.36e-05 |
| M_g3p_c | 0.00625 | 0.0179 |
| M_g6p_c | 0.00353 | 0.01 |
| M_lac__L_c | 0.0233 | 0.0607 |
| M_nadh_c | 0.251 | 0.47 |
| M_pep_c | 0.261 | 0.475 |
| M_pyr_c | 0.00106 | 0.0034 |

| pool | trajectories | reference inside the band | largest excursion outside it | at (s) | band width at the end (particles) |
|---|---|---|---|---|---|
| M_ptsg_P_c | 100 | 100.0% of handshakes | 0 | 1 | 5.6 |
| M_pep_c | 100 | 100.0% of handshakes | 0 | 1 | 10.7 |
| M_2pg_c | 100 | 100.0% of handshakes | 0 | 1 | 20.8 |

The excursion is the smallest relative observation noise below which phase 15 must exclude the pool (spec §3, amended 2026-09-26 B).

**The partition and the clamps**, per step size. A channel is noisy only while every species it moves holds at least 10 particles. A clamp creates mass, so each moiety's largest drift from the reference is shown beside what the clamps injected into it, both in particles at the end of the cycle, worst trajectory.

| h | trajectories | channel-steps noiseless | steps clamped | handshakes clamped by aₖ | redox drift / booked | adenylate drift / booked | guanylate drift / booked | phosphate drift / booked | carrier_ptsI drift / booked | carrier_ptsH drift / booked | carrier_Crr drift / booked | carrier_ptsG drift / booked | trna drift / booked |
|---|---|---|---|---|---|---|---|---|---|---|---|---|---|
| 0.01 | 100 | 45.94% | 6568 | 118084 | 59.5 / 69.3 | 0.00332 / 0 | 676 / 679 | 2.82e+03 / 1.13e+04 | 14.9 / 14.9 | 1.83 / 1.83 | 0.974 / 0.974 | 0.000201 / 0 | 0.00157 / 0 |
| 0.02 | 20 | 45.92% | 3897 | 23390 | 58.9 / 162 | 0.00221 / 0 | 332 / 333 | 2.04e+03 / 1e+05 | 2.79 / 2.79 | 0.123 / 0.124 | 0.000116 / 0 | 0.00015 / 0 | 0.00143 / 0 |

**The halving comparison.** Excursion at h = 0.02 (20 trajectories) against h = 0.01: `M_ptsg_P_c` 0.00699 against 0; `M_pep_c` 0 against 0; `M_2pg_c` 0 against 0.

**The mutation.** Enolase's catalytic constants halved, same seed, against the nominal band: `M_ptsg_P_c` 7.17, `M_pep_c` 7.14, `M_2pg_c` 4.55. The pool it moves most is `M_ptsg_P_c`.

## Check 8: the summation theorems on the frozen-expression variant

Section `check8` at commit ebe9e9d, Slurm job 17559265, on dave5, 2026-09-26T19:35:18.493, Julia 1.10.5.

The variant: the four ODE modules at their nominal enzyme counts, with `FrozenDemand` consuming charged tRNA and GTP at the published demand at the nominal pools. External lactate, the demand's accumulator and the held protein counts are held. 24 multipliers, one per reaction (spec §3, amended 2026-09-26).

Steady state in 19.4 s: max |f| 2.43e-16 mM/s against a largest reaction rate of 0.0569 mM/s, after 1 Newton iteration(s). Eigenvalues of the Jacobian within the conserved class: -2.09e+04 to -0.00597 /s, all negative.

**Conserved combinations: 10.** The eight of 14a's moieties the variant keeps (redox, adenylate, guanylate, phosphate and the four carriers) and the tRNA pair each lie in their span. The tenth is not one 14a asserts. With glucose and lactate both at the boundary it is a closed combination over the glycolytic intermediates, NAD/NADH, the carriers and the nucleotides, found numerically.

**Summation identities.** Worst concentration sum 3.61e-14 (should be 0), worst flux sum 2.58e-14 from 1. Tolerance 1e-6.
Zero-flux reactions, left out of the flux sums: `R_GK1`.

| multiplier | flux (mM/s, on) | Σⱼ C^J − 1 |
|---|---|---|
| R_PGI | -0.02843 (`M_g6p_c`) | 2.1e-14 |
| R_PFK | -0.02843 (`M_f6p_c`) | 2.6e-14 |
| R_FBA | -0.02843 (`M_fdp_c`) | -4.4e-15 |
| R_TPI | -0.02843 (`M_dhap_c`) | -8.0e-15 |
| R_GAPD | -0.05686 (`M_g3p_c`) | 1.6e-15 |
| R_PGK | -0.04698 (`M_13dpg_c`) | 9.3e-15 |
| R_PGM | -0.05686 (`M_3pg_c`) | 1.8e-15 |
| R_ENO | -0.05686 (`M_2pg_c`) | 8.9e-15 |
| R_PYK | -0.02239 (`M_pep_c`) | -7.8e-15 |
| R_LDH_L | -0.05686 (`M_pyr_c`) | -6.7e-16 |
| R_PGK3 | -0.009883 (`M_13dpg_c`) | 8.9e-15 |
| R_PYK3 | -0.006038 (`M_pep_c`) | -7.8e-15 |
| R_ADK1 | 0.04094 (`M_adp_c`) | 6.7e-16 |
| R_GK1 | 2.359e-16 (`M_atp_c`) | zero flux |
| R_PPA | 0.04094 (`M_pi_c`) | 1.3e-15 |
| R_GLCpts0 | -0.02843 (`M_pep_c`) | -6.2e-15 |
| R_GLCpts1 | 0.02843 (`M_ptsi_c`) | 1.2e-14 |
| R_GLCpts2 | 0.02843 (`M_ptsh_c`) | 1.2e-14 |
| R_GLCpts3 | 0.02843 (`M_crr_c`) | 8.4e-15 |
| R_GLCpts4 | 0.02843 (`M_g6p_c`) | 1.2e-14 |
| R_LACt | -0.05686 (`M_lac__L_c`) | 4.2e-15 |
| R_charging | -0.02047 (`M_atp_c`) | 6.7e-16 |
| frozen_translation | 0.02047 (`M_trna_c`) | 1.1e-15 |
| frozen_gtp | 0.01592 (`M_pi_c`) | 2.7e-15 |

Gene-level concentration control coefficients of the smallest steady pools, summing a gene's reactions (PGK and PYK each two):

| gene | `M_ptsi_P_c` | `M_13dpg_c` | `M_ppi_c` | `M_nadh_c` | `M_ptsh_P_c` |
|---|---|---|---|---|---|
| PGI | -8.17e-08 | -5.72e-08 | -1.52e-05 | -2.92e-08 | -1.39e-07 |
| PFK | -1.03e-06 | -7.2e-07 | -0.000227 | -3.68e-07 | -1.75e-06 |
| FBA | -1.16e-08 | -8.12e-09 | -0.00577 | -4.15e-09 | -1.97e-08 |
| TPI | -4.43e-09 | -3.1e-09 | -0.00257 | -1.58e-09 | -7.54e-09 |
| GAPD | 5.22e-09 | 3.65e-09 | -1.34 | 1.87e-09 | 8.87e-09 |
| PGK | -0.000445 | -1.01 | -0.000291 | -0.000441 | -0.000444 |
| PGM | 5e-05 | -0.014 | -0.000188 | -0.000943 | 4.99e-05 |
| ENO | -3.89e-06 | -0.00977 | -0.000211 | -0.000943 | -3.88e-06 |
| PYK | 0.613 | 1.18 | 1.25 | 0.605 | 0.611 |
| LDH_L | 1.02e-06 | 1.92e-06 | -0.000107 | -0.949 | 1.02e-06 |
| ADK1 | 0.000927 | 0.000762 | 0.00228 | 0.000916 | 0.000925 |
| GK1 | -2.7e-17 | -2.67e-17 | -4.21e-17 | -2.67e-17 | -2.7e-17 |
| PPA | 3.67e-12 | 2.57e-12 | -0.969 | 1.31e-12 | 6.24e-12 |

**The mutation.** Without `k_chg`'s multiplier the worst concentration sum is 0.807 and the worst flux sum is 0.384 from 1. The check reports: ArgumentError: Check 8 fails: the concentration control coefficients of :M_ppi_c sum to -0.8074879841980597, not 0 (tolerance 1.0e-6)

## Check 6: the nominal trajectory, and task 9.9 on the assembled model

Section `check6` at commit ebe9e9d, Slurm job 17559266, on dave5, 2026-09-26T19:35:16.792, Julia 1.10.5.

Seed 1410, one full cycle. Fractional growth **1.0382** (volume over its initial value).
- time to 1.02×: 2131 s
- time to 1.04×: not reached
- time to 1.06×: not reached
- `doubling_time` is refused in code: yes, reporting instead `fractional_growth` or `time_to_threshold`.

### Task 9.9: the charging stoichiometry, matched on cycle-mean flux

Over 5 paired seeds, the published AMP + PPi lumping charges 318.42 residues/s on the cycle mean. The two-ATP lumping's cycle-mean flux against `k2_scale`:

| k2_scale | 0.25 | 0.5 | 1 | 2 | 4 | 8 | 16 |
|---|---|---|---|---|---|---|---|
| flux (/s) | 114.64 | 177.30 | 273.65 | 324.07 | 322.66 | 308.10 | 324.27 |
| relative to published | -64.00% | -44.32% | -14.06% | +1.78% | +1.33% | -3.24% | +1.84% |

Matched at `k2_scale` = 1.85: flux 321.39, +0.93% from published.

| lumping | charging flux (/s) | ATP/ADP, cycle mean | ATP/ADP at the end | AMP (mM) | PPi (mM) | ADK1 net flux (/s) |
|---|---|---|---|---|---|---|
| published, AMP + PPi | 318.42 | 0.8301 | 0.6610 | 0.9579 | 0.0029 | -290.9 |
| two ATP → two ADP | 321.39 | 0.9391 | 0.7789 | 0.8143 | 0.00216 | -509.2 |

ATP/ADP differs by **13.13%** between the lumpings on the cycle mean, with live glycolysis, at matched flux. Phase 9's double pinned ADP and measured 0.035%.
ADK1's forward direction is AMP + ATP → 2 ADP. A negative net flux is the kinase making AMP.

## F5: the two external comparisons (task 14b.5)

Section `f5` at commit ebe9e9d, Slurm job 17559267, on dave5, 2026-09-26T19:35:16.335, Julia 1.10.5.

Assembled model, 20 seeds (1410–1429), one cycle each. Transcripts are time-averaged per gene after 600 s and over seeds. Fold change is each gene's protein count at the end over its count after the first handshake, with ptsG's cytoplasmic pool added back (phase 11's definition), averaged over seeds.

| comparison | threshold (§3) | measured | verdict |
|---|---|---|---|
| transcripts: Spearman | ≥ 0.7 | 0.850 | pass |
| transcripts: within 2× | ≥ 15 of 17 | 16 of 17 | pass |
| fold change: median | [1.7, 2.3] | 1.614 | **miss** |
| fold change: range | none < 1.5 or > 3.0 | 1.489 to 1.755 (1 below, 0 above) | **miss** |
| fold change against length | negative slope | 0.030 | **miss** |

| gene | predicted mRNA | measured | ratio | fold change | length (nt) |
|---|---|---|---|---|---|
| JCVISYN3A_0445 | 0.463 | 0.440 | 1.05 | 1.620 | 1284 |
| JCVISYN3A_0220 | 0.663 | 0.741 | 0.89 | 1.601 | 981 |
| JCVISYN3A_0131 | 1.084 | 0.347 | 3.13 | 1.503 | 894 |
| JCVISYN3A_0727 | 0.616 | 0.661 | 0.93 | 1.590 | 747 |
| JCVISYN3A_0607 | 2.304 | 2.178 | 1.06 | 1.658 | 1017 |
| JCVISYN3A_0606 | 0.650 | 0.656 | 0.99 | 1.640 | 1215 |
| JCVISYN3A_0729 | 0.426 | 0.519 | 0.82 | 1.539 | 1596 |
| JCVISYN3A_0213 | 1.439 | 1.678 | 0.86 | 1.600 | 1356 |
| JCVISYN3A_0221 | 0.910 | 0.957 | 0.95 | 1.615 | 1437 |
| JCVISYN3A_0475 | 1.641 | 1.894 | 0.87 | 1.623 | 957 |
| JCVISYN3A_0233 | 0.671 | 0.555 | 1.21 | 1.755 | 1722 |
| JCVISYN3A_0234 | 0.494 | 0.518 | 0.95 | 1.646 | 465 |
| JCVISYN3A_0694 | 0.437 | 0.507 | 0.86 | 1.489 | 270 |
| JCVISYN3A_0779 | 1.476 | 1.287 | 1.15 | 1.614 | 2238 |
| JCVISYN3A_0651 | 0.379 | 0.346 | 1.09 | 1.647 | 642 |
| JCVISYN3A_0344 | 0.336 | 0.292 | 1.15 | 1.577 | 561 |
| JCVISYN3A_0203 | 0.299 | 0.317 | 0.94 | 1.555 | 894 |

Per-seed fold-change medians: 1.549, 1.568, 1.554, 1.435, 1.529, 1.761, 1.549, 1.523, 1.545, 1.682, 1.531, 1.570, 1.670, 1.599, 1.546, 1.537, 1.622, 1.653, 1.554, 1.581.

## F3's mutation table, continued: checks 1b, 7 and 8

Checks 0 to 5 are in `corea_validation_result.md` (phase 14a).

| check | mutation | what it produced | named quantity |
|---|---|---|---|
| 1b | enolase's catalytic constants halved | excursion 7.17 outside the nominal band, against 0 nominal | `M_ptsg_P_c` |
| 7 | tRNA pool at 0.05 mM, a fifth of the asserted | translation's counter clips (asserted in `test/test_corea_validation_14b.jl`, 600 handshakes) | `tRNA_translat` |
| 8 | `k_chg`'s multiplier omitted | worst concentration sum 0.807, worst flux sum 0.384 from 1, against 1e-6 | ArgumentError: Check 8 fails: the concentration control coefficients of :M_ppi_c sum to -0.8074879841980597, not 0 (tolerance 1.0e-6) |
