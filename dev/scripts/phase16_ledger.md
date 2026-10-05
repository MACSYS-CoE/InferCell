# Phase 16 compute ledger

The cap is 50k CPU-h for the phase, and 5k for 16a (`spec/phases/16-recovery.md`
D16.6). CPU-h is wall time × cores, as reported by Slurm.

| Job | Task | What | Cores | Wall | CPU-h |
|---|---|---|---|---|---|
| 17754678 | 16a.1 | V2/V2b tests, first run | 1 | 0:03:39 | 0.06 |
| 17754679 | 16a.1 | Full-cycle V2, 3 cells, first run | 1 | 0:08:16 | 0.14 |
| 17755845 | 16a.1 | V2/V2b tests (of record) | 1 | 0:03:40 | 0.06 |
| 17755846 | 16a.1 | Full-cycle V2 and V2b (of record) | 1 | 0:08:25 | 0.14 |
| 17756205 | 16a.2 | 16a tests, identity failed (toy) | 1 | 0:04:21 | 0.07 |
| 17756874 | 16a.2 | 16a tests (of record) | 1 | 0:04:03 | 0.07 |

| 17758130 | 16a.3 | 16a tests, name clash | 1 | 0:04:07 | 0.07 |
| 17758372 | 16a.3/4 | 16a tests, cap bound failed | 1 | 0:04:29 | 0.07 |
| 17758553 | 16a.3/4 | 16a tests (of record) | 1 | 0:05:27 | 0.09 |
| 17800819 | 16a.8 | 16a tests with M0 (of record) | 1 | 0:05:28 | 0.09 |
| 17801360 | 16a.5 | 16a tests with block 1 | 1 | 0:06:43 | 0.11 |
| 17801198 | 16a.2 | V1 on M0, smoke (2 tasks) | 1 | ~0:04 each | 0.13 |
| 17801682 | 16a.2 | V1 on M0, 10⁴ paths (100 tasks) | 1 | ~0:28 each | 34.74 |
| 17801857 | 16a.6 | 16a tests with 16a.6 | 1 | 0:05:34 | 0.09 |
| 17801858, 17802079, 17802347 | 16a.6 | Core A′ clip scans (ENO) | 1 | | 0.29 |
| 17801684–17801688, 17802939–17802941, 17814376 | 16a.5 | V4 grids 1 to 3, merges, 20 chains | 1 | | 29.77 |
| 17814374, 17814375, 17830703 | 16a.6 | Core A′ clip scans (FBA, PGK3) | 1 | | 0.65 |
| 17814555, 17814556, 17814572, 17830702, 17830938 | 16a.5 | V4 grid 4, merge, finals, target check | 1 | | 4.77 |
| 17839207–17839209 | 16a.5 | V4 FBA-fine grid, merge, final (of record) | 1 | | 10.00 |
| 17839362, 17839571, 17840103, 17854107 | 16a.7 | 16a tests with CSMC (corrected after the #81 review: the two earlier rows overlapped) | 1 | | 0.41 |
| v5s_* (two smoke rounds) | 16a.7 | V5 smoke, both cases | 1 | | 0.67 |
| 17854476–17854480 | 16a.7 | V5, M0 (of record) | 1 | | 41.52 |
| 17854481–17854485 | 16a.7 | V5, toy (of record) | 1 | | 5.30 |
| 17877524–17877529 | 16a.7 | Variant smoke, N = 5, three windows | 1 | | 0.42 |
| 17877760, 17877761, 17877922, 17877923 | 16a.7 | Weight-collapse diagnostic (scratch) | 1 | | 0.26 |
| 17878615, 17878616, 17879029, 17879030 | 16a.7 | Collapse attribution, two passes (of record) | 1 | | 0.42 |
| 17879419, 17879420, 17880205, 17880206 | 16a.7 | Rank check, one failed run and the run of record | 1 | | 0.82 |
| 17878617–17878622 | 16a.7 | KP measurement at N = 50 (cell 4 failed under ancestor sampling) | 1 | | 11.00 |
| 17882299 | 16a.7 | 16a tests with the impossible-firing fix (of record) | 1 | 0:06:16 | 0.10 |
| 17882309–17882311 | 16a.7 | KP, cell 4 under ancestor sampling, rerun | 1 | | 2.30 |
| 17884470 | 16a.7 | KP merge (of record) | 1 | 0:01:22 | 0.02 |

| 17904315, 17904316 | 16a.7b | Rank check, five full-scale cells, and merge (of record) | 1 | | 2.56 |
| 17904432 | 16a.7a | 16a tests with annealing (of record) | 1 | 0:13:27 | 0.22 |
| 17905020–17905026 | 16a.7a | Annealed windows on M0, K ∈ {20, 50, 100}, N ∈ {2, 5} (of record) | 1 | | 22.01 |
| 17905037–17905039 | 16a.7a | V5 with annealing at K = 5, M0 (of record) | 1 | | 109.17 |
| 17905042–17905044 | 16a.7a | V5 with annealing at K = 5, toy (of record) | 1 | | 12.54 |
| 17905489 | 16a.7a | Inside one annealed window, K up to 1,000 (of record) | 1 | 1:20:57 | 1.35 |

| 17933661, 17933662 | 16a.7a | Collapse attribution at the cited detection floor (of record) | 1 | | 0.28 |

| 17934962 | 16a.7 | Bulk identifiability from 15.8's run (of record) | 1 | 0:03:31 | 0.06 |

| 17938700–17938702 | 16a.7 | F6 ENO cell: full-path information bound, 40 cells (of record) | 1 | | 5.07 |

| 18049605 | 16a.7c | 16a tests with bulk metabolites (of record) | 1 | 0:12:28 | 0.21 |
| 18052532–18052536 | 16a.7c | Bulk update rates on M0, N ∈ {5, 10, 20, 50} (of record) | 1 | | 35.36 |
| 18052537–18052542 | 16a.7c | V5 bulk half, M0 and toy (of record) | 1 | | (in the line above) |

| 18060223 | 15.7b | Bulk observations for 15.7, two cells resimulated (of record) | 1 | 0:03:22 | 0.06 |

| 18060205–18060213, 18061130–18061135 | 16a.7c | Fresh V5 references, M0 and toy, and the reference checks (of record) | 1 | | 25.01 |

| 18060214–18060217, 18061281, 18061355–18061363 | 16a.7c | V4 on block 1's bulk form: grids 1 to 3, FBA- and ENO-fine grids, chains, check, final (of record) | 1 | | 45.01 |

| 18098390 | 16a.7c | Per-cell update-rate merge (of record) | 1 | | 0.08 |

**16a total: about 403.6 CPU-h** (the rows' sum, recomputed after the #81 review,
which found the running total had drifted).
