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

**16a total: about 65.9 CPU-h** (V4's grid 4, final, and the FBA and PGK3 scans are still running).
