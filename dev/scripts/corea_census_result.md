# Check 7: the clipping census, scored for K5

Merged by `dev/scripts/corea_census_merge.jl` from 20 task files. Regenerate with `sbatch dev/scripts/corea_census.slurm`, then this script. The tasks' headers:

- commit 7e3c1c9, job 17549123, 2026-09-26T11:46:58.625, Julia 1.10.5
- commit 7e3c1c9, job 17549124, 2026-09-26T11:46:57.727, Julia 1.10.5
- commit 7e3c1c9, job 17549134, 2026-09-26T11:47:07.984, Julia 1.10.5
- commit 7e3c1c9, job 17549135, 2026-09-26T11:47:07.989, Julia 1.10.5
- commit 7e3c1c9, job 17549136, 2026-09-26T11:47:08.146, Julia 1.10.5
- commit 7e3c1c9, job 17549137, 2026-09-26T11:47:08.347, Julia 1.10.5
- commit 7e3c1c9, job 17549138, 2026-09-26T11:47:08.582, Julia 1.10.5
- commit 7e3c1c9, job 17549139, 2026-09-26T11:47:08.524, Julia 1.10.5
- commit 7e3c1c9, job 17549140, 2026-09-26T11:47:12.275, Julia 1.10.5
- commit 7e3c1c9, job 17549141, 2026-09-26T11:47:08.424, Julia 1.10.5
- commit 7e3c1c9, job 17549142, 2026-09-26T11:47:08.596, Julia 1.10.5
- commit 7e3c1c9, job 17549122, 2026-09-26T11:47:09.124, Julia 1.10.5
- commit 7e3c1c9, job 17549125, 2026-09-26T11:46:59.296, Julia 1.10.5
- commit 7e3c1c9, job 17549126, 2026-09-26T11:46:59.181, Julia 1.10.5
- commit 7e3c1c9, job 17549127, 2026-09-26T11:47:09.213, Julia 1.10.5
- commit 7e3c1c9, job 17549129, 2026-09-26T11:47:09.213, Julia 1.10.5
- commit 7e3c1c9, job 17549130, 2026-09-26T11:47:08.056, Julia 1.10.5
- commit 7e3c1c9, job 17549131, 2026-09-26T11:47:08.554, Julia 1.10.5
- commit 7e3c1c9, job 17549132, 2026-09-26T11:47:08.407, Julia 1.10.5
- commit 7e3c1c9, job 17549133, 2026-09-26T11:47:08.248, Julia 1.10.5

A drain clips when a consumer counter's pool pays less than was accrued, so a deficit is carried. The fraction is per drain, at the published 1 s drain.

## At published parameters

| seed | status | drains | drains clipped | fraction | counters that clipped (drains, first s) |
|---|---|---|---|---|---|
| 1310 | ok | 6300 | 253 | 0.04016 | `GTP_translat` (253, 953), `tRNA_translat` (20, 1603) |
| 1410 | ok | 6300 | 23 | 0.003651 | `GTP_translat` (23, 4914) |
| 14800 | ok | 6300 | 34 | 0.005397 | `GTP_translat` (33, 1580), `tRNA_translat` (4, 1927) |
| 14801 | ok | 6300 | 35 | 0.005556 | `GTP_translat` (34, 298), `tRNA_translat` (5, 298) |
| 14802 | ok | 6300 | 94 | 0.01492 | `GTP_translat` (94, 496) |
| 14803 | ok | 6300 | 21 | 0.003333 | `GTP_translat` (21, 1007), `tRNA_translat` (3, 3286) |
| 14804 | ok | 6300 | 61 | 0.009683 | `GTP_translat` (60, 1720), `tRNA_translat` (6, 2057) |
| 14805 | ok | 6300 | 32 | 0.005079 | `GTP_translat` (32, 2049), `tRNA_translat` (4, 3676) |
| 14806 | ok | 6300 | 133 | 0.02111 | `GTP_translat` (133, 1636) |
| 14807 | ok | 6300 | 253 | 0.04016 | `GTP_translat` (253, 1047), `tRNA_translat` (3, 3286) |

Seeds carrying a deficit: **10 of 10**. K5's threshold here is zero.

## Draw 0: freeing the constants changes nothing

Draw 0 at seed 14800: 34 drains clipped, against 34 for the published build at the same seed. Per-counter records **identical**.

## Across prior draws

200 draws; 200 completed, 0 failed.
Draws carrying any deficit: **146 of 200, 73.0%** (Wilson 95% interval 66.5–78.7%). K5's threshold is 5%.
Per-drain clipping fraction over completed draws: median 0.653, max 0.995.

| counter | draws in which it clipped | median drains clipped when it did | earliest first clip (s) |
|---|---|---|---|
| `ATP_mRNA` | 42 | 290 | 40 |
| `ATP_mRNAdeg` | 66 | 430 | 40 |
| `ATP_transloc` | 24 | 6 | 81 |
| `ATP_trsc` | 40 | 190 | 90 |
| `GTP_mRNA` | 29 | 3 | 65 |
| `GTP_translat` | 102 | 5824 | 35 |
| `tRNA_translat` | 93 | 3095 | 22 |

## K5 on the scoreboard (T3)

| criterion | threshold | measured | verdict |
|---|---|---|---|
| K5, published parameters | 0 seeds carrying a deficit | 10 of 10 | **fires** |
| K5, prior draws | ≤ 5% of draws clip | 73.0% of 200 | **fires** |

K5 **fires**. What follows from that — smoothing the drain or resizing a pool — is a separate decision (spec §12, 2026-09-25 C).
