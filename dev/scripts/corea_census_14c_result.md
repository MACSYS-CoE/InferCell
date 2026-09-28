# Phase 14c: check 7's census, rerun three ways

Merged by `dev/scripts/corea_census_14c_merge.jl`. Regenerate with `sbatch dev/scripts/corea_census_14c.slurm <mode>` for each of `d11`, `broad` and `gtp`, then this script. 14b's task files, which this merge also reads, are gitignored: in a fresh clone run `sbatch dev/scripts/corea_census.slurm` first. Every run is at the published clamped drain. The tasks' headers:

- commit 946bfaa, job 17623903, mode d11, 2026-09-28T16:17:50.492, Julia 1.10.5
- commit 946bfaa, job 17623904, mode d11, 2026-09-28T16:17:50.622, Julia 1.10.5
- commit 946bfaa, job 17623913, mode d11, 2026-09-28T16:17:50.674, Julia 1.10.5
- commit 946bfaa, job 17623914, mode d11, 2026-09-28T16:17:50.191, Julia 1.10.5
- commit 946bfaa, job 17623915, mode d11, 2026-09-28T16:17:50.228, Julia 1.10.5
- commit 946bfaa, job 17623916, mode d11, 2026-09-28T16:17:50.469, Julia 1.10.5
- commit 946bfaa, job 17623917, mode d11, 2026-09-28T16:17:50.558, Julia 1.10.5
- commit 946bfaa, job 17623918, mode d11, 2026-09-28T16:17:51.120, Julia 1.10.5
- commit 946bfaa, job 17623919, mode d11, 2026-09-28T16:17:51.130, Julia 1.10.5
- commit 946bfaa, job 17623920, mode d11, 2026-09-28T16:17:49.287, Julia 1.10.5
- commit 946bfaa, job 17623921, mode d11, 2026-09-28T16:17:49.788, Julia 1.10.5
- commit 946bfaa, job 17623900, mode d11, 2026-09-28T16:17:50.267, Julia 1.10.5
- commit 946bfaa, job 17623905, mode d11, 2026-09-28T16:17:50.352, Julia 1.10.5
- commit 946bfaa, job 17623906, mode d11, 2026-09-28T16:17:50.606, Julia 1.10.5
- commit 946bfaa, job 17623907, mode d11, 2026-09-28T16:17:51.402, Julia 1.10.5
- commit 946bfaa, job 17623908, mode d11, 2026-09-28T16:17:51.195, Julia 1.10.5
- commit 946bfaa, job 17623909, mode d11, 2026-09-28T16:17:50.842, Julia 1.10.5
- commit 946bfaa, job 17623910, mode d11, 2026-09-28T16:17:51.043, Julia 1.10.5
- commit 946bfaa, job 17623911, mode d11, 2026-09-28T16:17:50.972, Julia 1.10.5
- commit 946bfaa, job 17623912, mode d11, 2026-09-28T16:17:51.367, Julia 1.10.5
- commit 946bfaa, job 17623922, mode broad, 2026-09-28T16:17:50.680, Julia 1.10.5
- commit 946bfaa, job 17623923, mode broad, 2026-09-28T16:17:50.660, Julia 1.10.5
- commit 946bfaa, job 17623932, mode broad, 2026-09-28T16:18:03.619, Julia 1.10.5
- commit 946bfaa, job 17623933, mode broad, 2026-09-28T16:18:03.638, Julia 1.10.5
- commit 946bfaa, job 17623934, mode broad, 2026-09-28T16:18:03.568, Julia 1.10.5
- commit 946bfaa, job 17623935, mode broad, 2026-09-28T16:18:03.553, Julia 1.10.5
- commit 946bfaa, job 17623936, mode broad, 2026-09-28T16:18:04.140, Julia 1.10.5
- commit 946bfaa, job 17623937, mode broad, 2026-09-28T16:18:03.544, Julia 1.10.5
- commit 946bfaa, job 17623938, mode broad, 2026-09-28T16:18:04.377, Julia 1.10.5
- commit 946bfaa, job 17623939, mode broad, 2026-09-28T16:18:03.866, Julia 1.10.5
- commit 946bfaa, job 17623940, mode broad, 2026-09-28T16:18:03.039, Julia 1.10.5
- commit 946bfaa, job 17623901, mode broad, 2026-09-28T16:18:03.621, Julia 1.10.5
- commit 946bfaa, job 17623924, mode broad, 2026-09-28T16:17:50.816, Julia 1.10.5
- commit 946bfaa, job 17623925, mode broad, 2026-09-28T16:18:06.581, Julia 1.10.5
- commit 946bfaa, job 17623926, mode broad, 2026-09-28T16:18:06.581, Julia 1.10.5
- commit 946bfaa, job 17623927, mode broad, 2026-09-28T16:18:04.057, Julia 1.10.5
- commit 946bfaa, job 17623928, mode broad, 2026-09-28T16:18:03.334, Julia 1.10.5
- commit 946bfaa, job 17623929, mode broad, 2026-09-28T16:18:04.611, Julia 1.10.5
- commit 946bfaa, job 17623930, mode broad, 2026-09-28T16:18:03.391, Julia 1.10.5
- commit 946bfaa, job 17623931, mode broad, 2026-09-28T16:18:02.914, Julia 1.10.5
- commit 946bfaa, job 17623941, mode gtp, 2026-09-28T16:18:03.448, Julia 1.10.5
- commit 946bfaa, job 17623942, mode gtp, 2026-09-28T16:18:03.671, Julia 1.10.5
- commit 946bfaa, job 17623954, mode gtp, 2026-09-28T16:18:32.357, Julia 1.10.5
- commit 946bfaa, job 17623955, mode gtp, 2026-09-28T16:18:32.240, Julia 1.10.5
- commit 946bfaa, job 17623956, mode gtp, 2026-09-28T16:18:32.559, Julia 1.10.5
- commit 946bfaa, job 17623957, mode gtp, 2026-09-28T16:18:32.781, Julia 1.10.5
- commit 946bfaa, job 17623958, mode gtp, 2026-09-28T16:18:30.524, Julia 1.10.5
- commit 946bfaa, job 17623959, mode gtp, 2026-09-28T16:18:31.257, Julia 1.10.5
- commit 946bfaa, job 17623960, mode gtp, 2026-09-28T16:18:31.159, Julia 1.10.5
- commit 946bfaa, job 17623961, mode gtp, 2026-09-28T16:18:37.961, Julia 1.10.5
- commit 946bfaa, job 17623962, mode gtp, 2026-09-28T16:18:37.817, Julia 1.10.5
- commit 946bfaa, job 17623902, mode gtp, 2026-09-28T16:18:37.946, Julia 1.10.5
- commit 946bfaa, job 17623943, mode gtp, 2026-09-28T16:17:59.969, Julia 1.10.5
- commit 946bfaa, job 17623944, mode gtp, 2026-09-28T16:18:00.062, Julia 1.10.5
- commit 946bfaa, job 17623945, mode gtp, 2026-09-28T16:18:00.165, Julia 1.10.5
- commit 946bfaa, job 17623946, mode gtp, 2026-09-28T16:18:00.118, Julia 1.10.5
- commit 946bfaa, job 17623947, mode gtp, 2026-09-28T16:17:59.923, Julia 1.10.5
- commit 946bfaa, job 17623948, mode gtp, 2026-09-28T16:17:59.578, Julia 1.10.5
- commit 946bfaa, job 17623949, mode gtp, 2026-09-28T16:17:59.785, Julia 1.10.5
- commit 946bfaa, job 17623950, mode gtp, 2026-09-28T16:18:04.582, Julia 1.10.5

14b's census, set beside these, is `dev/scripts/corea_census_result.md` at `7e3c1c9`, recomputed here from its task files. A drain clips when a consumer counter's pool pays less than was accrued. A draw is starved when some counter clips at more than half its drains.

## 14c.2: K5's prior half, over D11's six targets

Each draw varies the three promoter strengths, `krnadeg`, and ENO's and FBA's forward constants from their priors, with the two reverse constants derived (task 14c.1). Everything else is at its published value.

Draw 0 (the freed build, no value written) at seed 14800: 34 drains clipped, against 34 for the published build at the same seed. Per-counter records **identical**.

200 draws; 200 completed, 0 failed.

| | draws clipping | Wilson 95% | starved draws | Wilson 95% | median per-drain fraction | starved draws' median | the others' median | max |
|---|---|---|---|---|---|---|---|---|
| 14c.2, D11's six, Haldane-consistent | 163 of 200, 81.5% | 75.5–86.3% | 8 of 200, 4.0% | 2.0–7.7% | 0.00611 | 0.622 | 0.00524 | 0.781 |
| 14b broad, independent draws | 146 of 200, 73.0% | 66.5–78.7% | 104 of 200, 52.0% | 45.1–58.8% | 0.653 | 0.969 | 0 | 0.995 |

The per-drain fraction is bimodal: a starved draw clips on most drains and the others on almost none. The median over all draws therefore falls in whichever mode holds more than half of them, and moves with the starved count rather than with how long a starved draw clips. Read the split columns.

| counter | draws in which it clipped | median drains clipped when it did | earliest first clip (s) |
|---|---|---|---|
| `GTP_translat` | 162 | 74 | 143 |
| `tRNA_translat` | 125 | 56 | 88 |

| criterion | threshold | measured | verdict |
|---|---|---|---|
| K5, prior draws over D11's targets (rescored, §12 2026-09-28 B) | ≤ 5% of draws clip | 81.5% of 200 (Wilson 75.5–86.3%) | **fires** |
| K5, prior draws over every informed constant (14b, kept as the stress result) | ≤ 5% of draws clip | 73.0% of 200 | **fires** |

## 14c.3: 14b's broad census, drawn Haldane-consistently

The same informed constants as 14b's census, at the same seed numbers, with every reverse constant derived rather than drawn. The draws are not paired with 14b's: a derived reverse constant consumes no random number, so the streams diverge after the first. Recorded, not gated: the difference from 14b is what the broken equilibrium constants did.

Draw 0 (the freed build, no value written) at seed 14800: 34 drains clipped, against 34 for the published build at the same seed. Per-counter records **identical**.

200 draws; 200 completed, 0 failed.

| | draws clipping | Wilson 95% | starved draws | Wilson 95% | median per-drain fraction | starved draws' median | the others' median | max |
|---|---|---|---|---|---|---|---|---|
| 14c.3, broad, Haldane-consistent | 141 of 200, 70.5% | 63.8–76.4% | 95 of 200, 47.5% | 40.7–54.4% | 0.271 | 0.957 | 0 | 0.994 |
| 14b broad, independent draws | 146 of 200, 73.0% | 66.5–78.7% | 104 of 200, 52.0% | 45.1–58.8% | 0.653 | 0.969 | 0 | 0.995 |

The per-drain fraction is bimodal: a starved draw clips on most drains and the others on almost none. The median over all draws therefore falls in whichever mode holds more than half of them, and moves with the starved count rather than with how long a starved draw clips. Read the split columns.

| counter | draws in which it clipped | median drains clipped when it did | earliest first clip (s) |
|---|---|---|---|
| `ATP_mRNA` | 39 | 214 | 147 |
| `ATP_mRNAdeg` | 66 | 73 | 148 |
| `ATP_transloc` | 18 | 6 | 487 |
| `ATP_trsc` | 38 | 152 | 172 |
| `GTP_mRNA` | 19 | 1 | 108 |
| `GTP_translat` | 102 | 5480 | 38 |
| `tRNA_translat` | 93 | 1644 | 41 |

## 14c.6: the published-parameter census with the GTP branch loosened

PGK3's and PYK3's catalytic constants, forward and reverse, are scaled by s, which scales each reaction's capacity by s at every handshake and keeps its equilibrium constant. Recorded, not gated.

| scale | seeds carrying a deficit | drains clipped, median (max) | `GTP_translat` drains, median (max) | `tRNA_translat` drains, median (max) | failed |
|---|---|---|---|---|---|
| 1 | 10 of 10 | 48 (253) | 47 (253) | 4 (20) | 0 |
| 2 | 10 of 10 | 10 (193) | 0 (4) | 10 (193) | 0 |
| 5 | 10 of 10 | 9 (429) | 0 (0) | 9 (429) | 0 |
| 10 | 9 of 10 | 22 (193) | 0 (0) | 22 (193) | 0 |
| 100 | 10 of 10 | 58 (195) | 0 (0) | 58 (195) | 0 |

At scale 1 the 10 seeds' records are **identical** to 14b's published-parameter census.

