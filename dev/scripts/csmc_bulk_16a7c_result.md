# 16a.7c: the path update under bulk metabolites

Code at `50e5af4`, with runs at `22d1f3e`.

**The tests pass,** 211 of 211 (job 18049605). They include:
- the bulk record recomputing from the latent and its seed;
- block 1 scoring the cells' mean;
- the CSMC bulk weight equalling the bulk likelihood with the particle in place;
- a sweep's returned latent equalling its path's replay.

## Update rates: KP no longer fires (jobs 18052532 to 18052536)

**The setup:**
- **The population:** M0 replicate 16085's 50 cells, observed in bulk at the
  replicate's σ = 0.0727 as σ_b.
- **The update:** cells 1 to 5 are updated in turn by PG, each conditioned on the
  other 49. The rest are held at their true paths.
- **The run:** 16 scans, so 80 cell-sweeps per N, from the true paths.
- **Cost** is per cell-sweep, in plain M0 cycles timed in the same job, with the
  first scan dropped.

| N | Cycles per cell-sweep | Lowest window rate | Windows below 10% | Per-window rate, windows 1 to 10 |
|---|---|---|---|---|
| 5 | 4.2 | 0.05 (window 1) | 3 of 10 | 0.05 0.06 0.07 0.12 0.19 0.25 0.24 0.29 0.46 0.56 |
| 10 | 8.4 | 0.28 (window 1) | 0 of 10 | 0.28 0.31 0.36 0.28 0.31 0.38 0.40 0.41 0.66 0.71 |
| 20 | 16.5 | 0.55 (window 7) | 0 of 10 | 0.61 0.65 0.70 0.69 0.70 0.70 0.55 0.57 0.74 0.78 |
| 50 | 40.3 | 0.59 (window 8) | 0 of 10 | 0.78 0.79 0.82 0.72 0.75 0.74 0.60 0.59 0.80 0.86 |

- **Every window clears the 10% threshold from N = 10.**
- **Cost is about 0.84N cycles per cell-sweep,** in line with D16.6's PG row.
- **Proposed: N = 10.** Its lowest window rate is 2.8 times the threshold, at half
  N = 20's cost. N = 20 is the fallback if V6 shows slow mixing in early windows.

## V5's bulk "on" half: passes, with a watch item (prep, kernel and merge jobs 18052537 to 18052542)

**The setup:** three cells share the data's transcripts. The reference is 200,000
random triples of kept runs, weighted by the bulk likelihood of their mean. The
kernel draws are 2,000 chains of five scans at N = 5. The functionals are averaged
over the three cells.

| Case | σ_b | Reference ESS | z: births w1, w2, w3 | z: translation w1, w2, w3 |
|---|---|---|---|---|
| Toy | 0.1 (the data's) | 23,353 | 2.51, 0.02, 0.33 | −0.19, −0.19, 1.60 |
| M0 | 0.8 (inflated 8×) | 11,588 | 0.89, **2.85**, **−2.81** | −0.81, 2.29, −0.79 |

**Both pass at 3 SE.**

**Watch item: births in window 2 on M0.** In every M0 V5 comparison so far, the
kernel puts more births in window 2 than the reference does.

| Run | Kernel | Reference | z |
|---|---|---|---|
| Per-cell, metabolites off | 0.018 | 0.013 | about 1.6 |
| Per-cell, metabolites off, annealed | 0.018 | 0.013 | about 1.6 |
| Per-cell, metabolites on | 0.015 | 0.011 | 1.38 |
| Per-cell, metabolites on, annealed | 0.013 | 0.009 | 1.24 |
| Bulk | 0.013 | 0.009 | 2.85 |

- **The comparisons are not independent.** Every one uses the same rejection
  reference set, so one low draw of that set would make them all agree.
- **The bulk run's −2.81 for births in window 3** is the mirror image: windows 2
  and 3 trade a birth.
- **What would decide it:** a fresh rejection reference with new seeds, about 30
  CPU-h.

## Still to do for 16a.7c

- **V4 on block 1's bulk form.** Block 1's bulk term is implemented and unit-tested,
  but the grid check has not been rerun.
- **The full-scale half of parent 15.7b.** 15.7's bulk observations need the latent
  record regenerated, about 1.8 CPU-h.
