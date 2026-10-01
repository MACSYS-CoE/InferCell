# 16a.7: the conditional SMC path update, and V5

## Unit tests (job 17854107 at `d701490`; `test/test_csmc.jl`, 40 of 40)

- **The weight's pieces.** The birth and death rates read from the driver, and
  each log transition probability, match a matrix exponential to 1e-12. The rate
  constant is set to k and then to 2k.
- **The reference particle** replays its window bitwise.
- **Per-particle streams** (sub-spec §12 2026-09-30): alike seeds give identical
  proposals, and different seeds differ. A proposal always hits the observed
  count.
- **V9:** restored particles share no driver state.
- **A sweep** (particle Gibbs, and ancestor sampling at lag 2) returns a path
  whose replay reproduces the data's transcripts.
- **M0's transcript map,** and one M0 window advanced.

## V5 (`dev/scripts/csmc_v5_16a7.jl`, at `644c140`)

The setup:
- **Rejection reference.** At fixed θ, the published model is simulated and the
  runs whose transcripts match the data at every window are kept. Each kept run is
  an exact draw of p(X | transcripts).
- **Kernel.** Particle Gibbs with N = 5. 2,000 chains each start from a distinct
  exact draw and run 5 sweeps, and only the final path is kept. If the kernel
  leaves the target invariant, those 2,000 paths are independent exact draws.
- **Functionals.** Each window's transcript births and the gene's translation
  count, 6 in all.
- **Metabolite likelihood off.** A two-sample χ² against the other half of the
  kept runs. 6 tests at p > 0.01 give a family-wise false-failure rate of about 6%.
- **Metabolite likelihood on.** The kept runs weighted by the panel's likelihood
  at the smallest σ, doubling from 0.1, whose ESS is at least 1,000. Expectations
  must agree within 3 SE.

### M0 cut to ptsG, three windows (jobs 17854476 to 17854480): passes

8,710 kept runs from 100,000 (acceptance 8.7%).

| Functional | Off: kernel / rejection mean | χ² (df), p | On (σ = 3.2, ESS 2,459): kernel / reference | z |
|---|---|---|---|---|
| births w1 | 0.018 / 0.017 | 0.05 (1), 0.826 | 0.015 / 0.011 | 1.38 |
| births w2 | 0.018 / 0.013 | 2.58 (1), 0.108 | 0.009 / 0.009 | −0.02 |
| births w3 | 1.008 / 1.011 | 1.17 (1), 0.280 | 1.012 / 1.010 | 0.39 |
| translation w1 | 3.446 / 3.447 | 9.22 (9), 0.417 | 3.696 / 3.722 | −0.51 |
| translation w2 | 3.519 / 3.478 | 13.32 (9), 0.149 | 2.401 / 2.408 | −0.14 |
| translation w3 | 5.276 / 5.304 | 20.88 (14), 0.105 | 5.502 / 5.578 | −1.02 |

Within 5 sweeps, windows 1, 2 and 3 changed in 97%, 100% and 100% of chains.
- **The births tests are weak in windows 1 and 2.** The data's transcripts force
  births there in fewer than 2% of runs, so those two χ² tests have little power.
  The translation counts carry most of the check.
- **The metabolite-on half needed σ = 3.2,** 32 times the data's, to reach an ESS of
  1,000 over the 17-pool panel. So it tests the weight at a weaker metabolite
  likelihood than the data's.

### The amplified toy (jobs 17854481 to 17854485): passes

10,924 kept runs from 500,000 (acceptance 2.2%).

| Functional | Off: kernel / rejection mean | χ² (df), p | On (σ = 0.1, ESS 2,089): kernel / reference | z |
|---|---|---|---|---|
| births w1 | 58.35 / 58.14 | 34.06 (40), 0.734 | 60.10 / 60.17 | −0.32 |
| births w2 | 34.07 / 34.04 | 34.91 (36), 0.520 | 32.90 / 32.81 | 0.49 |
| births w3 | 10.11 / 10.27 | 20.63 (21), 0.482 | 8.97 / 8.78 | 2.11 |
| translation w1 | 111.76 / 110.96 | 103.58 (102), 0.438 | 119.24 / 119.52 | −0.90 |
| translation w2 | 67.82 / 66.85 | 108.70 (87), 0.058 | 65.15 / 64.57 | 1.22 |
| translation w3 | 20.24 / 20.73 | 61.55 (53), 0.197 | 18.01 / 17.39 | 2.23 |

Within 5 sweeps, windows changed in 96%, 99% and 100% of chains.

**The amplified case's conditions, measured on matching runs.** The first constant
is 1/s.

| Constant from | Spread between runs: 90th/10th percentile, max/min | Each run's drop across that rebuild: median, 10th percentile |
|---|---|---|
| 59 s | 1.31, 2.32 | 1.69×, 1.49× |
| 119 s | 2.09, 9.52 | 3.82×, 2.98× |
| 179 s | 3.36, 36.13 | 5.75×, 4.45× |

- **Both 2× conditions hold from the second rebuild on.** At 119 s and 179 s the
  constants differ by at least 2× between runs, and each run's drops are 3 to 6×
  with transcription at order 1/s.
- **At the first rebuild they do not quite hold.** The spread is 1.31 at the
  90th/10th percentile, and the median drop is 1.69×.
- **Watch item.** With the metabolite likelihood on, window 3's births and
  translations both read high, at z = 2.11 and 2.23. That passes, but it is in the
  window where the constants differ most.

## What the smoke runs fixed before the full run

- The toy's pool, drained by translation cost, emptied and set k = 0 (job
  17853650). It now decays first order at a rate set by the live protein count.
- **A proposal that cannot reach the observed count threw** (job 17853653). Its
  target density is zero, so it now gets log weight −Inf and is not advanced. That
  is a change in `src/csmc.jl`.
- The merge read z = NaN for a functional with no spread on either side. It now
  passes only if the means agree.

## Not yet done in 16a.7

The variant choice is still to come: per-window update rates and cost per sweep,
on M0 and on one full-scale cell, for PG, PGAS and truncated PGAS at
N ∈ {5, 10, 20, 50}, recorded in §4.
