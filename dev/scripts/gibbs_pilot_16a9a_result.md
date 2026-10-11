# 16a.9a: the Gibbs chain on M0, and a pilot to size V6 and V7

Code at `226d983` (the chain), with the pilot at `6f8bcbc` and the tests of record
at `b030512`.

**The tests pass,** 28,658 of 28,658 (job 18109288 at `b030512`). The run before
it (job 18105666 at `226d983`) failed one test, Aqua's compat check on the new
Serialization test extra, which `b030512` bounds. Its 43 Gibbs tests passed. They
include:
- B2 being exact on the hybrid: its conditional's difference between two values of
  each promoter and of `krnadeg` equals the replay's path log-density difference
  plus the prior difference, to 1e-10;
- a base written at once replaying bitwise as one written in steps;
- the threaded B1 replay and the threaded chain being bitwise equal to serial at
  one seed (on 4 threads);
- a resumed chain continuing bitwise as the uninterrupted one;
- a sweep keeping the chain's invariants.

## The pilot (chains 18108335 to 18108338, merge 18108339, at `6f8bcbc`)

**The setup:**
- **The data:** M0 replicate 16085 (16a.7c's), 50 cells observed in bulk.
- **The chains:** 4, from θ and σ_b drawn from their priors with
  `Xoshiro(16_090_000 + k)`, and from paths drawn by transcript-weighted SMC at
  that θ. So the starts are overdispersed.
- **The sampler:** PG at N = 20, 200 sweeps per chain, the first 100 dropped.
  20 threads run a cell's particles in B3, and the cells in B1.

| Parameter (log) | R̂ (rank) | Bulk ESS | Tail ESS | Bulk ESS per sweep | Tail ESS per sweep | MCSE(q05)/SD | MCSE(q95)/SD | Posterior mean ± SD | Truth | Truth's quantile |
|---|---|---|---|---|---|---|---|---|---|---|
| S_JCVISYN3A_0607 | 1.090 | 43 | 185 | 0.108 | 0.463 | 0.180 | 0.128 | 1.528 ± 0.123 | 1.489 | 0.38 |
| S_JCVISYN3A_0779 | 0.999 | 222 | 145 | 0.555 | 0.363 | 0.208 | 0.106 | 0.700 ± 0.196 | 1.012 | 0.96 |
| krnadeg | 1.059 | 96 | 107 | 0.240 | 0.268 | 0.154 | 0.223 | 2.218 ± 0.079 | 2.243 | 0.63 |
| kcatF_R_ENO | 1.010 | 132 | 243 | 0.329 | 0.608 | 0.152 | 0.168 | 4.240 ± 0.043 | 4.182 | 0.09 |
| kcatF_R_FBA | 1.008 | 326 | 268 | 0.816 | 0.669 | 0.066 | 0.249 | 4.414 ± 0.206 | 4.207 | 0.17 |
| sigma_b | 1.018 | 280 | 205 | 0.699 | 0.513 | 0.201 | 0.080 | -2.623 ± 0.060 | -2.622 | 0.52 |

**Mixing.** The ptsG promoter, `S_JCVISYN3A_0607`, mixes slowest: bulk ESS 0.108
per sweep, with R̂ 1.090 over 100 kept sweeps per chain. `krnadeg` is next, at
R̂ 1.059. The pilot does not claim convergence by V6's standard. It measures the
rate that sizes V6.

**The rate is rough.** The slowest ESS per sweep rests on a bulk ESS of 43, so
the projections below could be off by tens of percent either way. The first few
hundred sweeps of 16a.9b refine it.

**Update rates.** Over 414 live cell-windows the lowest is 0.34, and none is below
D16.3's 0.10. The 86 cell-windows with empty true paths average 0.15.

**Truths.** Every truth's quantile lies between 0.09 and 0.96. That is one
replicate, so it is a sanity check, not calibration.

**Cost per sweep,** wall-clock at 20 threads: 199 s.
- **B3** takes 126 s and **B1** 73 s. B2 and the rebuild take under 0.01 s.
- That is 111.5 plain M0 cycles of 1.79 s each, or **1.11 CPU-h per sweep** at
  20 cores. Initialisation takes 141 s.
- **Core use is 62 to 63%** across the four chains (Slurm job reports), because
  B3 updates the cells in turn and only a cell's particles run in parallel.

## Projections

**V6** (4 chains; R̂ < 1.01, bulk and tail ESS ≥ 1,000, quantile MCSE ≤ 0.05 SD):
- ESS ≥ 1,000 at 0.108 per sweep needs 9,284 kept sweeps. The quantile MCSE needs
  9,927, so it binds.
- **N_ref = 20:** 10,327 sweeps in all, 2,582 per chain. That is **11,444 CPU-h**
  and **6.0 days** of wall-clock with the chains in parallel.
- **N_ref = 50:** the same sweeps, 22,322 CPU-h, and 11.6 days.

**V7** (one chain per replicate; burn-in 100 plus 100 effective draws, 1,028
sweeps):
- **R = 50: 56,978 CPU-h,** 2.4 days of wall-clock with every replicate in
  parallel.
- **R = 100:** 113,957 CPU-h.

**Against D16.6.** 16a.9's CPU-h are recorded but not held to 16a's 5k cap
(§12 2026-10-06). The phase as a whole has 50k, though. V6 at N_ref = 20 plus V7
at R = 50 comes to about 68k CPU-h, which is more than the whole phase cap. §12
2026-10-07 takes this up.

Regenerate: `dev/scripts/gibbs_pilot_16a9a.sh` from a worktree at `6f8bcbc`.
The per-sweep logs and checkpoints are in that worktree's
`dev/scripts/gibbs_pilot_16a9a/`.
