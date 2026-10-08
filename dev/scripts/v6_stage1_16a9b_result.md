# 16a.9b, V6 stage 1: the pilot's chains extended to 700 sweeps

Staged by §12 2026-10-07 (approved). Runs at `281d548`, from
`.worktrees/v6-16a9b` with the pilot's sweep-200 checkpoints copied in. The code
is the pilot's (`6f8bcbc`) apart from the test extras' compat.

**The jobs:** chains 18157561 to 18157564 resumed at sweep 200 and ran to 700,
each for 27.5 to 29.8 h. The merge, 18157565, took sweeps 351 to 700, the
script's default of dropping half. The script prints "Pilot of 16a.9a" in its
header because it is the pilot's script.

| Parameter (log) | R̂ (rank) | Bulk ESS | Tail ESS | Bulk ESS per sweep | Tail ESS per sweep | MCSE(q05)/SD | MCSE(q95)/SD | Posterior mean ± SD | Truth | Truth's quantile |
|---|---|---|---|---|---|---|---|---|---|---|
| S_JCVISYN3A_0607 | 1.024 | 230 | 446 | 0.164 | 0.319 | 0.097 | 0.093 | 1.546 ± 0.117 | 1.489 | 0.31 |
| S_JCVISYN3A_0779 | 1.006 | 958 | 883 | 0.684 | 0.631 | 0.090 | 0.048 | 0.699 ± 0.204 | 1.012 | 0.95 |
| krnadeg | 1.020 | 303 | 725 | 0.216 | 0.518 | 0.108 | 0.073 | 2.224 ± 0.076 | 2.243 | 0.60 |
| kcatF_R_ENO | 1.002 | 389 | 599 | 0.278 | 0.428 | 0.076 | 0.125 | 4.243 ± 0.042 | 4.182 | 0.07 |
| kcatF_R_FBA | 1.003 | 1210 | 1013 | 0.864 | 0.724 | 0.062 | 0.076 | 4.404 ± 0.189 | 4.207 | 0.13 |
| sigma_b | 1.004 | 1016 | 1149 | 0.726 | 0.821 | 0.067 | 0.057 | -2.625 ± 0.057 | -2.622 | 0.52 |

**Against V6.** V6 is not met yet.
- **R̂:** four of the six parameters are below 1.01. The ptsG promoter
  (1.024) and `krnadeg` (1.020) are not.
- **ESS:** two parameters reach 1,000 in both bulk and tail ESS, `kcatF_R_FBA`
  and `sigma_b`. `S_0779` is close, at 958 bulk and 883 tail.
- **Quantile MCSE:** none is within 0.05 SD on both quantiles. The worst is
  `kcatF_R_ENO`'s 95% quantile, at 0.125.

**Against the pilot.**
- **Mixing:** the slowest bulk ESS per sweep rose from 0.108 to 0.164, still on
  the ptsG promoter. The pilot's rate rested on a bulk ESS of 43, and this one
  rests on 230.
- **The binding criterion** is now the quantile MCSE, at 8,801 kept sweeps,
  rather than ESS ≥ 1,000, at 6,083.
- **Update rates:** the lowest over the 414 live cell-windows is 0.39, and none
  is below 0.10. The empty ones average 0.15.
- **Cost:** 205 s per sweep (B3 130, B1 75), or 1.14 CPU-h at 20 cores. Core use
  is 62 to 64%, as in the pilot.

**The merge's projection:** V6 needs 10,201 sweeps in all at N_ref = 20, which
is 11.6k CPU-h and 6.1 days. That figure drops the first half of each chain,
counting the sweeps already run.

**From here:** the four chains need 8,801 kept sweeps in all.
- **With 350 sweeps of burn-in per chain,** 1,400 are kept so far, so 7,401
  remain. That is about 1,850 more per chain: **8.4k CPU-h and 4.4 days.**
- **With 100 sweeps of burn-in,** 2,400 are kept so far, so about 1,600 more
  per chain: **7.3k CPU-h and 3.8 days.** The R̂ above supports that, but it has
  not been checked against the trace.

**V7's projection** at R = 50 is 54.6k CPU-h, using this merge's burn-in of 350.
With 100 sweeps of burn-in it is roughly 40k. *Superseded below:* the
burn-in-100 merge and the thread scaling bring it to 25.6k at 5 threads.

**Compute:** 2,308.3 CPU-h (Slurm CPUTimeRAW).

## The same chains with 100 sweeps of burn-in (merge job 18271545)

This keeps sweeps 101 to 700, which is 2,400 draws.

| Parameter (log) | R̂ (rank) | Bulk ESS | Tail ESS | Bulk ESS per sweep | MCSE(q05)/SD | MCSE(q95)/SD | Posterior mean ± SD |
|---|---|---|---|---|---|---|---|
| S_JCVISYN3A_0607 | 1.011 | 432 | 907 | 0.180 | 0.071 | 0.075 | 1.541 ± 0.118 |
| S_JCVISYN3A_0779 | 1.004 | 1459 | 1470 | 0.608 | 0.054 | 0.042 | 0.703 ± 0.199 |
| krnadeg | 1.012 | 560 | 1129 | 0.233 | 0.069 | 0.069 | 2.222 ± 0.077 |
| kcatF_R_ENO | 1.001 | 653 | 1015 | 0.272 | 0.078 | 0.068 | 4.243 ± 0.042 |
| kcatF_R_FBA | 1.002 | 2161 | 1731 | 0.900 | 0.047 | 0.072 | 4.405 ± 0.188 |
| sigma_b | 1.003 | 1726 | 1728 | 0.719 | 0.048 | 0.033 | -2.624 ± 0.059 |

**The posterior means move by at most 0.005,** which is 0.04 SD, against the
burn-in-350 merge. The ESS per sweep is about the same or higher; the largest drop is ENO's, from
0.278 to 0.272. Nothing indicates the chains
were still moving after sweep 100, so **100 sweeps of burn-in is taken as
adequate**. The pilot's chains started from prior draws, so this also bears on
V7's fresh chains.

**R̂ is 1.011 and 1.012** for the ptsG promoter and `krnadeg`, just above V6's
1.01. Every other parameter passes.

**The projection.** The quantile MCSE needs 5,807 kept sweeps, which is 6,207
in all, or 1,552 per chain. **So 852 more sweeps per chain, 3,408 in all.** The
cost of those depends on the thread count (`thread_scaling_16a9b_result.md`):

| Threads per chain | CPU-h | Wall-clock |
|---|---|---|
| 20 | 4.1k | 2.2 days |
| 10 | 3.4k | 3.5 days |
| 5 | 2.7k | 5.5 days |

**V7 at R = 50:** 656 sweeps per replicate, 100 of burn-in plus 100 effective
draws at 0.180 per sweep.
- **At 20 threads:** 39.9k CPU-h and 1.7 days.
- **At 5 threads:** 25.6k CPU-h and 4.3 days, on 250 cores.
