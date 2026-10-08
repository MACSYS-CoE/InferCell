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
With 100 sweeps of burn-in it is roughly 40k. Either way it is not fundable
within the phase cap alongside V6. Its sizing waits for the B3 profiling
(§12 2026-10-07).

**Compute:** 2,308.3 CPU-h (Slurm CPUTimeRAW).
