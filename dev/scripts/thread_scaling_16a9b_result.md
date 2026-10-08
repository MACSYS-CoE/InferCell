# 16a.9b: a Gibbs sweep's cost per thread count

§12 2026-10-07, option 2. Script `thread_scaling_16a9b.jl` at `2a200db`, from
chain 1's state at sweep 700 of V6 stage 1. Each job ran one warm-up sweep (701)
and two timed sweeps (702 and 703).

**The jobs:**
- 1 thread: job 18258719.
- 5 threads: job 18258720.
- 10 threads: job 18256436.
- 20 threads: job 18256437.

The first 1- and 5-thread jobs, 18256434 and 18256435, were cancelled by root
2 minutes in, on two different nodes at the same moment, before writing a log.
They were resubmitted unchanged.

**The chain is the same at every thread count.** All four runs end at the same
state, bitwise:
- cme = `[4.838625848831599, 2.112034916160309, 9.163265764941675]`;
- u = `[4.1882452276922955, 4.607664248000716]`;
- σ = `0.07227565235453373`.

So sweeps 702 and 703 are the same work in each run, and the comparison is
paired. Sweep 703's B1 is heavier than 702's in every run (about 1.9×). That is
the slice sampler's variable work, not noise between runs.

| Threads | Wall per sweep (s) | B3 core-seconds | B1 core-seconds | CPU-h per sweep | B3 efficiency | B1 efficiency | GC share |
|---|---|---|---|---|---|---|---|
| 1 | 2,890 | 1,724 | 1,167 | 0.803 | 1.00 | 1.00 | 1.4% |
| 5 | 562 | 1,719 | 1,090 | 0.780 | 1.00 | 1.07 | 3.6% |
| 10 | 359 | 2,266 | 1,324 | 0.998 | 0.76 | 0.88 | 4.2% |
| 20 | 219 | 2,752 | 1,629 | 1.216 | 0.63 | 0.72 | 7.5% |

The figures are the means of sweeps 702 and 703. Efficiency is the 1-thread
core-seconds over this row's core-seconds.

**5 threads costs 36% fewer CPU-h per sweep than 20,** at 2.6× the
wall-clock. At 5 threads B3 is as efficient as serial. Above 5, the per-window
barrier, the serial resampling copies and GC take over
(`profile_b3_16a9b_result.md`).

## What this does to 16a.9b

These use the burn-in-100 merge (`v6_stage1_16a9b_result.md`, job 18271545).

**V6:** the four chains, continued from sweep 700, need 852 more sweeps each,
3,408 in all.

| Threads per chain | CPU-h | Wall-clock |
|---|---|---|
| 20 | 4.1k | 2.2 days |
| 10 | 3.4k | 3.5 days |
| 5 | 2.7k | 5.5 days |

**V7 at R = 50:** one chain per replicate, 656 sweeps each.

| Threads per replicate | CPU-h | Wall-clock with every replicate in parallel |
|---|---|---|
| 20 | 39.9k | 1.7 days |
| 5 | 25.6k | 4.3 days, on 250 cores |

**Against the 50k phase cap:** about 3.6k has been spent. V7 at 5 threads plus
V6 at any thread count comes to 28.3k to 29.7k more, so 31.9k to 33.4k in all.
That leaves 16.6k to 18.1k for 16a.10, 16a.11 and 16b. 16b's budget is reopened
at 16a.11 in any case (D16.6). V7 at 20 threads would take the phase to about
47k, leaving almost nothing.

**Compute for these measurements:** 18.3 CPU-h, including the profile (job
18255035) and the two cancelled jobs.
