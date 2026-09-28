# 14c.5: the sampler model against the published one

Merged by `dev/scripts/corea_smoothing_14c_merge.jl` from 50 task files. Regenerate with `sbatch dev/scripts/corea_smoothing_14c.slurm`, then this script. The tasks' headers:

- commit b377dde, job 17625499, 2026-09-28T18:11:49.964, Julia 1.10.5
- commit b377dde, job 17625500, 2026-09-28T18:11:50.092, Julia 1.10.5
- commit b377dde, job 17625509, 2026-09-28T18:11:58.724, Julia 1.10.5
- commit b377dde, job 17625510, 2026-09-28T18:11:58.403, Julia 1.10.5
- commit b377dde, job 17625511, 2026-09-28T18:12:06.092, Julia 1.10.5
- commit b377dde, job 17625512, 2026-09-28T18:12:05.658, Julia 1.10.5
- commit b377dde, job 17625513, 2026-09-28T18:12:05.954, Julia 1.10.5
- commit b377dde, job 17625514, 2026-09-28T18:12:06.092, Julia 1.10.5
- commit b377dde, job 17625515, 2026-09-28T18:12:05.592, Julia 1.10.5
- commit b377dde, job 17625516, 2026-09-28T18:12:05.982, Julia 1.10.5
- commit b377dde, job 17625517, 2026-09-28T18:12:05.904, Julia 1.10.5
- commit b377dde, job 17625518, 2026-09-28T18:12:05.858, Julia 1.10.5
- commit b377dde, job 17625501, 2026-09-28T18:11:50.141, Julia 1.10.5
- commit b377dde, job 17625519, 2026-09-28T18:12:05.891, Julia 1.10.5
- commit b377dde, job 17625520, 2026-09-28T18:12:05.792, Julia 1.10.5
- commit b377dde, job 17625521, 2026-09-28T18:12:06.012, Julia 1.10.5
- commit b377dde, job 17625522, 2026-09-28T18:12:05.785, Julia 1.10.5
- commit b377dde, job 17625523, 2026-09-28T18:12:06.107, Julia 1.10.5
- commit b377dde, job 17625524, 2026-09-28T18:12:05.697, Julia 1.10.5
- commit b377dde, job 17625525, 2026-09-28T18:12:05.933, Julia 1.10.5
- commit b377dde, job 17625526, 2026-09-28T18:12:05.958, Julia 1.10.5
- commit b377dde, job 17625527, 2026-09-28T18:12:05.402, Julia 1.10.5
- commit b377dde, job 17625528, 2026-09-28T18:12:05.724, Julia 1.10.5
- commit b377dde, job 17625502, 2026-09-28T18:11:49.969, Julia 1.10.5
- commit b377dde, job 17625529, 2026-09-28T18:12:06.091, Julia 1.10.5
- commit b377dde, job 17625530, 2026-09-28T18:12:05.730, Julia 1.10.5
- commit b377dde, job 17625531, 2026-09-28T18:12:06.015, Julia 1.10.5
- commit b377dde, job 17625532, 2026-09-28T18:12:05.795, Julia 1.10.5
- commit b377dde, job 17625533, 2026-09-28T18:12:06.178, Julia 1.10.5
- commit b377dde, job 17625534, 2026-09-28T18:12:05.851, Julia 1.10.5
- commit b377dde, job 17625535, 2026-09-28T18:12:06.025, Julia 1.10.5
- commit b377dde, job 17625536, 2026-09-28T18:12:05.724, Julia 1.10.5
- commit b377dde, job 17625537, 2026-09-28T18:12:05.611, Julia 1.10.5
- commit b377dde, job 17625538, 2026-09-28T18:12:05.469, Julia 1.10.5
- commit b377dde, job 17625503, 2026-09-28T18:11:50.174, Julia 1.10.5
- commit b377dde, job 17625539, 2026-09-28T18:12:05.611, Julia 1.10.5
- commit b377dde, job 17625540, 2026-09-28T18:12:05.071, Julia 1.10.5
- commit b377dde, job 17625541, 2026-09-28T18:12:05.415, Julia 1.10.5
- commit b377dde, job 17625542, 2026-09-28T18:12:05.350, Julia 1.10.5
- commit b377dde, job 17625543, 2026-09-28T18:12:00.355, Julia 1.10.5
- commit b377dde, job 17625544, 2026-09-28T18:12:00.554, Julia 1.10.5
- commit b377dde, job 17625545, 2026-09-28T18:11:59.734, Julia 1.10.5
- commit b377dde, job 17625546, 2026-09-28T18:12:00.438, Julia 1.10.5
- commit b377dde, job 17625547, 2026-09-28T18:12:00.406, Julia 1.10.5
- commit b377dde, job 17625498, 2026-09-28T18:12:00.432, Julia 1.10.5
- commit b377dde, job 17625504, 2026-09-28T18:11:50.048, Julia 1.10.5
- commit b377dde, job 17625505, 2026-09-28T18:12:01.796, Julia 1.10.5
- commit b377dde, job 17625506, 2026-09-28T18:12:02.015, Julia 1.10.5
- commit b377dde, job 17625507, 2026-09-28T18:11:58.423, Julia 1.10.5
- commit b377dde, job 17625508, 2026-09-28T18:11:58.421, Julia 1.10.5

100 seeds, each a full cycle under the published model (clamped drain, fractional carry), the sampler model (drain smoothed at one particle, pools carried continuously), and the control (the published model with `krnadeg` scaled by 1 + 1e-5). 50 observables (every ODE state in particles with its carried remainder, every transcript, the volume) at 106 save points. Each difference is paired per seed and relative to max(|published|, floor): 500 particles for an ODE state, one copy for a transcript, none for the volume.

## The sampler model against the published model

**Coupling.** 47 of 100 pairs decouple within the cycle, at a median handshake of 3303 (earliest 179, latest 6281).

**Coupled half.** 406700 paired comparisons at save points before decoupling. The largest is 0.763% (`M_gtp_c`, seed 15059, t = 1920 s). Gate: every one within 1% — **passes**.

**Decoupled half, at t = 6300 s.** Mean over 100 seeds of the paired difference, per observable. 20 observables pass, 0 fail, and 30 are unresolved (SE above 1%).

The mean over seeds, its SE, and z = mean/SE, so a mean far from zero is visible even when the SE leaves it unresolved.

| observable | mean | SE | z | largest single seed | verdict |
|---|---|---|---|---|---|
| `M_gtp_c` | 37.8% | 16% | 2.33 | 866% | unresolved |
| `mRNA_JCVISYN3A_0131` | 16.5% | 7.6% | 2.18 | 400% | unresolved |
| `mRNA_JCVISYN3A_0607` | 15.9% | 10% | 1.59 | 500% | unresolved |
| `mRNA_JCVISYN3A_0213` | 15.8% | 6.8% | 2.33 | 400% | unresolved |
| `mRNA_JCVISYN3A_0475` | 15.4% | 7.3% | 2.12 | 400% | unresolved |
| `M_pi_c` | 14% | 8.8% | 1.59 | 622% | unresolved |
| `M_lac__L_c` | 13.5% | 7.4% | 1.82 | 537% | unresolved |
| `mRNA_JCVISYN3A_0727` | 11% | 6% | 1.84 | 300% | unresolved |
| `M_trna_c` | 9.09% | 6.1% | 1.48 | 390% | unresolved |
| `mRNA_JCVISYN3A_0220` | 8.5% | 5.7% | 1.48 | 300% | unresolved |
| `mRNA_JCVISYN3A_0203` | 7% | 4.3% | 1.62 | 300% | unresolved |
| `mRNA_JCVISYN3A_0221` | 3.92% | 6.6% | 0.59 | 300% | unresolved |

*Diagnostic, not the gate:* the difference of the ensemble means at t = 6300 s, over max(|published mean|, floor), which is not skewed by the gate's per-seed ratio. The largest |z| over 50 observables is 1.72 (`mRNA_JCVISYN3A_0233`, -8% ± 4.6%); 0 exceed 3, and 0 of the 20 with SE within 1% move by more than 1%.

## The control against the published model

**Coupling.** 3 of 100 pairs decouple within the cycle, at a median handshake of 1827 (earliest 1188, latest 4791).

**Coupled half.** 520650 paired comparisons at save points before decoupling. The largest is 0% (`none`, seed 0, t = 0 s). Gate: every one within 1% — **passes**.

**Decoupled half, at t = 6300 s.** Mean over 100 seeds of the paired difference, per observable. 48 observables pass, 0 fail, and 2 are unresolved (SE above 1%).

The mean over seeds, its SE, and z = mean/SE, so a mean far from zero is visible even when the SE leaves it unresolved.

| observable | mean | SE | z | largest single seed | verdict |
|---|---|---|---|---|---|
| `M_trna_c` | 1.54% | 1.5% | 1.00 | 154% | unresolved |
| `M_pi_c` | 1.39% | 1.4% | 1.00 | 139% | unresolved |
| `mRNA_JCVISYN3A_0131` | -1% | 1% | -1.00 | 100% | passes |
| `mRNA_JCVISYN3A_0213` | 1% | 1% | 1.00 | 100% | passes |
| `mRNA_JCVISYN3A_0220` | -1% | 1% | -1.00 | 100% | passes |
| `mRNA_JCVISYN3A_0221` | -1% | 1% | -1.00 | 100% | passes |
| `mRNA_JCVISYN3A_0694` | -1% | 1% | -1.00 | 100% | passes |
| `mRNA_JCVISYN3A_0727` | 1% | 1% | 1.00 | 100% | passes |
| `mRNA_JCVISYN3A_0779` | 1% | 1% | 1.00 | 100% | passes |
| `M_3pg_c` | 0.597% | 0.6% | 1.00 | 59.7% | passes |
| `M_lac__L_c` | 0.583% | 0.58% | 1.00 | 58.3% | passes |
| `mRNA_JCVISYN3A_0475` | 0.5% | 0.5% | 1.00 | 50% | passes |

*Diagnostic, not the gate:* the difference of the ensemble means at t = 6300 s, over max(|published mean|, floor), which is not skewed by the gate's per-seed ratio. The largest |z| over 50 observables is 1.01 (`M_f6p_c`, 0.288% ± 0.29%); 0 exceed 3, and 0 of the 49 with SE within 1% move by more than 1%.

## Verdict

| | pairs decoupled | coupled half | end-of-cycle observables failing | unresolved |
|---|---|---|---|---|
| sampler model | 47 of 100 | passes | 0 | 30 |
| control | 3 of 100 | passes | 0 | 2 |

14c.5's agreement gate: the coupled half **passes**. The end-of-cycle half resolves 20 of 50 observables, of which 0 fail, and it can neither pass nor fail the other 30: their SE exceeds 1%, because once a pair decouples the per-seed ratio is Monte Carlo noise, skewed upward. For those observables the evidence is the symmetric diagnostic above, which is not the gate.

The control was specified as matched to the sampler model's nudge in size. It is not matched in effect: it decouples 3 of 100 pairs against the sampler model's 47, so it shows that the model is deterministic enough to stay paired under a small nudge, and cannot serve as a null for the sampler model's decoupled half.

## Clipping under each run

The two models' counts use different predicates in effect. A smoothed debit leaves a residue `w·exp(−gap/w)`, and the census counts a clip while that residue exceeds 1e-6 of the accrual: whenever the pool is less than `w·ln(1e6/accrual)` particles above the accrual, at most about 14, although the debit paid in full to within 1e-6 of it. So the sampler row is not a like-for-like count, and it is recorded, not compared.

| | median drains clipped | seeds carrying a deficit |
|---|---|---|
| `clamped` | 64.0 | 93 of 100 |
| `smoothed` | 55.0 | 93 of 100 |
| `control` | 64.0 | 93 of 100 |

commit b377dde, job 17625499, task 0 of 50, 2026-09-28T18:11:49.964, Julia 1.10.5

## 14c.5 and 14c.7: the derivative across a clip, on the assembled model

Seed 14800: `GTP_translat` first clips at handshake 1580 (t = 1580 s). `kcatF_R_PGK3` is 140.8 at nominal, and that drain stops clipping at 188.543 (1.339 × nominal). The pool's sensitivity there, dP/d ln θ on the paying side, is 501.99 particles.

The slope, in ln θ, of the GTP state after the debit, either side of the clip, as a fraction of that sensitivity. Finite-difference step 1.99e-05 in ln θ. All three variants start from the published model's saved state.

| spacing either side (steps) | pool − accrual at the spacing (particles) | clamp, continuous pools | sampler model |
|---|---|---|---|
| 100 | ±1 | 0.9156 | 0.3923 |
| 30 | ±0.3 | 0.9153 | 0.125 |
| 10 | ±0.1 | 0.9151 | 0.0419 |
| 3 | ±0.03 | 0.9151 | 0.01258 |

A continuous derivative has a jump that falls with the spacing. The clamp's stays at the whole sensitivity, since one side pays and the other floors at zero.

The published model on the same 21 points, spaced 37 steps apart across the clip: the state is a whole number of particles at 18 of them, and the finite-difference slope is exactly zero at 21. Its derivative is zero wherever no one-particle step falls inside the stencil, so a gradient through the handshake sees no dependence on θ at all.

*Note added at merge.* The clamped jump reads 0.915, not one, because the sensitivity it is divided by is taken at 1.5 × the clip, where it is larger. The paragraph above, as the driver at the run's commit wrote it, says "the whole sensitivity"; the driver's wording is corrected for future runs.
