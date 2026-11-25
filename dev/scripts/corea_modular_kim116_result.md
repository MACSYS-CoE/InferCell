# KIM-116: the modular arm beside the joint

The reduced model was inferred one block at a time on KIM-92's dataset, and the result
sits beside KIM-92's joint posterior. The pre-registration is `prereg/kim-116`
(`spec/prereg/kim-116.md`, tag on 53e2861). There were two variants of the interface
plug-in (Decision log 2026-11-19).

**The factor cut is close to the joint, but costs more.** Every shift is at most 0.27 joint
SD and every contraction is within 0.006. It costs 1.44× the joint's CPU-h. Its stage 1
solves the same ODE over the same 21 coordinates as the joint, so it costs the same
(1.00×), and stage 2 comes on top.

**The data plug-in is cheaper, but not close.** It costs 0.80× the joint, and its critical
path is 56 h against 121 h. It fails the frozen "close" test: promoter 0607 shifts +0.65
joint SD (Monte Carlo SE 0.07). The bias is systematic: 16 of 17 promoters shift up, with
a mean shift of +0.20 SD. On D3 the plug-in matches the joint for dQ and for the families
at the later horizon. For Q and Q_kd its CRPS is about twice the joint's, because its
median sits about 1% high.

Only the plug-in has D3 scores (prereg §6). The close variant has none, and the scored
variant is not close.

## Run

- **Commits:** the stages ran from `dev/scripts/corea_modular_kim116.{jl,sh}` at 7b6c0b1
  (through 23 Nov) and 12c066a (after). The two commits differ only in the predict files
  (`git diff 7b6c0b1 12c066a`), so the modular driver is the same in every job.
- **D3:** the scores came from `dev/scripts/corea_modular_kim116_predict.{jl,sh}` at
  12c066a.
- **This note:** its tables are from `dev/scripts/corea_modular_kim116_compare.jl` at
  e129fcc (Slurm 21628001).
- **Archive:** `/fred/oz022/tkimpson/infercell-kim116-archive-2026-11-26/`, covered by
  `SHA256SUMS`. It holds:
  - `run/`: the modular run, both variants, including the pre-extension copies in
    `{cut,plugin}/pre_extension/`;
  - `predict/`: the plug-in's D3 run;
  - `profile/`: item 2;
  - `sacct.tsv`: every job in `run/jobs.tsv` and `predict/jobs.tsv`, plus the profile;
  - `compare/`: `compare.md` and its log.
- **The joint:** `/fred/oz022/tkimpson/infercell-kim92-archive-2026-11-17/`.
- **D3 beside it:** `/fred/oz022/tkimpson/infercell-kim95-archive-2026-11-19/`.
- **Job states:** `sacct` shows COMPLETED for all 31 jobs (88 allocations with array
  tasks).

## The gates (prereg §4)

Stage 1 had to reach R̂ ≤ 1.02 and bulk ESS ≥ 100 on the 18 φ coordinates, with one
extension allowed. Both variants needed the extension.

| variant | gate on 200 draws per chain | gate on 400 draws per chain | stage 2 |
|---|---|---|---|
| plug-in | fails 16/18: 0727 R̂ 1.0246, 0233 1.0218 (21315914) | **passes 18/18**: max R̂ 1.0099, min ESS 2264 (21375141) | 21375287: 200 inner chains, 0 numerical errors, mean acceptance 0.90 |
| cut | fails 17/18: 0445 R̂ 1.0299 (21457208) | **passes 18/18**: max R̂ 1.0114, min ESS 2022 (21596404) | 21597007: 200 inner chains, 0 numerical errors, mean acceptance 0.90 (lowest chain 0.69) |

- **Chains:** 21208676–79 for the cut and 21208681–84 for the plug-in.
- **Extensions:** 21462423–26 for the cut and 21316056–59 for the plug-in.
- **PPC** (21627293–94): 10 of 884 rows fall outside the central 99% band for each
  variant, against the joint's 11.
- **Merge:** 21627295.

## Cost

These tables are copied unchanged from `compare/compare.md`.

| arm | stage 1 CPU-h | stage 2 CPU-h | PPC + merge CPU-h | total CPU-h | ÷ joint | critical path, h |
|---|---|---|---|---|---|---|
| joint (KIM-92) | 452.5 | — | 2.0 | 454.6 | 1.00 | 121.0 |
| cut | 453.0 | 198.2 | 1.1 | 652.3 | 1.44 | 140.9 |
| plugin | 165.4 | 198.7 | 1.1 | 365.2 | 0.80 | 56.1 |

- **What is not in the table:** every arm shares KIM-92's precompile, prep (δ, S and the
  800-cell reference ensemble), probe and profile, 32.0 CPU-h paid once.
- **Reconciling with T1:** the joint's 454.6 plus the shared 32.0 is KIM-92's merge total,
  486.60. KIM-92's first pass adds 5.32, which gives T1's 491.9.
- **The merge's own total:** `run/modular_summary.md` says "1017.5 CPU-h beside the
  joint's 491.9". That adds both variants together against one joint, so it is not a
  comparison.
- **Stage 1:** the cut's stage 1 costs 1.00× the joint's chains, check and extension, as
  item 2's profile predicted (36.9 s against 37.9 s per gradient). The plug-in's costs
  0.37× (11.3 s per gradient, with no solve).
- **Stage 2:** this is 30% of the cut's cost and 54% of the plug-in's. Its size, K = 200
  stage-1 draws × NUTS(100, 0.8), was frozen in the pre-registration and never tuned for
  cost. So the plug-in's 0.80× belongs to this design. It is not a floor.
- **The critical path** adds the longest job in each group that ran side by side, and the
  sum of the rest. Time spent in the queue is left out. Queue time was large here: the
  partition went down on 19 Nov.

**D3 cost:** for the posterior arm alone (100 draws × 973 cells × 2 arms), the plug-in took
4,230.6 CPU-h with a 6.8 h critical path; the joint (KIM-95) took 4,272.3 CPU-h and 7.0 h.
These figures come from `sacct`. Item 8's `collect.md` gives 4,218.7 instead, because it
counts the chunks' timed seconds × threads.

## The posterior beside the joint

The full per-coordinate tables, with the Monte Carlo SE of each shift, are in
`compare/compare.md`. The merge's tables are in `run/modular_summary.md`.
`compare.md`'s shift and contraction columns match the merge's to 1e-10, and the script
checks this itself.

| | cut | plug-in |
|---|---|---|
| close (prereg §5) | **yes** | **no**: 0607 at +0.65 |
| max \|shift\| (joint SD) | 0.27 (0233; SE 0.08) | 0.65 (0607; SE 0.07) |
| mean shift over the 17 promoters | +0.01 (9 of 17 up) | +0.20 (16 of 17 up) |
| max \|Δ contraction\| | 0.006 | 0.009 |
| sd ratio, cut ÷ joint: median (range) | 1.00 (0.96–1.13) | 1.02 (0.90–1.12) |
| truth in the 90% interval | 19 of 21 (joint 20) | 19 of 21 (joint 20) |
| promoters' (mean − truth)/sd: mean, number negative | −0.79, 14 of 17 | −0.61, 14 of 17 |
| ENO/FBA–φ \|r\|: max, rms (joint 0.125, 0.053) | 0.074, 0.036 | 0.091, 0.031 |

**The Monte Carlo SE of a shift.** The cut's φ is constant within each stage-1 draw's
block of 10 inner draws, so its means rest on 200 distinct draws, not 2000. The joint's
rest on its ESS. With that, the plug-in's +0.65 on 0607 is about 2 SE past the 0.5 line,
and the cut's 0.27 is well inside it.

**The correlations the cut drops.** In the joint, ENO and FBA hardly correlate with φ: the
largest |r| over the 36 pairs is 0.125, and the rms is 0.053. In the cut these
correlations cannot exist, so its |r| is noise. Its floor is 1/√200 = 0.071, and its
largest |r| (0.074 cut, 0.091 plug-in) is at that floor. On the generous rung the
coupling the cut drops is weak, and that is why the cut sits close.

**The promoter offset is not from the bulk coupling.** KIM-92's joint put the promoters
about 0.8 SD below the truth (14 of 17 negative; Go entry, 2026-11-17). The cut's
stage 1 sees only the transcripts, and it gives the same offset: −0.79, with 14 of 17
negative. So the offset does not come from the bulk factor pulling on the promoters. It
is already in the transcript factor, on this dataset.

The plug-in's −0.61 is not "nearer the truth". It is the plug-in's upward bias (mean
+0.20 joint SD) laid over the same downward offset.

**Why the plug-in shifts.** Its stage 1 takes its transcription rates from the observed
panel, not from the mean-field solve. Item 2 measured those rates at the truth against
route (ii)'s: |Δ log k| has median 0.0145, 90th percentile 0.053, and maximum 0.112
(Slurm 21196497). The plug-in treats the noisy panel as if it were the true state, an
errors-in-variables cut, so a systematic bias in the rates moves the promoters. The upward
shift fits Q's median sitting about 1% high in D3. This mechanism is consistent with the
evidence but has not been tested separately.

## D3 scores: the plug-in beside the joint

This table is copied unchanged from `compare/compare.md`. It reads both archives'
`score.md`, written by KIM-95's scorer with no changes, on the same truth and prior files.

| target | truth ± SE | plug-in 5–50–95% | joint 5–50–95% | in 50/90: plug-in | joint | CRPS plug-in | joint | skill plug-in | joint |
|---|---|---|---|---|---|---|---|---|---|
| Q | 1.498e+07 ± 1.828e+04 | 1.486e+07 – 1.51e+07 – 1.532e+07 | 1.481e+07 – 1.501e+07 – 1.524e+07 | out/in | in/in | 6.631e+04 | 3.293e+04 | 0.9699 | 0.9851 |
| Q_kd | 1.504e+07 ± 1.755e+04 | 1.495e+07 – 1.519e+07 – 1.539e+07 | 1.49e+07 – 1.51e+07 – 1.532e+07 | out/in | in/in | 8.826e+04 | 4.155e+04 | 0.963 | 0.9826 |
| dQ | 6.423e+04 ± 2.32e+04 | -2.058e+04 – 8.035e+04 – 2.435e+05 | -2.714e+04 – 8.801e+04 – 2.183e+05 | in/in | in/in | 1.848e+04 | 1.857e+04 | 0.8696 | 0.869 |

| family | hit 50: plug-in | joint | hit 90: plug-in | joint | median species skill: plug-in | joint |
|---|---|---|---|---|---|---|
| pool | 39% | 34% | 89% | 92% | 0.9412 | 0.9452 |
| transcript | 76% | 80% | 100% | 100% | 0.9464 | 0.9514 |

- **Q and Q_kd:** the plug-in's median is about 1% above the truth. The truth falls outside
  its 50% band but inside its 90% band, and its CRPS is about twice the joint's. The skill
  against the prior stays above 0.96.
- **dQ, the knockdown's effect:** the plug-in scores the same as the joint. Both medians
  sit high by a similar amount, so the shared offset cancels in the difference.
- **The later horizon:** the plug-in is within a few points of the joint in both
  families, and median species skill differs by less than 0.005.
- **The D3 draws are item 7's:** `predict/prep.jls`'s `log_draws` equals
  `run/cut_draws_plugin.tsv` (Slurm 21628002). Item 8's scores therefore come from the
  posterior described above.

## What this means

- **For D7 and M10.** On the generous rung, the factor cut reproduces the joint but saves
  nothing, because its stage 1 still solves the ODE. The plug-in saves the solve but moves
  the promoters by up to 0.65 joint SD. Neither result gives a cheap and faithful modular
  arm on this model as designed. Before M10 runs the modular arm at every rung, these
  options are worth weighing:
  - a smaller stage 2 (K or the chain length);
  - correcting the plug-in for the panel's noise;
  - accepting the cut at the joint's price.

  M10 has left paper 1 (Decision log 2026-11-25), so this goes to the flagship paper.
- **For paper 1.**
  - The promoter offset is in the transcript factor, not the coupling. That bears on
    "promoter shift diagnosed or stated accurately".
  - The modular arm is a measured comparison: the cut is as good as the joint at the
    joint's price or more, and the plug-in is cheaper and biased.
  - The issue expected that "on the generous rung the exact transcripts pin the
    promoters, so the cut should sit close to the joint. The saving is unprofiled". That
    holds for the cut and does not hold for the plug-in. The saving is now measured.
