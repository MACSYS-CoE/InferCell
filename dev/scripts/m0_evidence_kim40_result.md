# KIM-40 and KIM-41: M0's crossed comparison and ptsG panel (M5)

**Run:** `dev/scripts/m0_evidence_kim40.sh` at `ed179c2` (main; the run's code is
`9310ddb`'s, checked before launch), into
`/fred/oz022/tkimpson/infercell-m5-kim40-2026-11-01/` (`SHA256SUMS`, `jobs.tsv`,
`sacct.tsv`, `audit.txt`, `score.md`, `score_datasets.tsv`). Archived with its inputs by
KIM-43, in `/fred/oz022/tkimpson/infercell-predict-archive-2026-11-06/runs/m0-kim40-kim41/`.
- **The design:** the M4 freeze (`93abbce`), KIM-69's scale (`5ec56d3`), and the
  scoring KIM-40 declared before the run (`9310ddb`), with KIM-70's docs pass
  (`ed179c2`).
- **Jobs:** precompile 19684097; prep 19684098; work 19684099 (1,150 tasks), and
  19827688 for task 587 (below); finish 19684100 (100 tasks); score 19684101.
- **Time and cost:** 95 h of wall-clock (1–5 November) and 12,978 CPU-h, at 150
  tasks at once. KIM-69 projected about 65 h and 9,807 CPU-h, from median per-cell
  times, and warned that the mean sets the total. Tasks' mean run times:

  | task | n | mean | CPU-h |
  |---|---|---|---|
  | δ | 100 | 12.5 h | 1,253 |
  | truth | 50 | 12.6 h | 630 |
  | hybrid block | 1,000 | 11.1 h | 11,088 |
  | finish | 100 | 0.06 h | 6 |

  No task came near its 24 h limit (the longest was 16.2 h). The largest MaxRSS was
  3.1 GiB, against 6 GB.

**Task 587 was resubmitted by hand.** Its first attempt started on gina12 at 19:53
on 2 November, and was `CANCELLED by 0` (root) 2 minutes later, with no log and no
output. No other job on gina12 was affected, and the cause is not known; the admin
logs are out of reach. `scontrol requeue` refused the task. At the user's request it
was resubmitted while the array ran, which the launcher's `missing` mode would have
refused until the queue cleared:
- job 19827688, on gina4, 11:27:34, completed long before `finish` started;
- the launcher's flags and the same tree, so the same task and seed, with gina12
  added to the excluded nodes;
- its id appended to `jobs.tsv` as a second `work` row, as `missing` would.

The block's rng depends only on (r, b), so this is the output the first attempt
would have written.

**The after-run checks** (spec §12, KIM-40) all pass:
- `progress`: work 1,150 of 1,150; finish 100 of 100; score done.
- `sacct.tsv` holds all six job ids. Every job completed but 587's first attempt, and
  none ran on a turin node.
- `dev/scripts/m0_audit.{jl,sh}` (job 20002599, on dave49), which writes `audit.txt`:
  all 2,402 records were read, each with `cpu_name` `znver3`, the run's commit and
  SMOKE off, and each work record names its own array task. They ran on 89 nodes, 13
  of them gina. Its command was
  `sbatch --chdir=<run tree> m0_audit.sh <outdir> <directory holding m0_audit.jl>`.
- `sha256sum -c inputs.sha256`: all 200 inputs unchanged.
- `SHA256SUMS` over the outputs: 3,662 files, checked straight after writing.
  Its own SHA256 is `29ca43c53b5eba87…`.

**The cells:** none failed, starved or fell below a floor. Clipped at least once:
0.77% of the hybrid's baseline cells under V7 and 0.75% under B7, and 0.41% and
0.40% under the knockdown. The remainder ran clean.

## How to read the coverage

The references are the ones KIM-40 declared before the run. A calibrated arm's
intervals, used as they are with 100 draws, cover **88.2% at nominal 90% and about
49% at nominal 50%**, with SEs over 50 datasets of about 4.6 and 7.1 points. V7's
thinned chain pulls these down slightly. The hybrid's MC term pushes its own
coverage up, by KIM-69's arithmetic to about 92% at a 10% MC share, before the
finite-sample effect.

## Q = `lactate_out` at 600 s (KIM-40, primary)

| arm | cover 50% | cover 90% | median 90% width | width / V7 hybrid | shift from V7 hybrid (SD) | z, mean ± SD |
|---|---|---|---|---|---|---|
| V7 hybrid | 62.0 ± 6.9 | 92.0 ± 3.8 | 36,900 | 1 | 0 | 0.03 ± 0.94 |
| V7 mean field | 58.0 ± 7.0 | 94.0 ± 3.4 | 36,700 | 0.988 [0.937, 1.02] | 0.136 | 0.17 ± 0.96 |
| V7 mean field + δ | 56.0 ± 7.0 | 90.0 ± 4.2 | 36,700 | 0.988 [0.937, 1.02] | 0.005 | 0.04 ± 0.99 |
| B7 hybrid | 50.0 ± 7.1 | 88.0 ± 4.6 | 31,100 | 0.823 [0.761, 0.93] | −0.110 | −0.09 ± 1.11 |
| B7 mean field | 54.0 ± 7.0 | 82.0 ± 5.4 | 28,700 | 0.772 [0.716, 0.892] | 0.017 | 0.08 ± 1.17 |
| B7 mean field + δ | 48.0 ± 7.1 | 84.0 ± 5.2 | 28,700 | 0.772 [0.716, 0.892] | −0.089 | −0.08 ± 1.21 |

- **The reference and the route.** V7 hybrid covers within 1 SE of its reference at
  90%, and is high at 50% (62%, 1.8 SE above 49%). B7 hybrid is on the reference at
  both, with intervals 18% narrower than V7's.
- **The forward-description effect** (V7 mean field against V7 hybrid): a shift of
  0.14 SD and the same width. δ removes the shift (0.005 SD).
- **The inference-route effect** (B7 hybrid against V7 hybrid): a shift of −0.11 SD
  and a width of 0.82.
- **The cheap deployed route** (B7 with the mean field) has 82–84% at 90%, about 1 SE
  below the reference, with a z SD of 1.17–1.21 and 8 of 50 truths in the lowest
  PIT bin (about 5 expected). Its shift without δ, 0.017 SD, is the two effects nearly cancelling, not
  an absence of either.
- **Fifty datasets cannot separate these arms at 90%:** all six are within 1.4 SE of
  88.2%.

## ΔQ under the ptsG knockdown, dose 0.3 (KIM-41)

| arm | cover 50% | cover 90% | median 90% width | width / V7 hybrid | shift from V7 hybrid (SD) | z, mean ± SD |
|---|---|---|---|---|---|---|
| V7 hybrid | 62.0 ± 6.9 | 98.0 ± 2.0 | 10,200 | 1 | 0 | 0.04 ± 0.74 |
| V7 mean field | 46.0 ± 7.0 | 90.0 ± 4.2 | 7,890 | 0.733 [0.69, 0.764] | −0.064 | −0.04 ± 1.00 |
| V7 mean field + δ | 48.0 ± 7.1 | 88.0 ± 4.6 | 7,890 | 0.733 [0.69, 0.764] | 0.027 | 0.09 ± 1.06 |
| B7 hybrid | 52.0 ± 7.1 | 96.0 ± 2.8 | 10,300 | 0.964 [0.858, 1.03] | 0.087 | 0.10 ± 0.78 |
| B7 mean field | 32.0 ± 6.6 | 86.0 ± 4.9 | 6,760 | 0.635 [0.597, 0.715] | 0.027 | 0.05 ± 1.15 |
| B7 mean field + δ | 26.0 ± 6.2 | 78.0 ± 5.9 | 6,760 | 0.635 [0.597, 0.715] | 0.147 | 0.16 ± 1.27 |

**The effect itself is resolved.** The true ΔQ has a median of −20,400 (range −79,500
to −3,500), and each truth is at least 11 SEs from zero (median 27). Every arm's 90%
interval excludes zero on all 50 datasets, with the right sign. *(This read-off from
`score_datasets.tsv` was not declared before the run.)*

**The rules are not met, so these coverages are not clean calibration numbers.** See
the next section. The intervals' spread and the truth's error are both dominated by
cell noise, in different ways per arm:
- **The hybrid arms are conservative.** At V7's median MC share of 0.503, half of the
  sample's variance is cell noise. Arithmetic from the table below, not a measurement:
  with the truth's SE at 0.295 of the parameter SD, a normal predictive gives a z SD of
  √((1 + 0.295²)(1 − 0.503)) = 0.735 and a coverage of about 97%. The run shows
  0.737 and 98%. B7's figures give 0.71 and about 98%; the run shows 0.78 and 96%. So
  98% and 96% say the intervals are too wide because of the noise, not that the
  posterior is too wide.
- **The width ratios are mostly that noise too.** With V7 hybrid's noise taken out
  (a factor of √(1 − 0.503) = 0.705), V7 mean field's 0.733 becomes about 1.04.
- **The mean-field arms carry no MC term, so the truth's error counts against them.**
  V7's cover at the reference. δ adds its own error to ΔQ (its SE is about 0.3 of the
  parameter SD), and its correction is small (the mean field's shift is 0.06 SD), so
  here it costs more than it gives.
- **B7 mean field** is the one arm clearly low: 32% and 26% at nominal 50%, 2.4 and 3.2
  SE below 49%, and 78% with δ at nominal 90%, 2.2 SE below 88.2%. Its PIT is U-shaped.
  The noise explains part of its z SD (about 1.05 without δ and 1.11 with, by the same
  arithmetic, against 1.15 and 1.27). The rest cannot be separated from the route's
  narrower intervals at these sizes.

## The 10% rules

KIM-69 set M0's sizes so that the MC term is at most 10% of the spread, and the
truth's and δ's SEs at most 10% of the predictive SD. It sized them on Q, on one
held-out posterior, and said they would first be seen across the 50 here.

| target | posterior | MC share, median (max; datasets > 10%) | truth SE / SD, median (max; > 10%) | δ SE / SD, median (max) |
|---|---|---|---|---|
| Q | V7 | 0.074 (0.098; 0) | 0.083 (0.099; 0) | 0.085 (0.099) |
| Q | B7 | 0.100 (0.202; 25) | 0.102 (0.156; 27) | 0.100 (0.151) |
| ΔQ | V7 | 0.503 (0.866; 50) | 0.295 (0.743; 50) | 0.301 (0.771) |
| ΔQ | B7 | 0.538 (0.708; 50) | 0.321 (0.529; 50) | 0.323 (0.464) |

- **Q under V7 meets the rules on all 50.** Under B7 the medians sit on the limit,
  since B7's posteriors are narrower than V7's.
- **ΔQ fails them on all 50, under both.** Pairing removes little of the cells' spread
  (√E[σ̂²] is 71,000 for ΔQ and 93,600 for Q), while the spread over θ, net of that
  noise, is about 5× smaller. The freeze's sizes were Q's.
- **The cost of meeting them on ΔQ,** at the median dataset: about 9× the cells per
  draw (n = (1 − f)/f · σ²/Var_θ gives about 8,100) and about 9× the truth and δ
  cells (about 85,000). At this run's rates that is on the order of 100,000 CPU-h.

**Not rerun: the user's decision, 2026-11-05** (spec §12). ΔQ's limit is stated in
prose beside its numbers.

## What this supports

- **On Q,** the exact posterior with the hybrid, and the mixed posterior with the
  hybrid, cover at the declared reference on M0, and B7's intervals are about 18%
  narrower. The mean field's shift is small and δ removes it. The B7–mean-field route
  sits about 1 SE low.
- **On ΔQ,** the panel shows that a model fitted to unperturbed data predicts the
  knockdown's effect with the right sign and scale on all 50 datasets. It does not
  give a clean calibration of those intervals: the hybrid's are widened by cell noise,
  and the truth's error is not small against the mean field's.
- **It does not** test the carry-over into the reduced model (a declared assumption, KIM-47).

The readings above are after the run, but no primary score was changed or re-scored.
