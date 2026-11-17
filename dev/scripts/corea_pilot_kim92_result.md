# KIM-92: the reduced-model pilot (M8's kill check). PASSES

**The kill check passes on the frozen gates** (`prereg/kim-92`, commit d2703cf). The
convergence gate failed on the first 200 draws per chain, and passed after the one
extension the pre-registration allows. The PPC passed with 11 of 884 rows outside the
band, against a limit of 44. Cost had no cap: it came to about 492 CPU-h, and the chains
took 67–81 h of wall-clock each. Whether that is affordable is KIM-93's call.

**Run:** `dev/scripts/corea_pilot_kim92.sh` at `de24b27`, from the detached worktree
`.worktrees/kim92-run-de24b27`. Each output records the commit it was made at:
- `dataset*.jls`: `a1be360`;
- `profile.jls` and `profile.md` (the re-profile): `4e52e37`;
- `prep.jls`: `16cc1a7`, which has the same prep code as `de24b27`;
- everything else: `de24b27`.

Those commits are on this PR's branch, which merges with a merge commit so they stay
reachable. After the run, the driver changed only in comments and the launcher's usage
message: `git diff de24b27 HEAD -- dev/scripts/corea_pilot_kim92.*`.

The output directory is `/fred/oz022/tkimpson/infercell-kim92-posterior-2026-11-10/`, holding
`pilot_summary.md`, `pilot_coordinates.tsv`, `pilot_ppc.tsv`, `pilot_draws.tsv`,
`gate_200.md`, `gate_400.md`, `metric.md`, `probe.md`, `profile.md`, `jobs.tsv`,
`sacct.tsv` and `SHA256SUMS`.

**Archives:**
- **The final run:** `/fred/oz022/tkimpson/infercell-kim92-archive-2026-11-17/` is a copy of
  that directory. Its `SHA256SUMS` covers 48 files and has the SHA-256 `9ec9675e31a30aed…`.
  All 48 pass `sha256sum -c` in both the original and the copy.
- **The first dataset and profile:** the archive's `pilot-2026-11-10/` subdirectory is a
  copy of `/fred/oz022/tkimpson/infercell-kim92-pilot-2026-11-10/`. It holds:
  - jobs 20454165–67 (precompile, dataset, first profile);
  - their logs and `sacct.tsv`;
  - the original dataset.

  Its own `SHA256SUMS` covers 9 files and has the SHA-256 `18828fbaf91d8aeb…`.
- **Before the extension:** the chains and checkpoints were copied to
  `/fred/oz022/tkimpson/infercell-kim92-archive-2026-11-15/`, with their own `SHA256SUMS`.

**Design:** the design is the pre-registration's, with one change: chains start from a
warm diagonal metric and checkpoint (Decision log 2026-11-12, made before any chain
started). The gates did not change.

**Where the numbers come from.**
- **The gates, recovery, contraction and cost table** come from `pilot_summary.md`, which
  merge 21052425 wrote.
- **Derived numbers** come from `dev/scripts/corea_pilot_kim92_reads.jl`, a read-only
  script with no model evaluation. That covers:
  - the per-chain means and halves;
  - the sampler statistics;
  - the conditioning numbers;
  - the Laplace comparison;
  - the standardised errors and θ_ref against the truth;
  - the seconds per gradient.

  It reads `chain_c.jls`, `extend_c.jls`, `metric.jls`, `prep.jls` and
  `pilot_coordinates.tsv`. It was run on the login node against the archive, so these
  numbers have no job id:
  `julia --project dev/scripts/corea_pilot_kim92_reads.jl /fred/oz022/tkimpson/infercell-kim92-archive-2026-11-17`.
  Each section of its output is numbered, and the note cites a section as (reads §n).
- **Single iterations' wall times** come from the progress lines in
  `logs/kim92_chain<c>_<job>.log`, which print one line every 5 iterations. The chain files
  do not store them.
- **Job states, CPU-h, wall times and MaxRSS** come from `sacct.tsv`. For jobs 20454165–67,
  they come from `pilot-2026-11-10/sacct.tsv` in the archive.

## The gates

| gate | criterion | result | verdict | job |
|---|---|---|---|---|
| convergence, first 200 draws per chain | R̂ ≤ 1.02 and bulk ESS ≥ 100 on all 21 | worst R̂ 1.0284 (log_S_0213), then 1.0256 (log_kcatF_R_ENO); smallest ESS 1190.7 | fail, so extend | check 20660937 |
| convergence after the one extension, 400 per chain | as above | worst R̂ 1.0102 (log_sigma_b), smallest ESS 2301.6 (log_krnadeg), 0 failing | **pass** | check 20968978 |
| PPC | at most 44 of 884 rows outside the central 99% band | **11 outside** | **pass** | PPC 21052424 |

**Why the first gate failed.** No chain is stuck (reads §1). On the two failing
coordinates, the four chains' means over draws 201–400 agree to within about 0.1
posterior SD:
- log_S_0213's means lie in 1.0158–1.0192, against a per-chain SD of about 0.03;
- ENO's means lie in 3.9474–3.9525, against an SD of about 0.045.

The 100-draw halves of each chain spread a little wider, at about 0.17–0.19 SD (0213:
1.0145–1.0197; ENO: 3.9472–3.9556). That is what half as many draws would give, and no
half drifts away from the others. The per-chain SDs differ by up to about 25%: ENO's run
from 0.037 to 0.048. Rank-normalised R̂ is the larger of the bulk and the folded R̂, so it
is sensitive to that difference in scale at 200 draws. After the extension, the two
coordinates have R̂ 1.0037 and 1.0041.

**Checklist item 5 (the identifiability diagnosis) was not triggered,** because the gate
passes after the extension. The quantities it asks for are reported anyway:
- **Contraction:** at least 0.919 on every coordinate (below).
- **Correlations:** one posterior correlation has |r| ≥ 0.5, between log_S_0607 and
  log_krnadeg (0.69).
- **Weak directions:** none. No marginal sits at its prior.

**The PPC's 11 rows** are rows 60, 169, 348, 430, 490, 495, 599, 634, 743, 828 and 855 of
`pilot_ppc.tsv`. They spread across the record and do not cluster. At 884 rows and a 99%
band, about 9 would be expected outside. As the pre-registration says, this tests whether
the fitted likelihood reproduces its own data. It does not test the mean-field
approximation.

## Recovery and contraction (reported, not gated)

From `pilot_summary.md` (merge 21052425), on 1,600 pooled kept draws: 4 chains × (200 + 200).

| coordinate | prior sd | posterior mean ± sd | truth | 90% interval | truth inside | (mean − truth)/sd | contraction |
|---|---|---|---|---|---|---|---|
| log_S_JCVISYN3A_0445 | 0.693 | 0.6383 ± 0.0349 | 0.6339 | [0.5803, 0.6980] | yes | +0.13 | 0.997 |
| log_S_JCVISYN3A_0220 | 0.693 | -0.4471 ± 0.0540 | -0.4240 | [-0.5358, -0.3584] | yes | -0.43 | 0.994 |
| log_S_JCVISYN3A_0131 | 0.693 | 1.5024 ± 0.0209 | 1.5208 | [1.4688, 1.5370] | yes | -0.88 | 0.999 |
| log_S_JCVISYN3A_0727 | 0.693 | 1.1106 ± 0.0232 | 1.1465 | [1.0727, 1.1487] | yes | -1.55 | 0.999 |
| log_S_JCVISYN3A_0607 | 0.693 | 4.3383 ± 0.0082 | 4.3459 | [4.3248, 4.3516] | yes | -0.93 | 1.000 |
| log_S_JCVISYN3A_0606 | 0.693 | 0.8095 ± 0.0322 | 0.7774 | [0.7573, 0.8609] | yes | +1.00 | 0.998 |
| log_S_JCVISYN3A_0729 | 0.693 | 1.1870 ± 0.0293 | 1.2156 | [1.1373, 1.2334] | yes | -0.97 | 0.998 |
| log_S_JCVISYN3A_0213 | 0.693 | 1.0169 ± 0.0300 | 1.0081 | [0.9693, 1.0656] | yes | +0.29 | 0.998 |
| log_S_JCVISYN3A_0221 | 0.693 | 0.7459 ± 0.0348 | 0.7606 | [0.6887, 0.8011] | yes | -0.42 | 0.997 |
| log_S_JCVISYN3A_0475 | 0.693 | 1.4902 ± 0.0214 | 1.5252 | [1.4548, 1.5278] | yes | -1.63 | 0.999 |
| log_S_JCVISYN3A_0233 | 0.693 | 1.3704 ± 0.0288 | 1.4056 | [1.3226, 1.4163] | yes | -1.23 | 0.998 |
| log_S_JCVISYN3A_0234 | 0.693 | 1.1717 ± 0.0201 | 1.2079 | [1.1380, 1.2043] | **no** | -1.80 | 0.999 |
| log_S_JCVISYN3A_0694 | 0.693 | 0.2139 ± 0.0246 | 0.2505 | [0.1753, 0.2548] | yes | -1.48 | 0.999 |
| log_S_JCVISYN3A_0779 | 0.693 | 1.4207 ± 0.0309 | 1.4698 | [1.3707, 1.4721] | yes | -1.59 | 0.998 |
| log_S_JCVISYN3A_0651 | 0.693 | -0.8888 ± 0.0545 | -0.8500 | [-0.9792, -0.7960] | yes | -0.71 | 0.994 |
| log_S_JCVISYN3A_0344 | 0.693 | -0.7997 ± 0.0519 | -0.7545 | [-0.8863, -0.7161] | yes | -0.87 | 0.994 |
| log_S_JCVISYN3A_0203 | 0.693 | 0.3076 ± 0.0346 | 0.3302 | [0.2517, 0.3667] | yes | -0.65 | 0.998 |
| log_krnadeg | 0.693 | 1.4534 ± 0.0067 | 1.4514 | [1.4424, 1.4644] | yes | +0.29 | 1.000 |
| log_kcatF_R_ENO | 0.157 | 3.9503 ± 0.0446 | 3.9502 | [3.8758, 4.0254] | yes | +0.00 | 0.919 |
| log_kcatF_R_FBA | 0.558 | 4.1110 ± 0.0891 | 4.2182 | [3.9671, 4.2554] | yes | -1.20 | 0.974 |
| log_sigma_b | 1.000 | -1.0100 ± 0.0234 | -1.0306 | [-1.0496, -0.9713] | yes | +0.88 | 0.999 |

- **Truth inside the 90% interval:** 20 of 21. There is one truth, so this is not a
  coverage rate (prereg §3).
- **The data reach** (contraction ≥ 0.5) **all 21 of 21.** The least-reached coordinates
  are kcatF_R_ENO (0.919) and kcatF_R_FBA (0.974), which are also the two with the
  tightest priors.
- **Standardised errors:** Σz² over the 21 is 22.6. That is unremarkable for 21
  coordinates taken one at a time.

**But the promoter strengths sit below the truth.**
- Of the 17 promoter strengths, 14 have posterior means below the truth. Their mean
  standardised error is −0.81, with an SD of 0.76.
- If the 17 errors were independent and calibrated, their mean would be 0 ± 0.24, so this
  is about 3.3 of those SEs away from zero.
- θ_ref, the MAP of the same likelihood, has the same shift: its mean (MAP − truth)/sd
  over the promoters is −0.86 (reads §5).
- So the shift comes from the likelihood on this dataset, not from the sampler.

It is not gated and does not change the verdict. One dataset cannot separate two causes:
- **noise in this dataset** shared by every promoter, for example through the common bulk
  scale;
- **a bias** from fitting the mean-field bulk likelihood to data from the stochastic
  simulator.

The second is exactly what KIM-95's scored predictions guard against (prereg §2.2,
§4). A second dataset seed would tell the two causes apart, but it is not part of this
issue. **For KIM-93:** read the promoters' predictions with this shift in mind.

**The posterior is close to Gaussian.** For each coordinate, the Laplace SD at θ_ref,
√diag((−H)⁻¹) from `metric.jls` (job 20651329), is between 0.960 and 1.083 of the
posterior SD of the 1,600 pooled draws (reads §4). The largest ratios are 1.083 for ENO
and 1.065 for FBA. The one strong correlation,
0607–krnadeg, is 0.669 in (−H)⁻¹ and 0.685 in the draws (reads §4).

ENO's and FBA's MAP values differ from their posterior means by 2.3 and 1.4 posterior SD.
Their posteriors are the least Gaussian, or the MAP stopped short along them; this run
cannot say which. For the other 19 coordinates, the posterior mean is within 0.35
posterior SD of θ_ref. The largest of those are log_S_0233 (0.35) and log_sigma_b (0.33)
(reads §5).

## The sampler

From the chain files (reads §2). All four chains ran on dave306, and the extensions on
dave305.

| chain | job | extension job | adaptation: mean steps per iteration, numerical errors, max depth | kept 201–400: mean steps, max depth | extension: mean steps, max depth | step size |
|---|---|---|---|---|---|---|
| 1 | 20660933 | 20968974 | 13.8, 9, 7 | 13.3, 4 | 13.9, 4 | 0.374 |
| 2 | 20660934 | 20968975 | 16.7, 11, 7 | 14.7, 4 | 14.8, 4 | 0.320 |
| 3 | 20660935 | 20968976 | 17.1, 12, 9 | 13.7, 4 | 14.0, 4 | 0.345 |
| 4 | 20660936 | 20968977 | 18.9, 6, 10 | 14.8, 4 | 14.8, 4 | 0.324 |

- **No numerical errors after adaptation.** All 1,600 kept iterations ran at depth 3 or 4.
- **Two pathological iterations**, both in adaptation under the warm metric (before
  iteration 100). The depth and steps come from the chain files, as each chain's deepest
  iteration (reads §2). The times come from the 5-iteration blocks of progress lines that
  hold each iteration:
  - chain 4, iteration 59: depth 10, 1,023 steps. Iterations 56–60 took 11.5 h (log 7.7 →
    19.2 h).
  - chain 3, iteration 87: depth 9, 731 steps. Iterations 86–90 took about 8.5 h.
- **The windowed metrics** (from iteration 101) shortened the step size.
  - Logged step sizes, one every 5 iterations, ran 0.19–1.23 (median about 0.6) over
    iterations 5–100, and 0.076–0.87 (median about 0.3) over 105–200.
  - Leapfrog steps per 25-iteration block (reads §2), across the 4 chains:
    - iterations 101–200: 295–588, mean 432;
    - iterations 1–100: 195–611, mean 301, leaving out the two blocks that hold the deep
      iterations (929 and 1,219);
    - the kept draws: 319–375.
  - So the second half of adaptation cost about 1.4 times the first. The re-estimate at
    iteration 150 brought no clear drop: chain 1's blocks fell, but chain 2's and chain 4's
    did not.
- **Against the probe.** The unit-metric probe (job 20545099, `probe.md`) took about 51
  steps per iteration at step sizes of 0.001–0.012. It ran 38 iterations in 24.22 h,
  including its setup, which projects to about 10–11 days per 400-iteration chain. The
  conditioning numbers (reads §3),
  with −H at θ_ref from `metric.jls`:
  - under the unit metric, κ(−H) = 528 (√κ 23);
  - under the warm diagonal metric, κ(D½ (−H) D½) = 5.98 (√κ 2.44), with D = diag(minv).
- **Gradient cost** (reads §6): each job's sampling seconds over its total leapfrog steps.
  That comes to 42.8–44.8 s per step for the chains on dave306, and 46.7–48.4 s for the
  extensions on dave305. The re-profile measured 40.3 s (job 20541310, `profile.md`). The
  per-step figure includes adaptation and Turing's per-iteration log-density evaluation.
  The chains' sampling timer may also include part of the 0.9 h suspend on dave306 (see
  Cost). Chain 1's log gives 66.2 h of sampling against 66.06 h Elapsed. That could raise
  their figure by up to about 1.3%, roughly 0.5 s.

## Cost

CPU-h come from `sacct` CPUTimeRAW. Every job is single-threaded except the PPC, which ran
on 16 threads.

| stage | jobs | CPU-h | wall-clock |
|---|---|---|---|
| precompile, first pass | 20454165 | 0.01 | 1 min |
| dataset (drawn at `a1be360`, copied by sha256; `DATASET_PROVENANCE.md`) | 20454166 | 4.51 | 4 h 31 min |
| first profile, before KIM-99 | 20454167 | 0.80 | 48 min |
| precompile, re-profile | 20541309, 20541310 | 0.43 | 26 min |
| prep: MAP 9,721 s and an 800-cell reference ensemble 16,332 s | 20545098 | 7.28 | 7 h 17 min |
| depth probe (unit metric, stopped at its 24 h budget) | 20545099 | 24.33 | 24 h 20 min |
| warm metric (−H by finite differences of the gradient) | 20651329 | 0.70 | 42 min |
| 4 chains, 200 + 200 | 20660933–36 | 299.09 | 67.0, 77.9, 76.8, 81.1 h (Start → End; see below) |
| check on 200 draws | 20660937 | 0.06 | 3 min |
| extension, 4 × 200 | 20968974–77 | 152.63 | 36.6, 38.5, 37.7, 39.8 h |
| check on 400 draws | 20968978 | 0.05 | 3 min |
| PPC, 200 draws × 10 replicates | 21052424 | 1.99 | 7.5 min |
| merge | 21052425 | 0.04 | 3 min |
| **total** | | **491.9** | |

- **Totals:** the merge's own total, over `jobs.tsv`, is 486.60 CPU-h. Jobs 20454165–67
  add 5.32, which gives 491.9.
- **The chains' wall-clock** is sacct's Start → End. Their Elapsed (66.1, 77.0, 75.9 and
  80.2 h) is about 0.9 h shorter for all four on dave306, which looks like a node-level
  suspend. The extensions' Elapsed equals Start → End. CPU-h follow Elapsed.
- **Peak memory** (sacct MaxRSS of the batch step): 5.0–5.5 GiB for the chains, and
  4.7–6.0 GiB for the extensions. The job report's "4.2 GB peak" for chain 1 is a
  different measure.
- **Wall-clock from the first submission to the last job's end:** 2026-11-10 01:47 to
  2026-11-17 14:24, which is 7 d 12.6 h. That includes the time spent on the probe and on
  the decision to use the warm metric.
- **Critical path of the design as it finally ran** (dataset, prep, metric, chains,
  extension, PPC and merge, without queue waits or the probe): about 133 h, or 5.5 days,
  on Elapsed. The chains and the extension make up 120 h of that.

## Regenerating the derived numbers

`dev/scripts/corea_pilot_kim92_reads.jl` prints every number above marked (reads §n). Run
it from the repo root against the archive:

    julia --project dev/scripts/corea_pilot_kim92_reads.jl /fred/oz022/tkimpson/infercell-kim92-archive-2026-11-17

Its header lists what each section computes. It deserialises the run's `.jls` files, so it
needs the Manifest those files were written with, which this PR does not change. On a login
node it takes about a minute, most of it loading packages.

The session also made one check that is not kept: an R̂ and ESS preview on 185 draws per
chain, made before check 20660937. It flagged the same two coordinates as the check.
