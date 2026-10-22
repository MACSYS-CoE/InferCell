# 16a.10, F10: a production-length chain is not resolved from the M0 reference

The rule was fixed before any divergence was computed (`spec/phases/16-recovery.md`
§12 2026-10-22). The run is approved in §12 2026-10-13.

**The runs.**
- **Production:** chain 5 on V6's dataset (M0 replicate 16085, 50 cells, bulk
  metabolites, PG at N = 20), started from the prior with `Xoshiro(16_090_005)`.
  - Job 18751635 at `a0060de`, from `.worktrees/f10-16a10`: 656 sweeps at 5
    threads, in 4 d 5.5 h, a mean of 557 s per sweep (median 550).
  - Smoke job 18751515 ran first.
  - The checkpoint and log are archived, with checksums, at
    `/fred/oz022/tkimpson/infercell-16a10-archive-2026-10-22/`.
- **Reference:** V6's chains 1 to 4 (`v6_16a9b_result.md`), at 1,650 sweeps each.
- **D16.5's case:** the PG case. Production and reference both run PG at N = 20, so
  the particle count costs mixing, not bias. F10 measures the Monte Carlo and
  mixing error of production's *length*.

**Reproduce:**
```
julia --project dev/scripts/f10_16a10_export.jl \
    /fred/oz022/tkimpson/infercell-16a9b-archive-2026-10-17/v6/gibbs_pilot_16a9a \
    /fred/oz022/tkimpson/infercell-16a10-archive-2026-10-22
/fred/oz022/tkimpson/miniconda3/bin/python3 dev/scripts/f10_16a10.py
```
- **Where the inputs are:** both are archives with checksums. V6's chains are
  byte-identical to the copies in `.worktrees/v6-16a9b`.
- **Python:** `/fred/oz022/tkimpson/miniconda3/bin/python3`, with numpy 2.2.6,
  scipy 1.15.3, matplotlib 3.10.5 and arviz 0.22.0. The system `python3` has no
  scipy.
- **What it ran on:** the login node. The export refuses unless all five chains
  share the replicate, the cells, N, the truth, σ_b and the parameter names.
- **What it writes:** `f10_16a10_draws.tsv` and `f10_16a10_truth.tsv`. The figures
  go to `dev/notes/figures/f10-16a10/`.
- **Check on the inputs:** the export reproduces V6's merge. The reference means
  and SDs below match `v6_16a9b_result.md` to every printed digit.

## The posteriors

Reference: sweeps 101 to 1,650 of chains 1 to 4, which is 6,200 draws.
Production: sweeps 101 to 656 of chain 5, which is 556 draws. All on the log scale.

| Parameter | Reference mean ± SD | Production mean ± SD | Δmean / SD_ref | SD ratio | Production ESS | Null chains' ESS |
|---|---|---|---|---|---|---|
| GAPD promoter, `S_JCVISYN3A_0607` | 1.545 ± 0.116 | 1.556 ± 0.137 | +0.09 | 1.19 | 63 | 111, 108, 32, 105 |
| ptsG promoter, `S_JCVISYN3A_0779` | 0.700 ± 0.199 | 0.692 ± 0.182 | −0.04 | 0.91 | 287 | 307, 382, 347, 325 |
| `krnadeg` | 2.225 ± 0.076 | 2.226 ± 0.077 | +0.02 | 1.01 | 75 | 122, 175, 63, 145 |
| `kcatF_R_ENO` | 4.242 ± 0.042 | 4.239 ± 0.041 | −0.09 | 0.98 | 55 | 175, 152, 135, 177 |
| `kcatF_R_FBA` | 4.408 ± 0.184 | 4.386 ± 0.190 | −0.12 | 1.03 | 491 | 419, 487, 404, 560 |
| `sigma_b` | −2.625 ± 0.060 | −2.621 ± 0.064 | +0.06 | 1.06 | 328 | 389, 406, 316, 364 |

**How the ESS is computed:** bulk ESS (arviz, rank-normalised and split, the kind
V6 and V7 report), for each single chain over sweeps 101 to 656.
- **Production on GAPD:** 63, below the 97 that §12 2026-10-13 projected from V6's
  pooled rate. The matched nulls give 32 to 111, so chain 5 is within their
  spread.
- **Production on ENO:** 55, against the nulls' 135 to 177.
  - Chain 5 wanders slowly on ENO: its two halves differ by 0.43 reference SD.
  - An unsplit Geyer estimate, used in this file's first version, read 165,
    because it misses that drift.
  - ENO's divergence still ranks 3 of 5.

## The divergence

W1 is the Wasserstein-1 distance against the reference, in reference SDs. The
matched nulls are each V6 chain's sweeps 101 to 656, against the other three
chains' sweeps 101 to 1,650 (4,650 draws). Production is read against all four
chains (6,200 draws).
- **The asymmetry.** The smaller reference adds a little noise to the nulls, which
  leans toward "not resolved".
- **Its size, checked:** production against each three-chain subset gives W1 of
  0.139 to 0.159 on GAPD, 0.115 to 0.127 on FBA, and 0.041 to 0.068 on `krnadeg`.
  FBA stays above every null under every subset, so the reading does not move.

| Parameter | Production | Matched nulls | V8 split-half floor | Production's rank of 5 | Reading |
|---|---|---|---|---|---|
| GAPD promoter | 0.152 | 0.084, 0.117, 0.120, 0.186 | 0.100 | 4 | not resolved |
| ptsG promoter | 0.109 | 0.044, 0.040, 0.052, 0.115 | 0.027 | 4 | not resolved |
| `krnadeg` | 0.041 | 0.110, 0.033, 0.161, 0.210 | 0.096 | 2 | not resolved |
| `kcatF_R_ENO` | 0.091 | 0.082, 0.048, 0.094, 0.118 | 0.058 | 3 | not resolved |
| `kcatF_R_FBA` | 0.120 | 0.043, 0.100, 0.056, 0.078 | 0.034 | 5 | above all four, by 0.020 |
| `sigma_b` | 0.075 | 0.042, 0.038, 0.079, 0.076 | 0.067 | 3 | not resolved |

**Reading, under the rule.**
- **Chain 5 ranks 5 of 5 on one parameter, `kcatF_R_FBA`.** About 1.2 of six are
  expected by chance.
- **Its margin over the largest null is 0.020 SD,** against nulls spread over 0.043
  to 0.100. The rule asks for a margin that is large against that spread, or for
  several parameters together, and neither holds.
- **So no parameter's divergence is a finding.** A production-length chain is
  **not resolved** from the reference on any marginal.
- **On FBA, the mean is 0.12 SD low.** That is FBA's posterior SD × 0.12, about
  0.022 on the log scale, which is about 2% in `kcat`.

**What V8's floor alone would have said.**
- **Chain 5** exceeds the split-half floor on five of the six parameters; only
  `krnadeg` is below it.
- **The nulls** are production-length chains by construction, and 19 of the 24
  exceed it too.
- **That is the failure the §12 entry anticipated:** the floor measures the noise
  of 3,100-draw halves, not of a 556-draw chain.

**What this does not show.**
- **Low resolving power.** One production chain is one draw of the divergence, and
  the nulls are four draws. A production-length chain carries about 0.04 to 0.2 SD
  of W1 from Monte Carlo noise alone, and a bias below that is invisible here.
- **The test is not independent of V6.** The matched nulls come from V6's own
  chains. They are exchangeable with chain 5 by construction, sharing the start
  rule, the code path and the burn-in, but they are not an independent
  reference.
- **One dataset and one truth.** V7 covers calibration across datasets, and this
  does not.
- **Mixing.** On the GAPD promoter, `krnadeg` and, in chain 5, ENO, a production
  chain can hold as few as 30 to 75 effective draws. The nulls range 32 to 177 on
  those three. That is the noise this comparison sees, and the noise 16b's
  single-chain fits would carry.

## Figures

In `dev/notes/figures/f10-16a10/`, as PDF and PNG:
- **`f10_marginals`:** each marginal as a KDE, with the reference filled, chain 5 in
  red, each null chain's sweeps 101 to 656 in grey, and the truth dashed.
- **`f10_contour`:** ptsG × ENO, with the reference's and chain 5's 50% and 90%
  highest-density contours.
  - Each KDE uses its own Scott bandwidth, and chain 5 has fewer draws, so its
    contours are drawn slightly wider: about 3% in SD.
- **`f10_divergence`:** each marginal's W1, beside the matched nulls and V8's floor.

## Compute

- Job 18751635: 5 CPUs × 4-05:30:05, which is 507.5 CPU-h.
- Smoke 18751515: 0.3 CPU-h.
- Phase 16 is at about 42.6k of 50k (`phase16_ledger.md`).
