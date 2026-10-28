# Phase 16 at M0: what the result licenses, and what it does not (16b.7)

Phase 16 closes at the M0 rung (`spec/phases/16-recovery.md` §12 "16a.11", call 2).
**The claim is only as large as the rung it was measured on** (16-recovery §1).
- **M0 is Core A′ cut to two genes over 600 s** (D16.4): ptsG (`JCVISYN3A_0779`)
  and GAPD (`JCVISYN3A_0607`).
- **The other 15 proteins are held at their proteomics counts** as fixed enzymes,
  a labelled reduction.
- **Each dataset has 50 cells,** a number set by the cap, not by D11's precision
  argument.
- **Five of D11's six targets are free,** plus σ_b. `S_PGI` is not in M0.

## Every result, with its rung and label

| Result | Task | Rung, cells | What it shows | Label | Source |
|---|---|---|---|---|---|
| V6, the reference posterior | 16a.9b | M0, 50 cells, one dataset | Converges at 4 × 1,650 sweeps; every truth inside its central 95% | M0 | `v6_16a9b_result.md` |
| V7, SBC | 16a.9b | M0, 50 replicates × 50 cells | 6 of 6 within their rank-ECDF bands | M0; 23 of 50 replicates below bulk ESS 100 | `v7_16a9b_result.md` |
| F10, production against the reference | 16a.10 | M0, one production chain | Not resolved from the reference on any marginal | M0; low resolving power (one chain, nulls from V6) | `f10_16a10_result.md` |
| Cost and K1's bound | 16a.11 | M0 measured; full scale projected | 507.5 to 678 CPU-h per M0 posterior; about 21k projected at full scale against K1's 1,000 | K1 fires at full scale: **budget-bound** | 16-recovery §12 "16a.11" |
| 16b.1, the moved truth | 16b.1 | M0, one dataset, one chain | ENO recovers: shrinkage 0.258, shift ratio 0.951 ± 0.021 | M0 | `m0_16b_result.md` (job 19378382) |
| 16b.2, the full set | 16b.2 | M0, V6 | Every truth inside its central 95%; all six shrink below 0.5 | M0, five of D11's six | `m0_16b_result.md` |
| 16b.3, coverage | 16b.3 | M0, R = 50 | 29 of 30 cells in band; ENO 46 of 50 at 95% | M0; K4's band restated at R = 50, from 100 | `m0_16b_result.md` (job 19399322) |
| 16b.4, F6 | 16b.4 | M0, one chain per stream | ptsG under metabolites only unmixed and not read; ENO under transcripts only not resolved from the prior | M0 | `m0_16b_result.md` (job 19395555) |
| 16b.5, K2 | 16b.5 | M0, 50 cells (K2 is defined at 1,000) | ENO and FBA under transcripts only below 0.05 nat | M0 read-off; K2 itself **budget-bound** | `m0_16b_result.md` (job 19395555) |
| 16b.6, K7 | 16b.6 | — | Not run | **Budget-bound** (not funded at 16a.11) | 16-recovery §11 |

## What it licenses

- **On M0, the three-block Gibbs sampler is calibrated at production settings.**
  - V7 passes SBC over 50 prior-drawn datasets.
  - Coverage is in band in 29 of 30 cells.
  - Truths are recovered in V6, and a truth moved by one prior SD is followed
    (16b.1).
- **On M0, the joint data resolve all five targets and σ_b from their priors**
  (F6's joint column).
- **On M0, each stream informs its own block.** Transcripts resolve the transcript
  parameters, and bulk metabolites resolve ENO, FBA and σ_b (F6).
- **The hybrid composition can be inferred at all.** One sampler runs a jump block
  coupled both ways to an ODE block, with every block updated, and is calibrated.
  That is evidence that the architecture works.

## What it does not license

- **Anything at full scale.** A Core A′ posterior is projected at about 21k CPU-h
  and was not run. Every number above is M0's, at 50 cells, never Core A′'s or
  syn3A's.
- **That information crosses between the blocks on M0.** Neither direction is
  measured on the exact sampler.
  - **Transcripts to metabolic parameters:** not resolved. ENO and FBA are at
    0.0041 and 0.0081 nat, below their floors, and below K2's 0.05 nat. At ENO
    this is the measured null (§12 2026-10-05).
  - **Metabolites to ptsG:** not measured. The transcript parameters did not mix
    under metabolites only (ESS 6 and 7), and the user chose to report them as
    unmixed (§12 2026-10-28). A redesign and rerun is kept as an option.
  - Until one of these is measured, the cross-block claim at M0 rests on track
    A's linearised estimate (ptsG at 0.47 of prior SD; `approximate-method.md`
    A1) and on the approximate route's fit. Neither is the exact posterior.
- **Precision beyond what single chains carry.** The 16b fits and F10 are one
  chain each. Some parameters hold 30 to 75 effective draws in a production chain
  (F10). F10 cannot see a bias below about 0.04 to 0.2 SD of W1.

## The scoping note's caveats

- **Core A′ is the best-measured region of the network.** M0 is a further cut of
  it, to two genes: ptsG, the headline, and GAPD, the best-determined promoter.
  Recovery here is not evidence for parts of syn3A with sparser data or wider
  priors.
- **Enzyme competition is partly removed.** Reactions that share a gene with
  Core A′'s enzymes but lie outside it are dropped. K7 names seven for PYK's gene,
  a count 16b.6 would have checked. On M0, 15 proteins are also fixed at their
  counts.
  - **K7's factor is not measured.** 16b.6 was not funded, so K7 is budget-bound.
    No conditioning factor is quoted beside the recovery claims, and none can be
    assumed to be 1.
- **Success is evidence that the architecture works, not that inference on the
  full model is well-posed.** A calibrated sampler on a two-gene cut says nothing
  about identifiability in the full model.
