# Response to the review of paper 1

**Date:** 2026-11-25
**Review:** [`manuscript_review.md`](manuscript_review.md) (2026-11-25), of the draft at `32b8496`.
**Item-by-item detail:** [`checklist.md`](checklist.md), 35 items with evidence, verdicts and tiers.

**Status of this note.** It records a discussion, not decisions. Tom is leaning towards one paper on
the full syn3A model. Nothing here is in the Decision log, Linear or the manuscript. Anything that
should bind becomes a Decision log entry.

Numbers marked *login-node* are quick reads of archived draws. Numbers marked *scaled* are
extrapolated from the reduced model's measured costs. Neither kind should be quoted in the paper.

---

## 1. Summary

- The review's facts mostly hold. Its counts and quoted numbers check out against the draft.
- Its strongest points are the uninformative knockdown, the single dataset, and the unexplained
  promoter shift.
- Its largest proposal is to drop the reduced model for the full cell. My first response was to
  disagree. Tom's reply moves towards it, so §4–§6 scope what it would cost.
- Where we now stand:
  - **Agreed:** no real data. *Nat. Comput. Sci.* is the target and PLOS Comput. Biol. the
    fallback, and neither needs it.
  - **Tom's lean:** one strong paper, not two, on the full syn3A model.
  - **My recommendation:** do the full cell, but keep the reduced model as a middle tier, and
    take three measurements before committing (§6).

## 2. What the review got right, and wrong

### Agreed with

| Point | Why | Checklist |
|---|---|---|
| The 30% ptsG knockdown tests nothing new | It moves Q by 0.4%, with an MC share of 1.09. The draft concedes it | 2b |
| One dataset and one truth | The draft's largest weakness | 2a |
| The promoter shift reads as unfinished | 3.3 SE, with "we have not run" a second dataset | 1b |
| Positioning against SBI, LNA and emulators is missing | 27 of 56 references are uncited | 4b |
| The first-order restriction is not listed as a limit | A referee from *E. coli* or yeast asks at once | 3b |
| Why NUTS, if Laplace matches it to 4–8% | The draft does not say | 4a |
| Why ENO and FBA | The draft never gives the reason | 3a |
| The §2 code listing shows an appendix use | Already open in `STATUS.md` | 6d |
| Figs. 1 and 2 are below the data figures | The caption of Fig. 2 warns the reader about its own arrows | F1, F2 |
| Reframe around targeted, per-module inference | It matches the framing note of 2026-10-26 and Spec D7 | A1 |

### Disagreed with, at the time of the checklist

| Point | My view | Checklist |
|---|---|---|
| Drop the reduced model for M0 plus the full cell | Disagreed. Now reopened by Tom: see §4 | B1 |
| Add a real-data case study | Disagree. A written non-goal. **Tom agrees** | 5c |
| Fit the reduced model to full-simulator data | Disagree for a reduced-model paper. It becomes natural if the full cell is in the paper | 5b |
| Promote the three graphs as a flagship result | Partly. The key direction has no exact result, and the saturation premise fails where 99.3% of cells clip GTP. Fix the contradiction with the Discussion; promote only if KIM-116 supports it | 6a |
| "M10 is mandatory" | Partly. It is the parked release-vs-M10 question. The snapshot scRNA-seq rung is a good addition | 5a |
| End the abstract on a forward-looking advance | Partly. The review splices two sentences, and its replacement claims things not shown | 6c |
| Delete "we did not read the supplements" | Only after someone reads them | 6b |
| Say neural SBI is prohibitive; integrate 27 references | Partly. Cite what the positioning needs | 4b |
| NUTS stress test with 10–30 metabolic constants free | Partly. A Hessian eigen-spectrum answers it for a few CPU-h | 3a |

### Stale or wrong in the review

- **KIM-90 is merged** (PR #189, 2026-11-23), with no cause named. The review treats it as pending.
  Its three todos can be filled once main is merged into the paper branch.
- **The disjoint-cohort test** needs a new exact run of 508–678 CPU-h. Tom declined it on
  2026-11-23.
- **`fig_model` already prints its clocks,** and `tab:reducedcost` already itemises the probe.
- **"GTP clipping is a pathology of the reduction"** is unverified. The clip rule is the published
  model's.

### A correction to my own first response

I said full-cell inference "is already PreCell's paper". That overstated it. PreCell infers 28
group multipliers from bulk measurements, with no single-cell data. The overlap with a full-cell
InferCell paper is the organism, not the method or the parameter set. The checklist is corrected.

## 3. Two findings from checking the review

Both are *login-node* reads of `infercell-kim92-archive-2026-11-17/pilot_draws.tsv`. Both need a
tracked script before they are quoted.

1. **The promoter shift is harder to call noise than the draft says.**
   - The 3.3 SE assumes independent errors. With the posterior's own correlations it is 3.1 SE.
   - 4 of 1,600 draws put the promoters' mean standardised error at or above zero.
   - The promoters' mean pairwise correlation is 0.01. Noise shared across promoters would show as
     correlation between them, and it does not.
   - So either this dataset is roughly a 1-in-400 draw, or the posterior is miscalibrated. "One
     dataset cannot separate the two possible causes" is too even-handed.
2. **The "MAP" in §5 is the MAP of a different target.**
   - The pilot computes θ_ref as "the bounded MAP at δ = S = 0" (`corea_pilot_kim92.jl:18`). The
     fitted likelihood has δ and S in it.
   - ENO's and FBA's marginals have skewness 0.05 and 0.06, so "skewed" does not explain their 2.3
     and 1.4 SD offsets. The likeliest reading is that δ and S move them. This is an inference.
   - "The MAP shares the shift, so it comes from the likelihood, not the sampler" needs rewording.
     What the numbers show is that the promoters sit 0.86 SD low before δ and S enter and 0.81
     after.

## 4. The direction Tom is leaning towards

Tom, after reading the checklist:

- **One paper, not two.** I read "two" as paper 1 on the reduced model and paper 2 on the full
  cell.
- **The full syn3A model, not the reduced one,** as the stronger paper.
- **SpeedyCell cannot be ported directly,** but the equations should carry over.
- **A forward validation comes with it:** run InferCell and SpeedyCell forward at fixed
  parameters and check they agree.
- **No real data.**

What this reopens in the checklist: B1 (was "disagree"), V1 (venue now named), 5b (the
misspecification test becomes natural), and RC1–RC3 (now the plan, not a side question).

## 5. Scope of a full-cell paper

### Size

| | Reduced model | Full syn3A |
|---|---|---|
| ODE species / reactions | 32 / 21 | 146 / 175 |
| Stochastic reactions | 52, all first order | 2,365, of which 363 second order |
| Genes | 17 | 455 |
| Mean-field state (what the gradient solves) | about 85, dense | about 2,800 |

Full-cell figures are from `SpeedyCell.jl/data/model/inventory.md`.

### Compute, line by line

All *scaled*. The anchors are measured: a gradient of 46 s at 21 parameters (KIM-99), of which the
solve is 12 s and the transcript factor 6.6 s per dual pass; chains of about 110 h each; 452 CPU-h
for the chains; 5,447 CPU-h for the scored predictions at 32 s a cell (KIM-95).

| Line | Estimate | What it rests on |
|---|---|---|
| Forward agreement with SpeedyCell, a few hundred cells each | under 100 CPU-h | Cheap; the cost is engineering |
| One gradient, about 20 free parameters, all transcripts observed | about 20 min | Solve 12 s → 400 s per pass if a sparse solver scales linearly; transcript factor 6.6 s → 180 s; two passes |
| One joint posterior on that block | about 11,000 CPU-h; about 115 days per chain, serial | ×25 on today's chains. Milan's limit is 7 days |
| Same, with the gradient threaded over genes and passes | 2–3 weeks per chain | Assumes 5–10× from threading; needs checkpointed chains |
| All 455 promoters free, joint | not feasible | Forward-mode cost is linear in parameters: about 8 h per gradient |
| Gene-expression block by the data plug-in | hundreds to a few thousand CPU-h, parallel by gene | 11.3 s per gradient now, with no ODE solve; genes decouple given observed pools |
| Scored predictions at KIM-95's sizes | 14,000–70,000 CPU-h | Low end assumes SpeedyCell's speed, about 80 s a cycle; high end assumes InferCell's present driver scaled up |
| Exact sampler on the full cell | out | 21,000 CPU-h was the projection for the reduced model |

What follows from the table:

- **Only targeted inference is feasible.** A joint posterior over everything is out. A block of
  about 20 parameters, or the plug-in over the gene-expression block, is in. This is the review's
  Point A, arrived at by cost.
- **A sparse stiff solve under dual numbers is a prerequisite.** A dense solve at 2,800 states is
  about four orders of magnitude worse than today's.
- **Wall-clock binds, not CPU-hours.** The CPU-hours are comparable to what the project has spent
  already (about 40,000 on the exact runs alone).
- **Scored prediction is the largest line,** and its bracket is the widest.

### What is not straightforward

The forward port is the easy part. SpeedyCell's model is neutral data files, rate strings and an
expression graph, that an InferCell module could read. Three things are harder.

1. **The mean-field closure is no longer exact.**
   - The 363 second-order reactions do not close.
   - Transcription depends on gene copy number, through 493 replication reactions.
   - The rate rules use `round`, `min`, `max` and thresholds.
   - The paper's claim that the means "close exactly, as linear ODEs" goes, and M0's exact check
     covers none of this. A new small exact check is needed, or the reduced model stays as the
     middle tier.
2. **Priors need their own import.** Kinetic constants are literals inside SpeedyCell's rate
   strings. Priors and provenance must come from the SBtab files, as they do now for glycolysis.
3. **The forward comparison is distributional, and needs a published-faithful mode.** InferCell
   corrects the base-count swap in transcription and takes the balanced column on 8 of 32
   glycolytic constants. SpeedyCell reproduces the reference as it is.

### Effect on the release

My guess is one to two months of engineering before the first full-cell posterior. That lengthens
the v0.3.0 hold by the same amount, unless the release is decoupled from the paper.

## 6. Recommendation

1. **Keep the reduced model as the middle tier.** It is finished, and it is the only place where
   a joint posterior can referee the modular route. The single paper is then three tiers:
   - M0: the exact check;
   - the reduced model: joint against modular;
   - the full cell: targeted inference and prediction, with the forward check against SpeedyCell.
2. **This depends on KIM-116.** If the modular arm lands close to the joint, the full-cell tier
   has a method with a measured error. If it does not, the full-cell tier has no validated route
   at an affordable cost.
3. **Take three measurements before committing.** None needs the port, and together they replace
   the two widest brackets in §5 with numbers. My guess is a few hours of work.
   - Time one full SpeedyCell cycle.
   - Build the first-order mean-field system from SpeedyCell's reaction table and time a sparse
     dual solve.
   - Time the transcript factor on 455 genes of synthetic counts.
4. **Do the text fixes now, whichever way this goes.** They carry over to any version of the
   paper: fill the KIM-90 todos, reword the promoter and MAP paragraphs, add the positioning
   paragraphs and missing citations, fix the code listing and the schematics.

## 7. What is still open

| Question | Whose | Depends on |
|---|---|---|
| Commit to the full cell, or stay with the reduced model | Tom | The three measurements; KIM-116 |
| Does the reduced model stay as a middle tier | Tom | — |
| Which block of about 20 parameters is inferred on the full cell | Tom | The uncertainty budget; what a prediction of interest depends on |
| Does v0.3.0 wait for the full-cell paper | Tom | The engineering estimate |
| M10, the degradation ladder: on which model, or dropped | Tom | The direction above |
| Replicate truths by MAP and Laplace (130–200 CPU-h) | Tom | Lifts the "one truth" non-goal |
| A second, resolvable intervention | Tom | On the full cell this may come free: a knockdown with a large effect |
| The contribution sentence | Tom | All of the above |

## 8. Not done

- The three measurements in §6 have not been started.
- No tier-A text fix has been applied to the manuscript.
- Nothing is recorded in the Decision log or Linear.
- `checklist.md`, this file and the pointer in `docs/handoff.md` are uncommitted.
