# InferCell framing: the goal, and where v0.3.0 sits in it

**Date:** 2026-10-26
**Status:** a statement of the project's goal, written to be checked against at
the next reassessment. It records a discussion. It changes nothing in
`spec/spec.md` or `spec/phases/approximate-method.md`. Anything here that should
bind becomes a §12 entry there.

**Related framing notes:**
- [`the-case-for-inference.md`](the-case-for-inference.md): why whole-cell models
  need inference at all.
- [`modular-bayesian-inference-heterogeneous-modules.md`](modular-bayesian-inference-heterogeneous-modules.md):
  the survey of modular methods, and its §7 (no shared parameters).
- [`model-uncertainty-and-selection.md`](model-uncertainty-and-selection.md):
  uncertainty over model structure, not just parameters.
- [`reduced-syn3a-scoping.md`](reduced-syn3a-scoping.md): the three-step ladder
  (reduced model, well-stirred syn3A, spatial 4D) and what step 1 keeps.
- [`2026-10-13-inference-paths.md`](2026-10-13-inference-paths.md): exact,
  modular and amortised paths after phase 16.
- [`2026-10-19-inference-framing-discussion.md`](2026-10-19-inference-framing-discussion.md):
  decisions A1–A7 and the open calls O1–O9.
- [`2026-10-19-paper-1-outline.md`](2026-10-19-paper-1-outline.md): paper 1's
  structure, the uncertainty budget (its §6) and the roadmap after v0.3.0.
- [`spec/spec.md`](../../spec/spec.md) §0–§1: the claim of record, not yet
  amended to match (O2).
- [`spec/phases/approximate-method.md`](../../spec/phases/approximate-method.md)
  §1: track B's claim.

## 1. The goal

The goal as first put: posteriors for every parameter of the full syn3A model,
then run the forward simulation probabilistically, or at least propagate the
priors forward.

**Restated:**

> Predictions from the full syn3A model with honest uncertainty, knowing where
> that uncertainty comes from and which measurements would shrink it.

The forward half stays. "Posteriors for every parameter" is dropped as the
target, for three reasons.

1. **The data cannot constrain most parameters.** The full model has thousands of
   them; the data are a proteome, some metabolomics, a doubling time and a few
   kinetic measurements. Models at this scale are sloppy (Gutenkunst et al. 2007).
   A joint posterior would be the prior in most directions and a poorly mixing
   ridge in a few. The reduced model shows it already: ENO under transcripts only
   stays at 0.999 of its prior SD.
2. **The cost does not close.** One exact posterior on the reduced model, about 8%
   of the cell, is projected at 60k–150k CPU-h. The cost is the number of
   simulations, not the speed of one, so a faster simulator does not rescue it.
   This is decision A1 (2026-10-19).
3. **It is a means, not the end.** A parameter's posterior matters only as far as
   it moves a prediction.

## 2. The stages

Each stage is useful on its own, and they come in order of difficulty.

| Stage | What it is | What it needs |
|---|---|---|
| **1. Prior forward UQ** | Every parameter drawn from its prior, run through the syn3a-ensemble port. Predictive distributions, e.g. whether the ~16% PEP collapse survives the model's own parameter uncertainty | Defensible prior widths (metabolic balancing first) and a parameter-override hook in the port. About 20 node-hours for 10³ draws × 100 cells |
| **2. Condition on known phenotypes** | Keep the draws that give a viable cell: doubling near 105 min, a plausible proteome. History matching or ABC on summaries | Stage 1's runs. Cheap, and probably the largest single gain: wide priors will make many cells that do not grow |
| **3. Goal-directed inference** | For a chosen prediction, infer only the parameters its variance budget says matter, and propagate the rest from priors | The uncertainty budget (paper-1 outline §6) and the modular and approximate methods. This is where the reduced-model work pays |
| **4. Real data** | Observation models for real assays; misspecification checks | Data and collaborators |

**Caveats.**
- Stage 1 is only as good as the priors. The reduced model alone has 68 asserted
  ones. Report every result under halved and doubled widths, and as "under these
  priors".
- Stages 1 and 2 use no new inference method. They run on the port, not on
  InferCell's re-implementation.

## 3. Where the current work sits

- **The reduced syn3A model is validation infrastructure, not the main road.** It
  is where inference across the ODE/stochastic boundary is shown to be correct
  (V6, V7, F10, track B) before it is trusted on the full cell. It feeds stage 3.
- **The full cell lives in syn3a-ensemble.** Stages 1 and 2 go through the port
  (decision of 2026-10-19: the full cell as a black box first).
- **v0.3.0 is the framework plus one formalism swap, validated on M0.** Gene
  expression has two descriptions, exact jumps and a mean-field ODE, each with its
  own inference (exact Gibbs; the mixed likelihood with NUTS). It is not a general
  demonstration of modular inference, and it is not the playground.

## 4. What drifted, as of 2026-10-26

- **The claim of record is stale.** `spec/spec.md` §1 claims calibration across
  the boundary on the full reduced model. Track B and the paper claim an
  approximate method that matches the exact one on M0. O2 is open; the parent §12
  has no entry after 2026-10-13.
- **The playground leaked into the release.** The playground's planned first study
  (exact vs modular vs amortised on M0) is in effect B7.
- **The claim softened.** On 2026-10-23 the Jensen gap became reported, not a
  stop, and the reference became "the most detailed composition". AV6 (B7) is now
  the method's only test, and what a failed AV6 means for the release is not
  written down.
- **"End of track B" is not one finish line.** B7 is about 2,000 CPU-h. B8 waits on
  O1. B9 is optional. 16b is closing separately.
- **No biological finding** (paper `STATUS.md` flag 1). The M0 uncertainty budget
  (about 50 CPU-h, needs O4) shows why inference matters, but it describes a
  reduced model on synthetic data under asserted priors. It is not a finding
  about a cell. Stage 1 on the full cell would be one.

## 5. Recommendations (not decided)

1. **Stop v0.3.0 at B7.** Its claim: the framework hosts a formalism swap, and the
   cheap route reproduces the exact posterior on M0 at a fraction of the cost.
   Decide before B7 runs what a failed AV6 means.
2. **Move B8 and the amortised and surrogate comparisons to the playground spec.**
3. **Scope the M0 uncertainty budget as a task,** with O4 decided first.
4. **Start stages 1 and 2 on the full cell soon, in parallel.** They are cheap, need
   no new methods, give a biological result, and show which parameters are worth
   expensive inference.
5. **Amend spec §1 (O2)** to whatever v0.3.0 is decided to claim.

## 6. The calls still open

| Call | What follows from it |
|---|---|
| O2: what the release claims | Rewrites spec §1; fixes most of the drift |
| B8 in or out of v0.3.0 (O1) | In: O1 now and a later release. Out: the release ends at B7 |
| O3: generic registry, or "demonstrated on one organism" | Whether "framework" in the title holds |
| O4 and the budget: does paper 1 need it | Methods paper or more; the venue |
| Stage 1 now or after the paper | Whether paper 1 or paper 2 carries the first biological result |
