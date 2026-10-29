# InferCell v0.3.0 and paper 1: completion plan (draft)

**Date:** 2026-10-29
**Status:** a working draft. It decides nothing. Anything here that should bind
goes into `spec/spec.md` §12 (or `spec/phases/approximate-method.md` §12) with a
date and a reason.

**Based on:**
- `dev/notes/2026-10-26-infercell-framing.md`;
- `dev/notes/2026-10-29-launch-readiness-and-predictive-demonstration.md` (the
  third-party assessment, called "the note" below);
- the discussion of 2026-10-29, which accepted most of the note and adjusted five
  points (§3).

**Repository state:**
- `origin/main` at `b58863d` (B7.8 merged);
- `paper-1-draft` at `ac4f4fd`, 111 commits behind main, with 46 placeholders;
- `issue-89-b8` at `a2416d7`, in progress in another session.

## 1. What the launch claims

Paper 1 is a methods paper. Its claim (O2, the note's wording):

> InferCell supports exact stochastic and mean-field descriptions across a coupled
> ODE–jump boundary. On M0, both inference routes pass simulation-based
> calibration, with small declared discrepancies. Posterior-updated parameters can
> be combined modularly with viability-conditioned uncertainty in the remaining
> parameters to produce forward and intervention predictions on the reduced model.

It does **not** claim:
- an all-parameter posterior;
- calibrated inference on the reduced model (Core A′);
- a biological finding.

Synthetic data and asserted priors throughout. In user-facing text say "reduced
syn3A model", not "Core A′".

## 2. Decisions to record before any work (step 1)

| Call | Proposed decision |
|---|---|
| O2: the claim | As in §1 |
| O3: generality | No generic registry. "A compositional inference framework demonstrated on reduced syn3A" |
| O4: the quantity | Lactate export, as `carbon_accounts`' `exported` (or integrated lactate out), after a sensitivity pilot. Minimum GTP, clipping and starvation are validity diagnostics on the unfiltered ensemble. Not fractional growth |
| B8 / O1 | Park B8 outside v0.3.0 and paper 1. Keep the branch and its evidence. Keep `corea_dataset` and bring it to main (its later use: a scored test on the reduced model). The parking is recorded from the B8 session, on its branch |
| δ for the predicted quantity | See §3.1. Proposed: report both, uncorrected and corrected |
| The viability rule | Fixed in writing before any prediction is read (§4.3) |
| The intervention | One ptsG promoter knockdown. Dose picked by a declared screen (§4.4) |
| The venue | Choose now: PLOS Comp Biol (a methods paper is enough) or Cell Systems (likely needs the full-cell stage 1 result) |

## 3. Adjustments to the note

### 3.1 The mean-field forward model has no δ for a new quantity

The mixed route is calibrated because the bulk factor is centred on m(θ) + δ, with
δ from an 800-cell hybrid ensemble at θ_ref for the **17 panel pools**. Lactate
export is not one of them. The cell "mixed posterior, mean-field forward" can
therefore run two ways:
- **uncorrected**, which shows the Jensen gap in prediction space;
- **corrected** by a hybrid ensemble at θ_ref, which costs the hybrid simulations
  the cheap route exists to avoid.

Report both. The uncorrected gap is the finding; the corrected one shows the fix
and what it costs.

### 3.2 Scored on M0, unscored on the reduced model

- **M0 is scored.** V7's 50 datasets come from M0 with every non-target parameter
  at its published value, and the fit holds them there. So on M0 there are no
  prior-held parameters, and coverage of the truth can be measured.
- **The reduced model is unscored.** The M0 posterior is carried into Core A′ and
  the remaining parameters come from their priors, filtered by viability. No truth
  matches those data. It shows that the method scales and what the uncertainty
  budget looks like there; it is not a calibration result.

The paper keeps the two apart in its text, figures and labels.

### 3.3 Viability on M0 conditions θ_I

On the reduced model, the cut is p(θ_I | data under M0) · p(θ_H | V). On M0 a
viability filter acts on θ_I alone and is extra data, so it changes the posterior
being scored. Either apply no filter on M0 (4% of D11's prior draws starve, per
14c.2) or condition jointly and report it.

### 3.4 Variance allocation by nested ordering

Two groups, inferred (θ_I) and held (θ_H). A fixed nested split, reported in both
orders, rather than Shapley effects. Plus the law of total variance for parameter
uncertainty against intrinsic noise, and measurement noise (σ_b) as a fourth term
only if the target is a future noisy observation.

### 3.5 Two prerequisites

- **The archived posteriors must hold draws.**
  - Exact: `/fred/oz022/tkimpson/infercell-16a9b-archive-2026-10-17/v7/sbc_v7_16a9b/replicateN.jls`
    (50 files of about 350 KB);
  - mixed: `/fred/oz022/tkimpson/infercell-89-b7-archive-2026-10-27/run/chain_N.jls`
    (116 files).

  Neither has been loaded to confirm this. If they hold only summaries, the
  crossed design needs refits.
- **The evidence package needs a costed pilot.** 50 datasets × 2 posteriors ×
  2 forward models × baseline and intervention × parameter draws × cells. The
  earlier ~50 CPU-h figure was for the budget alone.

## 4. The predictive demonstration

### 4.1 The interface

One public entry point, spelling not fixed:

```julia
predict(models, posterior;
        forward = :hybrid | :mean_field,
        completion = :none | :prior | :viability_conditioned_prior,
        quantity, estimand = :population_mean,
        intervention = nothing, delta = :none | :ensemble,
        rng, n_parameter_draws, n_cells)
```

- Reads B6's `Chains` and the V7 replicate files.
- Draws held parameters with `prior_draws.jl`'s Haldane-consistent sampler and
  skips the dead `enz_R_*` slots.
- Keeps provenance: parameter source per draw, seeds, reduction labels,
  viability decisions, failed trajectories.
- Forks the rng and restores the global state, as `infer` does (the route
  setups reseed the global generator per cell).
- Replaces or wraps the early all-ODE `posterior_predictive`.

### 4.2 The crossed comparison on M0 (primary, scored)

| Posterior | Hybrid forward | Mean-field forward |
|---|---|---|
| Exact Gibbs (V7) | Reference prediction | Effect of the forward description |
| Mixed NUTS (B7) | Effect of the inference route | The cheap deployed route (δ as §3.1) |

- Paired parameter draws and common cell seeds where possible.
- The estimand is a population mean (or an expected effect). The mean field has no
  single-cell distribution, so intrinsic cell-to-cell spread is reported
  separately from the hybrid.
- Report coverage of the truth and sharpness, over the 50 datasets.

### 4.3 The reduced-model budget (unscored, shows scale)

- The M0 posterior for θ_I, the priors for θ_H, and the viability filter V.
- V is declared in advance: its summaries, its thresholds, whether occasional
  clipping is allowed. "Clipped", "starved" and "viable" are not interchangeable.
- Report the acceptance rate, the effective sample size and the correlations the
  filter induces.
- Repeat under prior widths halved, default and doubled.
- The variance split as §3.4. The sentence to aim for: *x% of this prediction's
  uncertainty comes from parameters no measurement supports.*

### 4.4 The intervention (one panel, scored on M0)

- A ptsG promoter knockdown: an inferred parameter, crossing the
  gene-expression–metabolism boundary.
- Fit the unperturbed data only. Predict the paired effect
  ΔQ = Q_intervention − Q_baseline with common cell seeds.
- Across the 50 datasets: coverage of the true effect, sharpness, exact vs mixed
  posterior, hybrid vs mean-field forward.
- **A dose screen first.** It must keep the nucleotide pools away from K_d, since
  the one-way reduction rests on elasticities of 0.04–0.09. 20–30% is a starting
  hypothesis, not a decision.
- If it threatens the schedule, it is dropped. The baseline (§4.2, §4.3) is the
  required result.

## 5. The plan, in order

1. **Decide and record** §2's calls as a §12 entry on main. Park B8 from its own
   session; bring `corea_dataset` to main.
2. **Quick fixes**, each its own PR:
   - split `docs/src/api/corea-modules.md` under Documenter's 200 KiB limit (the
     DistributionsAD `@check_args` message is a non-fatal extension warning, not
     the failure);
   - README: calibration on M0 has passed; fix the 4D-model attribution
     (Thornburg et al.);
   - spec §1 to the new claim; STATUS flag 4's inconsistencies (21 vs 22
     reactions, the ENO-null reason);
   - write to the Luthey-Schulten lab about redistributing the parameter
     extracts. Fallback: fetch them during setup instead of vendoring them.
3. **Load the archived chains** (§3.5).
4. **Build `predict`** (§4.1) with unit tests, a parse-and-lower on the login node,
   and a smoke run of every stage from a worktree at a pushed commit.
5. **Pilots:** the cost of the evidence package, lactate's sensitivity to θ_I and
   θ_H, and the ptsG dose screen. Gate: project the full cost and confirm the
   scale before the run.
6. **The evidence run** via Slurm, archived with `SHA256SUMS`: §4.2, then §4.4,
   then §4.3.
7. **The paper:**
   - rebase `paper-1-draft` onto main;
   - fill §5.3 from B7 and B7.8, and §5.5 from step 6. Withdraw or reword the
     placeholders B8 would have filled (full-scale calibration, the abstract's
     injection count, Figure 5's third panel, Table 2's reduced-model row);
   - draw Figs 1, 3, 5, 6, 7 and fill Table 1 with sourced cells;
   - authorship, affiliations, venue material (Author Summary or Highlights);
   - `/review-scientific-manuscript`, then `/check-arxiv-llm-compliance` and
     `/check-refs`.
8. **The release:**
   - bump `Project.toml` to 0.3.0;
   - `CITATION.cff`, release notes, the quick starts rerun;
   - unit and docs CI green, and the suite by Slurm;
   - tag the commit the paper uses, archive it on Zenodo for a DOI, and point the
     preprint and the release at each other.

**In parallel, not blocking:** full-cell stage 1 (prior forward UQ through the
syn3a-ensemble port). It becomes paper 2 unless the venue makes it paper 1's
biological result.

## 6. Launch gates

v0.3.0 and paper 1 are ready only when:
- the quantity, the estimand, the viability rule and the intervention dose rule
  were fixed before results were read;
- the crossed comparison is reported with coverage and sharpness, and δ's
  handling is stated;
- the reduced-model budget reports acceptance rates and the prior-width
  sensitivities, and is labelled unscored;
- the variance split states its ordering and separates measurement noise;
- the paper has no result, figure, table, authorship or availability placeholders;
- unit and docs CI are green;
- README, docs, specs, paper and release notes make the same claim;
- redistribution of the upstream extracts has a defensible basis;
- the paper cites an immutable v0.3.0 archive with a DOI.

## 7. Open questions

- The venue (§2): it decides whether full-cell stage 1 blocks the launch.
- δ for the predicted quantity (§3.1): both, or one only.
- Whether the archived posteriors hold draws (§3.5).
- Whether 16b's option (b), a joint S_g-and-path move (~550 CPU-h), is needed. Only
  if the paper claims, on the exact sampler, that ptsG shrinks under metabolites
  only.
