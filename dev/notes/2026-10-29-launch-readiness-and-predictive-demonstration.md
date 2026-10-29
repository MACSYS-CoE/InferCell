# InferCell launch readiness and the predictive demonstration

**Date:** 2026-10-29  
**Status:** assessment and recommendation, written for discussion. This note is
not a decision of record and changes neither `spec/spec.md` nor
`spec/phases/approximate-method.md`. Decisions that should bind the release must
be recorded in the specifications and release issues.

**Repository basis:** `origin/main` at `b58863d` (B7.8 merged), the
`paper-1-draft` branch at `ac4f4fd`, and the in-progress `issue-89-b8` branch at
`a2416d7`. Live GitHub CI and pull-request state were checked on 2026-10-29.

**Related material:**

- [`2026-10-26-infercell-framing.md`](2026-10-26-infercell-framing.md): the
  project's goal, stages and open calls.
- [`2026-10-19-paper-1-outline.md`](2026-10-19-paper-1-outline.md): paper 1 and
  the uncertainty-budget demonstration.
- [`2026-10-19-inference-framing-discussion.md`](2026-10-19-inference-framing-discussion.md):
  decisions and open questions O1--O9.
- [`../scripts/av6_89_b7_result.md`](../scripts/av6_89_b7_result.md): B7 and
  B7.8's comparison of the mixed and exact routes.
- [`../scripts/corea_census_14c_result.md`](../scripts/corea_census_14c_result.md):
  the prior-draw clipping and starvation census.
- GitHub issues [#90](https://github.com/MACSYS-CoE/InferCell/issues/90) and
  [#91](https://github.com/MACSYS-CoE/InferCell/issues/91): track C and release
  documentation.

## 1. Executive judgement

The central recommendation is right: paper 1 should not attempt joint inference
over every parameter. It should demonstrate

> calibrated inference for the parameters the data can identify, explicit
> uncertainty distributions for the remainder, and honest forward prediction
> through both model descriptions.

This is stronger than a nominal all-parameter posterior that is expensive,
unvalidated and equal to the prior in most directions.

The inference core is close to ready. The paper and release are not. InferCell is
one major scientific deliverable -- the predictive demonstration -- plus a full
paper and release closure pass away from launch. Compute is not presently the
critical path.

## 2. How far from launch

| Area | Current state | Remaining launch gate |
|---|---|---|
| Exact and mixed inference on M0 | Substantially complete. The mixed route passes SBC for 6 of 6 parameters; median cost is 6.57 CPU-h per posterior. Small but resolved discrepancies remain. | Integrate the evidence into the paper and state the caveats plainly. |
| Probabilistic prediction | The required public feature is not built. The existing `posterior_predictive` is an early all-ODE implementation. | Hybrid prediction, chain ingestion, completion from priors, viability conditioning, interventions and reproducible seeds. |
| Paper | A substantial draft exists, but it branched before B6, is 111 `main` commits behind and contains 46 placeholder occurrences. | Evidence integration, the predictive result, five figures, Table 1, authorship and affiliations, venue material and review. |
| Release | Unit CI is green. Documentation CI is red. README and spec claims are stale, and `Project.toml` remains at 0.2.0. | Docs, claim audit, version bump, citation and DOI, release notes and archival release. |
| Scope | O2--O4 remain unresolved on `main`. The current B8 branch records the opposite of the recommendation below. | One explicit decision recorded consistently in the spec, paper, issues and README. |

### 2.1 Corrections to the earlier launch assessment

- The latest unit workflow on `b58863d` succeeded. The reported
  `g.fractional > 1.0` failure is stale.
- Documentation failed because `api/corea-modules.md` generated 201.95 KiB of
  HTML against Documenter's 200 KiB threshold. This is a real release blocker,
  but not a general failure to compile the documentation.
- B8 is not merely hypothetical. Its branch records "full B8, cap reopened with
  no ceiling," has a committed Core A-prime dataset generator, and had AV7 and
  profiling scripts in progress when inspected. Excluding B8 now requires an
  explicit superseding decision and a clean parking of that work.
- The paper is farther behind the evidence than "fill section 5.5 and draw the
  figures" implies. Its status still describes B6 and B7 as unbuilt.
- The vendored upstream parameter extracts are a release risk. The README says
  the upstream repository has no licence and therefore grants no redistribution
  rights. Before a high-profile release, obtain permission, stop redistributing
  the extracts and retrieve them during setup, or adopt another defensible
  solution.

## 3. What the predictive demonstration should claim

The scalable object is a modular, or cut, predictive distribution.

- The parameters inferred on M0 have a posterior under M0's forward map.
- Parameters completed from their original priors have not been inferred.
- Parameters retained after a viability screen are no longer at their original
  priors; they have a viability-conditioned or history-matched distribution.
- Injecting the M0 posterior into Core A-prime does not make it the joint
  posterior of Core A-prime. It is a deliberate modular transfer.

That distinction is a feature of the method, not a weakness. It is the honest
alternative to claiming a full posterior that the data and computation cannot
support. The code, Methods and figure labels should all use this language.

## 4. Compare inference and forward descriptions separately

Pushing posterior draws through both forward descriptions is necessary, but a
single comparison would mix two different approximations: the likelihood used
for inference and the description used for prediction. Use a crossed design:

| Posterior source | Hybrid forward model | Mean-field forward model |
|---|---|---|
| Exact Gibbs | Reference prediction | Forward-description effect |
| Mixed NUTS | Inference-route effect | Cheap deployed route |

This answers separately:

1. Did the approximate likelihood move the prediction?
2. Did the formalism swap move the forward prediction?

Use paired parameter draws where possible. The primary estimand should be a
population mean, an expected output, or an expected intervention effect. The
mean-field ODE cannot reproduce a single-cell predictive distribution, so it
must not be compared naively with raw single-cell hybrid variability. Report
intrinsic cell-to-cell variability separately.

## 5. Complete the posterior without pretending it is joint inference

Let

- `theta_I` be the inferred dynamical parameters;
- `theta_H` be the parameters completed from asserted priors or a
  viability-conditioned distribution;
- `xi` be intrinsic stochasticity;
- `V` be the predeclared viability event; and
- `Q` be the predicted quantity.

The demonstration samples a cut distribution of the form

```text
p_cut(theta_I, theta_H | data, V)
    = p(theta_I | data under M0) p(theta_H | V),
```

or a more explicit conditioned variant if viability depends on both groups.
This is not the full joint posterior under Core A-prime. Record the dependency
assumption and test its sensitivity.

The user-facing prediction entry point should support, at minimum:

```julia
predict(models, posterior;
        forward = :hybrid,
        completion = :viability_conditioned_prior,
        quantity,
        intervention = nothing,
        estimand = :population_mean,
        rng,
        n_parameter_draws,
        n_cells)
```

The exact spelling is not important. The important separation is between the
posterior source, forward description, completion policy, estimand and
intervention. Its result should retain parameter provenance, random seeds,
reduction labels, viability decisions and failed trajectories.

## 6. Viability is conditioning, not cleanup

The 14c census gives two distinct readings:

- with only D11's six targets varied Haldane-consistently, 8 of 200 draws (4.0%)
  are starved, although 163 of 200 (81.5%) clip at least once;
- with the broad set of informed constants varied Haldane-consistently, 95 of
  200 draws (47.5%) are starved.

Therefore "clipped," "starved" and "viable" must not be used interchangeably.
Before inspecting the prediction:

1. define the viability summaries and thresholds;
2. state whether occasional clipping is allowed;
3. record the acceptance rate and effective sample size;
4. preserve correlations induced by conditioning; and
5. repeat the conditioning separately under half, default and double prior
   widths.

If no starvation is part of `V`, starvation probability cannot also be the main
predicted quantity. Conditioning has already removed the event to be predicted.
Starvation remains useful as an unfiltered prior-predictive diagnostic and as a
model-validity panel.

## 7. Define the uncertainty split before running it

The unambiguous first decomposition is

```text
Var(Q | data, V)
  = Var_theta[E(Q | theta)]
  + E_theta[Var(Q | theta)].
```

It separates parameter uncertainty from intrinsic stochasticity. Viability
conditioning generally correlates `theta_I` and `theta_H`, so splitting the
first term between inferred and held parameters by an ordinary independent
ANOVA is not valid. Use either

- a fixed and reported nested ordering, with an ordering sensitivity; or
- Shapley effects or another variance allocation designed for correlated
  inputs.

M0's six inferred quantities are five dynamical parameters plus `sigma_b`.
`sigma_b` is observation noise, not a cell-dynamics parameter. It should not be
propagated into a latent biological quantity. If the target is a future noisy
observation, report four components:

1. inferred dynamical parameters;
2. prior-held or history-matched parameters;
3. intrinsic stochasticity; and
4. measurement noise.

## 8. Prediction quantity and intervention

### 8.1 Primary quantity

Use integrated lactate export, or the independent `exported` fraction returned
by `carbon_accounts`, as the primary quantity, subject to a short sensitivity
pilot confirming that both inferred and held parameters affect it.

Reasons:

- it is a smooth, continuous quantity;
- it is defined and metered in the reduced model;
- it depends on the coupled metabolic and gene-expression system;
- it survives viability conditioning; and
- `carbon_accounts` identifies `exported`, rather than the stoichiometric
  homolactic ratio, as the independent quantity.

Use minimum GTP, clipping and starvation as secondary robustness and validity
diagnostics on the unfiltered ensemble. Do not make fractional growth the
headline prediction: the reduced model explicitly does not reproduce syn3A's
doubling phenotype, and growth is partly a model closure rather than an
independent biological result.

### 8.2 Intervention panel

Include one mild ptsG promoter knockdown if it remains within the validated
regime. It is the best first intervention because it

- acts on an inferred parameter;
- crosses the gene-expression/metabolism boundary;
- exercises the discrete and continuous descriptions; and
- can reuse the 50 synthetic-data fits.

Fit only the unperturbed data. Predict the paired effect

```text
Delta Q = Q_intervention - Q_baseline
```

with common cell seeds where possible. Across the 50 datasets, report

- coverage of the true intervention effect;
- interval width or another sharpness score;
- exact-posterior versus mixed-posterior prediction; and
- hybrid versus mean-field forward prediction.

Coverage alone is insufficient because an arbitrarily broad interval covers.
A short, predeclared dose screen should choose a perturbation that does not push
the nucleotide pools into the regime where the mean-field reduction is known to
be at risk. A 20--30% knockdown is a starting hypothesis, not a decision to make
without the screen.

The baseline uncertainty budget remains the required result. The intervention
is one additional panel, not permission for a new intervention study to delay
the release.

## 9. Recommended open-call decisions

### O2: release and paper claim

Use a narrow, evidence-matched claim:

> InferCell supports exact stochastic and mean-field descriptions across a
> coupled ODE--jump boundary. On M0, both inference routes pass
> simulation-based calibration, with small declared discrepancies.
> Posterior-updated parameters can be combined modularly with
> viability-conditioned uncertainty in the remaining parameters to produce
> forward and intervention predictions on the reduced model.

Do not claim an all-parameter posterior, calibrated inference on Core A-prime,
or a biological result unless new evidence earns it.

### O3: framework generality

Do not build a generic registry for v0.3.0. Say "a compositional inference
framework demonstrated on reduced syn3A." Generality across organisms requires
another implementation, not more abstraction in this release.

### O4: quantity

Choose lactate export as the primary prediction, with GTP starvation as an
unfiltered validity diagnostic. Use the intervention effect on lactate export
as the strongest optional panel.

### O1/B8: release scope

If timely launch is the priority, park B8 outside v0.3.0 and paper 1. B8 is a
large full-scale calibration exercise without an exact Core A-prime reference.
The predictive demonstration completes the scientific argument more directly.

This recommendation supersedes neither the current B8 branch nor its recorded
decision by itself. If adopted, amend the spec explicitly, preserve the branch
and its evidence, and decide whether the already-built dataset generator has a
release use independent of AV7.

## 10. Methods paper versus biological paper

The reduced-model work uses synthetic data and asserted priors. Its result is a
methods demonstration, not a finding about a cell. This is acceptable if paper
1 is submitted as a methods paper.

The venue decision changes the critical path:

- for PLOS Computational Biology or a methods-oriented venue, M0 validation plus
  the reduced-model predictive and intervention demonstration can be sufficient;
- for Cell Systems, reviewers are more likely to require a biological finding,
  making prior forward UQ on the full syn3A port important.

Do not quietly turn the full-cell result into a launch blocker without choosing
the venue and paper type first. It can proceed in parallel and become paper 2 if
the first paper remains methodological.

## 11. Critical path to launch

1. **Freeze the claim and estimand.** Decide O2, O3, O4, B8, the primary
   quantity, the viability rule and the intervention dose-selection rule.
2. **Build one prediction interface.** It must drive both forward descriptions,
   read the existing chains, complete parameters, apply interventions and retain
   provenance and seeds.
3. **Run the evidence package.** The M0 crossed comparison is primary; the Core
   A-prime cut predictive distribution, viability and prior-width sensitivity
   demonstrate scale; the ptsG intervention is a panel of the same package.
4. **Integrate the paper.** Reconcile the paper branch with `main`, fill B7 and
   the predictive sections, build the figures and table, and audit every claim
   and scalar.
5. **Close release hygiene.** Fix docs CI, update README and the specs, resolve
   data redistribution, add `CITATION.cff`, choose authorship and venue, rerun
   tests and quick starts, and prepare release notes.
6. **Archive together.** Bump to 0.3.0, tag the code used by the paper, create a
   Zenodo DOI, and point the preprint and code release at one another.

## 12. Launch gates

Call v0.3.0 and paper 1 ready only when all of the following hold:

- the baseline predictive quantity and estimand are fixed before results are
  read;
- the exact/mixed and hybrid/mean-field predictive comparisons are reported;
- the viability rule, acceptance rates and prior-width sensitivities are
  reported;
- the variance allocation handles induced parameter dependence and separates
  measurement noise;
- intervention coverage is accompanied by sharpness, if the intervention panel
  is included;
- the paper contains no result, figure, table, authorship or availability
  placeholders;
- unit and documentation CI are green;
- README, docs, specs, paper and release notes make the same claim;
- redistribution of the upstream parameter extracts has a defensible basis; and
- the paper cites an immutable v0.3.0 archive with a DOI.

