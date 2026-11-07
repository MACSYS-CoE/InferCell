# Probabilistic trajectories: design

Status: approved by the author, 2026-11-07. Nothing here has run except the pilot in §6.

## Intent

Paper 1 presents InferCell as a framework. Its predictions are currently shown as end-of-run numbers
and scatter plots, and a reader never sees a prediction unfold in time with its uncertainty. We want a
time series with a proper shaded region.

The author added a second point: **this is a framework, so anyone should be able to release a
parameter to its prior and see what that does to the uncertainty of the prediction.** The figure has to
show that this works, and the capability has to be in the package and documented. It cannot live only
in a paper script.

## Goals

1. `predict` can return time series, in the public API, with its documented behaviour otherwise
   unchanged.
2. One new two-column figure, in the paper's story order: scored on M0, then unscored on the reduced
   model, then what happens when more parameters are released.
3. A documented recipe for releasing a parameter to its prior, tested in the suite.
4. The endpoints of the new runs reproduce the published numbers exactly, where the design replays a
   published run.

## Non-goals

- No new inference. No posterior changes, and no parameters fitted beyond M0's five plus σ_b.
- No new scoring. Coverage stays an endpoint result over the 50 M0 datasets. Time-resolved coverage is
  out of scope.
- No knockdown time series. KIM-42 could not resolve ΔQ's parameter variance at 4 cells per pair, and
  these runs have no more cells.
- No convenience wrapper over the module constructors (the builder's call in §5 explains why).

## Anti-goals: this fails even if it runs, if

- a reader could take panel A as validating panels B or C. The transfer from M0 to the reduced model is
  untested, and the caption must say so;
- the band reads as a prediction about a real cell. All of it rests on synthetic data;
- panel C reads as a posterior over 82 parameters. It is a sensitivity to priors, and 62 of those
  parameters were never fitted;
- the replicate or the released set is picked after looking at which looks best. Both are fixed in this
  document, before the runs.

## The uncertainties, and where each one goes

| Source | Count | Treatment |
|---|---|---|
| Inferred, θ_I | 5 | M0 posterior draws |
| Held promoters, θ_H | 15 | Prior, geometric SD 2. The width is our choice |
| Released kinetic constants (panel C only) | 62 | Prior at width 1, reverse constants derived by Haldane |
| Cell-to-cell noise | n/a | Its own band layer |
| Assay noise σ_b | 1, M0 only | Only on observed species in panel A, as a lighter outer band |
| Other constants | the rest of 260 | Not drawn: driver-written, derived or fixed. Named in the caption's source text |
| Model form | n/a | Not shaded. It stays in the model-form table |

The 260 and the group counts come from the priors tabulated on 2026-11-07 (§6).

## The figure

A new `fig_trajectories`, full width. It adds to `fig_predict_m0` and does not replace it.

**Panel A: M0, 0–600 s, against the truth.** Dataset: V7 replicate 1, exact posterior. It is chosen
because this posterior feeds the reduced model, and it is not chosen for how it looks. Three species, each with the posterior band and the true population mean:
- lactate exported, `Q` over time. **The fit did see this, noisily:** `M_lac__L_e` is in `M0_PANEL` and
  its only dynamics are the export flux over the medium's volume ratio. So the panel overlays those
  data, rescaled, and the caption says the quantity is observed;
- GTP, `M_gtp_c`, also in the panel, with its bulk data and a lighter outer band for a new bulk
  measurement (σ_b added);
- pyruvate, `M_pyr_c`, which is not in the panel. This is the one curve that is a prediction of
  something the fit did not measure, though it is tied to the measured pools by the model.

Q's predictive SD on this dataset is 1.8% of its mean (`score_datasets.tsv`), so a band on a raw
cumulative curve may be barely visible. The builder may plot lactate as a deviation from the posterior
mean, or add a residual axis.

Caption: one dataset, an illustration; calibration remains the endpoint coverage over 50 datasets. The
data are regenerated with `m0_dataset(16_200_001; ncells = 50)`, since the fit's file does not hold
them, and the caption or note says so. The regenerated truth must equal the fit's recorded truth, and
one small job at V7's commit `dfe570a` confirms the regenerated bulk data match.

**Panel B: reduced model, 0–6,300 s, unscored, width 1.** Three nested layers, for lactate export and
for GTP:
1. θ_I: quantiles over rows of the row means (θ_H averaged out, 80 cells per row);
2. θ_I + θ_H: quantiles over pairs of the pair means (4 cells per pair);
3. a new cell: quantiles over all 8,000 cells.

Bands are 90% (with a 50% inner band if it stays legible; the builder decides). GTP carries a reference
line at K_d, which makes the existing caveat visible. A small sub-panel plots the held share of Var(μ)
against time, from `variance_split` applied at each time, ending at the published 0.908.

Layer 2 includes the Monte Carlo term of 4 cells, about σ̂²/4. At 6,300 s that is about 7% of the
variance. The text says so, and the note gives the figure at the end of the run.

**Panel C: release to the prior, reduced model, unscored.** Lactate export and GTP, with the 20-parameter
band of layer 2 against the 82-parameter band. The filled bands use the viability rule's `per_theta_I`
weights, as the budget's headline does, with weighted quantiles. The same 82-parameter band without the
rule is drawn as an outline, because the rule may remove the starved tail that is the point of the
panel. `joint` goes in the result note. Released: every
`kcatF_*` and `km_*` of `CentralGlycolysis` and `NucleotideRecycling` that is currently fixed, 62
parameters. Reverse constants are derived from them, not drawn. Reported with it: the viability
acceptance, the starved fraction, and the median of the pair means at 6,300 s against panel B's, with
and without the rule. It is labelled a sensitivity and not a
posterior. The cut assumption is stressed harder than in panel B, because M0's posterior was fitted with
these constants at their nominal values.

The code to do it is short enough to show in the paper as a listing, which carries the framework point:

```julia
ms = d11_models()
ms[1] = CentralGlycolysis(enzymes = :translated,
                          free = [:kcatF_R_ENO, :kcatR_R_ENO, :kcatF_R_FBA, :kcatR_R_FBA,
                                  :kcatF_R_LDH_L, :kcatR_R_LDH_L])
p = predict(ms, posterior; horizon = COREA_CYCLE_S, completion = :prior, n_held = 20, n_cells = 4)
```

Without further arguments, `predict` draws every newly free parameter from its prior. Here that adds
`kcatF_R_LDH_L` to the held set, and `kcatR_R_LDH_L` is derived. The pilot ran this form.

## How it gets made

**Software.** A keyword on `predict`, `forward_cells` and `intervention_cells`: `trajectory = nothing`
by default, or a list of ODE-block species to keep, with a sampling grid (`every = 60.0` s). Each
`CellRun` keeps them. `_cell` already holds the `ValidationRun`, whose `ode` is recorded at every
handshake and is discarded today, so this adds no simulation. The default has to leave results bit for
bit identical: `test/fixtures/predict_reference.jl` is the guard. `CellRun` gains a field, so existing
`CellRun`s in archived `.jls` files will not open at the new commit. This run reads `w*.jls`, which
hold none, and the `.jls` files open only at their run's commit anyway. A docs page covers trajectories
and release.

**Replay, for the published endpoints.** Run from a pinned commit in a detached worktree.
- Panel B: `predict(d11_models(), pd; horizon = 6300, draws = rows, completion = :prior, n_held = 20,
  n_cells = 4, width = 1, rng = Xoshiro(19_542_000), trajectory = …)`. V's seeds are drawn after the cell
  seeds and the knockdown draws nothing, so the 8,000 baseline cells should match KIM-42's. At width 1
  V rejected nothing, so leaving V out changes no weight. **Check:** `p.mu == w1.mu` exactly, on Zen3
  nodes, from the archive's `runs/budget-kim42/`. A mismatch stops the work.
- Panel A: M0's `work` tasks for replicate 1 and the exact posterior, baseline only, 10 blocks of 10
  draws × 887 cells, plus the 9,849 truth cells. **Check:** each row's μ̂ and the truth's mean equal the
  archived values exactly. The builder reads the seeds from `dev/scripts/m0_evidence_kim40.jl`.
- Panel C: a new run with a new seed from M5's unused range, 40 draws × 50 held × (4 + 4) cells, under
  the frozen `ViabilityRule`. More columns than panel B, since 62 more parameters make 20 columns too
  few and V acts by column. A crossed matrix cannot be topped up, and starvation only shows over the
  full cycle, so a pre-check comes first: 20 pairs × 8 cells at 6,300 s, about 3 CPU-h. **Whatever the
  acceptance, the run goes ahead and the number is reported,** since heavy starvation is itself what
  release shows. Below 10% of pairs kept, the figure drops the filled conditioned band and shows the
  unconditioned one, labelled. The released set does not change without a §12 entry.

**Cost.** From earlier run rates, not measured for these runs:

| Run | Cells | CPU-h |
|---|---|---|
| A, replicate 1 | 88,700 + 9,849 | about 60 |
| B, width 1 | 8,000 | about 150 |
| C | 16,000 | about 290 |
| **Total** | | **about 500** |

Rates: M0 about 2.3 core-s per cell (`m0_evidence_kim40_result.md`), the reduced model about 65 core-s
(KIM-42's 1,295 h `TotalCPU` over about 72,000 cells). Storage: on the order of 1 GB at 60 s for the
two species plus the panel. Keep all ODE species only if the builder confirms it is small.

**Order and records.**
1. Software PR to `main`, gated with `/check-PR`, with a Linear issue.
2. The §12 entry fixing this design, then a smoke run at 630 s, then the Slurm runs, parse-checked first.
3. A result note with job ids, an archive outside `.worktrees/`, and `paper/figures/data` TSVs. The
   figure reads the TSVs, as the others do.
4. A paper PR: `make_figures.py`, `numbers.tex` entries by `\pn{}`, text in §3 and the abstract, and
   `paper/STATUS.md`.

## The builder's calls

- **No wrapper over the constructors.** `free =` plus the existing error message is short, and it works.
  A wrapper would have to reconstruct modules generically. We can add one if users ask.
- **Bands are empirical quantiles,** not mean ± SD, since Q is skewed once cells starve.
- **The 82 are fixed as above,** not chosen by effect size. Single-parameter release goes in the docs and
  the listing. The figure shows one group.
- **GTP as the second species in every panel,** because it is observed on M0 and ties to the K_d caveat.
- **Pyruvate as panel A's unobserved species,** fixed here before any run, because it is the end of
  glycolysis and is not in `M0_PANEL`.
- **Panel C's weights, its below-10% rule and the regenerated data,** as written above.

## §6: what the pilot showed

Throwaway runs, not results. Jobs 20231333, 20231949 and 20232368, at `d0a57c1`. 16 pairs, one cell each,
6,300 s, on 2 rows of V7 replicate 1.

| Drawn | Failed | Starved | Q, median (min–max) |
|---|---|---|---|
| 20 | 0 | 1 of 16 | 9.96e6 (5.69e6–1.18e7) |
| 60 (+ glycolysis) | 0 | 6 of 16 | 5.38e6 (7.39e5–9.87e6) |
| 82 (+ recycling) | 0 | 8 of 16 | 4.27e6 (8.17e5–1.46e7) |

So release works through the public API, and the framework's error message caught our own missing
`kcatR_*`. Starvation becomes common, and the rule that rejects it will do real work in panel C: at
`starved_fraction = 0.5`, acceptance could be low and panel C's cost per kept pair higher.

Priors on the reduced model: of 260 parameters, 203 are log-normal. Balanced `kcatF` have a median
geometric SD of 1.05 (`kcatF_R_LDH_L` 1.05), and `km_*` a median of 3.6 and a maximum of 9.4. **So
releasing one kcat barely moves a band, and the Michaelis constants are where the width is.** The listing
above is therefore a usage example, and not a claim that LDH matters.

## Decided by the author, 2026-11-07

**How the paper describes M0's score on Q.** The data include a noisy measure of Q, through
`M_lac__L_e`, so the Q coverage is not a test on an unmeasured quantity. Methods lists the pool, but the
results do not draw the connection. How far that measure drives the posterior has not been measured.
ΔQ is the quantity the data did not see. **The M0 results will say both of these plainly, and no run is
added for it.** The paper PR makes this change, with the abstract and discussion checked for the same
reading.
