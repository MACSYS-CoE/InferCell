# Story: InferCell, inference on a reduced whole-cell model of JCVI-syn3A

**Mode:** revise (restructure of the M0 draft on main at 24273c7)    **Venue:** none (preprint, two-column `arxiv.sty`)
**Source:** Linear project InferCell, description as of 2026-11-23 (contribution, section jobs, F1–F7, T1; F7 dropped with the ladder on 2026-11-25 and deleted on 2026-11-26, KIM-163); Decision log to 2026-11-23
**Signed off:** 2026-11-23 (Tom, by approving the KIM-98 drafting plan)
**Branch:** `tomkimpson/kim-98-plan-the-paper-1-rewrite-m6-around-reduced` (KIM-98), worktree `.worktrees/kim-98-paper`

## Contribution

InferCell infers a reduced whole-cell model of JCVI-syn3A across a coupled ODE–stochastic boundary and turns the posterior into predictions scored against held-out truth.

(The ladder clause was removed on 2026-11-25; see the amendment log.)

**The question:** can a whole-cell model with a real ODE–stochastic boundary be inferred from single-cell and bulk data, and how far do its posterior predictions carry?

## Section jobs

| Section | File | Job in the argument |
|---|---|---|
| Abstract | sections/00-abstract.tex | (via /write-abstract) contribution |
| 1 Introduction | sections/01-intro.tex | Whole-cell models give point predictions; mixed formalisms and latent paths block inference; InferCell carries a reduced whole-cell model from data to scored probabilistic prediction. |
| 2 InferCell | sections/02-infercell.tex | The framework: declared priors and reductions, two descriptions across the boundary, two inference routes, predict. |
| 3 The test system | sections/03-testsystem.tex | The reduced syn3A model, and M0 as the slice small enough for an exact check. |
| 4 The inference engine, checked on M0 | sections/04-m0.tex | Both routes pass SBC on M0; the mixed route is about 80× cheaper; the remaining gap is stated with its leads. Short. |
| 5 Inference on the reduced model | sections/05-inference.tex | One joint posterior over 21 parameters: the data reach all 21, the truth sits inside 20 of 21 intervals, the posterior reproduces its own data, and it costs about 492 CPU-h. The modular arm sits beside it (todo). |
| 6 From posterior to scored prediction | sections/06-prediction.tex | The later horizon and the ptsG knockdown, scored against the true simulator against the prior predictive; Q's uncertainty budget now comes from inference, not priors. |
| 7 Discussion | sections/08-discussion.tex | What carries (Q, transcripts) and what does not (pools, ΔQ); limits (one truth, synthetic data, well-specified model, the promoter shift); next: real data, a degradation ladder, the full cell. Conclusions resting on KIM-116 are red todos until KIM-162. |
| Availability | sections/09-availability.tex | as now |
| App. Methods | sections/app_methods.tex | Reproduction: driver, observation model, exact sampler, transcript factor, mean-field route, SBC, prediction, compute; add the reduced-model design (21 coordinates, fit window, warm metric). |
| App. M0 in detail | sections/app_m0.tex | M0's exact sampler, holding S fixed, the mean-field swap, what each stream informs, the crossed comparison, the gap leads (todo, KIM-90). |
| App. M0 predictions | sections/app_m0predict.tex | M0's scored predictions, the KIM-42 transferred budget, trajectories and releasing a constant. |
| App. Declared reductions | sections/app_reductions.tex | The model-form table. |
| App. Modular arm | sections/app_modular.tex | todo (KIM-116). |

## Floats

| ID | File | Source script | Claim (caption's first sentence) | Where | Width |
|---|---|---|---|---|---|
| F1 `fig:architecture` | figures/fig_architecture.pdf | figures/make_schematics.py | InferCell composes modules with declared priors and reductions. | main §2 | full |
| T `tab:frameworks` | inline | tables/frameworks_sources.md | (as now) current frameworks lack the four design features | main §1 | column |
| F2 `fig:model` | figures/fig_model.pdf | figures/make_schematics.py | The boundary in the reduced model is real and two-way. | main §3 | full |
| F3 `fig:approx` | figures/fig_mixed.pdf | figures/make_figures.py | Both routes calibrate on M0, and the mixed route costs about 6.6 CPU-h per posterior. | main §4 | full |
| T `tab:cost` (M0 cost) | inline | — | (as now; merges into T1 if it reads better) | main §4 | column |
| F4 `fig:reducedpost` | figures/fig_reduced_posterior.pdf (new) | figures/fig_reduced_posterior.py, data/kim92_*.tsv | The data reach all 21 parameters of the reduced model, and the posterior reproduces its own data. | main §5 | full |
| T1 `tab:reducedcost` | inline | numbers from KIM-92 note; M0 rows from numbers.tex | The mixed route makes a reduced-model posterior affordable at about 492 CPU-h. Rows: M0 exact, M0 mixed, reduced joint, reduced modular (todo, KIM-116). | main §5 | column |
| F5 `fig:reducedpredict` | figures/fig_reduced_predict.pdf (new) | figures/fig_reduced_predict.py, data/kim95_score_targets.tsv | The posterior's predictions cover the truth on Q and the transcripts and miss ten of the seventeen pools by 1–1.4 SD, against a prior predictive that is far wider. | main §6 | full |
| F6 `fig:reducedbudget` | figures/fig_reduced_budget.pdf (new) | figures/fig_reduced_budget.py, data/kim96_/kim42_budget_rows.tsv | Once every parameter is inferred, none stays prior-dominated; at the default prior width Q's parameter variance falls from 1.15e13 through KIM-42's 2.17e12 to 1.29e10 (wording: Decision log 2026-11-19, with the width and attribution qualifiers). | main §6 | column |
| appendix floats | fig_calibration, fig_graphs, fig_predict_m0, fig_crossed, fig_budget, fig_trajectories; tab:edges (main §2), tab:gain, tab:budget, tab:crossed, tab:modelform, tab:compute | existing | as now | appendices | as now |

## Claims and evidence

| # | Claim | Evidence | Section |
|---|---|---|---|
| C1 | Both routes calibrate on M0 (SBC) | F3, fig:calibration, existing `\pn{}` M0 keys | §4 |
| C2 | The mixed route is about 80× cheaper than the exact sampler on M0 | tab:cost, existing keys | §4 |
| C3 | Three M0 parameters sit 0.42–0.56 SD from the exact reference; the leads are double counting and the bulk factor, neither measured as a shift | KIM-90 note (unmerged) → todo until merged | §4, App. M0 |
| C4 | The joint posterior converges after one extension (R̂ ≤ 1.0102, ESS ≥ 2,302) | \pn{kim92.rhat}, \pn{kim92.essmin} | §5 |
| C5 | The data reach all 21 coordinates (contraction ≥ 0.919); truth in the 90% interval on 20 of 21 | F4, \pn{kim92.contrmin}, \pn{kim92.truthin} | §5 |
| C6 | The posterior reproduces its own data: 11 of 884 rows outside the 99% band (limit 44) | F4, \pn{kim92.ppcout} | §5 |
| C7 | The promoters sit low (mean z −0.81, 14 of 17), the MAP shares it, so it comes from the likelihood | \pn{kim92.promz}, \pn{kim92.mapz} | §5 |
| C8 | A reduced-model posterior costs about 492 CPU-h | T1, \pn{kim92.cpuh} | §5 |
| C9 | Modular arm vs joint | todo (KIM-116) | §5, §6, App. |
| C10 | Q and Q_kd are covered by ±1.4% intervals, skill 0.98; the prior is 24–26% low | F5, kim95.* | §6 |
| C11 | 900 of 901 transcript saves inside 90%; transcripts are mostly cell noise | F5, kim95.* | §6 |
| C12 | 10 of 17 pools shifted 1–1.4 SD (trioses high, energy charge split), θ-driven | F5, kim95.* | §6 |
| C13 | ΔQ is not resolved beyond "small next to cell noise" | kim95.dqmc | §6 |
| C14 | Q's parameter variance falls 1.15e13 → 2.17e12 → 1.29e10; ~170× (88–321) at default width; 44–1,140 over widths ½–2; not only the promoters | F6, kim96.* | §6 |

## Appendices: what is relegated, and why

- M0's SBC detail, the exact sampler, holding S fixed, the mean-field swap, what each data stream informs, the crossed comparison → App. M0 in detail. Why: §4 is a sanity check of the engine; the argument needs only that both routes calibrate and what they cost.
- M0's scored predictions, KIM-42's transferred budget, trajectories and releasing a constant → App. M0 predictions. Why: the main text predicts from the reduced model's own posterior; KIM-42 enters the main text only as F6's middle state.
- The declared-reductions table → App. Declared reductions. Why: reproduction and honesty about model form; one main-text sentence points to it.
- Modular arm detail → its own appendix (todo, KIM-162).

## Off-story (revise mode)

- 03-results "Probabilistic predictions on the reduced model, unscored" — moved to App. M0 predictions; superseded by §6. Shorten there.
- 03-results "What the mean-field swap costs" and "What each data stream informs" — moved to App. M0 in detail; §4 keeps one sentence each.
- 01-intro "This work" paragraph — rewrite: it describes the old three-step M0 story.
- 08-discussion — rewrite opening around the reduced model; keep M0 working rules only as they bear on the reduced model.
- Abstract — rewrite.

## Terms

| Use | Not |
|---|---|
| the reduced syn3A model, the reduced model | Core A′, corea |
| mixed route | hybrid route, approximate route |
| exact sampler | reference sampler |
| joint posterior | full posterior |
| modular arm | cut model (as a name for the arm) |
| later horizon | forecast window |
| contraction | shrinkage |
| detectable discrepancy, declared reduction | resolved discrepancy |
| the posterior-to-prediction step | same pipeline |
| true simulator | ground-truth model |
| degradation ladder, rung | ladder of realism, degradation level |
| factor cut, data plug-in (the two modular variants) | stage-1 cut, plug-in arm |
| the budget transferred from M0 (F6's middle state) | KIM-42 state (in prose) |
| MC share | noise fraction |
| CRPS skill | skill score |

## Notation

| Symbol | Meaning | First defined in |
|---|---|---|
| $\theta$ | the free parameters (21 on the reduced model) | §2 |
| $\theta^\ast$ | the true parameters that made the synthetic data | §3 |
| $Q$ | extracellular lactate: lactate a cell exports over the 6,300 s cycle, in particles, as a population mean over cells (`:lactate_out`) | §2 (sec:predictapi) |
| $Q_{\mathrm{kd}}$, $\Delta Q$ | $Q$ under the ptsG knockdown, and $Q_{\mathrm{kd}}-Q$ | §6 |
| $\sigma_{\mathrm{b}}$ (`\sigb`) | bulk assay noise | §2 / Methods |
| $\widehat{R}$ (`\Rhat`) | rank-normalised R-hat | §4 |
| $\mathrm{Var}_\theta(Q)$ | Q's variance across parameter draws (parameter variance) | §6 |
| $\theta_I$, $\theta_H$ | inferred and held sets of the cut | §2 (sec:predictapi) |
| $V$ | viability rule | §2 (as now) |

## Budgets

- main_words: 9500
- mean_max: 22
- long: 35
- long_max: 8
- abstract_max: 250
- term_sections: main
- allowed_numbers: 3, 17, 21, 50, 6300, 6,300

## Stages

- [x] 0 Orient (mode: revise/restructure, branch: tomkimpson/kim-98-plan-the-paper-1-rewrite-m6-around-reduced)
- [x] 1 Story signed off (2026-11-23)
- [x] 2 Scaffold (check_paper passes)
- [x] 3 Figures (F4, F5, F6 new; F1–F3 kept)
- [x] 4 Draft (waves 1–3)
- [x] 5 Edit passes (invariants ok)
- [ ] 5b Pangram — deferred (user, 2026-11-23: run stages 0–6 now)
- [x] 6 Joins pass
- [ ] 7 References and compliance — deferred
- [ ] 8 Cold read, final checks, report — deferred

## Pending results (red `\todo{}` until they land)

| What | Issue | Where in the paper |
|---|---|---|
| Modular arm: posterior vs joint, cost, plug-in D3 scores | KIM-116 | §5 paragraph, T1 row, §6 line, App. modular |

## Pangram

| Section | Words | Score after Stage 5 | Score after extra pass |
|---|---|---|---|
| (deferred) | | | |

## Figures review

| ID | Outcome | Note |
|---|---|---|

## Amendment log

### 2026-11-25 — the ladder leaves paper 1; the framing is decided
**Trigger:** Decision log 2026-11-25; reasoning in `framing_decision.md`.
**Change:** the contribution sentence loses "then shows how both degrade as the synthetic data become realistic". §7, F7, C15, App. Degradation ladder and the M10 row of "Pending results" go. The ⚠️ OPEN note under Contribution is closed. The framing sentence around the contribution waits for KIM-116. The tables above were edited to match on 2026-11-26 (KIM-163), when §7, F7 and the ladder appendix were deleted; the Discussion is now §7. The M0-gap row of "Pending results" went with KIM-160 (PR #192).
**Agreed by:** Tom, 2026-11-25.

### 2026-11-23 — §4's "gap is diagnosed"
**Trigger:** Decision log 2026-11-23 (KIM-90): the cross-covariance check names nothing.
**Change:** §4's job reads "the remaining gap is stated with its leads", not "diagnosed". The leads are red todo until KIM-90 merges.
**Agreed by:** Tom, 2026-11-23 (flagged as OPEN in the approved plan).

### 2026-11-23 — precision fixes from wave 1 (no claim changes)
**Trigger:** the figure and section agents read the archives directly.
**Change:**
- C3: the three M0 gap parameters are the ptsG promoter (S_0779), the decay rate (krnadeg) and ENO, per the KIM-90 note and `mix.gaplo`; S_0607 enters only the double-counting lead.
- C10: skill is 0.985 on Q and 0.983 on Q_kd, written "0.98–0.99".
- C12: the ten shifted pools run 0.94–1.40 SD ("about 0.9–1.4"); "θ-driven" holds for the trioses and FDP (MC share 0.02–0.11), not the nucleotides (0.31–0.62).
- C8: the 492 CPU-h is the whole run (dataset, profiles, the abandoned 24 CPU-h probe, PPC); T1's caption says so.
- C2: "about 80×" is 77× against the production chain and 103× against V7's mean; kept as "about 80".
**Agreed by:** within the signed-off claims; to be confirmed by Tom at review.
