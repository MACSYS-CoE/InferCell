# Manuscript Review: InferCell (Paper 1)

**Manuscript reviewed:** `paper/main.tex` (and `paper/sections/*.tex`, `paper/figures/`, `paper/references.bib`, `paper/STATUS.md`)  
**Target venues evaluated:** *Nature Computational Science* and *Cell Systems*  
**Date:** 2026-11-25

---

## 1. Executive Verdict & Journal Fit

**Bottom line:** The core problem is field-defining and the methodological rigor is unusually high, but **the manuscript is not yet positioned to clear editorial triage at *Cell Systems* or peer review at *Nature Computational Science***—even after completing the current `\todo{}` blocks (KIM-90, KIM-116, M10)—without addressing a small set of high-leverage scientific and framing gaps.

| Target Journal | Realistic Verdict (if completed as currently designed) | Rationale |
| :--- | :--- | :--- |
| ***Nature Computational Science*** | **Best natural fit.** Currently **Reject with encouragement to resubmit / Major Revision**; **genuinely competitive** after closing the gaps below. | *Nat. Comput. Sci.* actively publishes foundational inference frameworks, hybrid UQ methods, and scientific software even on synthetic benchmarks—*provided* the approximation's errors are diagnosed, scalability and generality are defended, and modern alternatives (SBI, Laplace/VI) are positioned. |
| ***Cell Systems*** | **High risk of desk rejection** (likely transfer offer to *Patterns*, *iScience*, or *PLOS Computational Biology*). | *Cell Systems* demands either new biological insight or immediate utility on real experimental data. Right now, the paper explicitly disclaims real data, disclaims biological claims about *JCVI-syn3A*, removes 92% of the cell, and infers 21 parameters on 1 synthetic dataset. |

**Recommendation:** Make ***Nature Computational Science*** the primary target unless a real-data application or biological discovery on *JCVI-syn3A* is added specifically for *Cell Systems*.

---

## 2. What Is Already Top-Tier (Protect and Amplify)

1. **A genuinely important, longstanding open problem.** Every landmark whole-cell model (*M. genitalium*, *E. coli*, *JCVI-syn3A* 2022/2026) remains a forward point-estimate simulator. Bridging the ODE–stochastic boundary for calibrated Bayesian inference directly answers the challenge posed since DREAM8 (Karr et al., 2015) and Babtie & Stumpf (2017).
2. **Exceptional statistical rigor and reproducibility.** Pre-registered kill criteria, Simulation-Based Calibration (SBC) with simultaneous bands, an exact Particle Gibbs + endpoint-conditioned CTMC reference sampler to benchmark the fast mixed NUTS route, bit-for-bit path replay, Slurm compute accounting, and programmatic `\pn{}` number provenance are far above the norm in computational biology.
3. **Sharp, memorable conceptual contributions.**
   - **Declared reductions** as first-class computable objects rather than buried supplement prose.
   - **Statistical cuts vs. dynamical reductions** across a boundary with *no shared parameters*.
   - **The "Three Graphs"** (coupling vs. likelihood vs. information graph, `sections/app_m0.tex` §`sec:info`)—a flagship-level systems-biology insight (see Major Issue #6).
4. **Publication-grade data figures.** Figures 3–6 (`fig_mixed.pdf`, `fig_reduced_posterior.pdf`, `fig_reduced_predict.pdf`, `fig_reduced_budget.pdf`) and `fig_trajectories.pdf` are crisp, legible, and clean.

---

## 3. Six Major Issues That Will Block Acceptance (Ranked by Impact)

### Issue 1: Unresolved Biases in the Core Approximation ("Mixed Route")
The flagship methodological contribution is the mixed route (~80× faster than exact PGAS), and on the reduced model the paper relies on it *without* an exact benchmark. Currently, three unresolved discrepancies sit side-by-side:

- **On M0 (`sections/04-m0.tex`, lines 63–78):** 3 of 6 parameters sit `0.42–0.56 SD` from the exact posterior, and the text states *"We have not found the gap's cause."*
- **On the reduced-model fit (`sections/05-inference.tex`, lines 47–58):** 14 of 17 promoter means sit below the truth (mean standardized error `−0.81 SD`, **3.3 SEs from zero**), and the text states *"One dataset cannot separate the two possible causes... A second dataset would tell them apart, and we have not run one."*
- **On reduced-model prediction (`sections/06-prediction.tex`, lines 56–72):** 10 of 17 metabolite pools miss by `0.9–1.4 SD` (G3P `+1.39`, DHAP `+1.23`, plus systematic energy-charge splitting across ATP/ADP/AMP and GDP/GMP).

**Why reviewers will reject on this:** Leaving an unexplained `3.3 SE` systematic bias in the headline fit—and explicitly stating *"a second dataset would tell them apart, and we have not run one"* when a fit costs 492 CPU-h on OzSTAR—reads as unfinished rather than thorough.

**High-value scientific opportunity:** In `sections/08-discussion.tex` (lines 67–73) and `sections/app_methods.tex` (lines 93–105), **99.3% of baseline cells clip GTP to zero** (crossing below $K_d$), whereas the mean-field ODE *does not clip*, evaluates fixed $\delta, S$ at a single MAP point (found at $\delta=0, S=0$), and ignores double-counting of the same cells across both likelihood factors. Diagnosing whether the `−0.81 SD` promoter shift and `10/17` pool shifts come from (a) a single unlucky noise draw, (b) fixed $\delta$ / GTP clipping below $K_d$, or (c) double-counting across streams would turn a vulnerability into **one of the strongest methodological findings in the paper**.

### Issue 2: The `N = 1` Reduced-Model Dataset & an Uninformative Knockdown
- **Single truth (`N = 1`):** Everything in Sections 5 and 6 rests on a single dataset and single truth $\theta^\ast$. Running even **3–5 replicate truths** on the reduced model (or using the modular/data-plug-in route if KIM-116 validates it as cheap and accurate) would separate random dataset noise from systematic bias and substantially strengthen Sections 5–6.
- **Unresolved intervention (`sections/06-prediction.tex`, lines 29–44):** The headline out-of-sample intervention (30% *ptsG* knockdown) changes $Q$ by only **0.4%** (`6.4 × 10^4` out of `1.5 × 10^7`), completely buried in cell noise (`MC share = 1.09`). As the text admits, *"covering $Q_{\mathrm{kd}}$ tests nearly the same thing as covering $Q$."*
  - **Fix:** Add an intervention that moves the cell in a biologically and statistically nontrivial way (e.g., 70–90% CRISPRi-style knockdown of *ptsG*, *eno*, or *gapd*, or a step change in extracellular glucose). Showing calibrated out-of-sample prediction of a *large, non-obvious* perturbation response (15–40% shift in flux/pools) is crucial for *Nat. Comput. Sci.* / *Cell Systems* impact.

### Issue 3: Parameter Scope (Only 2 Metabolic Constants) & First-Order Generality
Reviewers will flag two scope bottlenecks immediately:

1. **18 of 21 free parameters are constitutive gene-expression rates; only 2 (ENO, FBA $k_{\text{cat}}$) are nonlinear metabolic constants (`sections/03-testsystem.tex`, lines 52–55).** The ~62 other metabolic constants remain fixed at their nominal values. A referee will ask: *With exact 1-minute single-cell transcript counts, inferring 17 linear promoter rates and 1 decay rate is nearly a decoupled textbook birth–death problem. Does the mixed route still work when 10–30 sloppy, correlated metabolic parameters are free—where whole-cell non-identifiability (Gutenkunst 2007) actually lives?*
   - **Fix:** Explicitly justify why ENO and FBA were chosen, and ideally run a stress test freeing a broader block of metabolic constants (e.g., all glycolytic $k_{\text{cat}}$s) to demonstrate how InferCell handles sloppy metabolic manifolds (or clearly characterize where NUTS/warm-metric scaling breaks down).
2. **The first-order gene-expression restriction (`sections/02-infercell.tex`, lines 104–109, 189):** Exact closure of the mean-field ODEs and the exact birth–death transcript factor rely on gene expression being strictly first-order (constitutive, no feedback regulation, no Hill repressors/activators, no transcriptional bursting). That holds largely in minimal *JCVI-syn3A*, but not in *E. coli*, yeast, or mammalian cells.
   - **Fix:** Address this directly in Section 2 and the Discussion: show how the compositional interface accommodates nonlinear gene regulation (e.g., swapping the exact birth–death factor for FSP, LNA/moment closure, or a local neural likelihood surrogate inside the same typed-edge contract).

### Issue 4: Missing Positioning & Comparisons Against Modern Alternatives
A *Nature Computational Science* referee will immediately raise two sharp questions:

1. **Why pay 492 CPU-h (5.5 days wall-clock) for NUTS if Laplace matches it within 4–8%?** In `sections/05-inference.tex` (lines 71–75), Laplace SD / NUTS SD is `0.96–1.08` across all 21 parameters. Explicitly address why NUTS was still necessary—e.g., (a) the MAP optimizer missed the ENO and FBA posterior means by `2.3 SD` and `1.4 SD` (investigate whether this is skewness or optimizer stopping tolerance), and (b) along the degradation ladder (M10), sparser/noisier data will likely produce non-Gaussian ridges where Laplace fails and NUTS succeeds. Turning that into an explicit test makes the 492 CPU-h compute cost feel justified rather than excessive.
2. **Where does InferCell sit relative to neural Simulation-Based Inference (SBI), LNA, and linked emulators?** **27 of the 56 entries in `paper/references.bib` are never cited in the text**—including `cranmer2020frontier`, `geffner2023compositional`, `golightly2011bayesian`, `schnoerr2017approximation`, `lubitz2010parameter`, `lindsten2014particle`, `lindensantangeli2025increasing`, `kyzyurova2018coupling`, and `ming2021linked`.
   - **Fix:** Add a concise 2-paragraph positioning in the Introduction or Discussion contrasting InferCell's hybrid differentiable-likelihood approach against (i) black-box neural SBI (where 21+ parameters + 884 bulk rows + 200-cell trajectories at 32 s/cell makes amortised simulation budgets prohibitive without structural likelihood factorization), (ii) Bayesian parameter balancing (`lubitz2010parameter`, currently discussed in §1 without citing it), and (iii) recent Bayesian multimodel/linked-emulator methods (`lindensantangeli2025increasing`, `ming2021linked`). Also cite `lindsten2014particle` in `sections/app_methods.tex` (line 67, ancestor sampling) and Gneiting & Raftery (2007) for CRPS (`sections/app_methods.tex`, line 261).

### Issue 5: Bridging the "Purely Synthetic / Well-Specified" Gap
Currently, synthetic data come from the *exact same reduced model* being fitted, with exact single-cell transcript counts every 60 s (an assay that does not physically exist for 17 genes + 17 bulk pools simultaneously).

- **Completing M10 (the degradation ladder, `sections/07-ladder.tex`) is mandatory** before submission to either journal. Realistic rungs should include **destructive/snapshot scRNA-seq** (different cells at each time point, removing single-cell temporal autocorrelation that currently drives 98.9% of `krnadeg` shrinkage), **capture dropout / binomial thinning**, **sparser sampling intervals**, and **subsetted metabolite panels**. Identifying *which exact assay limitation kills cross-module identifiability first* gives experimentalists and systems biologists a concrete, citable result.
- **High-impact upgrade options:**
  - **Option A (structural misspecification test using companion paper):** Since a fast simulator of the **full well-stirred *JCVI-syn3A* model** (`kimpson2026whole`) exists, generate synthetic data from the *full* cell model and fit the *reduced* model. Testing InferCell under genuine structural misspecification—and showing how declared reductions and predictive checks detect it—would substantially elevate the paper.
  - **Option B (real-data case study for *Cell Systems*):** Fit a real experimental dataset (e.g., published *JCVI-syn3A* proteomics/metabolomics/growth constraints from Breuer 2019 / Thornburg 2022, or a well-characterized coupled metabolism–expression dataset in *E. coli* / yeast). Without *any* real data, *Cell Systems* remains a very steep climb.

### Issue 6: Buried Conceptual Gems & Over-Defensive Tone
The prose is admirably clear, concise, and free of hype, but two structural choices currently *undersell* the work to editors:

1. **Promote the "Three Graphs" and Elasticity Working Rule from Appendix B to the Main Text.** `sections/app_m0.tex` (lines 99–149, `fig_graphs.pdf`, Table `tab:gain`) contains arguably the paper's most generalizable scientific insight: *a two-way physical coupling graph collapses into a one-way information graph because saturated enzyme rate laws ($[\text{ATP}] \gg K_d$) kill the reverse elasticity, and a modeller can diagnose cross-module identifiability analytically before running a single MCMC chain.* Burying this in Appendix B under *"No claim in the main text rests on this subsection"*—while `sections/08-discussion.tex` (lines 29–33) explicitly builds its Discussion on it—is both a logical contradiction (already logged in `paper/STATUS.md`) and a missed opportunity for impact.
2. **Dial back self-undermining phrasing while keeping scientific honesty:**
   - **Table 1 footnote (`sections/01-intro.tex`, lines 93–101):** Remove *"We did not read the supplements"* (verify the supplements of those 4 papers so Table 1 is bulletproof) and reframe neutrally around design scope rather than apologizing for the comparison criteria. Update *"for the 2026 syn3A model we read the preprint"* once checked against the published *Cell* 2026 paper (`thornburg2026bringing`).
   - **Abstract (`sections/00-abstract.tex`, lines 18–23):** Ending on *"The data are synthetic, from one truth... show where it is wrong"* lands flat. State the synthetic scope cleanly earlier, and close on the forward-looking scientific advance (what compositional inference enables for experimental design, uncertainty-aware perturbation prediction, and scaling to full whole-cell models).
   - **Code snippet (`sections/02-infercell.tex`, lines 225–238):** Currently demonstrates releasing fixed constants for an appendix-only sensitivity check (`n_held = 20`) rather than the main paper's core `infer` → `predict` workflow. Replace with a clean 6-line snippet of the main paper's workflow, or move the release snippet to `sections/app_m0predict.tex`.

---

## 4. Presentation, Figures, and Technical Polish

- **Upgrade Schematic Figures 1 and 2 (`fig_architecture.pdf`, `fig_model.pdf`):** Unlike the sharp data figures, these two schematics read like draft matplotlib boxes:
  - In `fig_architecture.pdf`, the `predict` box crowds the short vertical arrows emerging from `Mean field` and `Mixed route`.
  - In `fig_model.pdf`, the four horizontal arrows sit evenly between columns rather than connecting to specific modules, forcing the caption to warn *"Arrows run between the blocks, not the boxes they touch."* Redrawing Figs. 1 and 2 with clean module/block enclosures, matched palette styling, and explicit 1 s / 60 s clocks will noticeably lift editorial first impressions.
- **Unfinished markers to clear (current counts):**
  - `22 \todo{}`, `2 \placeholder{}`, `1 \DRAFT{}` (M0 mixed-route wall-clock in `sections/05-inference.tex`, line 137), and `1 [CITE NEEDED]` (`sections/app_methods.tex`, line 261).
- **Housekeeping items:**
  - Resolve the open items already tracked in `paper/STATUS.md` (naming of $Q$ vs. "lactate exported"; clarifying that the prior state in Fig. 6 draws 20 parameters rather than 21; reconciling gradient timing `37.9 s` vs. `43–45 s`).
  - In `tab:reducedcost` (`sections/05-inference.tex`, line 126), separate the production posterior cost from the exploratory `24.3 CPU-h` unit-metric probe and PPC (or state cleanly in a breakdown column/footnote).
  - Remove redundant `\label{sec:results}` / `\label{sec:model}` pair in `sections/03-testsystem.tex` (lines 2 and 12).

---

## 5. Concrete Roadmap to Maximize Impact

Prioritized sequence for the strongest shot at ***Nature Computational Science***:

1. **Finish the planned triad (M10 + KIM-116 + KIM-90)**—ensuring KIM-90 *diagnoses* the M0 gap (e.g., splitting cells into disjoint transcript vs. bulk cohorts cleanly tests double-counting).
2. **Diagnose the reduced-model `−0.81 SD` promoter bias & `10/17` pool shifts** by running **2–4 additional synthetic truths** on the reduced model and testing whether GTP clipping below $K_d$ vs. fixed $\delta$ drives the bias.
3. **Add one strong, resolvable perturbation** in Section 6 (e.g., 80–90% knockdown of a rate-limiting glycolytic enzyme or an extracellular glucose shift) so the out-of-sample intervention test meaningfully tests predictive calibration.
4. **Add a compact scaling & method comparison:**
   - Free a larger block of metabolic constants (e.g., 10–15 glycolytic $k_{\text{cat}}$s) to test handling of sloppy metabolic subspaces.
   - Compare Laplace vs. NUTS along the M10 degradation ladder to empirically demonstrate where full MCMC outperforms Laplace.
5. **Restructure the narrative for flagship impact:**
   - Promote the **Three Graphs + elasticity working rule** (`fig_graphs.pdf`, `tab:gain`) into the Main Text.
   - Integrate the ~27 uncited literature references already in `paper/references.bib` (SBI, parameter balancing, LNA/stochastic kinetics, linked emulators).
   - Redraw Figs. 1 and 2 to journal-grade schematic quality and tighten defensive phrasing in the Abstract and Table 1 footnote.
6. *(Optional, but decisive if specifically targeting **Cell Systems**)*: Add either a **structural misspecification test** (fitting the reduced model to data from the companion paper's full *JCVI-syn3A* simulator) or a **real experimental data fit** (e.g., published *syn3A* omics/growth constraints).

---

## 6. Strategic Follow-Up: Full *JCVI-syn3A* (`PreCell` / `SpeedyCell.jl`) & Compositional / Targeted Inference

### Point A: Reframing Around Targeted / Per-Module Inference (Not Monolithic Joint Inference)

**Assessment:** Positioning InferCell around **compositional, targeted, and per-module inference**—rather than monolithic joint inference over all parameters—is a significantly stronger, more original, and more defensible thesis. Currently, the manuscript accidentally argues the opposite:

- The Abstract and Section 5 lead with *"one joint posterior over all 21 free parameters."*
- Section 6 (`fig_reduced_budget.pdf`) frames partial/modular inference (*"Five inferred, 15 at prior"*) as an inferior intermediate step *"superseded"* (`sections/app_m0predict.tex`, line 12) by full joint inference.
- The `free = [...]` / held-parameter ($\theta_H$) mechanism—which directly implements conditioning on *in vitro* / literature measurements (e.g., BRENDA enzyme kinetics, known replication rates) while inferring only remaining *in vivo* unknowns—is framed defensively (*"not a joint posterior... a cut... a sensitivity to priors"*).

**Why reframing around targeted / per-module inference elevates impact:**
1. **Matches real biological practice.** Monolithic joint inference over all ~5,000+ parameters of a whole-cell model is both computationally intractable and statistically ill-posed (sloppiness, non-identifiability, incompatible assay timescales). Biologists want to fix or prior-anchor parameters measured *in vitro*, infer only the unknown subset conditioned on those values, and quantify residual uncertainty.
2. **Unifies InferCell's existing components into a coherent design philosophy:**
   - **Targeted conditioning ($\theta_I$ inferred, $\theta_H$ fixed or held at prior):** A primary scientific workflow rather than a compromise.
   - **Uncertainty budgets (`predict` + variance decomposition):** Becomes an experimental-design tool identifying whether phenotypic uncertainty (e.g., doubling time or knockdown response) is governed by inferred parameters or unmeasured *in vitro* priors.
   - **The "Three Graphs" + coupling elasticities:** Provides the principled mathematical criterion for deciding *when* per-module inference (statistical cut) suffices vs. *when* strong cross-module feedback requires joint inference across a boundary.

---

### Point B: Dropping the Reduced Model (`Core A'`) in Favour of `M0` + Full *JCVI-syn3A* (`PreCell` / `SpeedyCell.jl`)

**Assessment:** If technically feasible in a reasonable timeframe, replacing the 3-tier hierarchy (`M0` → `reduced syn3A` → `full syn3A` in Paper 2) with a **2-tier architecture (`M0` + `Full JCVI-syn3A`)** would make the manuscript substantially cleaner and much higher impact.

**Why the reduced model currently creates editorial drag:**
- Requires 93 declared reductions (`S1_File_reductions.tsv`) just to document removing 92% of the cell.
- Forces disclaimers (*"we make no claim about JCVI-syn3A itself"*, refusing to report doubling time).
- Likely introduces artificial dynamical pathologies (e.g., `99.3%` of cells clipping GTP to zero).
- Causes flagship editors (*Cell Systems*, *Nat. Comput. Sci.*) to perceive the work as a toy-model demonstration.

**Advantages of the cleaner 2-tier architecture:**
1. **Tier 1 — Tractable Exact-Verification Benchmark (`M0`):** Retains the existing 40,000+ CPU-h exact PGAS vs. mixed NUTS validation, SBC, crossed predictions, and modular-cut vs. joint verification—with no need to apologize for its small size.
2. **Tier 2 — Full *JCVI-syn3A* Whole-Cell Model (`PreCell` / `SpeedyCell.jl`):** Locks/anchors well-characterized parameters at literature/BRENDA values and demonstrates compositional, targeted inference directly inside the full cell:
   - Infer a biologically meaningful subset of parameters (e.g., central carbon metabolism kinetics or a gene-expression block) via the mixed route and/or per-module inference.
   - Compare joint vs. per-module inference directly on that 15–25 parameter block *inside the full cell*—removing the need for a separate 92%-amputated model.
   - Predict **real whole-cell phenotypes** (doubling time, cell-cycle completion, essentiality/knockdowns with large resolvable effects, metabolic fluxes) and decompose predictive uncertainty across inferred parameters, held BRENDA priors, and intrinsic cell noise.
   - Opens the door to conditioning on or validating against **real published *JCVI-syn3A* data** (doubling time ~105 min, proteomics/transcriptomics from Breuer 2019 / Thornburg 2022), bringing ***Cell Systems*** genuinely into reach alongside ***Nature Computational Science***.

---

### Three Engineering Reality Checks Before Pivoting

1. **Differentiable mean-field swap vs. module-local AD on the full cell:**
   - Discrete full-cell processes (DNA replication initiation/fork progression, ribosome assembly, cell division) do not close as linear first-order ODEs.
   - Verify whether the full cell supports a smooth mean-field swap (treating discrete replication/division as fixed clocks or smooth approximations) **or** whether InferCell should perform **module-local inference** (e.g., smooth metabolic/expression module with data plug-in) and compose the resulting posterior back into the full hybrid `SpeedyCell.jl` simulator for `predict`. Either fits the compositional thesis cleanly.
2. **Per-gradient wall-clock cost on the full cell (~20 free parameters):**
   - On the reduced model, one ForwardDiff gradient costs ~38–45 s (~5.5 days/chain for NUTS). If full-cell gradients run in ≤60–90 s (or faster via module-local inference / adjoint / Enzyme), full-cell targeted NUTS is tractable on OzSTAR; if gradients take 10–30 min, **module-local inference + full-cell composition** is the clear path.
3. **Timeline / publication-strategy trade-off:**
   - **If `PreCell` / `SpeedyCell.jl` integrates cleanly in weeks:** Pivoting to `M0 + Full JCVI-syn3A` with compositional/targeted framing is a clear upgrade—very competitive for *Nature Computational Science* and viable for *Cell Systems* (especially with even one real *JCVI-syn3A* data constraint).
   - **If full-cell integration requires months of core solver engineering:** Finish the reduced-model manuscript for *Nature Computational Science* (adopting the compositional framing from Point A plus the diagnostic fixes in Sections 3–5), and reserve full-cell compositional inference + real data for a high-impact *Cell Systems* / *Nature Methods* follow-up.
