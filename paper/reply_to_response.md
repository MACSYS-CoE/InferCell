# Reply to the Review Response: Strategy & Pushback for InferCell (Paper 1)

**Date:** 2026-11-25  
**Target venue under consideration:** *Nature Computational Science* (primary), *PLOS Computational Biology* (fallback); *Cell Systems* only if a full-cell / real-data anchor is added.  
**Status:** Discussion & strategic evaluation (no edits made to the manuscript, Decision log, or Linear).

---

## Context & What We Are Replying To

This document evaluates the second agent's response to our initial manuscript review of InferCell Paper 1, with the goal of deciding where we agree, where we should push back, and how to assemble the strongest possible paper without falling into either editorial under-positioning or a multi-month engineering quagmire.

### Document Chain

1. **Manuscript under review:** [`main.tex`](main.tex) (draft at `32b8496` on the `KIM-98` branch), including [`00-abstract.tex`](sections/00-abstract.tex), [`01-intro.tex`](sections/01-intro.tex), [`02-infercell.tex`](sections/02-infercell.tex), [`03-testsystem.tex`](sections/03-testsystem.tex), [`04-m0.tex`](sections/04-m0.tex), [`05-inference.tex`](sections/05-inference.tex), [`06-prediction.tex`](sections/06-prediction.tex), [`07-ladder.tex`](sections/07-ladder.tex), [`08-discussion.tex`](sections/08-discussion.tex), [`app_methods.tex`](sections/app_methods.tex), and [`STATUS.md`](STATUS.md).
2. **Initial review:** [`manuscript_review.md`](manuscript_review.md) — identified six major blockers for flagship venues (unresolved mixed-route biases, `N = 1` reduced-model truth & uninformative 30% *ptsG* knockdown, 2-constant metabolic scope & first-order restriction, missing positioning vs. Laplace/SBI/emulators, synthetic well-specified gap / M10, and buried gems / over-defensive tone) plus two strategic reframings (Point A: targeted/per-module inference; Point B: dropping the reduced model for `M0` + full *JCVI-syn3A*).
3. **Responding agent's triage & counter-analysis (what we are replying to):**
   - [`checklist.md`](checklist.md) — 35-item repository-verified checklist assigning verdicts (`agree`, `partly`, `disagree`) and tiers (`A`, `B`, `C`), plus two new archive findings from `pilot_draws.tsv`.
   - [`response_to_review.md`](response_to_review.md) — synthesis following Tom's lean toward **one strong paper on the full *JCVI-syn3A* model (no real data)**, scoping full-cell compute and recommending a **3-tier paper** (`M0` + reduced model as middle tier + full cell) gated on KIM-116 and three timing measurements.

### Relevant Literature

Citations below match [`references.bib`](references.bib) (with one missing entry noted explicitly):

- **Whole-cell models & *JCVI-syn3A* biology / data:**
  - Karr et al. (2012), *A Whole-Cell Computational Model Predicts Phenotype from Genotype*, *Cell* 150(2):389–401 (`karr2012whole`, doi:10.1016/j.cell.2012.05.044) and WholeCellKB (`karr2012wholecellkb`).
  - Breuer et al. (2019), *Essential Metabolism for a Minimal Cell*, *eLife* 8:e36842 (`breuer2019essential`, doi:10.7554/eLife.36842).
  - Macklin et al. (2020), *Simultaneous Cross-Evaluation of Heterogeneous E. coli Datasets via Mechanistic Simulation*, *Science* 369(6502):eaav3751 (`macklin2020simultaneous`, doi:10.1126/science.aav3751).
  - Thornburg et al. (2022), *Fundamental Behaviors Emerge from Simulations of a Living Minimal Cell*, *Cell* 185(2):345–360.e28 (`thornburg2022fundamental`, doi:10.1016/j.cell.2021.12.025).
  - Thornburg et al. (2026), *Bringing the Genetically Minimal Cell to Life on a Computer in 4D*, *Cell* 189(9):2582–2597.e27 (`thornburg2026bringing`, doi:10.1016/j.cell.2026.02.009).
  - Agmon et al. (2022), *Vivarium: An Interface and Engine for Integrative Multiscale Modeling in Computational Biology*, *Bioinformatics* 38(7):1972–1979 (`agmon2022vivarium`, doi:10.1093/bioinformatics/btac049).
  - Kimpson (2026), *A Whole Minimal Cell in a Minute: A Validated, Fast Implementation of the Hybrid JCVI-syn3A Model*, in preparation (`kimpson2026whole`; `SpeedyCell.jl` / companion PreCell draft).
- **Whole-cell parameter inference, sloppiness & balancing:**
  - Gutenkunst et al. (2007), *Universally Sloppy Parameter Sensitivities in Systems Biology Models*, *PLoS Comput. Biol.* 3(10):e189 (`gutenkunst2007universally`, doi:10.1371/journal.pcbi.0030189).
  - Lubitz et al. (2010), *Parameter Balancing in Kinetic Models of Cell Metabolism*, *J. Phys. Chem. B* 114(49):16298–16303 (`lubitz2010parameter`, doi:10.1021/jp108764b) — currently uncited in text.
  - Karr et al. (2015), *Summary of the DREAM8 Parameter Estimation Challenge*, *PLoS Comput. Biol.* 11(5):e1004096 (`karr2015summary`, doi:10.1371/journal.pcbi.1004096).
  - Babtie & Stumpf (2017), *How to Deal with Parameters for Whole-Cell Modelling*, *J. R. Soc. Interface* 14(133):20170237 (`babtie2017how`, doi:10.1098/rsif.2017.0237).
  - Linden-Santangeli et al. (2025), *Increasing Certainty in Systems Biology Models Using Bayesian Multimodel Inference*, *Nature Communications* 16(1):7416 (`lindensantangeli2025increasing`, doi:10.1038/s41467-025-62415-4) — currently uncited in text.
- **Modular inference, cuts, discrepancy, stochastic kinetics & alternatives:**
  - Kennedy & O'Hagan (2001), *Bayesian Calibration of Computer Models*, *JRSS-B* 63(3):425–464 (`kennedy2001bayesian`, doi:10.1111/1467-9868.00294).
  - Andrieu et al. (2010), *Particle Markov Chain Monte Carlo Methods*, *JRSS-B* 72(3):269–342 (`andrieu2010particle`, doi:10.1111/j.1467-9868.2009.00736.x); Lindsten et al. (2014), *Particle Gibbs with Ancestor Sampling*, *JMLR* 15(63):2145–2184 (`lindsten2014particle`) — currently uncited at `app_methods.tex:67`.
  - Wood (2010), *Statistical Inference for Noisy Nonlinear Ecological Dynamic Systems*, *Nature* 466:1102–1104 (`wood2010statistical`); Price et al. (2018), *Bayesian Synthetic Likelihood*, *JCGS* 27(1):1–11 (`price2018bayesian`).
  - Plummer (2015), *Cuts in Bayesian Graphical Models*, *Stat. Comput.* 25(1):37–43 (`plummer2015cuts`); Jacob et al. (2017) (`jacob2017better`); Goudie et al. (2019) (`goudie2019joining`); Carmona & Nicholls (2020) (`carmona2020semimodular`).
  - Schnoerr et al. (2017), *Approximation and Inference Methods for Stochastic Biochemical Kinetics*, *J. Phys. A* 50(9):093001 (`schnoerr2017approximation`, doi:10.1088/1751-8121/aa54d9) — currently uncited.
  - Cranmer et al. (2020), *The Frontier of Simulation-Based Inference*, *PNAS* 117(48):30055–30062 (`cranmer2020frontier`, doi:10.1073/pnas.1912789117); Geffner et al. (2023), *Compositional Score Modeling for SBI*, *ICML* (`geffner2023compositional`) — currently uncited.
  - Kyzyurova et al. (2018) (`kyzyurova2018coupling`) and Ming & Guillas (2021) (`ming2021linked`) — linked emulators, currently uncited.
  - Hoffman & Gelman (2014) (`hoffman2014nouturn`), Talts et al. (2018) (`talts2018validating`), Vehtari et al. (2021) (`vehtari2021rank`), Säilynoja et al. (2022) (`sailynoja2022graphical`).
  - **Not yet in `references.bib`:** Gneiting & Raftery (2007), *Strictly Proper Scoring Rules, Prediction, and Estimation*, *JASA* 102(477):359–378 (needed for CRPS in [`06-prediction.tex`](sections/06-prediction.tex#L19-L21) and [`app_methods.tex`](sections/app_methods.tex#L261)).

---

## 1. Executive Summary: Agree vs. Pushback

The responding agent's triage in [`checklist.md`](checklist.md) and [`response_to_review.md`](response_to_review.md) is unusually disciplined and evidence-driven. It caught real factual/stale errors in [`manuscript_review.md`](manuscript_review.md) and discovered two genuine flaws in the current draft by checking `pilot_draws.tsv`.

| Topic | Agent's Position | Our Verdict | Actionable Takeaway |
|---|---|---|---|
| Fact-checks (KIM-90 merged, clocks on Fig. 2, Table 2 caption, spliced abstract quote) | Corrects review | **Agree** | Accept all corrections. |
| Archive findings (promoter shift = 3.1 SE w/ `r=0.01`; §5 "MAP" is $\theta_{\text{ref}}$ at $\delta=S=0$; ENO/FBA symmetric) | Flagged for Tier-A rewording | **Agree + Extend** | Immediately points to a 5-min **Jensen / GTP-clipping test** in the transcript factor (Pushback 3). |
| Cheap substitutes (10-truth MAP+Laplace; Hessian eigenspectrum for 10–15 $k_{\text{cat}}$s; 3 timing probes) | Recommended | **Agree** | Make the Hessian probe a real Tier-B item, not optional; run the true MAP (with $\delta, S$) first. |
| Three Graphs (`6a`) & PreCell/SBI positioning (`4b`) | Hold promotion pending KIM-116; nuance SBI | **Agree** | Fix Discussion–Appendix contradiction now; don't claim SBI is flatly prohibitive. |
| **3-tier paper** (`M0` + reduced middle tier + full cell) | Recommended if pivoting to full cell | **Push back** | Avoid a bloated 3-tier Frankenstein; if full cell works, make it the star and demote reduced model to a short validation bridge (Pushback 1). |
| **Full cell with zero real data** | Agreed ("no real data") | **Push back** | Paying a 1–2+ month full-cell engineering tax while refusing even one real *syn3A* anchor is a strategic mismatch (Pushback 2). |
| **Diagnosing the `-0.81 SD` promoter shift (`1d`)** | Pessimistic (citing M0 KIM-90) | **Push back** | Test the unadjusted mean-field pool ratio in the transcript likelihood immediately—takes minutes, 0 CPU-h (Pushback 3). |
| **M10 optionality & old non-goals (`5a`, `2a`, `2b`)** | M10 optional for methods venues; cites prior non-goals | **Push back** | For *Nat. Comput. Sci.*, a compact ladder with **destructive/snapshot scRNA-seq** is close to mandatory. Non-goals govern transparency, not paper quality (Pushback 4). |
| **Why NUTS (`4a`), Table 1 footnote (`6b`), Abstract (`6c`)** | Minimalist / defensive retention | **Push back** | Check Laplace at the true MAP + along M10; replace apologetic Table 1 wording with neutral framing; strengthen abstract closing without overclaiming (Pushback 5). |

---

## 2. Where We Strongly Agree with the Responding Agent

1. **Repository fact-checks.** KIM-90 is already merged (`ae731c6`) with no single cause named; the disjoint-cohort exact run (`508–678 CPU-h`) was already declined; `fig_model` already prints its 1 s / 60 s clocks; `tab:reducedcost`'s caption already itemises the 24.3 CPU-h probe; the abstract quote in [`manuscript_review.md`](manuscript_review.md#L81) spliced two sentences; and whether 99.3% GTP clipping is an artefact of the reduction is unverified because the clip rule is from `thornburg2022fundamental`.
2. **The two findings from `pilot_draws.tsv` ([`response_to_review.md`](response_to_review.md#L77-L97)).**
   - Finding 1 shows the promoter shift is *harder* to attribute to shared noise than [`05-inference.tex`](sections/05-inference.tex#L47-L58) suggests (`3.1 SE` accounting for correlations, `0.01` mean pairwise correlation, only `4/1,600` draws $\ge 0$).
   - Finding 2 catches a genuine misstatement in §5: $\theta_{\text{ref}}$ is the MAP at $\delta = S = 0$, *not* the MAP of the fitted likelihood, and ENO/FBA have skewness `0.05` and `0.06` (symmetric), so $\delta$ and $S$ are moving ENO/FBA rather than skew or optimizer early stopping.
3. **High-leverage, low-cost compute designs.**
   - Running **MAP (of the fitted likelihood with $\delta, S$) + Laplace on ~10 replicate truths** (`~130–200 CPU-h`) settles both `N = 1` and systematic bias at ~1/15th the cost of full NUTS replicates.
   - Running a **Hessian eigenspectrum probe** over 10–15 glycolytic $k_{\text{cat}}$s answers the sloppiness question (`gutenkunst2007universally`) for a few CPU-h rather than burning thousands of CPU-h on a wider NUTS run.
   - Taking **three timing measurements** before committing to the full cell prevents walking blindly into a wall-clock trap.
4. **Holding off on promoting the "Three Graphs" (`6a`).** Because the exact sampler did not mix on the metabolites-to-expression direction (`ESS = 39`) and GTP clips below $K_d$ in 99.3% of reduced-model cells, promoting it as a headline claim before KIM-116 validates the factor cut would invite a devastating referee critique. Resolve the contradiction between [`08-discussion.tex`](sections/08-discussion.tex#L29-L33) and `app_m0.tex:102` now; revisit promotion only with KIM-116 in hand.
5. **Nuancing neural SBI positioning (`4b`).** Because the companion PreCell draft trains a diffusion transformer on the full cell over 28 group multipliers from bulk data, claiming neural SBI is flatly "prohibitive" would contradict the group's own work. Contrast exact-verified mechanistic likelihoods over individual parameters + single-cell paths against amortised bulk-only inference (`cranmer2020frontier`, `geffner2023compositional`), and cite only the references needed for crisp positioning rather than forcing all 27 uncited `.bib` keys into the text.

---

## 3. Five Points of Pushback

### Pushback 1: The "3-tier paper" compromise risks a bloated, incoherent manuscript
In [`response_to_review.md`](response_to_review.md#L179-L186), the agent recommends doing the full cell *while keeping the reduced model as a full middle tier* (`M0` → `reduced model` → `full cell`).

- **The problem:** This reads like sunk-cost protection (~6,000 CPU-h spent on the reduced model). A 3-tier manuscript forces a *Nature Computational Science* reader through three distinct models where each tier breaks the previous tier's mathematical guarantee (`exact hybrid` → `exact first-order mean field` → `non-closing second-order + replication + threshold rules`), while still carrying the editorial drag of 93 declared reductions (`S1_File_reductions.tsv`) for an amputated cell right next to the full cell.
- **Alternative if committing to the full cell:** Make **Full *JCVI-syn3A*** the unmistakable centerpiece, keep `M0` as the exact-verification benchmark, and **demote the reduced model to a compact validation bridge (1 page or Supplement)** whose *sole* role is refereeing the modular cut against joint NUTS where joint NUTS is still affordable.
- **Check the full-cell route assumption, not just the 2,800-state joint ODE:** The agent's scary `~20 min/gradient` (`~115 days/chain` serial) number in [`response_to_review.md`](response_to_review.md#L135-L137) assumes a 2,800-state sparse ODE solve inside a joint 20-parameter fit *with all 455 transcripts observed*. If the full-cell route instead uses **module-local / data plug-in inference** (no 2,800-state ODE solve in the inner loop) composed into `SpeedyCell`-validated forward prediction for `predict`, the wall-clock bottleneck changes completely *and* sidesteps differentiating through 363 second-order reactions, 493 replication reactions, and `round`/`min`/`max` rules. The pre-commitment benchmarks should explicitly scope *both* the sparse-ODE targeted route and the module-local plug-in route.

### Pushback 2: Paying the full-cell engineering tax while ruling out *all* real data is a strategic mismatch
[`response_to_review.md`](response_to_review.md#L24-L26) locks in both *"one strong paper on the full syn3A model"* and *"Agreed: no real data."*

- **The problem:** On the reduced model, "no real data" is forced by construction (92% of the cell is removed). On the **full *JCVI-syn3A* model**, you pay 1–2+ months of engineering and 15,000–70,000 CPU-h to build a biologically complete cell—and then decline to confront it with a single experimental measurement.
- **What a referee will ask:** *"You implemented the entire JCVI-syn3A model—which already targets real doubling time (~105 min), proteomics, and metabolomics (`breuer2019essential`, `thornburg2022fundamental`)—yet you only ever fit synthetic exact-count trajectories. Why?"*
- **Takeaway:**
  - If you truly want a pure synthetic-data methodology paper (leaving real bulk data to PreCell), that argues **against** paying a 2-month full-cell engineering tax in Paper 1 unless the full-cell port is genuinely lightweight and yields dramatic biological phenotype tests (doubling time, essential-gene knockdowns) cleanly.
  - Conversely, if you *do* pay the full-cell engineering cost, adding **even one lightweight real-data predictive check or constraint** (e.g., published *syn3A* doubling time, proteome, or growth constraints) massively raises the ceiling—disarming the "purely synthetic" critique at *Nature Computational Science* and bringing *Cell Systems* back within reach.

### Pushback 3: Don't give up on diagnosing the `-0.81 SD` promoter shift — Findings 1 & 2 point to a 5-minute testable mechanism
In [`checklist.md`](checklist.md#L195-L205), the agent expresses pessimism about diagnosing the reduced-model bias (because KIM-90 ended with no single cause on M0) and notes that because promoters already sit `-0.86 SD` low at $\delta = S = 0$, the values of $\delta$ and $S$ are not what puts them low.

- **A concrete, untested mathematical suspect:** Combine Finding 1 (`14/17` promoters shifted negative, `0.01` pairwise correlation), Finding 2 (already `-0.86 SD` low before $\delta, S$ enter), and [`app_methods.tex`](sections/app_methods.tex#L74-L105):
  - The **bulk metabolite factor** includes the Kennedy–O'Hagan discrepancy offset $\delta$ (`kennedy2001bayesian`).
  - The **transcript birth–death factor** reads the **unadjusted mean-field pools** directly—with *no* $\delta$ and *no* Jensen correction.
  - Meanwhile, **99.3% of true cells clip GTP to zero** (crossing below $K_d$), and the transcription rate law is a **concave** Michaelis–Menten saturation function of nucleotide concentrations. By Jensen's inequality (compounded by zero-clipping), any systematic mismatch between the mean-field nucleotide multiplier $f([\text{ATP}]_{\text{MF}}, [\text{GTP}]_{\text{MF}})$ and the true cell-averaged multiplier $\mathbb{E}[f([\text{ATP}], [\text{GTP}])]$ acts as a **common deterministic scale factor across all 17 promoter strengths** inside the transcript likelihood.
  - This single mechanism simultaneously explains:
    1. why 14 of 17 promoters shift negative together,
    2. why their posterior pairwise correlation is `0.01` (it is a deterministic likelihood bias, not shared posterior covariance), and
    3. why the shift is already `-0.86 SD` at $\delta = S = 0$ and barely moves (`-0.81 SD`) when $\delta$ and $S$ are turned on in the bulk factor.
- **Immediate 5-minute check (0 CPU-h):** On the archived reference ensemble at $\theta_{\text{ref}}$, compute the time-resolved ratio $\mathbb{E}[f(\text{pools}_{\text{exact}})] / f(\text{pools}_{\text{MF}})$ over the fit window. If this ratio explains the shift, you have diagnosed the paper's biggest unexplained statistical anomaly—and can test whether a transcript-rate discrepancy term $\delta_{\text{tx}}$ (analogous to $\delta$) removes it. That turns a major vulnerability into a genuinely original methodological contribution.

### Pushback 4: M10 is not optional for *Nature Computational Science*, and prior non-goals shouldn't block a stronger paper
Throughout [`checklist.md`](checklist.md), the agent cites earlier internal decisions ("written non-goal", "reverses the 2026-11-09 decision", "parked release-vs-M10 question") as arguments against upgrades—especially treating M10 as optional (`5a`).

- **Why push back:** Pre-registration and Decision logs exist for **scientific transparency**, not to lock you into shipping a weaker manuscript. Referees evaluate the paper, not Linear milestones.
  - Keep the frozen 30% *ptsG* knockdown reported and explicitly label a stronger, resolvable knockdown as *post hoc* (here the agent's protocol in [`checklist.md`](checklist.md#L241-L247) is spot-on).
  - For ***Nature Computational Science***, a paper whose only data design is exact 1-minute live-cell transcript counts across 200 cells (an assay that does not physically exist) will face Major Revision or Rejection. A **compact M10 degradation ladder**—specifically incorporating the **destructive/snapshot scRNA-seq rung** (which the agent rightly praised as the review's best addition to Spec D6) plus dropout/bulk-only—should be treated as required for *Nat. Comput. Sci.*

### Pushback 5: Sharpen `4a` (Why NUTS vs. Laplace?), `6b` (Table 1 footnote), and `6c` (Abstract)
- **Why NUTS when Laplace matches widths to 4–8% (`4a`):** Answering *"we only know Laplace works because we ran NUTS"* is honest, but leaves a loose end. First, compute the **true MAP + Laplace of the fitted likelihood (with $\delta, S$)** on the current dataset to test whether Laplace also matches NUTS on ENO/FBA *means*. If it does, report Laplace as a validated fast route on generous data—and use the **free Laplace-vs-NUTS comparison along the M10 degradation ladder** (the Hessian is already computed for the warm metric) to show *where* sparser/noisier data break Gaussianity and make NUTS genuinely necessary.
- **Table 1 footnote (`6b`):** Agree you must read the four supplements (`karr2012whole`, `macklin2020simultaneous`, `thornburg2022fundamental`, `thornburg2026bringing`) before removing *"We did not read the supplements."* **Disagree** with keeping *"These columns are InferCell's own design aims, so the comparison favours it"* ([`01-intro.tex`](sections/01-intro.tex#L93)). Keep the honesty without sounding apologetic:
  > *"Columns reflect capabilities required for compositional Bayesian inference on dynamics; earlier frameworks were built primarily as forward simulators calibrated to real experimental data rather than synthetic inference benchmarks."*
- **Abstract closing (`6c`):** Agree the review spliced two sentences and proposed claims the paper does not prove. **Disagree** that the current closing ([`00-abstract.tex`](sections/00-abstract.tex#L21-L23)) is optimal—ending on *"where it is wrong"* and *"what it does not cover"* undersells the work. Move the synthetic-data scope sentence cleanly into the body, and close on what the paper *actually achieves*: turning whole-cell models from point-estimate simulators into auditable probabilistic frameworks whose parameter uncertainty, approximation error, and declared reductions propagate explicitly to predictions.

---

## 4. Recommended 4-Gate Decision Protocol

Do **not** commit to Full-Cell Paper 1 vs. Reduced-Model Paper 1 until **four fast, cheap gates** complete (days, not months; `<250 CPU-h` total):

1. **Gate 1 — Jensen / GTP-clipping ratio check (minutes, 0 CPU-h):** Compute $\mathbb{E}[f(\text{pools}_{\text{exact}})] / f(\text{pools}_{\text{MF}})$ on the archived reduced-model reference ensemble to test Pushback 3.
2. **Gate 2 — True MAP + 10-truth Laplace replicate check (`~150–200 CPU-h`):** Run the true fitted-likelihood MAP (with $\delta, S$) on the current dataset (settling Finding 2 on ENO/FBA and `4a`) plus ~10 new synthetic truths via MAP + Laplace (settling `1b`, `1d`, and `2a`).
3. **Gate 3 — KIM-116 results:** Wait for the modular factor cut & data plug-in comparison against the joint posterior. This governs both **Point A** (targeted/compositional framing) and whether full-cell modular inference has a validated approximation error.
4. **Gate 4 — Full-cell feasibility benchmarks (a few hours):** Run the agent's three measurements (one `SpeedyCell` cycle, a ~2,800-state sparse dual ODE solve, and the 455-gene transcript factor) *and* verify whether a **module-local / data plug-in route** avoids differentiating through the full cell's 363 second-order and 493 replication reactions.

### The Fork After Gates 1–4

- **Branch A — One Flagship Full-Cell Paper** (*choose only if* KIM-116 validates the modular/plug-in route AND Gate 4 shows targeted full-cell inference + scored predictions fit within ~1–2 months and a realistic OzSTAR budget):
  - Make **Full *JCVI-syn3A*** the star (real whole-cell phenotypes such as doubling time + a strong resolvable knockdown), keep `M0` as the exact-verification benchmark, and compress the reduced model into a tight validation bridge for joint-vs-modular calibration.
  - Seriously revisit adding **at least one lightweight real-data anchor** (`breuer2019essential` / `thornburg2022fundamental`) so the full-cell engineering investment pays off at *Nature Computational Science* (and potentially *Cell Systems*).
- **Branch B — Upgraded Reduced-Model Paper 1 + `v0.3.0` Release** (*choose if* Gate 4 shows full-cell non-closure or `15k–70k CPU-h` scoring will turn into a 3–6 month engineering slog):
  - Target *Nature Computational Science* (*PLOS Computational Biology* fallback) with:
    1. **Point A framing** (targeted/compositional inference, joint posterior as yardstick), conditional on KIM-116.
    2. **Resolved/diagnosed promoter shift** + 10-truth Laplace replicate check.
    3. **One strong, resolvable knockdown** (post hoc, clearly labeled alongside the frozen 30% result).
    4. **Compact M10 degradation ladder** including **destructive/snapshot scRNA-seq** and the free **Laplace-vs-NUTS** comparison along the ladder.
    5. **Hessian eigenspectrum probe** on 10–15 glycolytic $k_{\text{cat}}$s (`3a`).
    6. **All Tier-A polish** (KIM-90 todos, positioning paragraphs + missing citations including Gneiting & Raftery 2007 and `lindsten2014particle`, Table 1 supplement checks & neutral footnote, tightened abstract, redrawn Figs. 1–2, and main-workflow §2 code listing).
  - Ship `v0.3.0`, then execute **Full *JCVI-syn3A* + Real Experimental Data** cleanly as Paper 2.
