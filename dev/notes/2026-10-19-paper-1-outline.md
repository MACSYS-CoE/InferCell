# InferCell paper 1 (v0.3.0): outline

**Date:** 2026-10-19
**Status:** temporary working outline. The reasoning and the decisions behind it
are in [`2026-10-19-inference-framing-discussion.md`](2026-10-19-inference-framing-discussion.md)
(A4 to A6). Nothing here is a claim yet: every result marked *not built* or
*running* must exist before its sentence is written.

**Superseded in part (2026-10-23).** The draft now lives in `paper/` (branch
`paper-1-draft`); see `paper/STATUS.md`. Since this outline was written:
- F10 is done (`dev/scripts/f10_16a10_result.md`), not running.
- §4.3 predates the 2026-10-23 decision (`spec/phases/approximate-method.md` §12):
  the reference is the most detailed composition, not the truth, and the Jensen gap
  is a declared reduction corrected by δ, with KG reported, not gated.
- The measured R (64 on M0, 70 on the reduced model, at σ_b 10%) replaces the
  estimates in the discussion note's §7.
- The paper is journal-agnostic (Cell Systems or PLOS Computational Biology), and its
  lead is framework plus problem: *InferCell: inference-first whole-cell modelling*.

## Title, tagline, message

**Working titles:**
- *InferCell: composable whole-cell models with calibrated uncertainty*
- *Calibrated Bayesian inference across the stochastic–deterministic boundary of a whole-cell model*
- *What crosses the boundary? Scalable, calibrated inference in a hybrid whole-cell model*

Avoid "propagation" (readers hear forward Monte Carlo), "fast" (an M0 posterior
costs 500–700 CPU-h; speed is syn3a-ensemble's claim), and "whole-cell models"
in the plural (one reduced model).

**Tagline candidates** (decision O8):
- *Composable whole-cell models with calibrated uncertainty.* Preferred: V7 earns "calibrated".
- *Whole-cell models that say how sure they are.*
- *Inference-first whole-cell modelling: typed couplings, traceable priors, calibrated posteriors.*

**The message in one paragraph.** Whole-cell models fix their parameters, so they
cannot say how sure they are. Inference is hard because the models mix
formalisms, the stochastic paths are latent, and exact inference costs millions
of simulations. InferCell composes heterogeneous modules through typed coupling
edges, carries a prior and its provenance on every parameter, and matches the
inference to each module's description and each data stream. On a model reduced
from the published JCVI-syn3A whole-cell model, an exact sampler gives a
calibrated reference. An approximate method, which exploits the weak reverse
coupling and averages out the stochastic paths for bulk data, matches the
reference at a fraction of the cost, and then scales to the full reduced model,
where metabolite data inform a gene-expression parameter they never measure.

**The scalars for the abstract** (spec §6, updated): cost per posterior, exact
against approximate; the reverse-coupling gain range (0.04–0.09); the
replication counts (V7's 50; 5.3's injections); 95% coverage per target; ptsG
shrinkage under metabolites only, with ENO's null under transcripts only.

## 1. Introduction: why inference, and why nobody does it

**Content:**
- Whole-cell models simulate a cell from its parts (Karr 2012, Macklin 2020,
  Thornburg 2022, the 4D model of 2026). They are forward simulators with fixed
  parameters. A prediction is one sample with no error bar.
- Every parameter is uncertain; that uncertainty propagates; a framework with no
  place for a distribution cannot say whether a prediction is robust or fragile
  (`the-case-for-inference.md`).
- **Why it is not done, four reasons:**
  1. No place for uncertainty: parameters are stored as single numbers.
  2. Mixed formalisms: no single likelihood or sampler covers ODE and jump blocks,
     and the stochastic paths are latent.
  3. Cost: the number of simulations, not the cost of one.
  4. No parameter is shared across the seam (D13), so the standard modular
     methods (cut, melding, EP over shared quantities) have nothing to pass.
- Prior art: Babtie and Stumpf (2017) argue for this and do not do it.
  Thornburg et al.'s parameter balancing is Bayesian but metabolism-only, from
  literature values, with no inference from dynamics through the coupled model.
- **Contributions** (four bullets matching sections 2–6).

**Status:** written in substance in `the-case-for-inference.md`.

## 2. InferCell: design

**Content:** three principles.
1. **Parameters are distributions with provenance.** Every parameter carries a
   prior, a width, an informedness label and a source file (T1).
2. **Modules couple through typed edges.** Seven edge kinds on three clocks; the
   resolver checks completeness before anything runs; every reduction is
   declared, with its measured cost (T2).
3. **Inference matches the description.** A module declares its formalism; a data
   stream declares its observation model; the dispatcher picks the inference.

**Status:** 1 and 2 done (phases 0–15). 3 not built (the dispatcher).
**Open:** O3, a generic registry or "demonstrated on one organism".

## 3. The reduced syn3A model

**Content:**
- Twenty-two reactions in four ODE modules; seventeen genes in three jump
  modules; exchange every 1 s, rate constants every 60 s. About 92% of the cell
  removed, keeping the best-measured region.
- Why reduced: small enough for an exact reference, and still a genuine hybrid
  boundary from a published model.
- No parameter crosses the seam: the blocks couple through state (D13).
- M0, the two-gene rung, for the exact comparison.
- What the model does not claim: no syn3A phenotype comparison (the code refuses
  it).

**Figure 1:** the coupling figure (`dev/notes/figures/corea-coupling/`), with the
edge kinds and clocks.
**Status:** done.

## 4. How it works

### 4.1 Composition and the hybrid driver
Typed edges, the resolver, the 1 s handshake, the 60 s rebuild, exact path
replay. **Status:** done.

### 4.2 How information crosses the seam
- Forward, a distribution over the interface trajectory; backward, a likelihood
  over it. A cut drops the backward factor.
- The three graphs: coupling, likelihood dependence, information.
- Statistical cuts against dynamical reductions.
- **The gain table (F1):** forward order one; reverse 0.044–0.051 (transcription)
  and 0.091 (translation), because the published rate laws are saturated.

**Figure 2:** the three graphs for Core A′, side by side.
**Figure 3 or Table 1:** F1.
**Status:** F1 done. The saturation-versus-rebuild test is not run (discussion §6).

### 4.3 Two descriptions of gene expression, and the inference each one gets

| Description | Inference | Data stream it serves |
|---|---|---|
| Exact stochastic (jump) | Three-block Gibbs: θ_ODE by gradient-free MH, θ_CME by exact 1-D conditionals, paths by conditional SMC with exact transcript bridges | Reference for everything |
| Mean field (linear ODEs, exact means) | Mixed likelihood: exact transcript factor × Gaussian bulk factor; NUTS | Bulk metabolites |
| Exact discrete chain per gene | Closed-form transition probabilities | Single-cell transcripts, in both routes |

The choice is per data stream. The reductions the mean field needs are declared:
the Jensen gap, continuous drains, mean-field feedback. **R** decides the
regime (discussion §7).

**Status:** the exact route is done. The mean-field modules, the ODE path running
catalytic, drain and rebuild edges, the transcript factor and the dispatcher are
**not built**.

## 5. Validation

| | What it shows | Evidence | Status |
|---|---|---|---|
| **5.1 The reference** | The exact sampler is correct and calibrated | Block tests; V6 (4 chains × 1,650 sweeps, R̂ ≤ 1.007, ESS ≥ 1,081, 6.8k CPU-h); V7 (R = 50, all six within bands at 1 − 0.05/6; 23 of 50 below ESS 100); F10 (production against reference, ~13× cheaper) | V6, V7 done; F10 running |
| **5.2 Approximate vs exact on M0** | Same posterior at a fraction of the cost. ptsG still shrinks under metabolites only | V7's 50 datasets refitted; marginals, the ptsG–ENO contour, rank-ECDF, cost | Not built. **Gate:** M0 must carry a metabolites → ptsG signal (linearised check first) |
| **5.3 The full reduced model, D11's six** | Calibrated over ≥ 200 injections. Widths reach the information bound. F6 (transcripts only, metabolites only, joint) at full scale | Injections from the prior, data from the full two-way model; bound from `ensemble_jacobian` | Not built. Needs the cap reopened (O1) |
| **5.4 All ~110 parameters (optional)** | Which parameters the data constrain | Shrinkage map, checked by calibration only | Not built. Needs gradients (route (ii)) |

**Figure 4:** approximate against exact on M0 (overlaid marginals, the ptsG–ENO
contour, a cost bar).
**Figure 5:** rank-ECDF with simultaneous bands: exact on M0, approximate on M0,
approximate on Core A′.
**Figure 6:** F6 at Core A′ scale, with the information bound drawn on.
**Table 2:** cost per posterior: exact reference, exact production, approximate;
M0 and Core A′ (the Core A′ exact entry projected).

## 6. Demonstration: the uncertainty budget of a prediction

**Content:** infer what the data constrain, and propagate the rest.
- Posterior draws for the inferred parameters, prior draws for every other one,
  run through the coupled simulator: the predictive distribution of a quantity
  that was not observed.
- Split its variance three ways: the inferred parameters, the parameters left at
  their priors, and intrinsic noise. Variance-based sensitivity over the rest.
- The sentence to aim for: *"x% of the uncertainty in this prediction comes from
  parameters that no measurement supports."* (68 asserted priors in the full-cycle
  report, some of them dead slots.)
- Report the split with the asserted widths halved and doubled. Exclude dead slots
  (`enz_R_*`, overwritten every handshake). Draw reverse constants
  Haldane-consistently.
- Quantity (O4): fractional growth, lactate export, or how often GTP empties (the
  port finds 0.71% of full-model cells). Must stay inside the reduced model.

**Figure 7:** the budget.
**Cost:** ~50 CPU-h on M0; ~90–450 CPU-h for a prior-only budget on Core A′.
**Status:** not built.

## 7. Discussion

### 7.1 What the results mean

- What scales and what does not: exact inference stops at M0; the approximate
  method's limits are R and the Jensen gap.
- The three graphs as a principle: the observation model sets the inferential
  interface. Single-cell streams stay local; cross-module information travels on
  bulk observables.
- Limits: synthetic data only; one organism; the asserted priors.

### 7.2 Future work

One short paragraph in the paper, in this order. Do not list anything cheap
enough that a reviewer would ask "why not now?". The generic registry (O3) is
one: do it before the paper, or narrow the "framework" claim.

1. **The method's own limits,** one sentence each:
   - condition the metabolite factor on the observed transcripts, the coherent
     refinement, which matters once R ≳ 1;
   - pools the Gaussian approximation cannot handle (Pi, 3PG, 2PG at the
     one-particle floor or the clip), by sampling paths only for those pools;
   - regimes where the reverse feedback is strong: the one-way reduction rests
     on elasticities of 0.04–0.09, and a perturbation that drains the NTP pools
     toward K_d (0.1 mM) removes that basis.
2. **The full well-stirred syn3A model (Step 2) through the syn3a-ensemble port,
   then real data.** The obvious next paper, using both projects.
3. **Comparing model structure, and comparing formalisms,** which the formalism
   swap makes possible (`model-uncertainty-and-selection.md`).
4. **Goal-directed inference and experimental design,** framed by the closing
   question: *which parts of a coupled simulator must be retained, and which
   latent variables must be inferred, to preserve the information a prediction
   needs?*

The spatial model (Step 3) and the theory get one line each at most. Spelling
them out promises more than we want to be held to.

## After v0.3.0: the roadmap (planning, not paper text)

### What paper 1 leaves open

| Item | Where it is decided or done |
|---|---|
| Condition the bulk factor on the observed transcripts (R ≳ 1) | v0.4 |
| Non-Gaussian pools near the floor or a clip: path sampling for those pools only | v0.4 |
| Strong-feedback regimes, where the one-way reduction fails | Step 1b; perturbation studies |
| Per-module calibration (phase 17; the modular note's open problem 5) | Phase 17 |
| The reverse direction for FBA (ENO is a measured null) | K2 at 1,000 cells (16b.5) |
| An exact posterior on Core A′, if not funded for the paper (O5); exact coverage at that scale stays out of reach | 16a.11 |

### The longer programme

**Up the ladder of models** (spec §2)
- **Step 1b:** Core A′ plus the RNDR/dNTP branch, which makes the NTP pools live.
  The rung where the reverse coupling could start to matter.
- **Step 2:** the full well-stirred syn3A model, with the syn3a-ensemble port as
  the simulator. Open: the full gene-expression network has assembly, competition
  and degradation steps. If it is not first order, the means are no longer exact,
  and the mean-field route needs moment closure or a learned surrogate.
- **Step 3:** the spatial 4D model. Spatial resolution enters as a testable
  hypothesis against the well-mixed one (`spatial-dynamics-and-rdme.md`).
- **A second organism,** which needs the generic registry (spec §7 non-goal).

**Methods** (the modular note's §6 open problems become reachable)
- Message passing across many modules. With the paths integrated out, parameters
  appear in several likelihood factors, so EP and semi-modular inference apply.
- Per-module amortisation for stochastic modules with no closed-form means
  (bursty promoters, feedback inside a block): a neural likelihood or a learned
  message operator.
- Unbiased gradients through the stochastic block (StochasticAD), for
  structured variational inference.
- Theory for the three graphs: when modular inference recovers the joint
  posterior, and whether R and the Jensen gap can be computed in advance.

**Real data** (Step 2's business)
- Observation models for real assays: bulk metabolomics with unknown cell numbers,
  and single-cell RNA counts with capture efficiency.
- Misspecification: posterior predictive checks, and semi-modular inference
  weights per module. Separate the uncertainty within a model from the
  uncertainty about whether the model is right.

**Model structure** (`model-uncertainty-and-selection.md`)
- Structural comparison, such as constitutive against bursty promoters, by
  predictive score (ELPD, leave-one-stream-out), not Bayes factors.
- Comparison across formalisms: does the ODE or the stochastic description
  explain the data better? The normalisation hooks are already tracked.

**Goal-directed inference and experimental design**
- Generalise paper 1's uncertainty budget: spend inference on the parameters
  that drive a chosen phenotype.
- Experimental design: which measurement carries the most information across the
  boundary.
- Predictions under intervention. Recompute the information graph for each one,
  because a weak channel can become strong.

**Performance**
- The ODE step: FBDF + KLU, one process per core, GPU batching.
- The exact sampler: joint parameter–path moves, correlated pseudo-marginal
  methods, particle-parallel sweeps.
- Surrogates: the spec's two proposals, a learned promoter prior and a
  differentiable metabolic surrogate.

### Likely order

1. **v0.4:** the coherent refinement and the non-Gaussian pools, plus the generic
   registry if it is not already done.
2. **Paper 2:** Step 2 on the full cell through the port. Forward UQ of a
   phenotype first (is the ~16% PEP collapse robust to the model's own parameter
   uncertainty?), then inference.
3. **Then:** real data.

## Supplementary material

- The exact sampler in full: blocks, bridges, conditional SMC, V1–V9.
- The kill-criteria scoreboard (T3) and the reductions with measured costs (T2).
- The provenance table (T1), machine-generated.
- The ODE solver benchmark (FBDF + KLU at 3.4×, measured, not adopted).
- The cost arithmetic of `2026-10-13-inference-paths.md`.
- The `iterative_infer` correction, as a cautionary note (discussion §8).

## Reviewer questions to answer in the text

| Question | Answer |
|---|---|
| "It only works on a two-gene toy." | 5.3: the approximate method on the full reduced model, with injections |
| "Why not mean field everywhere?" | Per-stream choice: transcripts are low-copy single-cell counts and keep the exact likelihood |
| "Calibration does not show you kept the information." | 5.3's efficiency check against the information bound |
| "Can I use it for my model?" | O3: the generic registry, or say "demonstrated on one organism" |
| "Is the exact reference itself right?" | 5.1: V6, V7, F10 |

## Literature

From the repo's notes, already used there. Check every entry against the source
before citing. Entries marked † are not yet cited anywhere in the repo.

**Whole-cell models**
- Karr et al. (2012), *Cell*. The first whole-cell model.
- Macklin et al. (2020), *Science*. The E. coli whole-cell model.
- Thornburg et al. (2022), *Cell* 185:345. The syn3A model this one is reduced from.
- Luthey-Schulten et al. (2026), *Cell*. The 4D minimal cell.
- Agmon et al. (2022), *Bioinformatics*. Vivarium.
- Hutchison et al. (2016), *Science*; Breuer et al. (2019), *eLife*. syn3A itself.
- The syn3a-ensemble preprint (forthcoming). The fast port.

**Parameters in whole-cell and systems-biology models**
- Babtie and Stumpf (2017), *J. R. Soc. Interface*. The case for inference in whole-cell models.
- †Lubitz et al. (2010), *J. Phys. Chem. B*. Parameter balancing, the method behind the metabolic priors.
- †Gutenkunst et al. (2007), *PLoS Comput. Biol.* Sloppiness.
- †*Increasing certainty in systems biology models using Bayesian multimodel inference*, *Nat. Commun.* (2025). Authors to check.

**Modular Bayesian inference**
- Liu, Bayarri and Berger (2009). Modularization.
- Plummer (2015), *Stat. Comput.* Cuts.
- Jacob, Murray, Holmes and Robert (2017). *Better together?*
- Carmona and Nicholls (2020). Semi-modular inference.
- Yu, Nott and Smith (2023), *Statist. Sci.* Variational cut.
- Chakraborty, Nott, Drovandi, Frazier and Sisson (2023), *Stat. Comput.* Likelihood-free cut.
- Liu and Goudie (2025), *JRSS-B* 87:1171. Cutting feedback on a general DAG.
- Frazier and Nott (2025), *Bayesian Analysis*. Generalized-Bayes cut.
- Goudie, Presanis, Lunn and De Angelis (2019), *JASA*. Markov melding.
- Manderson and Goudie. Chained melding.

**Message passing**
- Minka (2001). Expectation propagation.
- Barthelmé and Chopin (2014), *JASA*. EP-ABC.
- Vehtari et al. (2020), *JMLR*. EP as a way of life.

**Exact inference for stochastic kinetics**
- Andrieu, Doucet and Holenstein (2010), *JRSS-B*. Particle MCMC.
- Andrieu and Roberts (2009), *Ann. Statist.* Pseudo-marginal.
- †Doucet, Pitt, Deligiannidis and Kohn (2015), *Biometrika*. The log-likelihood variance a pseudo-marginal chain tolerates (R's threshold).
- Golightly and Wilkinson. Particle MCMC for stochastic kinetic models.
- Hobolth and Stone (2009), *Ann. Appl. Stat.* Endpoint-conditioned CTMC bridges.

**Approximations and simulation-based inference**
- Schnoerr, Sanguinetti and Grima (2017), *J. Phys. A*. LNA, moment closure, CLE.
- Wood (2010), *Nature*. Synthetic likelihood.
- †Price, Drovandi, Lee and Nott (2018), *JCGS*. Bayesian synthetic likelihood.
- Cranmer, Brehmer and Louppe (2020), *PNAS*. Simulation-based inference.
- Geffner, Papamakarios and Mnih (2023), *ICML*. Compositional score modelling.
- Kennedy and O'Hagan (2001). Calibration with discrepancy.
- Kyzyurova, Berger and Wolpert (2018); Ming and Guillas (2021). Linked emulators.

**Validation**
- Talts et al. (2018). Simulation-based calibration.
- Säilynoja, Bürkner and Vehtari (2022), *Stat. Comput.* Rank-ECDF simultaneous bands.
- †Saltelli et al. (2010), *Comput. Phys. Commun.*; †Sobol' (2001). Variance-based sensitivity, for section 6.

**Related inference on gene expression**
- †Lin and Buchler (2019), *J. Chem. Phys.* Exact hybrid Monte Carlo for gene-expression models from single-cell snapshots.

## Before writing

1. F10 finishes.
2. The feasibility checks and the M0 ptsG linearised check pass (discussion §14).
3. Route (ii), the transcript factor and the dispatcher are built.
4. 5.2 passes; the cap is reopened for 5.3 (O1).
5. A §12 entry amends spec §1 (O2).
6. T3 is filled.
