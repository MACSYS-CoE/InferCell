# Spec: Core A′ — inference across a whole-cell ODE/stochastic boundary

**Status:** draft
**Created:** 2026-09-03  ·  **Last amended:** —

This is the authoritative document for the Core A′ work. It supersedes
`openspec/`, which moves to `dev/archive/openspec/` and is retained only so its
frozen interface requirements and four drafted module designs stay citable. It
also supersedes the wave ordering in `dev/plans/reduced-syn3a-wave-plan.md`; see
§4 D0 for why that ordering changed.

It does **not** supersede `dev/notes/reduced-syn3a-scoping.md`. That note decides
*what* Core A′ is and has survived two rounds of external review. Where this spec
and the scoping note disagree on the model's content, the note wins and this
spec is stale. Where they disagree on sequencing or on inference design, this
spec wins, because the note does not address either.

---

## 1. The claim

**Parameter uncertainty propagates correctly across a genuine ODE/stochastic
boundary in a model derived from a published whole-cell simulation.** Concretely:
in Core A′ — 21 metabolic reactions and 17 genes carved out of Thornburg et al.
2022's well-stirred JCVI-syn3A model, with three deterministic modules and four
stochastic ones exchanging state every simulated second and rate constants every
sixty — posteriors over parameters that cross the boundary are *calibrated*, not
merely centred: coverage on synthetic data is nominal at stated levels, rank
statistics are uniform within test, and what the data does and does not constrain
is characterised rather than assumed.

The enabling result, and not the claim itself, is that a composable
typed-interface architecture can host heterogeneous formalisms from a published
whole-cell model, execute the coupling between them, and support inference over
the whole with every reduction and every provenance choice enumerable from the
composed model.

**What this claim deliberately excludes.** Core A′ removes roughly 92% of the
cell. It will not reproduce syn3A's 105-minute doubling time, its metabolite
concentrations, or its proteome, and no result here is to be compared to them.
That is a consequence of the reduction, not a defect in it. Growth is reported as
fractional growth or time-to-threshold, never as a doubling-time prediction. The
bar this work is held to is in the table below, and the left column is not it.

| Not the bar | The bar |
|---|---|
| Reproduces the 105 min doubling time | Parameters recover from synthetic data with nominal coverage |
| Matches syn3A metabolite concentrations | Posteriors are calibrated, checked per module and jointly |
| Predicts the observed proteome | Identifiability is characterised: we know what the data constrains |
| Forward trajectories look biological | Uncertainty crosses the boundary correctly, and the gain in each direction is measured |

---

## 2. Background

### The ladder, and which rung this is

Step 1 of three toward running InferCell on JCVI-syn3A. Step 2 is the full
well-stirred model of Thornburg et al. 2022; Step 3 is the spatially resolved
4DWCM. This spec is Step 1a: the reduction defined in
`dev/notes/reduced-syn3a-scoping.md` as **Core A′**, built and run to synthetic
parameter recovery. Step 1b, which adds the ribonucleotide-reductase branch and
makes the nucleotide-triphosphate pools genuinely live, is out of scope (§7).

### What Core A′ is, in one table

Read `dev/notes/reduced-syn3a-scoping.md` for the derivation and
`dev/notes/figures/reduced-syn3a-coupling/fig1r_state_graph_reduced.png` for the
coupling graph. The summary:

| | |
|---|---|
| ODE block | 21 reactions in 3 modules: glycolysis through lactate (10), the phosphotransferase cascade plus lactate export (6), nucleotide recycling (5) |
| Stochastic block | 17 genes × 3 reactions plus one translocation = 52 reactions in 4 modules: transcription, translation, transcript decay, lumped tRNA charging |
| State | 32 dynamic species and 5 chemostats, named and ordered once in the registry |
| Copy-number span | mRNA 0–2, protein 266–1355, metabolites ~200–74,000 particles |
| Coupling | 4 channels on 2 timescales: 3 at 1 s, 1 at 60 s, plus volume-mediated dilution |
| Growth | Surface-area accounting; reaches ~1.07× initial volume, so the published 2× cap is never approached |

The reduction is honest about what it preserves: a genuine ODE/stochastic split,
bidirectional coupling, low-copy discrete species, and autocatalytic closure. It
is equally explicit that Core A′ is the *best-measured* region of the network,
that enzyme competition is partly removed, and that several treatments — the
lumped charging step, the chemostatted CTP/UTP/amino-acid pools, exogenous
membrane growth — are ours and not the published model's. Every one of those
caveats must reach any reported result (§6, T2).

### What wave 0 delivered, and is not rebuilt here

The interface contract is merged, archived, and frozen. Its requirements live at
`dev/archive/openspec/specs/corea-interface/`. In code:

| File | What it holds |
|---|---|
| `src/organisms/coreA/registry.jl` | `COREA_SPECIES` — 32 dynamic + 5 chemostat entries, each with group, treatment, copy-number regime, initial value, geometric standard deviation, source file and informedness. Plus `species_index`, `held_value`, `species_in_group` and eleven other helpers |
| `src/edges.jl` | The seven `CouplingEdge` kinds, each validating in an inner constructor: mass, currency, deferred counter, catalytic, rate constant, volume, clamped |
| `src/resolver.jl` | `resolve_coupling` — throws on eight classes of boundary inconsistency, reports unowned states, dead ends, chemostat exemptions, gradient obstructions and deviations |
| `src/loader.jl` | `read_source_table`, `load_parameter`, `ambiguity_report` — SBtab import that refuses to choose silently between two source files and records what it rejected |
| `src/labels.jl` | `reduction_declarations`, `reduction_report` — enumerates every departure from the published model |
| `src/parameters.jl` | `InferParameter` with a provenance field, `ParameterSource`, the informedness vocabulary |

829 tests pass. **No Core A′ sub-model exists.** The only thing exercising the
contract is the `CoreAStub` double in `test/corea_test_models.jl`, which
implements neither `dynamics` nor `reactions`.

### The four framework gaps, verified in source

These reorder the work, so they are stated with references rather than
described. Each was read in the source this session.

**G1 — mixed-formalism composition is refused in one line.**
`build_problem` routes on `_determine_formalism`, which returns `:mixed` for a
heterogeneous set, and the mixed branch is `error("Mixed formalism composition
not yet supported")` at `src/orchestrator.jl:27`. Core A′ is three ODE modules
and four jump modules. There is no hybrid path, no stub, and the refusal is
pinned by a test at `test/test_stochastic_ge.jl:58`. The same branch swallows a
uniform `:sde` composition, which is why the copy-number table's "chemical
Langevin sits naturally here" for the 266–1355-copy enzymes is unimplemented too.

**G2 — there is no execution layer for coupling.** Six of the seven edge kinds
are declare-and-validate only. Searching `src/` and `test/` for periodic
callbacks, discrete callbacks, callback sets, `tstops` or operator splitting
returns nothing outside the `affect!` closures inside individual jump
definitions. `solve` is called in five places, always as a single one-shot call.
So the 1 s handshake, the 60 s rebuild and the volume rescaling are postulated by
the contract and written nowhere. The resolver's own docstring says as much: a
deferred counter "is exempt — the hook debits it, the RHS never reads it."
`DiffEqCallbacks` is already present in the Julia depot, so adding it needs no
network fetch from a compute node.

Note that `src/boundary.jl` is **not** this layer despite its name. It is
inference-level message passing between two independently solved problems. It
does not couple trajectories.

**G3 — jump composition mis-indexes state and parameters.**
`_build_jump_problem` computes the per-module `SubModelContext` index maps,
passes them to validation, and then concatenates the raw jumps and discards the
contexts. But `reactions(m)` receives only the model and its closures index
locally, so two jump modules both write position one. The one test that composes
two jump models passes because it *expects* a throw, and the throw comes from the
duplicate-state-name check rather than from any index machinery — give the
modules disjoint names and state is silently corrupted. Also, `reactions`
receives no input vector at all, so a jump module can neither read nor write a
peer's state. Translation cannot read the transcripts transcription owns, and
decay cannot decrement them.

One of the drafted designs asserts the opposite: that `build_problem` "already
dispatches on `formalism = :jump`; the existing stochastic gene-expression models
exercise that path." True for one module, false for a composition. That is the
most load-bearing wrong assumption in the four drafts, and it is corrected here.

**G4 — the observation model is one rigid form.** The only likelihood in the
codebase is `_infercell_model` at `src/inference.jl:1-24`. It observes *every*
state through `MvNormal(sol[:, i], sigma)` with one scalar noise term read from a
fixed slot. Core A′ spans five orders of magnitude of copy number, so a single
additive term makes the small pools — 13DPG at ~198 particles, 2PG at ~549, GMP
at ~236 — contribute nothing to the likelihood while phosphate at 17.8 mM
dominates it. Those small pools are exactly the informative ones (§3, check 8).
Separately, `observe` on a trajectory ensemble *averages across replicates* and
`compute_summary_stats(...; times=…)` returns per-time means only, which is the
path the simulation-based sampler actually uses. Both destroy the distributional
information the 0-to-2-copy regime carries, which the scoping note identifies as
the whole point of that regime.

There is also no packaged calibration or coverage tooling. What exists is
script-level confidence-interval checking in two Slurm drivers, with an in-file
comment conceding it is not research-grade.

### Prior work this rests on

`dev/notes/well-stirred-minimal-cell.md` maps the source model, including the
`Mode`-column trap that governs every imported parameter.
`dev/notes/modular-bayesian-inference-heterogeneous-modules.md` supplies the
tractability rungs, the prior-double-counting hazard, and the instruction that a
reduced exact reference must exist before any approximation is trusted.
`dev/notes/model-uncertainty-and-selection.md` supplies the reason model
comparison is predictive rather than evidence-based here.

---

## 3. Scientific validity

### Model

- **ODE block.** 21 reactions. Every intracellular reaction takes the published
  simulator's modular rate law, which is `Rxns.Enzymatic` in the source and *not*
  the SBtab `KineticLaw` column, which is inert for these reactions:

  ```
  v = E · ( kcatF·Π(Sᵢ/KmSᵢ) − kcatR·Π(Pⱼ/KmPⱼ) )
        / ( Π(1+Sᵢ/KmSᵢ) + Π(1+Pⱼ/KmPⱼ) − 1 )
  ```

  One rate-law term per unit of stoichiometry, so a coefficient of two appears
  squared. The phosphotransferase cascade is mass-action forward/reverse pairs
  instead, and lactate export is passive diffusion, `P·(lac_c − lac_e)·3/r_cell`.
- **Stochastic block.** Constant-rate jumps. Transcription is constitutive: one
  reaction per gene firing on gene presence, with no promoter switching. Decay
  is first order in the transcript at one global catalytic constant over
  transcript length. Translation is first order in the transcript. The lumped
  charging step is `M_trna_c + ATP → M_trna_chg_c + AMP + PPi`.
- **Coupling.** Protein counts overwrite enzyme concentrations every 1 s;
  accrued expression costs are debited against the metabolite pools every 1 s;
  the pools set the stochastic block's rate constants at a 60 s rebuild;
  membrane-protein counts set surface area, which sets volume, which rescales
  every count-to-concentration conversion.
- **Discretisation.** Stiff implicit integration for the ODE block with pinned
  tolerances (§3, tolerance principle). Exact stochastic simulation for the jump
  block. The two advance by a 1 s split-operator exchange, matching the published
  model's own handshake.

### Assumptions and approximations, each with its regime

| Assumption | Holds when | What breaks outside it |
|---|---|---|
| Metabolites are continuous and deterministic | Pools are large. True for ATP at ~74,000 particles | Marginal at ~198 particles (13DPG), where relative Poisson noise is ~7% — and those are the pools with the largest concentration control coefficients. Bounded by check 1b, not assumed away |
| Enzyme counts enter the ODE deterministically | 266–1355 copies gives 3–6% relative Poisson noise | Core A′ therefore supports a claim about *parameter* uncertainty crossing the boundary, and honestly cannot support one about *stochasticity* crossing it. That is Step 1b's claim |
| One lumped charged-tRNA pool replaces 20 per-amino-acid chains | Untested. **Ours, not the model's** | The pool size is a quantity we chose, and it may set the strength of the strongest candidate reverse channel. See R8 |
| CTP, UTP and the amino-acid pool are chemostatted | Their sources are in modules Core A′ cuts. **Ours, not the model's** | Stops being defensible the moment Step 1b makes the pools live |
| Non-ptsG membrane growth is exogenous | ptsG is Core A′'s only membrane protein. **Ours, not the model's** | Growth reaches only ~1.07×, so the published 2× cap is never tested |
| The 60 s rebuild is piecewise-constant | It is what the published model does | A continuous variant is a different and arguably better model. Declarable, defaults to published, and the deviation is measured rather than argued. See R11 |
| Michaelis constants come from the balanced table, not the hand-patch layer | Gives every constant a prior and one provenance | Eight of 32 differ from what the published simulator runs, one by 227×. **Ours**, labelled. See D3 |
| The base-to-nucleotide mapping is corrected | The published code permutes three of four | A departure from the published model, labelled, with the published mapping retained behind a keyword. See D4 |

### The tolerance principle

Adopted from the drafted module designs, which derive conservation tolerances
from the integrator rather than hand-picking them, and sharpened in one way:

> **The assertion is not "residual below tolerance". It is "residual scales with
> the integrator tolerance".** Run every conservation check at `(abstol, reltol)`
> and again at `(abstol/10, reltol/10)`, and require the residual to fall by at
> least fivefold.

A structural leak does not shrink when the solver is tightened, whatever its
magnitude; a numerical residual does. A fixed threshold can be passed by
loosening it, and this cannot. Report the single-run bound as a secondary
number, derived as

```
tol_C = N_restarts · Σᵢ |nᵢ| · max(abstol, reltol · maxₜ|xᵢ(t)|)
```

The `N_restarts` factor is why check 4 must run a full cycle: at 6,300 handshakes
a per-run bound understates the accumulated bound by nearly four orders of
magnitude.

Pin `abstol = 1e-10` mM, `reltol = 1e-8`, a stiff solver, and 60 s save points.
**Every check also carries a mutation test** that perturbs a stoichiometry and
asserts the check fails naming the right quantity. A conservation test that
cannot fail is not evidence.

### Check 0, which precedes all six in the note

**The count-to-concentration round-trip policy.** Across ~6,300 handshakes, each
writing particle counts back as integers, deterministic rounding gives a bias
that accumulates linearly — up to 0.5 particle × 6,300 = 3,150 particles, about
0.156 mM of drift per species, roughly 140× the integrator bound. Stochastic
rounding accumulates as the square root, about 0.002 mM. Carrying the fractional
remainder is exact.

Fix one policy and assert its signature: exactly zero under fractional carry,
square-root scaling under stochastic rounding, linear under deterministic
rounding. **Reject deterministic rounding.** Recommend fractional carry, because
it is the only policy under which the tolerance principle above tests the
integrator rather than the rounding. Without this check first, a rounding
artefact and a genuine leak are indistinguishable, and every downstream
conservation number is uninterpretable.

### The correctness checks

| # | Check | What is asserted | Scope |
|---|---|---|---|
| 0 | Round-trip policy | Exact-zero, square-root or linear residual scaling per policy; deterministic rejected | Handshake phase, re-asserted on assembly |
| 1 | Non-negativity | Every state above the negative of its own integrator bound at every save point, naming the first state and time to violate | Per-module as a smoke test; evidence only assembled |
| 1b | Particle-floor honesty | Report every state's minimum in particles; flag any below 500. Cross-check the three smallest pools against a chemical-Langevin ensemble; exclude from the likelihood any observable whose ODE trajectory leaves the ensemble's 90% band by more than the assumed observation noise | Assembled |
| 2 | Carbon balance | Glucose in against lactate out plus intermediates plus biomass; and the homolactic ratio, which is **analytically exactly 2.000** | Assembled only — spans transport, glycolysis and export |
| 3 | Redox | NAD⁺ + NADH invariant | Per-module (glycolysis); exactly invariant there, so assembly adds nothing |
| 4 | Adenylate and guanylate, **over a full 6,300 s cycle** | Each moiety separately, corrected for declared drain and inbound flux. Three configurations: all five recycling reactions (conserved); adenylate kinase removed (**pool exhausted at ~141 s**); pyrophosphatase removed (pyrophosphate unbounded) | Assembled only. Standalone it is a tautology, which is why the drain double is the acceptance criterion and not scaffolding |
| 4b | Phosphate closure | Free phosphate plus every phosphorylated species with pyrophosphate counted twice, minus flux across the two inbound mass edges. Two forms: exact with the GTP-branch reactions inactive, flux-corrected with them active. **Subtract the flux; do not relax the tolerance** | Assembled only |
| 5 | Carrier conservation | **Four** independent sums, one per phosphotransferase carrier, each named separately so a failure localises | Per-module before translation exists; restated assembled as "conserved up to what translation adds", again by subtraction rather than by widening |
| 6 | Nominal trajectory | All of the above at published parameters, plus check 7's census | Assembled |
| 7 | Clipping census | Count handshakes at which any deferred counter carries a deficit, at published parameters and across 200 prior draws. Require **zero at published parameters**. If more than 5% of prior draws clip, the non-smooth drain is in the operating regime rather than at a measure-zero point, gradient-based sampling of the ODE block is invalid as posed, and smoothing becomes mandatory | Assembled |
| 8 | Metabolic control analysis identities | On the frozen-expression chemostatted variant at steady state, the summation theorems hold exactly: flux control coefficients sum to one per flux, concentration control coefficients sum to zero per metabolite. Require both within 1e-6, computed by automatic differentiation through the steady-state solve | Assembled, and independent of every balance check because it tests derivatives rather than balances |
| 9 | Observation scale recovery | Generate at a known noise scale, infer it, assert recovery inside the calibration band | Inference phases |

Check 8 is the one genuinely analytic result the implementation must reproduce,
and the same Jacobian is the object the observable design needs (§4 D9).
Check 9 exists because no test anywhere in `test/` recovers a known noise scale,
and a variance-versus-standard-deviation misreading of the multivariate normal
constructor would rescale every credible interval in the project while passing
every other check.

### External validation, and what its absence costs

**Two external comparisons exist. Both are aggregate scalars. Neither is a
phenotype.**

1. **Predicted transcript steady states against measured counts.** Predicted
   0.44–2.87 copies against measured means 0.29–2.18, landing within about 1.5×.
   This is genuinely out of sample: the parameters come from proteomics, the
   comparison data from transcriptomics. Assert Spearman correlation at least
   0.7 across the 17 genes, and agreement within a factor of two for at least 15
   of 17.
2. **Protein fold change over one cycle against the published histogram**
   (median just below 2, tail to ~6×, long genes underproducing). This is *not*
   circular: initial protein count and transcript steady state are both
   proportional to the proteomics count, so the fold change is independent of it.
   It therefore tests the global constants and the promoter normalisation, not
   the individual promoters. Assert a median across the 17 genes in [1.7, 2.3],
   no gene below 1.5 or above 3.0, and the sign of the fold-versus-length slope
   matching the published direction.

Analytic rather than external: the homolactic ratio of exactly 2.000, the net
two ATP per glucose, and check 8's identities.

**Everything else is an internal invariant.** No metabolomics comparison, no
flux comparison, no doubling time.

**The cost, stated plainly.** Core A′ can pass every check and still be wrong
about the biology in ways no check would catch, so every quantitative statement
about the *cell* must be labelled inherited-from-the-published-model rather than
validated. More sharply: both external checks reduce to a single scalar each, so
a systematic error common to all 17 genes — a wrong ribosome turnover constant, a
wrong charging stoichiometry, a wrong tRNA pool size — has nowhere to show up.
**There is no external check in Core A′ with enough degrees of freedom to
localise such an error.** That is why the kill criteria in §8 are
internal-consistency criteria rather than goodness-of-fit criteria.

---

## 4. Approach

### D0. Framework before modules, and the kill risk on a toy first

**Decision.** Two protocol-changing framework pieces land before any module, the
1 s handshake is proven on a two-module toy with wall-clock measured, and only
then are the seven modules built.

**Why.** `dev/plans/reduced-syn3a-wave-plan.md` places all seven modules first
and all three interfaces afterwards, and it makes the framework files read-only
to a module branch. Those two rules are individually right and jointly
impossible, because two of the framework changes alter the contract the modules
are written against: letting a module contribute to a state it does not own
changes what `dynamics` returns, and fixing jump composition changes what
`reactions` receives. Discovering either after the modules land re-opens seven
merged pull requests and spends the seven-way parallelism twice.

The handshake is separate and more urgent. It has zero lines of code (G2), a
mixed composition is refused in one line (G1), and it is the only genuine kill
risk in the project. Proving it needs two toy doubles, not seven real modules.
`dev/notes/modular-bayesian-inference-heterogeneous-modules.md` prescribes
exactly this: "Two modules, one ODE and one CTMC, sharing one parameter… This is
your reference." The scoping note's own principle is pilot before scaling.

**Alternatives considered.** Keeping the wave plan's order, accepting that the
three ODE and four stochastic modules would be revisited — rejected because the
cost is seven pull requests against one. Landing the protocol changes first but
building the handshake directly on the real modules — rejected because it defers
the kill measurement past the point where acting on it is cheap.

**Consequence to record.** The four drafted module designs each open by stating
that the interface is frozen and read-only. That was true when they were written
and is now false in two specific respects. Phases 6 to 12 consume the *amended*
protocol, and each amendment is a phase deliverable rather than a module's
concern.

### D1. `Mode` is the column that runs

**Decision.** Every imported kinetic constant and initial concentration comes
from the balanced table's `Mode` column, with its geometric standard deviation,
and the source file recorded.

**Why.** The published simulator reads `Mode`. The alternatives are wrong by
large factors: the unconstrained geometric mean and, worse, the "catalytic rate
constant geometric mean", which is the geometric mean of forward and reverse and
so is orders of magnitude below the forward rate for a near-irreversible
reaction. Reading the wrong column gives a model that silently starves or has
imaginary headroom. `dev/notes/well-stirred-minimal-cell.md` records getting this
wrong twice before landing on `Mode`.

### D2. The nucleotide file governs the nucleotide species, and the trap is wider than recorded

**Decision.** Every value that appears in more than one balanced file declares
which file governs, through the loader's existing machinery, and loading fails
rather than choosing silently.

**Why, with the numbers.** The scoping note flags two GTP-branch reactions as
appearing in both the central and nucleotide files with different modes. **All
five nucleotide-module reactions do, and three disagree by more than the two on
record:**

| Reaction | nucleotide `Mode` | central `Mode` | ratio | central gstd |
|---|---|---|---|---|
| PGK3 forward | 140.8 | 319.46 | 2.3× | 2.1e+63 |
| PYK3 forward | 672.84 | 1873.98 | 2.8× | 2.1e+63 |
| ADK1 forward | 319.16 | 228.06 | 1.4× | 7.0e+35 |
| **GK1 forward** | **410.23** | **7.60** | **54×** | 2.1e+63 |
| **PPA forward** | **646.73** | **583611.61** | **902×** | 1.7e+33 |

The nucleotide file governs all five. The central file's values are not a rival
estimate but an artefact of balancing with no data, and the geometric standard
deviations say so mechanically: 1e33 to 1e63 against 1.05 in the nucleotide file.
**Taking the central file's pyrophosphatase constant would have run the enzyme
902× too fast**, which is the class of error that hides because it makes a
conservation check pass.

The same trap runs in the concentration table and cuts the other way. AMP, GTP,
GDP and GMP read 0.0832, 1.6627, 0.2981 and 0.0117 mM in the nucleotide file and
sit at the 0.1 mM prior default in the central file. Reading the guanylate
species from the central file would understate that pool by an order of
magnitude — about 4,000 particles rather than 39,800.

**Consequence.** Provenance is not defensive paperwork. On the recycling module
it is the difference between a working model and a plausible-looking broken one.

### D3. Michaelis constants come from the balanced table, and that is a labelled deviation

**Decision.** The 32 Michaelis constants of the glycolytic module are imported
from the balanced `Parameter` table, not from the `Quantity` table the published
simulator reads. This is registered as a departure from the published model.

**Why.** These constants exist twice in the source file and **eight of 32
disagree**:

| Identifier | `Quantity` (what runs) | balanced `Mode` (ported) | gstd |
|---|---|---|---|
| `kmc_R_PGI_M_g6p_c` | 0.28 | 22.9419 | 1.0512 |
| `kmc_R_PGI_M_f6p_c` | 0.15 | 3.3488 | 1.2986 |
| `kmc_R_PFK_M_atp_c` | 0.117 | 0.0255 | 4.8229 |
| `kmc_R_FBA_M_fdp_c` | 0.005 | 0.2447 | 3.5730 |
| `kmc_R_FBA_M_dhap_c` | 0.095 | 0.0064 | 4.5985 |
| `kmc_R_FBA_M_g3p_c` | 1.0 | 0.0044 | 8.4265 |
| `kmc_R_GAPD_M_nad_c` | 1.3 | 5.6969 | 1.1808 |
| `kmc_R_PGM_M_2pg_c` | 1.47 | 0.0278 | 6.0168 |

The balanced column gives every constant one provenance and a real inherited
prior, which is what the inference problem needs; the `Quantity` column is a
hand-patch layer carrying no uncertainty at all. The cost is that Core A′ is not,
on these eight numbers, running what the published simulator runs, and two of
them move by 82× and 227×. So the choice is **ours**, registered through
`reduction_notes` so `reduction_declarations` enumerates it.

This is a third instance of the column trap the scoping note records twice. The
note's own summary — "~35 Michaelis constants from the balanced file (with
priors)" — assumes the balanced column is what runs. It is not, and the note
gains a correction-of-record row.

**Alternatives considered.** The `Quantity` values, which run but carry no
uncertainty, so eight parameters — the loose ones, exactly the ones inference
would target — would become asserted priors. Or `Quantity` values with balanced
widths, which runs what runs and keeps the widths but leaves eight priors whose
median is not their point value, a subtler thing to explain in every downstream
result than a column choice.

**A pleasant asymmetry worth recording.** All 17 Michaelis constants in the
nucleotide file agree between the two tables. So the two ODE modules end up with
different provenance stories through no inconsistency of ours, and the recycling
module labels no deviation on this axis while its sibling labels one.

### D4. The base-to-nucleotide mapping is corrected, and the correction is labelled

**Decision.** Transcription charges each base against the nucleotide it is
actually polymerised from. The published permutation is retained behind a keyword
so the difference can be measured.

**Why.** The published rate function assigns three of four base counts to the
wrong nucleotide pool, because the local variable names run A, U, C, G while the
dictionary keys run A, C, G, U:

| Assigned | Gets the count of | Charged against | Should be |
|---|---|---|---|
| A | A | ATP | A — correct |
| U | **C** | UTP | U |
| C | **G** | CTP | C |
| G | **U** | GTP | G |

It is present on both code paths, including the production 60 s rebuild.

The effect on the rate constant is about one percent, because the denominator is
dominated by transcript length. **The effect on the coupling is not.** Over the
17 genes the base counts are A 7236, C 2078, G 3094, U 5868, so the bug weights
the GTP term by the uracil count rather than the guanine count and makes every
transcription rate constant **1.9× more sensitive to the live GTP pool than it
should be.** That is the only reverse channel in the model, so porting the bug
would mean demonstrating bidirectional coupling with a gain a typo nearly
doubled.

**This is the unusual case where the default is the labelled one.** Correcting is
a departure, so the correction appears in `reduction_notes` while the published
mapping does not.

### D5. The promoter proxy is inherited, and doubly flagged

**Decision.** Per-gene promoter strength is the published proteomics proxy,
protein copy number over 180, inherited verbatim. Two separate declarations are
registered: that it is a proxy, and that conditioning on proteomics would be
circular.

**Why.** Replacing this proxy is the point of the later surrogate work, and a
baseline you have not implemented cannot be improved on. But it is not a measured
parameter — it is back-calculated from the steady-state protein abundance the
model is supposed to predict — and the same protein counts also set every enzyme
concentration in three other modules. So almost any protein-level observable
closes a loop through it. The constraint is recorded in code at the point the
proxy is defined, so the observable design in phase 15 meets it when choosing
rather than when interpreting a posterior. See D8 for the operational rule.

### D6. Lactate export is not optional, and external lactate carries a volume ratio that is ours

**Decision.** Lactate export is implemented as passive diffusion with the
permeability and cell radius carried as separate named quantities, and external
lactate is integrated as a dynamic state with its accumulation scaled by a
medium-to-cell volume ratio of 1e5.

**Why export is required.** Glycolysis produces two lactate per glucose. At Core
A′'s corrected demand of ~1,106 glucose/s that is ~13.9 M lactate over a cycle:
**~691 mM at initial volume and ~345 mM even at doubled volume**, exceeding the
entire measured phosphate pool by a factor of forty. Without export the core is
not a sustained pathway, it is a pathway that poisons itself in minutes, and
carbon balance is open so check 2 cannot even be posed. The earlier Core A spec
deleted this reaction; that was an error and it is recorded so it is not
repeated. Correcting the glucose demand for tRNA charging made the case for the
exporter twice as strong, not weaker.

**Why the ratio, and why 1e5.** The published model clamps every external
species, so published export drains into a pool pinned at zero and never
saturates. The registry instead lists external lactate among its 32 dynamic
states, and the resolver rejects a clamp on a state a module integrates.
Integrating it in a single shared volume would make export saturate as the pools
equilibrate — the 691 mM would not leave, it would move outside — which is a
different model, not a faithful one. So the ratio reproduces the published
unsaturated efflux through a genuinely dynamic state. It is chosen against a
stated criterion rather than picked: external lactate after a full cycle must
stay below one percent of steady cytosolic lactate, which 1e5 meets at 0.47% and
1e4 misses at 4.7%. **The ratio is ours**, registered, and a test must show that
a ratio of one visibly saturates so the default's purpose is demonstrated.

**Why the radius is not folded in.** `3P/r` collapses to a single constant at 200
nm, and folding it would be tempting and wrong, because the growth phase makes
the radius grow and the export rate constant then changes with it.

### D7. Phosphotransferase initial conditions come from the model's own proteomics fractions

**Decision.** The eight carrier phospho-states are initialised from published
copy number times published proteomics fraction, and both are recorded as the
published model's rather than asserted by us.

**Why.** The registry marks all eight `not_imported`, and the source offers two
candidates that disagree. The transport table gives a uniform 0.024 mM to all
eight, a placeholder putting 969 copies behind every protein and contradicting
the proteomics counts for three of the four carriers. The model's own data files
give the real answer: unphosphorylated/phosphorylated fractions of 0.05/0.95 for
three carriers and 0.15/0.85 for ptsG, in a `proteomics_fraction` column read at
startup rather than a disabled setting. Applying those to the copy numbers
reproduces the totals exactly — 353, 314, 290 and 831 — so both the totals and
the split are the published model's.

**The ten cascade constants and the membrane permeability carry no priors
anywhere in the source.** The transport file has no `Parameter` table at all.
They are point values, so every prior on the eleven is **asserted by this
project**, marked as such, enumerable from the composed model, and fixed by
default so nothing is inferred through an invented prior without a decision.

### D8. Observable design: what is safe, and the one rule that matters most

**Decision.** Condition on per-gene transcript counts as counts, and on
metabolite concentrations, jointly. Rule out fluxes as a likelihood term and
protein counts from the headline.

**Why each.**

| Observable | Verdict | Reason |
|---|---|---|
| Per-gene transcript counts, per replicate cell, 60 s | **Mandatory** | Not used to set any Core A′ parameter. The only stream that separates the 17 promoters from each other and identifies the decay constant orthogonally to promoter amplitude |
| Metabolite concentrations, 60 s — the small glycolytic intermediates and the adenylate/guanylate pools | **Mandatory** | Independent of every parameter's source. Carries both the concentration-control signal and the expression-drain signal. This is the *clean* route to promoter strengths, and the route that survives to Step 2 |
| Cell volume / fractional growth | Include, low weight | The only second route to the ptsG promoter. Weak, but structurally distinct |
| Reaction fluxes | **Rule out as a likelihood term** | At 3.6% of the tightest enzyme's capacity, glycolytic flux is demand-limited, so flux control coefficients are near zero and flux is nearly blind to enzyme abundance. Concentrations carry the same information without the degeneracy. Keep fluxes as a reported diagnostic |
| Per-gene protein counts | **Rule out of the headline** | Two reasons. It closes the promoter-proxy loop (D5). And it will be unavailable non-circularly at Step 2, so a demonstration whose power comes from it demonstrates something that does not transfer |
| Protein initial conditions as free parameters | **Rule out; hold fixed** | The circularity in its sharpest form: freeing them lets the likelihood re-estimate the exact quantity the promoter prior was built from |
| Doubling time or any phenotype | Rule out | Out of scope by construction |

**The rule that matters most, and it is easy to get wrong.** Circularity lives in
the prior-and-likelihood *pair*, not in the observable alone. For synthetic
recovery the data comes from a chosen truth, not from the proteomics table, so no
dataset is used twice. What bites instead is a tautology:

> **For any coverage or calibration result, the ground truth must be drawn from
> the prior. Never from the nominal published values.** If the prior is centred
> on the proxy, the truth is *set* to the proxy, and the observable is
> protein-side, then prior and likelihood agree by construction and the posterior
> comes out tight and covering for no reason at all. A fixed-truth-at-nominal run
> is a smoke test and must be labelled as one.

**Two existing-code defects this forces.** `observe` on a trajectory ensemble
averages replicates, and the summary statistics on the path the sampler uses keep
per-time means only. Transcript counts must be carried as a replicates × genes ×
times array, and the summaries must gain distributional terms — a per-gene count
histogram over replicates plus a lag-one autocovariance for the decay constant —
or be superseded entirely by an exact likelihood (see D10).

### D9. The identifiability structure, and the figure that settles the observable choice

**Decision.** Compute the concentration control coefficient matrix — 17 enzymes
against roughly 20 metabolites, at nominal parameters, by automatic
differentiation through the steady-state solve, with check 8's summation
identities as the correctness guard. Choose the metabolite panel from its largest
rows.

**Why this and not a sensitivity scan.** The scoping note lists "which observable
is most informative about the boundary-crossing parameters" as an open question
and says to answer it by simulation. The control coefficient matrix *is* the
answer, it is cheap, and it comes with an analytic self-check.

**The structure it will reveal, and one hard limit.** The deferred cost counters
are **aggregate**, so the expression-drain channel delivers one identified linear
combination of the 17 promoter strengths, weighted by gene length — not 17
parameters. The enzyme-concentration channel delivers gene-specific control
coefficients, but the observable pools are shared. So metabolite-side data
identifies roughly one aggregate direction plus a small number of control
directions in a 17-dimensional space, while transcript counts identify all 17.
Both streams are in the likelihood, and the reason is not that more data is
better — it is that they identify structurally different subspaces.

**The gain asymmetry is a finding, not an obstacle.** The two forward channels
are order one: the ATP pool turns over in about 109 s and the GTP pool in about
30 s, both driven by expression cost, so a twofold change in aggregate promoter
strength is unmissable in the adenylate and guanylate pools. The reverse channel
has a measured gain of 0.044 to 0.051. **The published coupling is therefore
nearly unidirectional in gain even though it is bidirectional in topology**, and
that is a substantive quotable statement about whole-cell model architecture
rather than a shortcoming of the reduction. It is worth more than a symmetric
result would have been.

### D10. Three inference rungs, and the reference is mandatory

**Decision.** Build a reduced exact reference first, then the production sampler,
then the protocol comparison.

**M0, the reference.** A deliberately small Core A′: two genes, 600 s horizon,
small enough that a very long run is a defensible ground truth. Data-augmented
blocked sampling: augment the stochastic path, so conditional on it the cost
drains are known deterministic forcing and the ODE parameters go to a
gradient-based sampler at tightened tolerances; conditional on the ODE
trajectory, and hence on the rate-constant sequence, the stochastic block's
likelihood for observed counts is closed-form. This targets the true joint
posterior and is the only unambiguous ground truth the project will get.
**Without it, a calibration pass is unattributable** — it could be a good
approximation or a converged wrong answer.

**M1, production, on the full 17 genes.** Same blocked structure. The honest
difficulty is that the 1 s drain means the ODE sees the stochastic *path*, not
just its 60 s endpoints, so a closed-form kernel alone does not give a joint
likelihood. Force the choice with a measurement rather than an argument:

> Compare the nominal trajectory under 1 s, 5 s and 60 s drain granularity. If
> the largest relative difference in any observed pool is below one percent, use
> the coarser aggregate drain. Given the ATP pool's 109 s and the GTP pool's
> 30 s turnover, 60 s is expected to fail this and 5 s to pass.

Whichever passes becomes a labelled reduction with a *measured* cost. If path
augmentation is genuinely needed, the fallback is pseudo-marginal sampling with a
bootstrap particle filter, which must be built — the single largest missing piece
if the closed-form route does not survive contact.

**M2, the protocol comparison.** Run the same problem, the same data and the same
targets through the joint posterior and through the explicit cut in
`src/boundary.jl`, and compare their calibration. One model, two interfaces, two
calibrations, and the typed interface contract is what makes the swap a one-line
change — which is the architecture claim's evidence.

**Alternatives considered.** Including the iterative exchange protocol as a third
arm, with its predicted prior-double-counting signature. Rejected: it is a
proof-of-concept rather than a production protocol, and measuring its calibration
is investment in something the project has already decided not to carry forward.
Its known limitation is stated in §7 and no compute is spent on it.

**Why the reference matters to the budget.** If the stochastic block's likelihood
is closed-form, a posterior costs a few CPU-hours and a calibration study a
couple of CPU-weeks. On the simulation-based path as currently coded, the same
study is CPU-years. That is a two-order-of-magnitude difference resting on one
structural claim, which is why it is verified in a named phase rather than
assumed (§9).

### D11. The first inference target set

**Decision.** Six free parameters, plus per-modality observation noise.
Everything else fixed at published values, then freed selectively.

| Target | Nominal | Why this one |
|---|---|---|
| `S_GAPD` | 7.528 | Strongest promoter, ~52 transcription events per cycle, ~14% relative standard error per cell. The best-determined promoter: if this does not recover, nothing will |
| `S_PGI` | 1.478 | Near-weakest, ~7.9 events per cycle, ~36% per cell. Deliberately the hard case, so the set spans the achievable-precision range |
| `S_ptsG` | 4.617 | **The headline parameter.** ptsG is Core A′'s only membrane protein, so it is the only promoter whose uncertainty crosses the boundary by *two* independent routes: the expression drain and volume-mediated dilution of every concentration |
| `krnadeg` | 3.5044 /s | Identified by the transcript autocorrelation, orthogonally to promoter amplitude |
| `kcat_ENO` | 62.18 | **Positive control.** Tightest prior in the core and the capacity-tightest enzyme, so its concentration control coefficients are the largest available. Must recover tightly |
| `kcat_FBA` | 59.7 | **Stress control.** Loosest prior in the core (gstd 1.747). The diagnostic is shrinkage: its posterior must be visibly wider than ENO's |

**One deviation from the scoping note, with the reason.** The note recommends a
shared gene-expression global, naming the polymerase turnover constant. That
constant is **exactly degenerate** with the promoter strengths: only their
product enters the rate law, and the cap on that product never binds — the
largest promoter-scaled turnover in the core is 8.85 nt/s against a ceiling of
180. So the pair admits a multiplicative ridge, and the aggregate metabolite
channel is a function of the product only. Fixing 14 promoters at published
values would anchor the scale, but then the anchor's correctness becomes
load-bearing for a headline parameter, which is a bad trade. The decay constant
is the honest substitute: it sets transcript lifetime, which the autocorrelation
identifies independently of amplitude.

Free the polymerase constant in a later round *as a deliberate demonstration of
the ridge*. That is a result, not a failure.

**Two hard constraints.** Hold initial protein counts fixed at proteomics values
(D8). And replace the single scalar observation noise: the states span 0.0098 to
17.8 mM, so one additive term makes the small informative pools contribute
nothing. Require per-modality noise on a log scale — lognormal for metabolites,
a count model or exact counts for transcripts. This is a genuine defect for this
application (G4), not a preference.

**Replicate count, derived rather than chosen.** Relative standard error on a
promoter strength from one cell is the inverse square root of its event count, so
36% at the weakest and 14% at the strongest. Reaching 10% on the weakest needs
about 13 cells and 5% about 51. For a metabolic parameter detected through the
reverse channel, the end-to-end gain is the product of a control coefficient and
0.045, so a three-sigma detection of an e-fold change needs on the order of 11
cells optimistically and 850 pessimistically. **The reverse channel being weak is
not fatal, it is expensive, and the price is quotable.** Default the synthetic
ensemble to 200 cells with 850 as the contingency.

### D12. What per-module calibration means operationally

The modules are not independently observable in the assembled model. In the
synthetic setting the interface variable can be generated and held, which is
exactly what makes a conditional calibration check well defined:

1. **ODE block given the stochastic path.** Draw ODE parameters from the prior,
   draw a path from the prior predictive, simulate the ODE with the path as known
   forcing, infer with the path held at its simulated truth, take ranks.
2. **Stochastic block given the rate-constant sequence.** The mirror image.
3. **Joint**, for the shared decay constant — the only place the composition can
   go wrong without either conditional failing.
4. **Message calibration.** Record each module's outgoing message at each
   replication and rank the truth within it. Prior double-counting shows up in
   the message before it shows up in the final posterior, which is what makes
   this diagnostic worth building. It is directly the survey's open problem 5.

**The attribution rule, stated once.** Calibration checks detect sampler error
and approximation error together. Run the two conditionals; if both pass and the
joint fails, the fault is the composition. That is the only route by which a
calibration failure means "architecture" rather than "sampler", and it is why M0
is mandatory rather than nice to have.

Two mechanics: ranks for a weighted posterior must be computed with the weights,
and use the rank-ECDF test with simultaneous bands rather than a bare histogram,
because a ten-bin histogram has no usable power at the replication counts this
project can afford.

---

## 5. Data and reproducibility

### Inputs and provenance

Everything Core A′ imports comes from `Luthey-Schulten-Lab/Minimal_Cell` at
commit `db048ac`. Nothing is fetched at run time.

| Input | Upstream | How it arrives |
|---|---|---|
| Glycolytic kinetics and initial conditions | central balanced SBtab | Vendored derived extract, ~65 rows, reshaped to the loader's single-identifier convention including the `conc_<species>` form so the registry-agreement check fires |
| Recycling kinetics and initial conditions | nucleotide balanced SBtab, **and** the central file's rival values for the same identifiers | Two extracts, so the ambiguity is real rather than described |
| Cascade constants and permeability | transport TSV | Vendored extract with no uncertainty column, so every row loads as asserted |
| Carrier initial conditions | proteomics counts × proteomics fractions | Vendored, with the derivation recorded |
| Per-gene sequence data | genome record, proteomics table, measured transcript counts | One extract, three upstream sources, one column each, so no cross-file ambiguity |

Extracts live under `src/organisms/coreA/data/` with a `README.md` recording the
upstream filename, the commit, which table each row class came from, that the
reshape renames identifiers and changes no value, and the exact command that
regenerates it. Regeneration scripts live under `dev/scripts/`. **Every extract
must round-trip:** re-running the command from the README produces an unchanged
file.

The raw upstream tables are not vendored, because the loader's single-identifier
convention does not fit the balanced `Parameter` table's three-part key and
`src/loader.jl` is not amended for it.

### Environment

Julia 1.10.5, available on compute nodes only through a module load. Dependencies
are pinned in `Project.toml`; `DiffEqCallbacks` and the stiff-solver subpackages
are already in the depot, so no new dependency requires a network fetch. **If a
new dependency is genuinely needed, the depot must be populated from the login
node first**, because compute nodes have no network and the failure mode is a
job dying in seconds on a host-resolution error.

Tests run through `sbatch test/run_tests.slurm`, which sets the offline flag that
skips the one quality check needing the package registry. The integration gate
is a separate environment variable and is off by default, so the Slurm job runs
the fast suite only. Long-running work gets its own driver and Slurm wrapper
alongside the existing three, each self-asserting.

### Seeds, cost and what is regenerable

Every synthetic dataset stores its own ground-truth parameters, its seed, and the
`reduction_report` of the model that produced it. Two runs at one seed must be
identical. Wall-clock per trajectory is unmeasured today and is the binding
practical number, because the inverse problem needs thousands; it is measured
twice, on the toy in phase 3 and at full scale in phase 13, and it drives the
first kill criterion.

Regenerable: every extract, every synthetic dataset, every figure. Irreplaceable:
nothing, provided the extracts' regeneration scripts and the upstream commit
reference stay correct — which is what the round-trip requirement protects.

---

## 6. Outputs

The claim decomposes into five sub-claims, and every output hangs off exactly
one: **C1** the boundary is genuine; **C2** the forward model is correct; **C3**
uncertainty propagates; **C4** it propagates *correctly*; **C5** the architecture
is the enabler.

| ID | Output | Claim | The number it delivers |
|---|---|---|---|
| **F1** | **Channel gain table.** Four rows — enzyme concentration, expression drain, pools to rate constants, volume — each with its analytic form and its measured value | C1 | The reverse channel at 0.044–0.051; the forward channels as pool turnover times of ~109 s and ~30 s; volume's contribution. **The most important table in the work, and it appears early** |
| **F2** | **Concentration control coefficient heatmap**, 17 enzymes × ~20 metabolites, with the summation identities as guard | C1, and the observable choice | Answers the scoping note's open question. Its largest rows *are* the metabolite panel |
| F2b | Flux control coefficients, same layout | C1 | Shows them near zero at 3.6% utilisation — the quantitative reason fluxes are ruled out of the likelihood |
| **F3** | **Invariant residual against integrator tolerance**, log-log, one line per invariant, slopes required positive. Beside it a **mutation table**: per check, the injected error, the residual it produced, the bound it exceeded | C2 | The tolerance principle, and the evidence that every check can fail |
| F4 | Conservation residual against handshake count under the three rounding policies | C2 | Exact-zero, square-root, linear. Justifies the policy and forestalls a rounding artefact being read as a leak |
| **F5** | **The two external comparisons.** Predicted against measured transcript steady states, 17 points, log-log with a twofold band; and the protein fold-change histogram with the published median marked | C2, external | The only two places Core A′ touches data it did not consume |
| F5b | Analytic drain-noise calculation against simulated | C1 | Closed-form variance of the cumulative expression drain, corroborating the reduction's noise estimate from an independent route |
| **F6** | **The shrinkage table.** All six targets × three data configurations — transcripts only, metabolites only, joint — reporting posterior over prior standard deviation and the prior-to-posterior divergence | **C3. This is the figure the whole work is for** | Two load-bearing cells: the ptsG promoter under *metabolites only* (shrinks ⟹ uncertainty crossed toward the ODE block) and the ENO constant under *transcripts only* (shrinks ⟹ crossed the other way) |
| F7 | Leave-one-stream-out: refit dropping each transcript stream and each metabolite stream, reporting the change per target | C3 | Identifies load-bearing single points of failure, and streams that never informed anything |
| **F8** | **Rank-ECDF difference plots with simultaneous bands**, per target, at the largest affordable replication count | C4 | Primary calibration evidence. Histograms are secondary |
| F9 | Coverage curve: nominal against empirical at 50, 80, 90, 95 and 99 percent, per target, with binomial intervals | C4 | The two controls are the diagnostic rows, with the kill-criterion band drawn on |
| **F10** | **Reference comparison** on the two-gene system: exact blocked posterior against production posterior, overlaid marginals plus a two-dimensional contour, plus the divergence between them | C4 — **without this a calibration pass is unattributable** | The number that makes every other calibration claim mean something |
| **F11** | **The architecture figure.** Two boundary protocols — joint and explicit cut — on the same data and targets, with calibration and coverage panels for each | **C5** | One model, two interfaces, two calibrations. Prediction: the joint is uniform; the cut is biased and over-wide on whichever parameter drew its identifiability from the severed feedback |
| **F12** | **Message-calibration histograms** for the outgoing boundary message under each protocol | C5, and the methodological contribution | Shows a site miscalibration before it reaches the posterior |
| F13 | Jacobian singular-value spectrum for the target set, with and without the polymerase constant freed | C3 | The multiplicative ridge as a singular value dropping by orders. Uses the existing `check_identifiability` |
| F14 | Wall-clock per trajectory against horizon and gene count, the budget per rung, and the achieved replication count | K1, reproducibility | The scoping note's binding practical number |
| **T1** | **Provenance table**, generated by the loader and `reduction_declarations` — every parameter with its value, prior, width, informedness, source file and rejected alternatives | Honesty, C2 | Machine-generated, never typed. The audit trail for the cross-file errors of record |
| **T2** | **Reduction declarations with measured costs** — the lumped charging step, the chemostatted pools, exogenous membrane growth, the corrected mapping, the rounding policy, the drain granularity, the lactate volume ratio, any smoothed counter | Scope honesty | Each row carries the *measured* cost where measured. This table is what stops a result depending on our choice unremarked |
| **T3** | **Kill-criteria scoreboard** — each criterion, its threshold, its measured value, pass or fail | Falsifiability | **Published whether or not everything passed.** A scoreboard with a fail on it is more credible than one without |
| T4 | Target-set summary: truth, posterior mean, 90% interval, coverage, shrinkage | C3, C4 | The conventional recovery table |

**Scalars that must appear in any abstract:** wall-clock per trajectory; the
reverse-channel gain range; the achieved replication count; 95% coverage per
target; the shrinkage of the ptsG promoter under metabolite-only data; and the
coverage gap between the joint posterior and the cut.

---

## 7. Non-goals

- **Reproducing syn3A's phenotype.** No comparison to the 105-minute doubling
  time, the measured metabolome, or the proteome. Growth is fractional growth or
  time-to-threshold, and the code refuses the comparison rather than merely the
  prose.
- **Step 1b.** The ribonucleotide-reductase and dNTP branch, which is what makes
  the nucleotide-triphosphate pools genuinely live and adds a genuinely
  stochastic boundary. Out of scope, and the five-reaction nucleotide module here
  stops being defensible the moment it lands.
- **Real data.** Step 2's business. Everything here is synthetic, which is the
  only route that admits a ground truth to check against.
- **The surrogate work.** Replacing the promoter proxy with a learned
  sequence-to-strength prior, and replacing the metabolic integrator with a
  differentiable surrogate. This spec implements the baseline both are measured
  against and no more.
- **Model comparison and structural uncertainty.** Every comparison quantity is
  a functional of a per-model posterior and cannot be computed before parameter
  inference works. What is committed now is that structural identity and
  normalisation semantics are tracked, so the layer can be added without
  rebuilding the interfaces.
- **Adding a two-state promoter.** A defensible extension and a real
  model-comparison question, but a deliberate departure from the published model
  rather than a port of it.
- **The chemical-Langevin formalism**, despite being a valid declared value and a
  natural fit for the 266-to-1355-copy enzymes. Unimplemented, and not needed for
  a claim about parameter uncertainty.
- **Fixing the iterative exchange protocol.** It is a proof-of-concept of
  cross-module exchange, not a production protocol. Its known limitation — naive
  iterated conditioning re-counts the prior each pass — is stated and no compute
  is spent measuring it.
- **Decoupling the framework from Core A′.** The resolver and loader still reach
  for the global registry. Parameterising them over a registry object is what a
  second organism costs, and it is deliberately unpaid.

---

## 8. Kill criteria

Compute is unconstrained, so each of these is a threshold on what can be
*claimed*, not on what can be afforded. All are recorded on the scoreboard (T3)
whether they fire or not.

**K1 — the calibration claim, through feasibility.**
*Fires when* measured wall-clock forces fewer than 50 calibration replications on
the assembled model. Concretely, above roughly 10 s per full-cycle trajectory on
the simulation-based path.
*Why 50:* below it the rank-ECDF test cannot distinguish a 95%-to-85% coverage
error from noise, so the posterior is uncheckable — and an uncheckable posterior
is not a result. This is a scientific criterion, not a budgetary one.
*Pre-committed fallback ladder, in order:* shorten the horizon; reduce to five
genes; run calibration on the reduced reference system only. If the ladder
bottoms out, the headline is restated from "is calibrated" to "runs, with
coverage at a single truth", **and that restatement appears in the abstract.**

**K2 — the bidirectional claim.**
*Fires when*, with 1,000 synthetic cells and transcript-only data, the
prior-to-posterior divergence is below 0.05 nats for **every** metabolic
parameter. That threshold is already the convergence tolerance in the existing
exchange protocol, so it is not invented for the occasion, and `kl_divergence` in
`src/boundary.jl` computes it today.
*If it fires:* the bidirectional claim is dead as stated. Report the forward
direction as the headline, since it is strong, and the reverse direction as a
*measured bound with the replicate count needed to detect it* — which is a real
result about published whole-cell architecture, not a failure — and promote Step
1b, where the pools are genuinely rate-limiting.
*Explicitly rejected:* the version of this criterion phrased on the coupling gain
alone. A gain of 0.045 is not fatal; D11 shows it costs tens to hundreds of
cells, which is affordable. **The kill has to be on the posterior, not the gain.**

**K3 — a calibration failure attributable to the architecture.**
*Fires when* both conditional checks pass and the joint fails its simultaneous
band at 95%.
*Not a kill — a redirection.* It says the composition, not the samplers, is
wrong, and it is the most scientifically interesting outcome available. Requires
the two conditionals to have been run, which is why D12's attribution rule is not
optional.

**K4 — the controls.**
*Fires when*, at 100 replications and nominal 95%, empirical coverage for the
tight-prior control falls outside 90 to 99 out of 100 (the binomial band).
*Second and sharper trigger, prior regurgitation:* posterior over prior standard
deviation **above 0.9** for the tight-prior control. A parameter with the core's
tightest prior and a directly informative observable must shrink; if it does not,
"recovery" is the prior reappearing and the coverage pass is empty. Require
shrinkage below 0.5 for at least one target.
*Third trigger:* if the loose-prior and tight-prior controls come out with the
*same* posterior width, the likelihood is not contributing what it appears to.
Investigate before reporting either.
*On power:* at 100 replications the 95% level detects overconfidence better than
the 50% level. Report both; the 95% level carries the test.

**K5 — the non-smooth drain is in the operating regime.**
*Fires when* any deficit is carried at published parameters, or when more than 5%
of prior draws clip (check 7).
*If it fires:* the ODE block is not differentiable where it matters,
gradient-based sampling is invalid as posed, and the smoothed model must be built
and validated against the clipped one before any posterior is reported.
**This is the criterion most likely to actually fire, and it fires in phase 3,
which is the good news.**

**K6 — structural non-identifiability of the target set.**
*Fires when* `check_identifiability` reports a rank below six for the chosen set,
or a condition number above 1e6. Computed *before* any sampling, at the cost of
one Jacobian.
*Predicted instance:* freeing the polymerase constant alongside the promoters
drops one singular value by orders (D11). Reparameterise to the product plus a
scale anchor rather than sample through the ridge.

**K7 — the reduction is doing the statistical work.**
*Fires when* the posterior width on the pyruvate-kinase constant changes by more
than twofold once the seven dropped reactions sharing its gene are stubbed in as
a competing sink.
*If it fires:* Core A′ is better-conditioned than the same subnetwork inside the
full model by a measurable factor, and that factor is quoted alongside every
recovery claim. The scoping note names this as something to watch; this makes it
a number.

**Rejected as criteria:** wall-clock in itself, since there is no compute
constraint and it only binds through K1; the coupling gain in itself, per K2; and
"forward trajectories look biological", which is out of scope by design.

---

## 9. Open questions

- **[NEEDS CLARIFICATION: is the stochastic block's likelihood closed-form?]**
  Conditional on the 60 s rate constants, transcription, decay and translation
  are zeroth or first order and the nucleotide pools enter as rate constants
  rather than reactants — so the network looks monomolecular, which would make
  the transition kernel closed-form and factorised over genes. **The lumped
  charging step has two reactants and may break this.** Settled in phase 12. The
  inference budget changes by roughly two orders of magnitude on the answer
  (D10), so this is the highest-value open question in the spec.
- **[NEEDS CLARIFICATION: what handshake granularity do the pools require?]**
  Settled by the 1 s / 5 s / 60 s comparison in D10 — in miniature in phase 4,
  at full scale in phase 13. Expected to rule out 60 s, given the GTP pool turns
  over in half the rebuild interval.
- **[NEEDS CLARIFICATION: is the charged-tRNA channel stronger or weaker than
  transcription's?]** The scoping note hopes it is stronger and would rescue the
  reverse direction. A first pass at the published rate law suggests the opposite
  — that both denominators are dominated by polymer length, bounding the whole
  reverse direction at a few percent structurally rather than accidentally.
  Settled in phase 10, and it changes how the result is framed.
- **[NEEDS CLARIFICATION: does protein fold change over a cycle reach two?]** A
  first pass suggests it may fall well short, which would indicate the lumped
  charging step throttling translation. Settled in phase 10, and far cheaper to
  diagnose there than to explain at the end.
- **Does the lumped charging step's stoichiometry change the answer?** The
  published AMP-and-pyrophosphate form puts 553 events per second through the
  adenylate kinase; a two-ATP-to-two-ADP form puts zero. That changes how much
  the kinase's equilibrium constant can move the ATP/ADP ratio the rebuild reads,
  and so whether the kinase belongs in the target set. One cheap comparison run,
  in phase 12.
- **Does the 60 s rebuild stay piecewise-constant?** Declarable, defaults to
  published, and the deviation is measured rather than argued (R11).
- **Is the promoter proxy inherited verbatim at Step 2?** Barely matters for
  synthetic recovery, where it only sets the truth we recover. Leaning toward
  inheriting it as the baseline the surrogate is measured against.
- **Whether the full-cycle checks belong in the default test suite.** They are
  the strongest evidence and the slowest thing added. Deciding needs the measured
  wall-clock. **Shortening the interval is not an option** — it is precisely the
  error of record.

---

## 10. Risks

Each risk is paired with the signal that arrives *before* it costs weeks. Several
arrive in phase 3 or in the module phases rather than at the end, which is the
point of D0's reordering.

| # | Risk | Early warning | Pre-committed response |
|---|---|---|---|
| R1 | The reverse channel is too weak for the bidirectional claim | **Measure the charged-tRNA elasticity the week the translation module lands** (phase 10). Below 0.02 and the reverse direction is bounded at a few percent by the structure of both published rate laws | K2's response: reframe as a measured bound plus the replicate count needed, and promote Step 1b |
| R2 | The non-smooth drain operates in-regime, so the ODE block is not differentiable where it matters | **Any non-zero deficit carried on the first full-cycle run at published parameters** (check 7, phase 3). Do not wait to meet it as sampler divergences, which look like a step-size problem | Smooth the drain; validate the smoothed model against the clipped one. K5 |
| R3 | Integer rounding drift swamps the integrator and is read as a conservation leak | **The adenylate residual failing to shrink when tolerances tighten tenfold.** One extra run distinguishes the two | Fractional-carry rounding (check 0). The magnitudes differ by ~140×, so this is not hypothetical |
| R4 | The reference system is skipped because production looks fine | **A calibration pass reported with no reference comparison beside it.** This is a process risk, so the warning is procedural | F10 is a phase deliverable, not an optional extra. A calibration claim without it is unattributable |
| R5 | The message family's finite support truncates every hand-off | **Acceptance rate collapsing across populations**, or more than 1% of posterior mass within one bandwidth of the density estimate's grid edge. The existing prior returns zero outside its grid, and the sampler's support gate then rejects everything beyond the previous stage's range | A fat-tailed message family or a moment-matched Gaussian site; log the edge-mass fraction every iteration |
| R6 | The six targets are not identifiable as a set | **A rank below six or a condition number above 1e6, computed before any sampling.** Costs one Jacobian | Reparameterise to the product plus an anchor. K6 |
| R7 | Calibration is unaffordable, so the claim cannot be made | **The phase 13 wall-clock exceeding 10 s per trajectory** | K1's ladder, pre-committed so it is not negotiated under pressure |
| R8 | The lumped charging step — ours, not the model's — determines the answer | Two signals. The stoichiometry comparison shifting any target's posterior mean by more than half a posterior standard deviation. And, more likely, **the charged-tRNA pool size appearing in the top three by sensitivity in F2** — a quantity we chose setting the strength of the strongest candidate reverse channel | Report the channel gain as a *function* of the pool size rather than a point value; label the lumping in T2 with its measured posterior sensitivity |
| R9 | Core A′ is better-conditioned than the same subnetwork inside the full model, so success over-claims | **Measurable now, without waiting:** the posterior width on the pyruvate-kinase constant with the seven dropped sibling reactions stubbed in as a competing sink | Quote the factor alongside every recovery claim. K7 |
| R10 | Everything rests on synthetic data, so misspecification is untested by construction | **No internal early warning exists, and that is why this is a stated scope limit rather than a mitigated risk.** The nearest signal is F5: if either external check fails, the model is already misspecified against the data that does exist | State the limit in the abstract. Do not claim robustness that was not tested |
| R11 | The 60 s rebuild is an implementation artefact the posterior depends on | **Nominal trajectory at 60 s against 6 s rebuild cadence; any observed pool differing by more than one percent.** With the GTP pool turning over in half the interval, expect this to fire | A labelled decision in T2 with the posterior sensitivity reported. Turns an open question into a measurement |
| R12 | Protein fold change misses two, so the reduction throttles translation | **The arithmetic, run the week the translation rate law exists.** A first pass puts one gene near 1.3 against a published median near 2 | Diagnose in phase 10 — the charging step or the ribosome constant — rather than explain it in phase 16 |
| R13 | The observation model's scale is wrong and every credible interval is misreported | **No noise-scale recovery test exists anywhere in `test/`.** A variance-versus-standard-deviation misreading of the multivariate normal constructor rescales every interval and passes every other check | Check 9, written before the first posterior |
| R14 | Summary statistics discard the low-copy information the reduction exists to exercise | **The ensemble observer averages replicates and the summaries keep per-time means only** — both true in the code today, both contradicting the scoping note's own recommendation | Carry transcripts as counts per replicate; add distributional summaries and a lag-one autocovariance, or supersede summaries with the exact likelihood |
| R15 | The seven module phases are written against a protocol that then changes again | **A framework need surfacing during a module phase.** That is an interface bug, and the wave plan's rule stands: it is a conversation on the main branch, not an edit on a module branch | Amend §12 with the change and its date, then land it as its own phase before the affected modules continue |

---

## 11. Task list

Seventeen phases plus a phase 0, each one reviewable pull request. Ordering is
D0's: the two protocol changes, then the kill phase on a toy, then the drivers,
then the modules, then assembly, validation and inference.

**Parallelism.** Phase 1 runs alone. Then a four-way fan-out: phases 6, 7 and 8
(the ODE modules) run concurrently with the framework track of phases 2 to 5.
Then phase 9 before phases 10 and 11, because both read the transcripts phase 9
owns; phase 12 is parallel to them. Phases 13 to 17 are strictly serial, because
each validation check catches errors the next would mask.

### Phase 0 — Retire OpenSpec, keep its findings

**Goal:** make this spec the only live spec system without losing a derived
number.
**Done when:** every finding in the four drafted changes appears in this spec or
in a note, and `openspec/` is read-only history under `dev/archive/`.
**PR:** _all six tasks done in the working tree alongside this spec, not yet
committed. `git add -A` is required before committing, because the four drafted
changes were untracked and are otherwise lost._

- [x] 0.1 Move `openspec/` to `dev/archive/openspec/` — verify by `git log
  --follow` on `dev/archive/openspec/specs/corea-interface/edge-kinds/spec.md`
  reaching the wave-0 commit, and by the four currently-untracked changes being
  `git add`ed in the same commit, since untracked files are not in history and
  would otherwise be lost outright.
- [x] 0.2 Diff the absorbed findings against the originals — verify by a checklist
  confirming that §4 D2, D3, D4, D6 and D7 each carry every number from the
  corresponding proposal and design: the eight differing Michaelis constants, the
  five cross-file ratios, the base-count totals and the 1.9× sensitivity, the
  lactate volume-ratio table, and the four proteomics fractions.
- [x] 0.3 Add the correction-of-record rows the four drafted changes each promised
  to `dev/notes/reduced-syn3a-scoping.md` — verify by the note's inline
  correction table gaining a row for the Michaelis-constant column, the
  base-mapping bug, the wider cross-file trap, the carrier phospho-split, and the
  transcription pyrophosphate the note's phosphate accounting omits; and by each
  row's numbers matching this spec's.
- [x] 0.4 Add a superseded header to `dev/plans/reduced-syn3a-wave-plan.md`
  naming this spec and recording that D0 changed its ordering and why — verify by
  the header naming the two protocol changes that move earlier, so a reader of
  the wave plan alone cannot follow the stale order.
- [x] 0.5 Drop the OpenSpec pointers from `CLAUDE.md` and `README.md` and add the
  spec workflow — verify by grepping both for `openspec` and finding only
  archive-path references.
- [x] 0.6 Confirm nothing executable changed — verify by `sbatch
  test/run_tests.slurm` passing at 829, unchanged. **Done:** job 16167788,
  829/829 in 1m12s, and `git status` shows `src/` and `test/` untouched.

### Phase 1 — A module contributes to a state it does not own

**Goal:** make an outbound mass or currency edge execute, so a shared pool can
receive terms from every module that produces it.
**Done when:** a two-module ODE composition in which module A produces a species
module B owns shows B's state rising at A's declared rate, and every outbound
mass or currency edge in the composition maps one-to-one onto an executed
contribution.
**PR:** _not started_

- [ ] 1.1 Choose and record the mechanism — a second return value from
  `dynamics`, a separate protocol function, or a global-index write buffer —
  verify by the pull request recording what each costs at 32 states and by the
  chosen one adding no dependency that is absent from the Julia depot.
- [ ] 1.2 Extend the protocol with the chosen channel, defaulting to empty —
  verify by `build_problem(TranscriptionTranslation())` producing byte-identical
  initial conditions, parameters and trajectory to the pre-change run, and by all
  829 existing tests passing unchanged.
- [ ] 1.3 Accumulate contributions into the global derivative by
  registry-resolved index in `_build_rhs` — verify by a two-module test asserting
  the owner's derivative equals its own term plus the contributor's exactly, and
  that the contributor's own slice is unchanged.
- [ ] 1.4 Reject a contribution to a chemostatted species, to a species no module
  owns, and to a name outside the registry — verify by three tests each throwing
  and naming both the species and the contributing module.
- [ ] 1.5 Cross-check contributions against declared edges in both directions —
  verify by a test that deleting a contribution while keeping its outbound edge
  throws naming both, and that adding a contribution with no edge throws
  likewise; the drift check must be live, as the untyped-inputs check already is.
- [ ] 1.6 Re-take the state-container decision at 32 states — the problem is
  currently built out-of-place over a static vector — verify by benchmarking
  static-out-of-place against mutable-in-place at 32 states and asserting the
  chosen path allocates nothing per right-hand-side call.
- [ ] 1.7 Run the suite and hand off — verify by `sbatch test/run_tests.slurm`
  passing above 829 by the tests added, and by `docs/handoff.md` recording that
  outbound edges now execute rather than merely resolve.

### Phase 2 — Several jump modules compose

**Goal:** make `_build_jump_problem` honour the index contexts it already
computes, so a jump module writes only its own global slice and can read and
write a declared peer's state.
**Done when:** two distinct jump modules — one owning a species, the other
incrementing and decrementing it — compose and produce a trajectory whose
statistics match a hand-written single-module equivalent within Monte Carlo
error.
**PR:** _not started_

- [ ] 2.1 Pin the bug before fixing it — verify by a test composing two doubles
  with disjoint state names and asserting that *today* the second module's effect
  mutates the first module's state, so the fix has a witness rather than a claim.
- [ ] 2.2 Give `reactions` the module's index context — verify by instrumenting
  each double's rate and effect functions and asserting every read and write
  falls inside that module's own state range.
- [ ] 2.3 Use the contexts in `_build_jump_problem`, which currently computes and
  discards them — verify by asserting the global initial-condition vector's
  length equals the sum of per-module state counts and that each module's counts
  sit at its own offsets.
- [ ] 2.4 Let a jump module read a declared peer's state in a propensity — verify
  by a propensity proportional to the peer's count and a test where doubling the
  peer's initial count doubles the measured firing rate over many replicates.
- [ ] 2.5 Let a jump module write a declared peer's state, gated on an outbound
  edge — verify by a decay-shaped double decrementing a transcript it does not
  own, and by the same write throwing once the edge is removed.
- [ ] 2.6 Preserve single-module behaviour exactly — verify by the existing
  stochastic and bursty gene-expression tests passing unchanged including the
  Fano-factor assertion, and by rewriting the "composition fails" test to assert
  what now actually fails (duplicate ownership) rather than what no longer does.
- [ ] 2.7 Suite and handoff — verify by `sbatch test/run_tests.slurm` passing and
  by the handoff recording that the four gene-expression modules may now assume
  composability, correcting the drafted design that assumed it already held.

### Phase 3 — The 1 s handshake, on a two-module toy. **The kill phase.**

**Goal:** prove a jump block and an ODE block can be advanced by a 1 s
split-operator exchange — counts to initial conditions, integrate, counts back,
plus the deferred debits — on the smallest pair that exercises all of it, and
measure what it costs.
**Done when:** a toy of one jump gene-expression module and one ODE metabolite
module runs 600 s of handshake in which the protein count sets the ODE rate law
and the ODE pool pays a deferred cost under the published clamped policy, and the
wall-clock per simulated second is recorded from a compute node.
**Kill criterion, stated in the phase:** if the handshake cannot be expressed
inside the sub-model protocol without a per-module special case, or if the
measured cost extrapolates above 10 s per 6,300 s trajectory, the project stops
here and the architecture is revisited rather than scaled. See K1.
**PR:** _not started_

- [ ] 3.1 Choose and record the composition mechanism — an outer split-operator
  loop over stepped integrators, a discrete callback on a single problem, or a
  jump process over an ODE problem — verify by the pull request recording, for
  each, whether the published piecewise-constant 1 s handshake is faithfully
  representable, and by the chosen one adding no dependency absent from the
  depot.
- [ ] 3.2 Replace the mixed-formalism refusal at `src/orchestrator.jl:27` with a
  hybrid build path — verify by a mixed composition returning a driver rather
  than throwing, and by the same composition *still* throwing with a named error
  when a required exchange declaration is absent.
- [ ] 3.3 Implement the count-to-concentration conversion through the registry at
  fixed volume with an explicit rounding policy — verify by check 0: exact-zero
  round-trip residual under fractional carry, and by the policy being a stated
  field on the driver rather than implicit in the arithmetic.
- [ ] 3.4 Isolate exchange error from integration error — verify by a toy whose
  ODE has zero derivative: after 600 handshakes every count is unchanged
  *exactly*, so any drift is attributable to the exchange alone.
- [ ] 3.5 Debit the deferred counters under the published clamped policy — verify
  by a test where the cost exceeds the pool: the pool floors at zero, the carried
  deficit equals the shortfall exactly, and the next step debits it; and by the
  existing gradient report naming the edge when a differentiable sub-model is
  present.
- [ ] 3.6 Execute the enzyme-concentration channel — protein counts overwriting
  the ODE's enzyme concentration through a catalytic edge's parameter slot —
  verify by a test that a step change in the toy's protein count changes the ODE
  flux by the ratio of counts and changes nothing else.
- [ ] 3.7 Run check 7's clipping census on the toy — verify by zero carried
  deficits at nominal parameters and by the fraction of 200 prior draws that clip
  being recorded, since K5 fires here or nowhere.
- [ ] 3.8 Measure and record wall-clock — verify by a Slurm run reporting seconds
  per simulated second and the 6,300 s extrapolation, an explicit pass or fail
  against the phase's stated budget, and the number written into `docs/handoff.md`
  and into the scoping note's open-questions list, which carries it as unmeasured
  today.
- [ ] 3.9 Suite, count and verdict — verify by `sbatch test/run_tests.slurm`
  passing and by the kill-criterion verdict written down as a sentence rather
  than implied by the pull request merging.

### Phase 4 — The 60 s rebuild

**Goal:** execute the only reverse channel — live pools to recomputed jump rate
constants, piecewise-constant at the interval the edge declares.
**Done when:** the toy's propensities change only at interval boundaries, the
interval is read from the declared edge rather than hard-coded, and the measured
elasticity of a rate constant to its upstream pool is reported.
**PR:** _not started_

- [ ] 4.1 Decide where mutable rate constants live — in the flat parameter vector
  the builder already returns as mutable, or as mutable state on the sub-model as
  the drafted transcription design chose — verify by the pull request recording
  the consequence for the parameter-substitution call every inference path uses,
  since a constant held on the struct is invisible to it.
- [ ] 4.2 Schedule the rebuild from the edge's declared interval — verify by a
  test that a 60 s edge refreshes ten times in 600 s and a 30 s edge twenty
  times, with the count read from the trajectory rather than asserted.
- [ ] 4.3 Assert piecewise-constancy — verify by sampling propensities between
  refreshes and asserting they are bitwise unchanged, so the coupling is the
  published piecewise-constant one and not an accidental continuous one.
- [ ] 4.4 Handle a continuous cadence explicitly — verify by either implementing
  it and asserting `reduction_declarations` labels it a deviation, or rejecting it
  with a named error; silently discarding the cadence field is the failure to
  avoid.
- [ ] 4.5 Measure the channel's gain on the toy — verify by an elasticity
  diagnostic and by the number recorded where phase 9 compares it against the
  0.044–0.051 the real transcription module predicts.
- [ ] 4.6 Run D10's granularity comparison in miniature — verify by nominal
  trajectories at 1 s, 5 s and 60 s drain granularity with the largest relative
  difference in any pool reported, so phase 13 inherits a measurement rather than
  an assumption.
- [ ] 4.7 Suite and handoff — verify by `sbatch test/run_tests.slurm` passing and
  by the handoff recording that bidirectional coupling is now executed rather
  than declared.

### Phase 5 — Growth and volume

**Goal:** make the count-to-concentration conversion volume-aware, so growth
dilutes the ODE block as the published model's does.
**Done when:** a volume edge on the toy's membrane-protein state drives radius
and volume, both conversion directions read the live volume, and the published
initial surface area returns a 200 nm radius exactly.
**PR:** _not started_

- [ ] 5.1 Promote a membrane-protein declaration to the protocol with an empty
  default, settling the open question the drafted transport design leaves for
  this phase — verify by a composition sweep finding the flagged states without
  naming any module type, and by a module declaring none composing unchanged.
- [ ] 5.2 Implement the surface-area, radius and volume chain — verify by
  asserting the published initial 502,831 nm² returns exactly 200.0 nm, and by
  the 28.0 nm² footprint carrying its calibrated-not-measured label, with the
  upstream docstring's contradictory 35 nm² recorded and the code's value taken.
- [ ] 5.3 Route both conversion directions through live volume — verify by a test
  that a volume increase dilutes every ODE concentration by exactly the volume
  ratio with no change in any count.
- [ ] 5.4 Cap growth at exactly twice initial volume as the published model does —
  verify by a test driving counts past the cap and asserting volume stops rather
  than growing.
- [ ] 5.5 Reproduce the scoping note's ptsG arithmetic — verify by a test
  asserting that doubling 831 copies moves surface area from 502,831 to 526,099
  nm², radius from 200.0 to 204.6 nm, and volume to about 1.07×.
- [ ] 5.6 Refuse the doubling-time comparison in code — verify by the composed
  model exposing the reporting constraint that this is fractional growth or
  time-to-threshold, and by a test asserting the constraint is retrievable rather
  than living only in prose.
- [ ] 5.7 Suite and handoff — verify by `sbatch test/run_tests.slurm` passing.

### Phase 6 — Central glycolysis

**Goal:** the ten reactions from glucose-6-phosphate through lactate, owning the
eleven glycolytic intermediates and the two redox species.
**Done when:** the thirteen owned states integrate non-negatively over a full
cycle from the registry's initial conditions at published parameters, the redox
pair is conserved to a bound that shrinks with the solver tolerance, and a
mutation to the dehydrogenase's stoichiometry makes that check fail.
**PR:** _not started_

The drafted task list at
`dev/archive/openspec/changes/add-central-glycolysis/tasks.md` holds the
fine-grained sub-checks and the spot values to assert; it is the reference for
this phase's detail. Three amendments, from D0: the outbound currency edges now
**execute**, so add an assertion that each moves the mass its rate law says; drop
the "declared but unexecuted" clause from its handoff note; and use the real
recycling module in place of its held-energy double wherever phase 8 has landed,
keeping the double only for standalone runs.

- [ ] 6.1 Vendor the derived extract of ~65 rows and its regeneration script —
  verify by the row classes and counts matching (10 forward and 10 reverse
  catalytic constants, 32 Michaelis constants, 13 concentrations), by four spot
  values matching upstream byte for byte, and by the README's command
  round-tripping to an unchanged file.
- [ ] 6.2 Implement the sub-model, its thirteen states and the modular rate law
  generic over substrate and product counts — verify by the state set equalling
  the registry's glycolytic and redox groups with strictly increasing indices, and
  by unit tests that equal forward and reverse terms give zero net rate and that
  the per-reaction term counts sum to 32.
- [ ] 6.3 Import every value through the loader with log-normal priors from mode
  and geometric standard deviation, fixed by default — verify by every parameter
  reporting the central file as its source, by the two prior-default species
  keeping their width, and by the registry-agreement check firing when a
  concentration row is mutated in a temporary copy.
- [ ] 6.4 Set the ten enzyme concentrations from copy number at the registry's
  volume, marked nominal and overridable, and declare the ten protein counts as
  inputs so translation later supersedes them — verify by each value equalling
  copies over 20,180 to six decimals, by none being the published no-rule default
  of 0.001 mM, and by overriding one scaling exactly the rates that enzyme
  catalyses.
- [ ] 6.5 Declare the boundary: currency edges on the energy species and mass
  edges on the shared intermediates, every peer unnamed — verify by the edge
  count, kinds and directions matching, by no edge naming a redox species, and by
  standalone resolution succeeding while listing the energy species as unowned.
- [ ] 6.6 Register the reduction notes: the dropped oxidase with its three
  reasons, the Michaelis-constant column choice naming all eight differing
  constants with both values, the held currencies, and the nominal enzyme
  concentrations — verify by `reduction_declarations` returning all four as
  sentences naming the affected reactions.
- [ ] 6.7 Integrate and check redox balance — verify by non-negativity at every
  save point, by the redox residual satisfying the tolerance principle across two
  tolerance settings, and by the mutation test failing and naming the conserved
  pool and the drift.

### Phase 7 — Phosphotransferase transport and lactate export

**Goal:** the five-step cascade that imports glucose against phosphoenolpyruvate,
plus passive lactate export, owning the eight carrier phospho-states and external
lactate.
**Done when:** the four carrier sums are independently conserved to
tolerance-derived bounds, each carrier's two forms sum to its published copy
number exactly, and the export rate matches the permeability law at the
registry's radius.
**PR:** _not started_

Reference detail at
`dev/archive/openspec/changes/add-pts-transport/tasks.md`. One amendment: the
membrane-protein declaration is the protocol function phase 5 promoted, not a
plain function on this module, which settles that design's own open question.

- [ ] 7.1 Vendor two extracts — eleven rate constants with no uncertainty column
  and nine initial conditions with their derivation — verify by every constant
  loading with informedness `asserted` and appearing in the composed model's
  asserted-prior enumeration, and by the balanced glycolytic parameters not
  appearing there.
- [ ] 7.2 Implement the cascade as mass-action forward and reverse pairs and the
  export law with permeability and radius carried separately — verify by a
  hand-checked derivative at one state vector, and by the export rate equalling
  the permeability law rather than a folded constant, so phase 5's growing radius
  changes it.
- [ ] 7.3 Set the eight initial conditions from copy number times proteomics
  fraction (D7) — verify by each carrier's forms summing to 353, 314, 290 and 831
  copies exactly, and by both the totals and the split recorded as the published
  model's rather than asserted.
- [ ] 7.4 Integrate external lactate with the volume ratio (D6) — verify by
  external lactate after a full cycle staying below one percent of steady
  cytosolic lactate at the default ratio, by a ratio of one visibly saturating
  export, and by the ratio registered as ours.
- [ ] 7.5 Declare the boundary, including external glucose clamped at 40 mM with
  published origin — verify by that clamp *not* appearing in the
  what-is-ours enumeration while the chemostatted nucleotide and amino-acid pools
  do, since the published model clamps external species too.
- [ ] 7.6 Flag ptsG through phase 5's protocol function with its 28.0 nm²
  footprint — verify by a composition sweep finding both phospho-forms without
  naming this module's type.
- [ ] 7.7 Check carrier conservation as four separate assertions — verify by four
  independent bounds each naming its own carrier, by a mutation to one cascade
  step failing that carrier's check while the other three still pass, and by the
  bounds satisfying the tolerance principle.

### Phase 8 — Nucleotide recycling

**Goal:** the five reactions that make GTP and close the adenylate, guanylate and
phosphate moieties, owning the four adenylate, three guanylate and one
pyrophosphate species — the pools every other module routes energy through.
**Done when:** over a full 6,300 s cycle against a charging drain, adenylate and
guanylate are each conserved and ATP stays positive; removing the adenylate
kinase exhausts the pool at about 141 s and removing the pyrophosphatase leaves
pyrophosphate unbounded.
**PR:** _not started_

Reference detail at
`dev/archive/openspec/changes/add-nucleotide-recycling/tasks.md`. Its charging
drain double is the acceptance criterion and not scaffolding; it gains one task,
recording its stoichiometry so phase 12 can assert the real module reproduces it.

- [ ] 8.1 Vendor two extracts — the governing nucleotide file and the central
  file's rival values for the same identifiers — verify by every one of the 35
  identifiers being held by both, so a truncated rival file fails loudly rather
  than quietly disarming the ambiguity check.
- [ ] 8.2 Declare a governing file for all five reactions and all eight initial
  conditions (D2) — verify by loading failing when a governing declaration is
  removed, by `governing_choices` returning each value with the file chosen and
  the file rejected, and by the pyrophosphatase constant being the nucleotide
  file's 646.73 and not the central file's 583,611.
- [ ] 8.3 Implement the five reactions with repeated stoichiometry expanded as the
  published builder does — verify by the two reactions with a doubled product
  each naming that species twice in the rate law, giving 17 Michaelis constants
  across five reactions rather than 19, and by the pyrophosphatase carried
  reversibly.
- [ ] 8.4 Set the three new enzyme concentrations, and assert the two shared with
  glycolysis agree — verify by a cross-module test asserting one enzyme is not
  run at two concentrations, since two modules carry a nominal for each of the
  two shared genes.
- [ ] 8.5 Declare the boundary: mass edges for species not owned, currency edges
  only where this module is the sole producer or sole consumer of a pool it owns —
  verify by the edge count and kinds matching, and by a test composing all three
  ODE modules' declarations asserting no species is described as both mass and
  currency in one direction.
- [ ] 8.6 Build the charging drain double and run three full-cycle configurations
  — verify by the all-reactions run conserving both moieties with ATP positive and
  pyrophosphate settling bounded, the kinase-removed run exhausting adenylate
  within an order of magnitude of 141 s, and the pyrophosphatase-removed run
  rising without bound.
- [ ] 8.7 Check phosphate closure in both forms — verify by exact invariance with
  the GTP-branch reactions inactive, and by the flux-corrected form with them
  active *subtracting* the inbound flux rather than relaxing the bound.
- [ ] 8.8 Record the drain's stoichiometry for phase 12 — verify by the recorded
  rate being 553.1 per second and by phase 12's acceptance test citing it.

### Phase 9 — Transcription

**Goal:** seventeen genes, one constitutive transcription reaction each, with the
corrected base mapping, the five cost counters, and rate constants that are
recomputable but never self-recomputed.
**Done when:** seventeen transcripts simulate over a full cycle as non-negative
integers, each gene's time-averaged count is within a factor of two of its
measured mean, the seventeen rate constants fall in 1.26e-3 to 8.29e-3 per second,
and the GTP elasticity reproduces 0.0079 to 0.0117 under the corrected mapping.
**PR:** _not started_

Reference detail at
`dev/archive/openspec/changes/add-corea-transcription/tasks.md`. One amendment:
its assertion that simulating past the declared cadence leaves the constants
unchanged becomes "phase 4's driver refreshes them, and the module still never
refreshes itself".

- [ ] 9.1 Vendor the per-gene extract from three upstream sources — verify by
  seventeen rows whose four base counts sum to the transcript length, by the
  totals matching A 7236, C 2078, G 3094, U 5868, and by the generator failing
  loudly on a missing locus rather than defaulting, since a silently absent gene
  shows up only as a model with sixteen transcripts.
- [ ] 9.2 Implement the seventeen reactions, genes carried as fixed quantities
  rather than states — verify by the state count being 17 transcripts plus 5
  counters, by none being a registry species, and by a recorded note that adding
  replication later must promote the genes to states.
- [ ] 9.3 Implement the rate constant with the corrected base mapping as default
  and the published permutation behind a keyword (D4) — verify by the correction
  appearing in `reduction_notes` while the published mapping does not, by the
  turnover cap never binding (the largest is 8.85 against a ceiling of 180), and
  by both mappings computable so the difference is one argument away.
- [ ] 9.4 Take all four nucleotide concentrations from balanced tables, never the
  first-minute setup constants — verify by the four values and their widths
  matching the balanced files, and by a recorded note that the nucleotide file
  repeats two of the setup constants, which makes the wrong choice look
  corroborated.
- [ ] 9.5 Expose one recomputation entry point and confirm the module never calls
  it — verify by simulating past the declared cadence with the constants
  unchanged, and by phase 4's driver changing them when it is present.
- [ ] 9.6 Declare the five deferred counters and the two clamped nucleotide pools
  with origin ours — verify by the counter table matching what each drains into,
  by all five taking the published clamped policy so none is a labelled
  deviation, and by three kinds coexisting on ATP inbound as the contract requires.
- [ ] 9.7 Register both promoter-proxy declarations (D5) — verify by
  `reduction_declarations` returning one entry saying it is a proxy and a second
  naming the circularity and the seventeen parameters affected, and by a test
  asserting this module's copy numbers agree with the metabolic modules', since
  one number has three consumers.
- [ ] 9.8 Simulate a full cycle and report the elasticities — verify by
  non-negative integer counts throughout, by each gene's mean within a factor of
  two of measured, by the elasticity to all four pools falling in 0.044 to 0.051,
  and by the GTP-alone elasticity being reported under both mappings so the 1.9×
  effect is measured rather than argued.

### Phase 10 — Translation

**Goal:** one translation reaction per transcript plus ptsG translocation, with
the rate constant built from the live lumped charged-tRNA pool, publishing the
seventeen protein counts the ODE block and the growth channel read.
**Done when:** seventeen protein counts simulate over a full cycle as
non-negative integers, the energy counter equals exactly twice the residues
translated, and the counts feed the metabolic modules' enzyme concentrations
through executed catalytic edges rather than nominal stand-ins.
**PR:** _not started_

- [ ] 10.1 Extract per-gene amino-acid counts and residue totals for the
  seventeen loci from the genome record — verify by the residue counts summing to
  3,484,518 for a full proteome doubling, matching the scoping note's own figure,
  and by four spot lengths matching the recorded 746, 574, 155 and 90 residues.
- [ ] 10.2 Implement seventeen jumps catalytic in the transcript count, reading
  transcripts as a phase-2 peer state — verify by firing one gene's reaction
  raising only that gene's protein by one, leaving the transcript unchanged, and
  raising the energy counter by exactly twice that gene's residue count.
- [ ] 10.3 Implement the rate constant with the lumped pool substituted for the
  twenty per-amino-acid pools — verify by the lumping registered as ours, and by
  the three polymerase-capacity constants carried with the note that all are
  computed at initial volume and are among the genuinely frozen quantities.
- [ ] 10.4 Add the ptsG translocation reaction — verify by only that locus
  carrying it, by translocation being what increments the state phase 5's volume
  edge reads, and by no cytosolic protein having it.
- [ ] 10.5 Declare the boundary: an inbound rate-constant edge on the charged
  pool at the declared cadence, a deferred counter on GTP under the published
  policy, an outbound mass edge returning uncharged tRNA, seventeen outbound
  catalytic edges naming each rate law's parameter slot, and a volume edge on
  ptsG — verify by the exact edge count, kinds and directions, and by the
  catalytic edges carrying no mass, with an attempt to include one in a
  conservation check rejected rather than counted as zero.
- [ ] 10.6 **This module's own check: residue-to-energy closure.** Assert the
  accumulated counter equals exactly twice the residues translated at every write
  point — verify by the assertion passing and by a mutation halving one gene's
  residue count failing and naming that gene.
- [ ] 10.7 Add protein degradation at the published rate — verify by the count
  decaying with the expected mean lifetime absent translation, and by degradation
  declared as returning nothing to the chemostatted amino-acid pool with the
  exemption recorded, so a later change making that pool live reinstates the
  check.
- [ ] 10.8 **Settle two open questions from §9 here.** Measure the elasticity of
  the translation rate constant to the charged pool, and the protein fold change
  over a full cycle — verify by both numbers recorded against the predictions
  (a reverse channel bounded at a few percent, and a fold change that may fall
  short of two), and by either shortfall diagnosed to the charging step or the
  ribosome constant *in this phase*, since R1 and R12 both say diagnosing here is
  far cheaper than explaining at the end.
- [ ] 10.9 Suite and handoff — verify by `sbatch test/run_tests.slurm` passing and
  by the handoff recording that the metabolic modules' nominal enzyme
  concentrations are now superseded by live counts.

### Phase 11 — Transcript decay

**Goal:** one decay reaction per transcript at the published length-dependent
global rate, returning the four monomers to the pools that can accept them.
**Done when:** each decay returns exactly its own base counts, the two
recyclable monomers reach the pools the recycling module closes, the two
chemostatted monomers are exempt with the exemption recorded, and the guanylate
return over a cycle matches the scoping note's ~24,000 particles.
**PR:** _not started_

- [ ] 11.1 Implement seventeen decay jumps with the published single global
  constant over transcript length — verify by exactly one rate parameter across
  all seventeen genes, by each half-life being a function of length alone, and by
  the upstream hand-tuning comment recorded rather than silently inherited.
- [ ] 11.2 Decrement the transcript this module does not own, through phase 2's
  declared peer write — verify by firing decay lowering only that transcript, by
  it being unable to fire at zero copies, and by the write throwing once the
  outbound edge is removed.
- [ ] 11.3 Return the four monomers per firing — verify by the returned counts
  equalling that gene's four base counts exactly, cross-checked against phase 9's
  extract so the two modules cannot disagree on a base count.
- [ ] 11.4 Declare outbound counters on the two recyclable monomers and take the
  chemostat exemption on the other two — verify by `resolve_coupling` reporting
  the exempt pair among its chemostat exemptions with the exemption recorded.
- [ ] 11.5 Add the decay energy cost counter the published hook drains — verify
  by it accruing per firing and by the composed model reporting which registry
  species it debits.
- [ ] 11.6 **This module's own check: monomer closure.** Assert that over a
  composed transcription-and-decay run the total monomers returned equals the
  total bases polymerised, per moiety — verify by the assertion passing and by a
  mutation to one gene's returned counts failing and naming the moiety.
- [ ] 11.7 Record the leak this module closes — verify by the composed model
  reporting the guanylate return over a cycle against the guanylate pool, about
  60%, as the recorded reason the guanylate kinase is in the core.
- [ ] 11.8 Suite and handoff — verify by `sbatch test/run_tests.slurm` passing.

### Phase 12 — Lumped tRNA charging

**Goal:** the one reaction that is ours rather than the published model's, owning
the charged and uncharged tRNA pools, labelled everywhere it can reach a result.
**Done when:** the tRNA pair is conserved over a full cycle, the charging flux
reproduces the published demand, phase 8's drain double is retired against the
real module, and `reduction_declarations` reports the lumping.
**PR:** _not started_

- [ ] 12.1 Implement the single reaction owning both tRNA states — verify by the
  state set equalling the registry's tRNA group with strictly increasing indices.
- [ ] 12.2 Assert the initial total pool, which the registry records as not
  imported for both forms — verify by the value loading with informedness
  `asserted`, appearing in the composed model's asserted-prior enumeration, and
  carrying the note that the published chain's rates are hard-coded inline with
  no priors to inherit.
- [ ] 12.3 Declare the three currency edges — verify by the kinds and directions
  matching, and by composing with phase 8 asserting no species is described as
  both mass and currency in one direction.
- [ ] 12.4 Register the lumping — verify by `reduction_declarations` returning it
  under `lumping` as a sentence naming what it replaces (twenty chains of five
  reactions) and why the published product stoichiometry was preferred over a
  two-ATP form that would close the same moieties by fiat.
- [ ] 12.5 **This module's own check: tRNA conservation.** Assert the pair's sum
  equals its initial value at every save point over a full cycle — verify by the
  check satisfying the tolerance principle, and by a mutation that creates tRNA
  rather than transferring it making it fail.
- [ ] 12.6 Retire phase 8's double — verify by composing the real module with the
  recycling module and asserting the charging flux matches the double's 553.1 per
  second within the tolerance the double used, so the stand-in is shown to have
  stood in correctly.
- [ ] 12.7 **Settle §9's highest-value open question here.** Determine whether the
  stochastic block's likelihood is closed-form given this reaction's two
  reactants — verify by a written determination with its consequence for the
  inference budget stated, and by phase 16's method choice citing it. If it is
  not closed-form, record the fallback (a particle filter) as work phase 16
  inherits.
- [ ] 12.8 Answer the stoichiometry question — verify by a comparison run under
  both lumpings recording the kinase traffic (553 per second against zero) and
  the resulting ATP-to-ADP ratio the rebuild reads, and by the scoping note's
  open-questions entry being replaced with the measured answer.
- [ ] 12.9 Suite and handoff — verify by `sbatch test/run_tests.slurm` passing.

### Phase 13 — Assemble Core A′ and assert structural completeness

**Goal:** compose all seven modules into one hybrid model and make the assertion
the wave-0 contract deliberately withholds.
**Done when:** every registry dynamic state is owned, no dead ends remain, every
declared edge is executed, the model runs a full 6,300 s cycle, and the
per-trajectory wall-clock is recorded.
**PR:** _not started_

- [ ] 13.1 Add a completeness mode — every dynamic state owned and no dead ends,
  as a *failure* rather than a report — verify by the assembled composition
  passing it and by a composition with one module removed failing it and naming
  the unowned states and stranded moieties. The contract is explicit that a
  successful resolve is not evidence of closure, so **this assertion must be
  shown to fail** or it is not an assertion.
- [ ] 13.2 Assert every declared edge is executed — verify by a test
  cross-checking the resolved graph against the driver's registered
  contributions, debits, rebuilds and volume reads, so a declared-but-inert edge
  fails the build rather than passing silently.
- [ ] 13.3 Fix inference dispatch for a hybrid composition — verify by a test that
  a mixed composition either dispatches on the composition or throws a named
  error, since dispatch currently reads only the first module in the vector and
  would silently take whichever path that module declares.
- [ ] 13.4 Settle the drain granularity with D10's measurement at full scale —
  verify by nominal trajectories at 1 s, 5 s and 60 s granularity, by the largest
  relative difference in any observed pool reported, and by the chosen
  granularity registered as a labelled reduction with its *measured* cost.
- [ ] 13.5 Run a full cycle at published parameters with a stiff solver and pinned
  tolerances — verify by a trajectory over the full interval with no solver
  failure, and by the solver and tolerances recorded on the model rather than in
  a test file.
- [ ] 13.6 Enumerate what is ours across the whole composition — verify by
  `reduction_report` returning the lumped charging step, the chemostatted
  nucleotide and amino-acid pools, the corrected base mapping, the eleven asserted
  transport priors, the Michaelis-constant column choice, the lactate volume
  ratio, the rounding policy, the drain granularity and exogenous membrane
  growth, and by the count matching the phases that registered them.
- [ ] 13.7 Measure and record per-trajectory wall-clock at full scale — verify by
  a Slurm run reporting seconds per cycle, comparison against phase 3's
  extrapolation, and an explicit statement of how many trajectories the inference
  budget affords. **K1 and K7 are decided by this number.**
- [ ] 13.8 Suite and handoff — verify by `sbatch test/run_tests.slurm` passing,
  with the full-cycle run marked and Slurm-gated if it is too slow for the default
  suite. Shortening the interval is forbidden; it is precisely the error of
  record.

### Phase 14 — The validation checks, in the scoping note's order

**Goal:** run the balance checks on the assembled model in the order in which
each catches errors the next would mask.
**Done when:** checks 0 through 8 all pass over a full cycle, and each has a
mutation test showing it can fail.
**PR:** _not started_

- [ ] 14.1 Re-assert check 0 on the assembled model — verify by the round-trip
  residual showing the signature its policy predicts across 6,300 handshakes, so
  no downstream residual can be a rounding artefact.
- [ ] 14.2 Check 1 and 1b — verify by non-negativity naming the first state and
  time of any violation rather than reporting a global failure, and by the
  particle-floor report flagging every state below 500 particles and
  cross-checking the three smallest against a chemical-Langevin ensemble, with
  any observable outside the band excluded from the likelihood.
- [ ] 14.3 Check 2, carbon balance — verify by closure satisfying the tolerance
  principle, by the homolactic ratio equalling 2.000 to integrator tolerance
  since it is analytically exact, and by a mutation removing lactate export
  failing with cytosolic lactate heading for the recorded ~691 mM.
- [ ] 14.4 Check 3, redox — verify by the module-local check from phase 6
  restated at composition scope rather than reimplemented, so there is one
  implementation and one bound.
- [ ] 14.5 Check 4, adenylate and guanylate **over a full cycle** — verify by both
  moieties closing over 6,300 s, by each being asserted separately so a failure
  names the moiety, and by an explicit assertion that the interval is a full
  cycle, since 144 s of it looked fine before the dead end was found.
- [ ] 14.6 Check 4b, phosphate closure — verify by closure after accounting for
  transcription's pyrophosphate, which the scoping note's own accounting omits,
  and by the flux correction being a subtraction rather than a relaxed bound.
- [ ] 14.7 Check 5, carrier conservation — verify by four independent bounds
  naming the carrier that drifts, restated as conserved up to what translation
  adds, again by subtraction.
- [ ] 14.8 Check 7 and check 8 — verify by the clipping census reporting zero
  carried deficits at published parameters and the clipping fraction across 200
  prior draws (K5), and by the metabolic control analysis identities holding
  within 1e-6, which is the one analytic result the implementation must reproduce.
- [ ] 14.9 Check 6, the nominal trajectory — verify by the trajectory produced and
  reported as fractional growth or time-to-threshold, with the doubling-time
  comparison refused in code as phase 5 established, not merely in prose.
- [ ] 14.10 Show each check can fail — verify by one mutation test per check, each
  failing the intended check and naming the intended quantity, collected into the
  mutation table of output F3.

### Phase 15 — Observables and synthetic data

**Goal:** choose what to condition on, generate synthetic data at known
parameters, and close the promoter-proxy circularity by construction rather than
by memory.
**Done when:** a synthetic dataset exists at recorded ground-truth parameters
drawn from the prior, the observable set is justified against D8, and any
observable that closes the proxy loop is rejected by a check.
**PR:** _not started_

- [ ] 15.1 Implement the circularity guard — verify by an observable set
  including protein counts, or freeing initial protein counts, being rejected and
  naming the seventeen promoter parameters it would double-count, and by
  metabolite time courses and transcript counts passing.
- [ ] 15.2 Implement the truth-drawing rule from D8 — verify by a coverage or
  calibration run refusing a truth set to the nominal published values, and by a
  fixed-truth-at-nominal run being permitted only when labelled a smoke test.
- [ ] 15.3 Emit the candidate observables at 60 s — verify by a trajectory writing
  per-species counts, per-reaction fluxes and volume at the published cadence, and
  by transcripts carried as a replicates-by-genes-by-times array of counts rather
  than a replicate mean, which the current ensemble observer would produce.
- [ ] 15.4 Replace the single scalar observation noise with per-modality noise on a
  log scale — verify by check 9 recovering a known noise scale inside the
  calibration band, and by a test that the smallest pools now contribute to the
  likelihood, which under one additive term they did not.
- [ ] 15.5 Compute output F2, the concentration control coefficient matrix, with
  check 8 as its guard — verify by the summation identities holding and by the
  metabolite panel being chosen from its largest rows rather than by argument,
  which settles the scoping note's open question.
- [ ] 15.6 Run the influence audit — verify by a report giving, per candidate
  observable and per target parameter, how much the observable moves it, so the
  data model is chosen by measurement.
- [ ] 15.7 Fix the target set and generate the dataset — verify by the six
  parameters of D11 being free with everything else fixed, by the ensemble sized
  at 200 cells with the derivation recorded, by two runs at one seed being
  identical, and by the stored dataset carrying its truth, its seed and the
  `reduction_report` of the model that produced it.
- [ ] 15.8 Check the target set is identifiable before sampling — verify by
  `check_identifiability` reporting full rank and a condition number below 1e6
  for the six, and by the same check on the seven-parameter set including the
  polymerase constant showing the predicted ridge (K6, output F13).

### Phase 16 — Recovery, coverage, and the reference

**Goal:** infer the target set, check the posterior covers truth at nominal
rates, and validate the production sampler against an exact reference.
**Done when:** coverage over repeated synthetic datasets is within Monte Carlo
error of nominal for every target, the tight-prior control shrinks and recovers,
the loose-prior control is visibly wider, and the reference comparison exists.
**PR:** _not started_

- [ ] 16.1 Build M0, the two-gene reference, and run it long — verify by the
  blocked sampler targeting the joint posterior on a system small enough that a
  very long run is defensible, and by the run's length justified rather than
  chosen. **This is a prerequisite, not an appendix:** a calibration pass without
  it is unattributable (D10, R4).
- [ ] 16.2 Choose the production route and record it — verify by the choice
  justified against phase 12's determination on the closed-form likelihood, and
  against what the gradient report says about the assembled model's non-smooth
  boundary edges.
- [ ] 16.3 Recover the tight-prior control alone — verify by the posterior
  concentrating on truth well inside its prior, by a deliberate perturbation of
  the truth moving the posterior with it, and by shrinkage below 0.9. **If the
  control does not recover, the machinery is wrong and this phase stops** (K4).
- [ ] 16.4 Recover the full six-parameter set — verify by every marginal
  containing truth and by the run's cost reported against phase 13's measured
  per-trajectory budget.
- [ ] 16.5 Compute coverage over repeated datasets — verify by nominal 50, 80, 90,
  95 and 99 percent intervals covering truth at those rates within binomial
  error, per parameter, with the replicate count stated, and by the tight-prior
  control falling inside K4's band.
- [ ] 16.6 Produce the shrinkage table, output F6 — verify by all six targets
  reported under transcripts-only, metabolites-only and joint data, and by the two
  load-bearing cells present: the ptsG promoter under metabolites only and the
  tight-prior control under transcripts only. **These two cells are what "crosses
  the boundary" means operationally.**
- [ ] 16.7 Compare against the reference and against the cut — verify by output
  F10 reporting the divergence between the production and reference posteriors on
  the two-gene system, and by output F11 running the same targets through the
  joint posterior and the explicit cut with both calibrations reported.
- [ ] 16.8 Decide K2 — verify by the prior-to-posterior divergence under
  transcript-only data reported for every metabolic parameter against the 0.05-nat
  threshold, and by the verdict written as a sentence with its consequence for how
  the result is framed.
- [ ] 16.9 Label what the result does and does not license — verify by the report
  carrying the scoping note's own caveats: that Core A′ is the best-measured
  region of the network, that enzyme competition is partly removed, that success
  is evidence the architecture works rather than that inference on the full model
  is well-posed, and K7's measured factor if it fired.
- [ ] 16.10 Fill the kill-criteria scoreboard, output T3 — verify by all seven
  criteria carrying a threshold, a measured value and a verdict, **published
  whether or not everything passed.**

### Phase 17 — Per-module calibration

**Goal:** check that each module's posterior is calibrated rather than merely
centred, and that a failure can be attributed to the composition rather than the
sampler.
**Done when:** rank statistics are uniform within test for the ODE block given
the stochastic path and for the stochastic block given the rate-constant
sequence, the joint is checked for the shared parameter, and coverage on the
boundary-crossing quantities is reported separately.
**PR:** _not started_

**Tasks are written when phase 16 has run.** This is deliberate rather than
lazy. Calibration for a modular, message-passing system has no established
protocol — the survey states that one will likely have to be built — and the
natural formulation, calibrating each site's outgoing message, depends on the
coupling structure phases 3 to 5 actually produce and on which sampler phase 16
settles on. What is committed now is the shape in D12: the four checks, the
attribution rule, the weighted-rank mechanic, and the rank-ECDF test with
simultaneous bands rather than a bare histogram. Anything more detailed written
now would be invented, and an executor would believe it.

### Beyond the horizon

Named so the spec's edges are visible, each a single honest line rather than a
fabricated task list.

- **Identifiability characterisation of the boundary-crossing parameters.**
  Needs phase 16's posteriors. The existing `check_identifiability` is a local
  Jacobian rank on an ODE-only problem and will not apply unchanged to a hybrid
  model.
- **Freeing the remaining parameters selectively.** Fixed at published values
  first, then freed to see whether the inference notices. Which ones and in what
  order is decided by what phase 16 finds.
- **Step 1b**, Core A′ plus the ribonucleotide-reductase and dNTP branch: the
  rung that adds a genuinely stochastic boundary and makes the nucleotide pools
  live. It is also where the five-reaction nucleotide module here stops being
  defensible.
- **The two surrogate proposals, and Step 2 against real data.** Neither is
  specifiable from here.

---

## 12. Amendment log

_No amendments yet._

Reality will contradict this spec. When it does, the change is recorded here with
its date and its reasoning, never applied silently. An amendment names what
changed, what evidence forced it, and which phase found it. R15 is the most
likely source: a framework need surfacing during a module phase is an interface
bug, and the wave plan's rule stands — that is a conversation on the main branch,
not an edit on a module branch.
