# Spec: Core A′ — inference across a whole-cell ODE/stochastic boundary

**Status:** in progress — phases 0 to 5b done (PRs #41 to #46, #48), with
phase 6 (#50), phase 7 (#49), phase 8 (#52), phase 9 (#59), phase 10 (#51),
phase 10b (#57), phase 11a (#62), phase 11 (#64) and phase 12 (#60); the
fan-out is complete. Phase 13 is split (§12, 2026-09-24): 13a, the framework
fixes assembly needs, is done (#66), and so is 13b, assembly (#68). Phase 14
is split (§12, 2026-09-25): 14a, the balance checks, is done (#70), and 14b,
the derivative, ensemble and census checks, is next
**Created:** 2026-09-03  ·  **Last amended:** 2026-09-26

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

## 0. Executive summary

- **What this is.** A specification for building and testing Core A′, a small
  organism carved out of the published JCVI-syn3A whole-cell model of Thornburg
  et al. 2022, inside InferCell.jl. It is the first of three steps toward running
  InferCell on the full syn3A model. Steps two and three are the full
  well-stirred model and the spatially resolved one, and neither is addressed
  here.

- **The question.** Whole-cell models mix formalisms: metabolism is a system of
  ordinary differential equations, gene expression is a stochastic jump process,
  and the two exchange state on a clock. Can parameter uncertainty be propagated
  *correctly* across that seam? Not just "does inference run", but: when we fit
  the model to data, do the resulting posteriors mean what they claim to mean?

- **What Core A′ contains.** Twenty-two metabolic reactions in four deterministic
  modules (glycolysis to lactate, sugar transport and lactate export, nucleotide
  recycling, and a lumped tRNA charging step), and seventeen genes in three
  stochastic modules (transcription, translation, transcript decay). The blocks
  swap state every simulated second and rate constants every sixty. About 92% of
  the cell is removed. It is the best-measured region of the network, which is
  why it was chosen.

- **How the two blocks talk.** Through state, not through parameters. Enzyme
  copy numbers set metabolic rates in one direction; nucleotide and charged-tRNA
  pools set transcription and translation rates in the other. No single rate
  constant appears in both blocks. This was discovered during specification,
  turned out to be a property of the published model rather than of our
  reduction, and reshaped the inference design (§4 D13, §12 amendment 2).

- **What "correct" means here, and how it is tested.** Everything is run on
  synthetic data generated from known parameter values, because that is the only
  way to have a ground truth. The posteriors must then cover the truth at the
  stated rate: a 90% interval should contain the true value about 90% of the
  time over many repeated fits. Rank statistics must be uniform. And an exact
  reference posterior on a two-gene toy must match the production one, so that a
  calibration pass can be attributed to the method and not to luck (§6 F8 to
  F10).

- **The result the work exists for.** A table of how much each of six target
  parameters' uncertainty shrinks under three data regimes: transcripts only,
  metabolites only, and both. Two cells carry the claim. A gene-expression
  parameter (the ptsG promoter strength) must shrink when only metabolite data
  is supplied, and a metabolic parameter (an enolase catalytic constant) must
  shrink when only transcript data is supplied. Each shows information crossing
  the seam in one direction (§6 F6).

- **What it is not.** It does not reproduce syn3A's 105-minute doubling time, its
  metabolite concentrations or its proteome, and no result is to be compared to
  them. Removing 92% of a cell guarantees that, and the code refuses the
  comparison rather than relying on the prose to. It does not fit real data:
  two published measurements are used as sanity checks only (§6 F5). It does not
  compare alternative model structures or build surrogates. Those are later
  steps (§7).

- **Where it could fail.** The largest risk is cost: if a full trajectory takes
  more than roughly ten seconds, the repeated fits that a calibration claim
  needs become unaffordable and the claim cannot be made. This is tested on a
  two-module toy before any real module is built, so the project can stop early
  and cheaply if it must (§11 phase 3, §8 K1). Every kill criterion is scored
  and published whether or not it fired (§6 T3).

- **How the work is organised.** Seventeen phases plus a phase 0 and a phase
  5b, each one
  reviewable pull request. Framework fixes come first, then the kill test on the
  toy, then the seven modules, then assembly, validation and inference. Phases
  fan out where they are independent and run serially where each check would
  mask the next (§11).

- **How the document is kept honest.** The decisions are recorded in §4, each
  with its reasoning and the alternatives rejected. Every departure from the
  published model is declared and, where possible, its cost measured (§6 T2).
  When reality contradicts the spec, the change is dated and logged in §12 with
  the evidence that forced it, never applied silently.

- **Where to read next.** `dev/notes/reduced-syn3a-scoping.md` for what Core A′
  is and why; `docs/handoff.md` for where the work currently stands; §11 for
  which phases are ticked; §12 for what has already changed and why.

---

## 1. The claim

**Parameter uncertainty propagates correctly across a genuine ODE/stochastic
boundary in a model derived from a published whole-cell simulation.** Concretely:
in Core A′ — 21 metabolic reactions and 17 genes carved out of Thornburg et al.
2022's well-stirred JCVI-syn3A model, with four deterministic modules and three
stochastic ones exchanging state every simulated second and rate constants every
sixty — posteriors over parameters whose information crosses the boundary are
*calibrated*, not merely centred: coverage on synthetic data is nominal at stated levels, rank
statistics are uniform within test, and what the data does and does not constrain
is characterised rather than assumed.

The enabling result, and not the claim itself, is that a composable
typed-interface architecture can host heterogeneous formalisms from a published
whole-cell model, execute the coupling between them, and support inference over
the whole with every reduction and every provenance choice enumerable from the
composed model.

**What "crosses the boundary" means here, precisely.** No parameter is *shared*
between the two blocks; §4 D13 establishes that this is inherited from the
published model rather than chosen by us. So the claim is not that one quantity
appears in both blocks' rate laws and is informed by both data streams. It is
that a parameter local to one block has its posterior moved by the *other*
block's data, through the state coupling. That is a harder claim and a narrower
one, and §6 F6 is where it is cashed out: the ptsG promoter strength shrinking
under metabolite-only data, and the tight-prior catalytic constant shrinking
under transcript-only data.

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
`dev/notes/figures/corea-coupling/fig1c_state_graph_corea.png` for the coupling
graph, which draws the charging step in the ODE block. Its predecessor,
`dev/notes/figures/reduced-syn3a-coupling/fig1r_state_graph_reduced.png`, draws
it in the stochastic block as the wave plan had it and is kept as a record. The
summary:

| | |
|---|---|
| ODE block | 22 reactions in 4 modules: glycolysis through lactate (10), the phosphotransferase cascade plus lactate export (6), nucleotide recycling (5), the lumped tRNA charging step (1) |
| Stochastic block | 52 reactions in 3 modules: 17 transcription, 17 translation, 17 decay, plus one translocation for ptsG. Matches the scoping note's own count. ~~Protein degradation would be seventeen more and is **ours if kept** — see phase 11~~ **Amended 2026-09-23:** there is none to keep — upstream defines `ptnDegRate` and never builds a reaction from it (task 11.9, §12) |
| State | 32 dynamic species and 5 chemostats, named and ordered once in the registry |
| Copy-number span | mRNA 0–2, protein 266–1355, metabolites ~200–74,000 particles |
| Coupling | 6 channels: 3 at 1 s, 2 at 60 s, and the volume channel, which is on neither clock. It runs in both directions after phase 5b — counts set the volume, and the volume re-enters an ODE rate law — and counts as one channel because one geometry carries both halves. Enumerated in §6 F1 |
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
not yet supported")` at `src/orchestrator.jl:27`. Core A′ is four ODE modules
and three jump modules. There is no hybrid path, no stub, and the refusal is
pinned by a test at `test/test_stochastic_ge.jl:58`. The same branch swallows a
uniform `:sde` composition, which is why the copy-number table's "chemical
Langevin sits naturally here" for the 266–1355-copy enzymes is unimplemented too.

**G2 — there is no execution layer for coupling.** Six of the seven edge kinds
are declare-and-validate only. **(Read as of 2026-09-03. Phase 1 made mass and
currency execute, phase 3 the catalytic and deferred-counter channels,
phase 4 the rate-constant one and phase 5 the volume chain, so **one of the
seven — clamped — remains declare-only**, its held value still travelling as
a fixed parameter beside the declaration, and **no phase schedules it**. The
*inbound* half of the volume channel — volume re-entering an ODE rate law,
which task 7.2 wants for the lactate exporter — was not built either, and
phase 5 refused such an edge by name rather than letting it resolve and never
run; **phase 5b built it**, so an inbound edge now carries a `param_slot` and a
`quantity` and executes at every handshake. G1 and G3 are
likewise closed. §2 is kept as the dated snapshot the phase ordering was
derived from, not refreshed.)** Searching `src/` and `test/` for periodic
callbacks, discrete callbacks, callback sets, `tstops` or operator splitting
returns nothing outside the `affect!` closures inside individual jump
definitions. `solve` is called in five places, always as a single one-shot call.
So the 1 s handshake, the 60 s rebuild and the volume rescaling are postulated by
the contract and written nowhere. The resolver's own docstring says as much: a
deferred counter "is exempt — the hook debits it, the RHS never reads it."
`DiffEqCallbacks` is already present in the Julia depot, so adding it needs no
network fetch from a compute node.

Note that `src/boundary.jl` is **not** this layer despite its name. It is
inference-level message passing between two independently solved problems, and it
does not couple trajectories. For Core A′ it in fact passes **nothing**: it
decides what to hand over by matching parameter names across the blocks
(`boundary_condition`, `:111-135`), and the two namespaces are disjoint. §4 D13
records why, and §7 records the consequence.

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
tractability rungs and the instruction that a reduced exact reference must exist
before any approximation is trusted. Its prior-double-counting hazard, and most
of its method families, presuppose a quantity shared between modules; §4 D13
records that Core A′ has none, and that note now carries a section on what would
change if one appeared.
`dev/notes/model-uncertainty-and-selection.md` supplies the reason model
comparison is predictive rather than evidence-based here.

---

## 3. Scientific validity

### Model

- **ODE block.** 22 reactions: the 21 published metabolic reactions plus the
  lumped tRNA charging step. Every intracellular metabolic reaction takes the published
  simulator's modular rate law, which is `Rxns.Enzymatic` in the source and *not*
  the SBtab `KineticLaw` column, which is inert for these reactions:

  ```
  v = E · ( kcatF·Π(Sᵢ/KmSᵢ) − kcatR·Π(Pⱼ/KmPⱼ) )
        / ( Π(1+Sᵢ/KmSᵢ) + Π(1+Pⱼ/KmPⱼ) − 1 )
  ```

  One rate-law term per unit of stoichiometry, so a coefficient of two appears
  squared. Three reactions take a different form: the phosphotransferase cascade
  is mass-action forward/reverse pairs, lactate export is passive diffusion,
  `P·(lac_c − lac_e)·3/r_cell`, and the lumped charging step is mass-action
  bimolecular:

  ```
  M_trna_c + ATP → M_trna_chg_c + AMP + PPi
  v = k_chg · [M_trna_c] · [M_atp_c]
  ```

  `k_chg` is **asserted and calibrated by us**, not imported: the published
  chain's four constants are hard-coded inline with no priors to inherit, so any
  constant here is asserted whichever form is chosen. See D13 for why this
  reaction is in the ODE block and D14 for how `k_chg` is fixed.
- **Stochastic block.** Constant-rate jumps. Transcription is constitutive: one
  reaction per gene firing on gene presence, with no promoter switching. Decay
  is first order in the transcript at one global catalytic constant over
  transcript length. Translation is first order in the transcript. **Every
  reaction in this block is zeroth or first order**, which is what makes the
  conditional likelihood of D10 closed form.
- **Coupling.** Six channels. Every 1 s: protein counts overwrite enzyme
  concentrations; accrued expression costs are debited against the metabolite
  pools; and translation's charged-tRNA consumption is debited against the
  charged pool with the uncharged pool credited. At the 60 s rebuild: the
  nucleotide pools set transcription's rate constants, and the charged-tRNA pool
  sets translation's. On neither clock: membrane-protein counts set surface area,
  which sets volume, which rescales every count-to-concentration conversion —
  and, after phase 5b, that same geometry re-enters an ODE rate law as a
  parameter, which is this channel's inbound half rather than a seventh channel.
  Enumerated with measured gains in §6 F1.
- **Discretisation.** Stiff implicit integration for the ODE block with pinned
  tolerances (§3, tolerance principle). Exact stochastic simulation for the jump
  block, which after D13 carries about 9,900 events per cycle, or 1.6 per
  handshake interval, rather than the 3.49 million it would carry with charging
  in it. The two advance by a 1 s split-operator exchange, matching the published
  model's own handshake.

### Assumptions and approximations, each with its regime

| Assumption | Holds when | What breaks outside it |
|---|---|---|
| Metabolites are continuous and deterministic | Pools are large. True for ATP at ~74,000 particles | Marginal at ~198 particles (13DPG), where relative Poisson noise is ~7% — and those are the pools with the largest concentration control coefficients. Bounded by check 1b, not assumed away |
| Enzyme counts enter the ODE deterministically | 266–1355 copies gives 3–6% relative Poisson noise | Core A′ therefore supports a claim about *parameter* uncertainty crossing the boundary, and honestly cannot support one about *stochasticity* crossing it. That is Step 1b's claim. D13 strengthens this caveat: the single most-fired process in the published stochastic block is now deterministic **by our choice** |
| One lumped charged-tRNA pool replaces 20 per-amino-acid chains | Untested. **Ours, not the model's** | The pool size is a quantity we chose. It sets the strength of the reverse channel into translation **and**, after D13, the low-pass filter on the dominant forward channel. See R8, which this upgrades |
| The lumped step is integrated deterministically in the ODE block | The registry marks both tRNA species metabolite-scale, and the step fires ~553/s. The scoping note already placed it here | The *formalism* is ours even though the placement is the note's. See D13 |
| `k_chg` is calibrated so the steady flux reproduces 553.1/s at nominal pools | Holds at nominal. **Amended 2026-09-23:** "nominal pools" is a total of 0.25 mM at 80% charged, both asserted by us, with ATP at the registry's 3.6529 mM | Off-nominal the flux is state-dependent rather than fixed, so the published demand figure becomes an emergent quantity rather than an input. See D14 |
| CTP, UTP and the amino-acid pool are chemostatted | Their sources are in modules Core A′ cuts. **Ours, not the model's** | Stops being defensible the moment Step 1b makes the pools live |
| Non-ptsG membrane growth is exogenous | ptsG is Core A′'s only membrane protein. **Ours, not the model's** | Growth reaches only ~1.07×, so the published 2× cap is never tested |
| The 60 s rebuild is piecewise-constant | It is what the published model does | A continuous variant is a different and arguably better model. Declarable, defaults to published. ~~The deviation is measured rather than argued.~~ **Amended 2026-09-05:** phase 4 refuses to *execute* a continuous cadence — the outer-loop mechanism of task 3.1 cannot represent it — so the deviation that is measured is a shorter piecewise-constant interval, which is what R11 varies. A genuinely continuous cadence remains declarable and labelled, and unmeasured. See R11 and §12 |
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
loosening it, and this cannot.

Report the single-run bound as a secondary number, derived as

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

**One exception, added 2026-09-10, and it is fenced so it cannot become a
loophole.** Some invariants are preserved by the *derivative* exactly, in
floating point, before any integration happens: the composed right-hand side
returns their weighted sum of derivative terms as literally `0.0`. Such a
residual contains no local truncation error, so it cannot fall when the solver
is tightened, and demanding that it fall fails a check that is *stronger* than
the demand was written for.

**The criterion is mechanical, and the phase must assert it rather than argue
it.** A check may use the exception only where a test evaluates the composed
right-hand side at sampled states and asserts `sum(nᵢ · duᵢ) === 0.0` — bitwise,
not `≈`. That is one assertion, it is specific to a *composition* rather than to
a moiety, and it cannot be talked into by a phase that would prefer an easier
check. **Crossing a module boundary is not the criterion.** It is the usual
*reason* the criterion fails, and it is neither necessary nor sufficient: a
moiety wholly inside one module fails it as soon as any reaction touches it with
unmirrored stoichiometry, and mirrored contributions *across* a currency edge
can satisfy it. Where the assertion does not hold, the fivefold fall stands —
unless the second branch below applies.

**A second branch, added 2026-09-11 from phase 8: exactly conserved by the
right-hand side, but not bitwise.** The branch above is written for a
composition whose derivative terms are exact IEEE negatives, which is what
happens when both terms come from one reaction in one module. A moiety can be
exactly conserved — `nᵀf ≡ 0` in real arithmetic, every reaction's contribution
cancelling term by term — and still not evaluate to `0.0` in floating point,
because the terms arrive from *several* modules and are summed in an order that
does not cancel bitwise. Phase 8's adenylate is the case: `(−adk1 − gk1) +
(2·adk1 + gk1) + (−adk1)`, with two further modules adding into the same slots.
There is no local truncation error in such a quantity either — the true
derivative of the invariant is exactly zero — so it cannot fall when the solver
is tightened, and the fivefold fall fails it for the same reason it fails the
first branch.

**Its criterion is mechanical too, and tighter than the first branch's
integrated bound.** A check may use this branch only where a test evaluates the
composed right-hand side at sampled states and asserts

```
|Σᵢ nᵢ · duᵢ|  ≤  1 ulp of the conserved sum,  at every sampled state
```

which is a statement about a *single evaluation* and so cannot absorb an
accumulating leak. It is not a softened `≈`: one ulp is the tightest non-zero
bound expressible, and a moiety with any unmirrored stoichiometry misses it by
orders of magnitude rather than by a factor. Phase 8 measured, on its own
composition over 106 trajectory states, adenylate at exactly **1.000 ulp** of
3.9539 mM, guanylate **0.125**, flux-corrected phosphate **0.115** — while a
mutated stoichiometry puts adenylate at 3.70 mM, which is 1.7e16 ulps. **Then
the ladder above applies over at least six decades, with one substitution.**
The first branch's flat 100-ulp bound is not the right yardstick here: phase 8
measured the integrated residual wandering between 128 and 1,119 ulps of the
sum over nine decades with no trend, while the evaluation count grew 36-fold —
so it neither falls nor accumulates, and a flat ulp count is a magic number
fitted to one composition. The bound is instead **`tol_C`, which §3 already
defines**: require every rung at least **three orders of magnitude below its
own `tol_C`**, and the largest rung within **100×** of the smallest. That
tightens automatically as the solver tightens, so unlike a fixed threshold it
cannot be passed by loosening; phase 8's nine rungs clear it by 3.9 orders at
the tightest and 12.1 at the loosest, with a spread of 8.7×. Report the ladder,
not a single pair.

**Amended 2026-09-25 (phase 14a): the per-evaluation bound is `n` ulps, where
`n` is the number of weighted terms in the moiety.** One ulp was fixed on
phase 8's three-term adenylate. At composition scope the carbon sum has 13
terms and the phosphate sum 21, and they arrive from several modules, so they
round more. On the assembled model carbon measures 3.94 ulps per evaluation
and phosphate 4.50. Against the floating-point summation bound
`n·ε·Σ|wᵢ·duᵢ|` that is 0.04 and 0.02, which is roundoff. Their integrated
residuals are 2.3e-9 and 8.7e-9 of `tol_C`. A mutated stoichiometry still
misses the widened gate by 11 to 18 orders, measured on the six
stoichiometry mutants (job 17539616). So it cannot absorb a leak. A moiety
that passes at one ulp still reports `:one_ulp`. The integrated ladder
condition is unchanged. See §12.

**The two branches are one exception with two gates, not a licence to argue.**
Neither may be claimed by reasoning about linearity, about which module a
moiety sits in, or about whether the trajectory is restarted. A phase asserts
the gate its composition actually passes, or it keeps the fall.

**What the exception then asserts, with numbers, so it is as checkable as the
fall it replaces.** Run the ladder over at least six decades of `(abstol,
reltol)`. Require every rung's residual within **100 ulps** of the conserved
sum, and the largest rung within **100×** of the smallest. Report the ladder,
not a single pair. **Amended 2026-09-23:** a composition that passes the first
gate's bitwise test may instead assert the second gate's `tol_C`-relative
ladder: every rung at least three orders below its own `tol_C`, the largest
within 100× of the smallest. Two conditions keep that from becoming a way to
pass by choosing rungs. The asserted ladder must **end at the pinned
tolerances**, and **every rung run must be reported**, asserted or not, from a
Slurm run of record, since the residual is roundoff and differs between
machines. The flat 100-ulp bound is measured in ulps of the conserved sum, but
the residual is roundoff set by the composition's *largest* states. A small
moiety beside large pools therefore fails the flat bound without leaking
anything. Phase 9's 0.25 mM tRNA pair is bitwise zero at every evaluation and
drifts 29 to 2,043 ulps of its sum (job 17044321). See §12.

**And it is not literally flat, which matters.** The integrated residual is
linear-algebra roundoff in the implicit solver's stage solves, accumulated over
the step count, so it *rises* mildly as tolerances tighten and the solver takes
more steps. Phase 6 measured 2 ulps at `(1e-5, 1e-3)` and 60 at `(1e-12,
1e-10)`, wandering non-monotonically in between — a 30× spread over eight
decades, against the ≥5× *fall per decade* the principle demands. The claim is
"bounded at the floor and not falling", never "constant".

**Scope: one continuous solve.** The exception is about the right-hand side, and
a handshake is not the right-hand side. At composition scope the driver rewrites
every diluting ODE concentration at each of ~6,300 handshakes, so a concentration
sum is not conserved at all once the cell grows, and even the particle sum
accumulates per-handshake rounding. **A check inheriting a module-local flatness
bound into the assembled model must restate it in particles and carry the
`N_restarts` factor above.** Task 14.4 is the one this binds.

Either way the check carries **the same mutation test**, and which rule applies
is stated per check rather than chosen per run. See §12.

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
| 0 | Round-trip policy | Exact-zero, square-root or linear residual scaling per policy; deterministic rejected. **Amended 2026-09-25:** on the assembly the growth laws do not separate, since the pool fractions are not systematically biased there. What is asserted is fractional carry at roundoff, within 10⁻⁶ of `tol_C`, and deterministic rounding at least 10⁶ times it. The carry's roundoff does accumulate with handshake count, which is what `N_restarts` in `tol_C` is for. Stochastic rounding and the growth ratios are recorded by the driver, not asserted (§12) | Handshake phase, re-asserted on assembly |
| 1 | Non-negativity | Every state above the negative of its own integrator bound at every save point, naming the first state and time to violate | Per-module as a smoke test; evidence only assembled |
| 1b | Particle-floor honesty | Report every state's minimum in particles; flag any below 500. Cross-check the three smallest pools against a chemical-Langevin ensemble; exclude from the likelihood any observable whose ODE trajectory leaves the ensemble's 90% band by more than the assumed observation noise. **The tRNA pair is in scope and its particle count is ours** (D14), so this check is also what bounds the pool size we assert. **Amended 2026-09-26 (phase 14b):** the observation noise is not fixed until task 15.4, so 14b reports, per pool, the largest relative excursion of the ODE trajectory outside the ensemble's 90% band. That is the smallest noise at which the pool would be excluded. Phase 15 applies the exclusion (§12). **Amended again 2026-09-26:** on the assembly the smallest pools sit at about one particle for most of the cycle, where neither the ODE nor a diffusion describes them. A pool whose cycle median is under 10 particles is excluded outright, with no ensemble. The cross-check runs on the three pools with the smallest medians of at least 10 (§12) | Assembled |
| 2 | Carbon balance | Glucose in against lactate out plus intermediates plus biomass; and the homolactic ratio, which is **analytically exactly 2.000** | Assembled only — spans transport, glycolysis and export |
| 3 | Redox | NAD⁺ + NADH invariant | Per-module (glycolysis); exactly invariant there, so assembly adds nothing. **Exact means exact**, and phase 6 measured it: the module-local residual is bounded at 2 to 60 ulps across eight decades of solver tolerance, so what this row asserts *on the module* is the tolerance principle's exact-invariant branch and not the fivefold fall. **"Assembly adds nothing" holds only for the derivative.** Growth dilution rewrites every concentration at each handshake, so the assembled restatement is in particles and carries the `N_restarts` factor; task 14.4 owns it (amended 2026-09-10; see §12) |
| 4 | Adenylate and guanylate, **over a full 6,300 s cycle** | Each moiety separately. After D13 charging is inside the ODE block and its ATP→AMP+PPi transfer conserves adenylate internally, so there is no declared drain to correct for and the check needs only the inbound mass flux. Three configurations: all five recycling reactions (conserved); adenylate kinase removed; pyrophosphatase removed (**amended 2026-09-10:** pyrophosphate strands the phosphate moiety and stalls the pathway; it cannot diverge in a model whose phosphate is closed, which this one is — see §12). **The kinase-removed assertion is a threshold crossing, not exhaustion** — under mass action the drain is proportional to ATP, so ATP decays exponentially and never reaches zero. Assert the time at which ATP falls below 1% of its initial value, and reconcile against the 144 s the scoping note computes for a constant drain | Assembled only. Standalone the recycling module conserves both moieties trivially, so the charging module or a drain double must be composed with it |
| 4b | Phosphate closure | Free phosphate plus every phosphorylated species with pyrophosphate counted twice, minus flux across the two inbound mass edges. Two forms: exact with the GTP-branch reactions inactive, flux-corrected with them active. **Subtract the flux; do not relax the tolerance.** **Amended 2026-09-11:** standalone, nothing consumes GTP, so the inbound flux is exactly `u[GTP] − u0[GTP]` — a state difference, which keeps the corrected closure a linear functional and inside §3's exception. Assembled, translation consumes GTP and that proxy is wrong; the restatement needs an accumulated flux, and a quadrature does carry tolerance-scaling error, so the fall binds there (task 14.4, meaning 14.6). See §12. **Amended 2026-09-25 (phase 14a):** no quadrature is needed. The assembled closure is whole-cell. It counts transcripts by sequence length, the CTP and UTP chemostat rows by what each counter carries, and the accrual not yet applied. Every boundary term is then an exact integer, and the closure passes the exception's second gate at `n` ulps. Charging's contribution is internal and exactly closed after D13 — ATP's three phosphates become AMP's one plus pyrophosphate's two — so it needs no correction, unlike transcription's pyrophosphate, which does | Assembled only |
| 5 | Carrier conservation | **Four** independent sums, one per phosphotransferase carrier, each named separately so a failure localises | Per-module before translation exists; restated assembled as "conserved up to what translation adds", again by subtraction rather than by widening |
| 6 | Nominal trajectory | All of the above at published parameters, plus check 7's census | Assembled |
| 7 | Clipping census | Count handshakes at which any deferred counter carries a deficit, at published parameters and across 200 prior draws. Require **zero at published parameters** (scored for K5 rather than gated, amended 2026-09-25; see task 14.8). If more than 5% of prior draws clip, the non-smooth drain is in the operating regime rather than at a measure-zero point, gradient-based sampling of the ODE block is invalid as posed, and smoothing becomes mandatory. **The charged-tRNA counter is the one to watch and D13 made it so**: it debits ~553 residues per second against a pool of order 10³ particles, a buffer of seconds, so it can clip on ordinary Poisson fluctuation in translation firing. Either the pool is sized so it provably cannot clip, or this census is what catches it. **Amended 2026-09-05:** the census is scored at the **published 1 s drain**. Phase 4's `drain_interval` changes what it counts in both directions — each debit is `steps_per_drain` times larger against the same pool, and there are that many fewer of them — so a coarse drain would fire or dilute this check on our choice rather than on the model. `clipping_census`'s fraction is per drain, not per handshake, for the same reason | Assembled |
| 8 | Metabolic control analysis identities | On the frozen-expression chemostatted variant at steady state, the summation theorems hold exactly: flux control coefficients sum to one per flux, concentration control coefficients sum to zero per metabolite. Require both within 1e-6, computed by automatic differentiation through the steady-state solve. **The perturbation set is 18 quantities, not 17**: the lumping deleted charging's synthetase, so `k_chg` is the perturbable parameter that reaction contributes and the identity fails without it. The frozen-expression variant must also clamp the tRNA pair or specify the transfer as constant forcing, because with expression frozen nothing consumes charged tRNA, the pool saturates and the charging flux goes to zero. **Amended 2026-09-26 (phase 14b):** the 18 quantities cannot close the sums. Only thirteen genes multiply a rate, and PGK and PYK each multiply two. The four PTS genes are carrier totals, and the steps are bilinear in two carriers. Export and the demand forcing have no enzyme. So the perturbation set is **one rate multiplier per reaction** of the frozen variant: the 22 ODE reactions plus its tRNA and GTP demands. Gene-level coefficients are reported by summing a gene's reactions, and carrier-total responses are reported apart. The coefficients come from the implicit function theorem on automatic-differentiation Jacobians at the solved steady state. A reaction with zero steady flux is named and left out of the flux sum (§12) | Assembled, and independent of every balance check because it tests derivatives rather than balances |
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
   0.332–2.406 copies against measured means 0.292–2.178, with sixteen of the
   seventeen within a factor of two and FBA (`JCVISYN3A_0131`) the one outlier.
   *An earlier draft of this row read 0.44–2.87 "within about 1.5×"; that band
   was computed from `rnaDegRate`, which the published model defines and never
   executes — see §12's 2026-09-10 phase 10 entry, amendment B.* This is
   genuinely out of sample: the parameters come from proteomics, the
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

**A supporting argument D13 supplies.** Charging, as an ODE module, contributes
to ATP, AMP and pyrophosphate that nucleotide recycling owns. Executing that
contribution is exactly what phase 1 builds. So the framework-first ordering is
not only about avoiding rework: the charging relocation is unimplementable
without phase 1, which is one more reason it cannot sit at the end.

**Consequence to record.** The four drafted module designs each open by stating
that the interface is frozen and read-only. That was true when they were written
and is now false in two specific respects. Phases 6 to 12 consume the *amended*
protocol, and each amendment is a phase deliverable rather than a module's
concern.

### D1. `Mode` is the column that runs

**Decision.** Every imported kinetic constant and initial concentration comes
from the balanced table's `Mode` column, with ~~its~~ **the only** geometric
standard deviation ~~,~~ **the table offers,** and the source file recorded
(amended 2026-09-10; see §12).

**What "its" cannot mean, measured in phase 6.** The balanced `Parameter` table
has exactly one width column and it is `UnconstrainedGeometricStd` — there is no
balanced width to pair with the balanced mode. For most rows that costs nothing,
because balancing did not move the value: across phase 6's extract **all 13
concentrations and all 32 Michaelis constants have `Mode` equal to
`UnconstrainedGeometricMean` exactly**, so the width belongs to the median it is
given. **Ten of the twenty catalytic constants do not**, which is expected —
they are what the thermodynamic balancing constrains — and there the prior's
width is inherited from a distribution not centred on the prior's median. The
worst is `R_TPI`'s reverse constant, mode 4 against an unconstrained geometric
mean of 65,341.7 carrying a 1.05 width. **This is ours, labelled, and it is the
population inference targets are drawn from**, so any recovery claim on a
catalytic constant states it.

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

**What eight input deltas do not tell you, measured in phase 6.** The
declaration as first written enumerated the eight constants and their factors,
which understates it. At the registry's initial concentrations the imported
constants **reverse the sign of the pathway's entry reaction**: R_PGI runs at
**−0.3325 mM/s** on the balanced column against **+4.5963 mM/s** on the
`Quantity` column the published simulator runs, and R_FBA is slower by **253×**
(0.00901 against 2.27985 mM/s). Both recover once the pools move — the module
integrates non-negatively over a full cycle and drains its carbon to lactate
either way — but **no statement about an initial flux in this module is the
published model's**, and `reduction_declarations` now says so with these
numbers rather than only with the input deltas. See §12 (2026-09-10).

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
should be.** That is the only reverse channel **into transcription**, so porting
the bug would mean demonstrating bidirectional coupling with a gain a typo nearly
doubled. It is not the only reverse channel in the model: after D13 the
charged-tRNA pool is an ODE state feeding translation's rate constants at the
same rebuild, which is a second one. An earlier draft of this decision said
"the only reverse channel in the model", contradicting the assumptions table and
R1; that is corrected here.

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
A′'s corrected demand of ~1,106 glucose/s that is ~13.9 M lactate over a cycle
— and note after D14 that this demand is an *emergent* quantity contingent on
`k_chg` and the tRNA pool rather than an input, so the figures below are the
nominal case rather than a bound:
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
the radius grow and the export rate constant then changes with it. Phase 5b
builds the channel that carries it: an inbound `VolumeEdge` fills the radius
slot at every handshake.

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

Two more join that set after D13 and D14, making thirteen: the charging rate
constant `k_chg`, and the total tRNA pool size, for which the registry records no
value at all. They come from a different module, so the enumeration is no longer
transport-only. **Amended 2026-09-23: fourteen.** D14's derivation also needs
the nominal charged fraction, a third asserted charging quantity (§12). The
pool size and the fraction are fixed at construction: varying them means
rebuilding `TrnaCharging`, not freeing a parameter.

**Thirteen counts the asserted *rate* constants, and phase 7 reports a larger
number for a different set.** `asserted_prior_params` selects on informedness,
so it also returns every initial condition whose source carries a point value
with no quantified uncertainty, plus the two geometry-and-ratio scalars this
project asserts and the external-glucose clamp — 23 from `PtsTransport` alone.
The thirteen here are the ones an inference phase might *free*; the 23 are every
prior in the module that is ours rather than inherited. Both are real; a result
that quotes one should say which.

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
30 s, so a twofold change in aggregate promoter strength is unmissable in the
adenylate and guanylate pools. The two are not the same shape, though, and D13
is why. The guanylate side is driven directly by translation's cost counter. The
adenylate side is ~81% charging, which after D13 is an ODE reaction, so its
forward channel is **two-stage and buffered by the tRNA pool**: translation
consumes charged tRNA, the uncharged pool refills, charging then draws ATP. The
pool size sets the lag, and we assert it. The reverse channel
has a measured gain of 0.044 to 0.051. **The published coupling is therefore
nearly unidirectional in gain even though it is bidirectional in topology**, and
that is a substantive quotable statement about whole-cell model architecture
rather than a shortcoming of the reduction. It is worth more than a symmetric
result would have been.

### D10. The conditional scheme, and the reference is mandatory

**Decision.** Exploit the conditional structure the boundary actually has, build a
reduced exact reference first, then the production sampler. There is no protocol
comparison; D13 explains why one is not available.

Because no parameter is shared across the boundary (D13), the coupling is a fixed
point in state space rather than a common entry in θ. That makes the following
factorisation the natural one, and it is what the spec commits to:

```
theta_ODE | path, data   ->  gradient-based sampler
theta_CME | path         ->  closed form
path      | theta, data  ->  mechanism decided in phase 16
```

**Two things about this that are easy to get wrong.**

*The middle block's conjugacy is real but narrow.* Given the path, each gene's
firing count has a likelihood proportional to `k^n exp(-k ∫…)`, which is gamma in
`k`. But the free parameters enter `k` as the product of polymerase turnover and
promoter strength, so **the conjugate parameterisation is exactly the degenerate
direction** D11 rejects the polymerase constant for and K6 tests. That is a
convergence rather than a coincidence: the identifiable parameterisation and the
tractable one are the same, and both point at the product.

*Conditioning on the path does **not** make the first block smooth.* An earlier
draft of this decision claimed the cost drains become known deterministic
forcing once the path is fixed. They do not: the drain is clipped at zero against
the pool, and the pool depends on `theta_ODE`, so the non-smoothness survives
conditioning. K5 governs inside the conditional step exactly as it does outside
it, and check 7's census is what decides whether it bites.

**M0, the reference.** A deliberately small Core A′: two genes, 600 s horizon,
small enough that a very long run is a defensible ground truth. The three-block
scheme above, run to convergence with a deliberately expensive path update. This
targets the true joint posterior and is the only unambiguous ground truth the
project will get.
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
>
> **Amended 2026-09-05, from phase 4's miniature run: "the nominal trajectory"
> is an ensemble, paired per seed, and sampled at drain-aligned instants.**
> Coarsening the drain changes how often the hook clears the counters and the
> driver therefore rebuilds the SSA's propensity aggregation, which consumes
> randomness — so one seed gives three different paths and a single-path
> difference is Monte Carlo noise. And between drains a coarse configuration
> holds up to `drain − interval` seconds of unpaid cost, so its pools sit high
> by exactly that: an unaligned sampling instant measures a deterministic
> bookkeeping offset, `(drain − interval)·cost_rate/pool`, and not the dynamics.
> Report the difference at a pre-chosen instant, the maximum with the
> multiplicity it was selected from, and the sawtooth separately. See §12.

Whichever passes becomes a labelled reduction with a *measured* cost.

**The path update is deliberately left open, and phase 16 closes it.** After D13
the stochastic block is purely first order and factorises over genes, and
promoter strength reads the fixed proteomics count rather than the live protein
count, so there is no feedback inside the block. That may admit exact forward
filtering over a truncated state space, or no augmentation at all. Committing to
a particle method now would foreclose the cheaper route D13 just opened. Phase 16
chooses among three options and records which:

| Option | What it needs |
|---|---|
| Particle Gibbs with ancestor sampling | The stochastic block written against the `SSMProblems` interface |
| Exact forward filtering | A defensible truncation, hardest for the protein counts, which reach 1355 |
| No augmentation | The 60 s aggregate drain passing the granularity measurement above. **Annotated 2026-09-25:** it fails, at +288% (task 13.4), so this option is closed unless the model changes |

**Two installed facts that constrain the first option**, recorded here so nobody
meets them the hard way. Ancestor sampling is implemented only for the
`SSMProblems` state-space interface, not for Turing's model macro, so a `@model`
cannot reach it. And Turing's compositional Gibbs does accept a particle sampler
on a latent block, but its task-copying layer deep-copies only arrays and
references, so a solver object held in a particle's state is **silently shared
across forked particles** rather than copied. That is a correctness trap, not a
performance note.

**Alternatives considered.** A protocol comparison against the explicit cut and
the iterative exchange, which earlier drafts made the architecture claim's
headline evidence. Both are unavailable: they decide what to pass by matching
parameter names across the blocks, and D13 shows that set is empty. §7 records
the consequence and the survey note carries what would be needed if a shared
parameter ever appeared.

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
| `kcat_ENO` | 62.18 | **Positive control.** ~~Tightest prior in the core~~ **Tight prior (gstd 1.170)** and the capacity-tightest enzyme, so its concentration control coefficients are the largest available. Must recover tightly. *Amended 2026-09-10:* phase 6's extract shows it is not the tightest — eight forward constants are, `kcatF_R_LDH_L` at 1.0512 tightest of all. The pairing still works, because what the control needs is a tight prior on the capacity-tightest enzyme, and only ENO is both |
| `kcat_FBA` | 59.7 | **Stress control.** ~~Loosest prior in the core (gstd 1.747)~~ **Loosest *forward* prior in the core (gstd 1.747)**. The diagnostic is shrinkage: its posterior must be visibly wider than ENO's. *Amended 2026-09-10:* three **reverse** constants are looser still — `kcatR_R_PFK` at 33.03, `kcatR_R_PYK` at 11.34, `kcatR_R_FBA` at 8.36 — and two of them exceed the prior-default width, so they enumerate as uninformed. The 1.747-against-1.170 contrast the diagnostic rests on is unaffected |

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

**The ridge is also the conjugate direction**, which is worth knowing before any
Jacobian is computed. D10's closed-form conditional for the stochastic block is
gamma in the product of polymerase turnover and promoter strength, so the
tractable parameterisation and the identifiable one coincide. K6 tests this
numerically; the analytic form says the answer in advance.

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
3. **Joint**, over the whole target set. Note what this is *not*: there is no
   shared parameter to check, because D13 establishes there is none. What the
   joint catches is the composition itself — the path update, and whether the two
   conditionals are consistent with each other — which is the only place the
   scheme can go wrong without either conditional failing.

There is no fourth check. An earlier draft had message calibration, ranking the
truth within each module's outgoing message to catch prior double-counting
before it reached the posterior. With no shared parameter there is no outgoing
parameter message to rank, so the diagnostic is vacuous here. It is worth
building for a system that does have one, and the survey note records that.

**The attribution rule, stated once.** Calibration checks detect sampler error
and approximation error together. Run the two conditionals; if both pass and the
joint fails, the fault is the composition. That is the only route by which a
calibration failure means "architecture" rather than "sampler", and it is why M0
is mandatory rather than nice to have.

Two mechanics: ranks for a weighted posterior must be computed with the weights,
and use the rank-ECDF test with simultaneous bands rather than a bare histogram,
because a ten-bin histogram has no usable power at the replication counts this
project can afford.

### D13. No parameters are shared across the boundary, and that is inherited

**Decision.** Accept that Core A′ has no parameter appearing in both blocks, treat
it as a first-stage assumption rather than a defect to fix, and move the lumped
tRNA charging step to the ODE block, which the scoping note already required.

**The finding, with its evidence.** The two blocks' parameter namespaces are
disjoint. The ODE side carries `kcatF_R_*`, `kcatR_R_*`, `km_R_*_M_*`,
`conc_M_*`, the eleven asserted transport constants and now `k_chg`; the
stochastic side carries 17 promoter strengths and a handful of gene-expression
globals. Nothing is in both. Four independent lines of evidence say this is the
published model's structure and not ours:

- The derived parameter-coupling figure in
  `dev/notes/figures/minimal-cell-coupling/` draws metabolism and gene expression
  as **separate panels with no edge between them**, and its Panel C is titled
  "Gene expression: no factor structure at all".
- Metabolism went through SBtab parameter balancing, with priors, balanced modes
  and posterior widths. Gene expression did not: its globals are written into the
  source with, in that figure's words, "no prior, no posterior, no provenance".
  Two blocks parameterised by two different processes.
- Upstream, 22 of 23 ODE synthetase reactions are commented out precisely because
  charging moved to the stochastic block. Where one reaction could have carried
  parameters on both sides, one copy was disabled deliberately.
- `src/boundary.jl` decides what to pass by matching parameter names
  (`boundary_condition`, `:111-135`). Run on Core A′ it passes nothing.

**Where sharing does live: inside each block, not across.** The published
stochastic block is a star — roughly 19 scalars set all 455 of its rate
constants, so perturbing one moves all of them together. Core A′ inherits the
hub of that star at 52 reactions rather than 455, which is exactly the ridge D11
rejects the polymerase constant for. On the metabolic side, Wegscheider
constraints tie equilibrium constants to standard chemical potentials within each
balanced file, while **32 compounds carry a different potential in two
separately balanced files** (`mu0_cross_file.csv`), so the assembled model is not
thermodynamically consistent as a whole. That is the same class of hazard D2
documents, and it is the one place a genuine shared-parameter problem exists.

**Why the charging step moves, and why that is a correction rather than a
departure.** The conditional scheme of D10 is only feasible if the stochastic
block's path is small enough to represent. It was not:

| | events per cycle | per 1 s handshake |
|---|---|---|
| charging, with it in the stochastic block | 3,494,391 | 555 |
| everything else in that block | 9,873 | 1.6 |

Charging is 353 times every other stochastic event combined, which is a
forward-cost problem before it is an inference one, against K1's 10-second
budget. But the decisive argument is that **the scoping note already places it on
the ODE side and always did**: its CME-side count is "17 genes × 3 reactions
plus one `TranslocRate` = 52 reactions", with no charging reaction; both tRNA
species appear in its list of dynamic ODE states; and the charged pool sits in
its energy-interface table beside ATP, ADP, Pi, AMP, PPi, GTP, GDP and GMP.
`dev/plans/reduced-syn3a-wave-plan.md` assigned charging to the stochastic block,
this spec inherited that, and the precedence rule at the top of this document
says the note wins on model content. The frozen registry agrees: it marks both
tRNA species `:metabolite` regime, the field that records which formalism is
defensible.

Independently, the lumping deleted the structure that justified the upstream
placement. The published model puts charging in the stochastic block with 20
chains, an explicit synthetase protein reactant drawn from its own translation,
and per-amino-acid resolution. Core A′ lumps all of that away, leaving a bulk
flux between two metabolite-scale pools with none of the properties that
warranted a discrete treatment.

**Two consequences that are costs, not benefits.** First, the dominant forward
channel becomes two-stage. ATP turns over in about 109 s against ~74,000
particles, so about 679 per second, of which charging was 553 — roughly 81%.
With charging inside the ODE block its flux contains no stochastic-block
quantity, so the adenylate forward channel survives only through the tRNA
transfer: translation consumes charged tRNA, the uncharged pool refills, charging
then draws ATP. **The tRNA pool size therefore sets a low-pass filter on the
dominant forward channel, and that pool size is a quantity we assert.** This
upgrades R8. It is not an artefact of the move — the published model has the same
two-stage structure — but the pool size is ours. Second, the non-smooth boundary
moves onto that pool rather than disappearing; see check 7.

**What the move fixes.** As a stochastic module, charging declared currency edges
on ATP, AMP and pyrophosphate, which a jump process cannot execute against ODE
states. It would have needed three more clipped deferred counters. In the ODE
block it contributes to those pools directly through the mechanism phase 1
builds, which is what that phase exists for.

**No frozen contract is touched.** The registry, the seven edge kinds, the
resolver and the loader are unchanged, and so are the 32 dynamic states. Only
which module integrates the tRNA pair changes, and the registry never said.

**The extension path.** Nothing here forecloses a shared-parameter treatment
later. Unfreezing cell volume on the stochastic side would create one, since the
published model computes its polymerase capacities from a frozen initial radius
while the ODE side reads volume live. Step 1b's dNTP branch is the other
candidate. What that would require, and what the cut and iterative protocols need
in order to run at all, is written up in
`dev/notes/modular-bayesian-inference-heterogeneous-modules.md` rather than here.

### D14. `k_chg` is calibrated, and the calibration is not the test

**Decision.** Fix `k_chg` so the charging flux reproduces the published residue
demand at nominal pools, mark it `:asserted`, hold it fixed by default, and
report the total tRNA pool size as a separate asserted quantity.
**Amended 2026-09-23:** and the nominal charged fraction too, which the
derivation cannot do without — `k_chg = demand / ([M_trna_c]·[M_atp_c])` needs
the uncharged pool, not the total. Defaults 0.25 mM total and 0.8 charged, both
`:asserted`; see §12.

**Why a rate constant rather than a fixed flux.** A fixed flux would make
charging a constant drain, which is what the published demand figure of 553.1 per
second describes. But a constant drain drives ATP to zero in finite time and is
insensitive to the pool it draws from, so it cannot participate in the two-stage
forward channel D13 describes. Mass action is the minimal form that responds to
both pools.

**The calibration must not become the acceptance test.** 553.1 per second is
derived: 3,484,518 residues over a 6,300 s cycle. If `k_chg` is set so the steady
flux equals 553.1 and the module is then checked against 553.1, the test verifies
arithmetic. So the flux target is derived independently, from published residue
demand and the asserted pool size, and `k_chg` is the **derived** quantity that
gets reported. The module's own acceptance test is conservation and the
threshold-crossing behaviour of check 4, not flux agreement.

**The pool size is the load-bearing assertion.** The registry records no initial
value for either tRNA species. The pool size sets the charging flux at fixed
`k_chg`, the buffer that decides whether check 7's counter clips, and the lag on
the dominant forward channel. It is one number with three consequences, so it is
reported as an asserted quantity in its own right rather than folded into
`k_chg`, and check 1b is what bounds it.


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
| Per-gene sequence data | genome record, annotation compilation, second genome record, proteomics table, measured transcript counts | One extract, five upstream sources, one column each, so no cross-file ambiguity. Three is the count for the sequence columns alone; the protein copy number needs two more, because the reduced genome record carries no protein ids (§12, 2026-09-10 phase 10, amendment C) |
| Charging rate constant, tRNA pool size and charged fraction | nothing — **asserted by us** | No upstream value exists. `k_chg` is derived from published residue demand, the asserted pool size and the asserted nominal charged fraction (D14, amended 2026-09-23); all three are recorded as asserted, not vendored |

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
are already in the depot, so no new dependency requires a network fetch. So are
the packages D10's path update might need: AdvancedPS 0.6.2, SSMProblems 0.5.2,
Turing 0.39.6 and Libtask 0.8.8, all pinned in `Manifest.toml`.

**Two facts about those packages belong here rather than in a phase**, because
one of them is a correctness trap. Ancestor sampling is implemented only for the
`SSMProblems` state-space interface, not for Turing's model macro, so a `@model`
cannot reach it. And Turing's compositional Gibbs does accept a particle sampler
on a latent block, but its task-copying layer deep-copies only arrays and
references, so a solver object held in a particle's state is **silently shared
across every forked particle** rather than copied. Nothing errors; the results
are simply wrong. **If a
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
twice, on the toy in phase 3 and at full scale in phase 13b, and it feeds the
first kill criterion, which is scored per posterior rather than per trajectory
(§8 K1, amended 2026-09-24).

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
| **F1** | **Channel gain table.** Six rows, matching §3's coupling bullet: enzyme concentration, the expression-cost drain, the tRNA transfer, nucleotide pools into transcription's rate constants, the charged pool into translation's, and volume. Each with its analytic form and its measured value. Volume's row carries **both** directions after phase 5b: the dilution gain, and the gain of whatever rate law reads the geometry — for the lactate exporter's `3P/r` that is `∂ ln rate / ∂ ln r = −1`, analytic and needing no measurement | C1 | The reverse channels at 0.044–0.051 and the charged-tRNA figure R1 measures; the guanylate forward channel as a ~30 s turnover; the adenylate forward channel as a two-stage gain with the tRNA pool's lag stated separately, per D13. **The most important table in the work, and it appears early** |
| **F2** | **Concentration control coefficient heatmap**, 17 enzymes × ~20 metabolites, with the summation identities as guard | C1, and the observable choice | Answers the scoping note's open question. Its largest rows *are* the metabolite panel |
| F2b | Flux control coefficients, same layout | C1 | Shows them near zero at 3.6% utilisation — the quantitative reason fluxes are ruled out of the likelihood |
| **F3** | **Invariant residual against integrator tolerance**, log-log, one line per invariant, slopes required positive — except for an invariant that passes either gate of §3's exact-conservation exception, whose line is bounded near the floating-point floor with a slope of zero, positive or negative, drawn with the ulp band marked and labelled as such rather than counted as a failure (amended 2026-09-10 and 2026-09-11; see §12). **The label must say which gate**: bitwise-zero derivative, or one ulp per evaluation, or `n` ulps for an `n`-term sum (amended 2026-09-25). Phase 8's adenylate is the second and wanders non-monotonically over nine decades, so "slightly negative" is too narrow a description of what a passing flat line looks like. Beside it a **mutation table**: per check, the injected error, the residual it produced, the bound it exceeded | C2 | The tolerance principle, and the evidence that every check can fail |
| F4 | Conservation residual against handshake count under the three rounding policies | C2 | Exact-zero, square-root, linear. Justifies the policy and forestalls a rounding artefact being read as a leak |
| **F5** | **The two external comparisons.** Predicted against measured transcript steady states, 17 points, log-log with a twofold band; and the protein fold-change histogram with the published median marked | C2, external | The only two places Core A′ touches data it did not consume |
| F5b | Analytic drain-noise calculation against simulated | C1 | Closed-form variance of the cumulative expression drain. **Must be computed after D13, not before:** removing 3.49 M near-smooth charging events and leaving ~10⁴ translation events carrying residue-weighted increments makes the drain noise relatively *larger*, not smaller. That recomputation is the measured cost for T2's formalism row |
| **F6** | **The shrinkage table.** All six targets × three data configurations — transcripts only, metabolites only, joint — reporting posterior over prior standard deviation and the prior-to-posterior divergence | **C3. This is the figure the whole work is for** | Two load-bearing cells: the ptsG promoter under *metabolites only* (shrinks ⟹ uncertainty crossed toward the ODE block) and the ENO constant under *transcripts only* (shrinks ⟹ crossed the other way) |
| F7 | Leave-one-stream-out: refit dropping each transcript stream and each metabolite stream, reporting the change per target | C3 | Identifies load-bearing single points of failure, and streams that never informed anything |
| **F8** | **Rank-ECDF difference plots with simultaneous bands**, per target, at the largest affordable replication count | C4 | Primary calibration evidence. Histograms are secondary |
| F9 | Coverage curve: nominal against empirical at 50, 80, 90, 95 and 99 percent, per target, with binomial intervals | C4 | The two controls are the diagnostic rows, with the kill-criterion band drawn on |
| **F10** | **Reference comparison** on the two-gene system: exact blocked posterior against production posterior, overlaid marginals plus a two-dimensional contour, plus the divergence between them | C4 — **without this a calibration pass is unattributable** | The number that makes every other calibration claim mean something |
| **F11** | **The architecture figure**, rebuilt after D13 removed the protocol comparison. Three panels: the hybrid composition running a full cycle at all, which no prior version of this framework could do; the machine-generated enumeration of every departure from the published model; and the completeness assertion of phase 13 shown both passing on the assembled model and failing on a deliberately incomplete one | **C5** | C5's evidence is now that the composition executes and audits itself, not that two interfaces calibrate differently. Weaker as a comparison, and honest: with no shared parameter there is no second protocol to compare against |
| F13 | Jacobian singular-value spectrum for the target set, with and without the polymerase constant freed | C3 | The multiplicative ridge as a singular value dropping by orders. Uses the existing `check_identifiability` |
| F14 | Wall-clock per trajectory against horizon and gene count, the budget per rung, and the achieved replication count | K1, reproducibility | The scoping note's binding practical number |
| **T1** | **Provenance table**, generated by the loader and `reduction_declarations` — every parameter with its value, prior, width, informedness, source file and rejected alternatives | Honesty, C2 | Machine-generated, never typed. The audit trail for the cross-file errors of record |
| **T2** | **Reduction declarations with measured costs** — the lumped charging step and its now-deterministic formalism, `k_chg`, the tRNA pool size and its charged fraction, the chemostatted pools, exogenous membrane growth, the corrected mapping, the rounding policy, the drain granularity, the lactate volume ratio, the capped rate-law geometry, any smoothed counter | Scope honesty | Each row carries the *measured* cost where measured. This table is what stops a result depending on our choice unremarked |
| **T3** | **Kill-criteria scoreboard** — each criterion, its threshold, its measured value, pass or fail | Falsifiability | **Published whether or not everything passed.** A scoreboard with a fail on it is more credible than one without |
| T4 | Target-set summary: truth, posterior mean, 90% interval, coverage, shrinkage | C3, C4 | The conventional recovery table |

**Scalars that must appear in any abstract:** wall-clock per trajectory; the
reverse-channel gain range; the achieved replication count; 95% coverage per
target; and the shrinkage of the ptsG promoter under metabolite-only data
together with the tight-prior control's shrinkage under transcript-only data.
Those last two are what "uncertainty crosses the boundary" reduces to once D13
establishes there is no shared parameter to point at.

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
- **The cut and iterative exchange protocols, and any comparison against them.**
  Both decide what to hand between blocks by matching parameter names, and D13
  shows that set is empty for Core A′, so neither passes anything: the cut
  degenerates into inferring each block alone and the iterative loop converges at
  the first pass having moved nothing. This is a structural mismatch, not a
  tuning problem. The iterative protocol is additionally a proof-of-concept
  rather than a production route, with a known prior-double-counting limitation.
  What either would require if a shared parameter ever appeared is written up in
  `dev/notes/modular-bayesian-inference-heterogeneous-modules.md`.
- **A shared-parameter treatment of the boundary.** Accepted as absent for this
  stage (D13). Unfreezing cell volume on the stochastic side would create one and
  is the named extension; it is not attempted here.
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
the assembled model. ~~Concretely, above roughly 10 s per full-cycle trajectory on
the simulation-based path.~~ **Amended 2026-09-24:** the 10 s proxy is retired.
It was never derived from an evaluation count, and replications are independent,
so a cluster runs them in parallel; a 60 s trajectory is affordable. What binds is
the **cost of one posterior**: sequential evaluations per posterior times
per-trajectory wall-clock, since D10's alternating updates do not parallelise
within a chain. K1 is scored in phase 16, once task 16.2 fixes the path-update
mechanism and so the evaluation count, from task 13.7's measured wall-clock.
**No numeric bound is set yet, deliberately** (2026-09-24): the first inference
on the minimal model comes first, and the bound is set from where that lands,
before phase 16 scores K1. T3 records that it was set after the measurement.
**D13 loosened this considerably** by removing about
350 times the stochastic block's event count, which is the largest single change
to the budget in this spec.
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
of prior draws clip (check 7). **Scored at the published 1 s drain** (amended
2026-09-05): the fraction is per drain rather than per handshake, and a coarse
`drain_interval` is a labelled reduction whose census is about our choice, not
the published model's.
*If it fires:* the ODE block is not differentiable where it matters,
gradient-based sampling is invalid as posed, and the smoothed model must be built
and validated against the clipped one before any posterior is reported.
~~**This is the criterion most likely to actually fire, and it fires in phase 3,
which is the good news.**~~ **Amended 2026-09-05:** it is scored on the
assembled model in phase 14, not on phase 3's toy, whose drain-to-pool ratio is
three orders of magnitude gentler. See §12. D13 did not reduce this risk, it relocated it: the
charged-tRNA counter debits ~553 residues per second against a pool of order 10³
particles, so it can clip on ordinary Poisson fluctuation in translation firing.
The non-smooth boundary now sits on a pool whose size we assert (D14), which is
either an opportunity — size it so it provably cannot clip — or a new way to
fail.

**K6 — structural non-identifiability of the target set.**
*Fires when* `check_identifiability` reports a rank below six for the chosen set,
or a condition number above 1e6. Computed *before* any sampling, at the cost of
one Jacobian.
*Predicted instance:* freeing the polymerase constant alongside the promoters
drops one singular value by orders (D11). Reparameterise to the product plus a
scale anchor rather than sample through the ridge. The prediction is analytic
rather than empirical: D10's closed-form conditional is conjugate in exactly that
product, so the ridge is visible before the Jacobian is computed.

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

- **RESOLVED by D13, and recorded rather than deleted because it was the
  spec's highest-value open question.** The obstruction was the lumped charging
  step's two reactants. With that step in the ODE block, every remaining reaction
  in the stochastic block — transcription, translation, decay, translocation and
  ~~protein degradation if phase 11 keeps it~~ (**amended 2026-09-23:** there
  is none; task 11.9, §12) — is zeroth or first order conditional
  on the 60 s rate constants, so the transition kernel is closed form and factorises over genes,
  and the promoter conditional is conjugate with the polymerase constant fixed
  (D10, D11). **Two residual caveats, neither a tractability problem.** The cost
  counters are deterministic accumulators entering no propensity, so they do not
  break the kernel's factorisation; what they do is restrict the ODE block to
  aggregate linear combinations of the 17 promoter strengths, which is D9's
  identifiability finding. And the factorisation is over genes *given the path*:
  the ODE conditional does not factorise. Verified by construction in the
  transcription and translation phases rather than by a determination task.
- ~~**[NEEDS CLARIFICATION: what handshake granularity do the pools require?]**
  Settled by the 1 s / 5 s / 60 s comparison in D10 — in miniature in phase 4,
  at full scale in phase 13. Expected to rule out 60 s, given the GTP pool turns
  over in half the rebuild interval.~~ — **resolved 2026-09-25 by task 13.4:
  the published 1 s drain, and nothing coarser.** Over 50 paired seeds, the
  largest final paired difference is **−52.9% ± 9.6%** for a 5 s drain
  (`M_3pg_c`) and **+288% ± 29%** for 60 s (`M_pi_c`), against D10's one
  percent (job 17298253). D10 expected 5 s to pass. It does not.
- ~~**[NEEDS CLARIFICATION: is the charged-tRNA channel stronger or weaker than
  transcription's?]** The scoping note hopes it is stronger and would rescue the
  reverse direction. A first pass at the published rate law suggests the opposite
  — that both denominators are dominated by polymer length, bounding the whole
  reverse direction at a few percent structurally rather than accidentally.
  Settled in the translation phase, and it changes how the result is framed.
  **D13 made this question well posed.** With charging in the stochastic block
  the charged pool was an intra-block read, not a channel at all; now it is a
  genuine second reverse channel alongside the nucleotide pools. Measuring it
  needs the charging module composed, or a double — see R1.~~ — **resolved
  2026-09-24 by task 11.10: stronger.** The translation constant's elasticity to
  the charged pool is **0.0909 to 0.0911**, about 1.9× transcription's four-pool
  0.044 to 0.051 (job 17131232). The first pass's "few percent" does not hold.
  That number rests on the lumped pool's per-amino-acid share (amendment
  2026-09-24 A). Reading the whole pool instead would give 0.005.
- ~~**[NEEDS CLARIFICATION: does protein fold change over a cycle reach two?]** A
  first pass suggests it may fall well short, which would indicate the lumped
  charging step throttling translation. Settled in phase 11, and far cheaper to
  diagnose there than to explain at the end.~~ — **measured 2026-09-24 by task
  11.10: no, and whether charging is the cause depends on the composition.**
  - With the constants held at the nominal pool (jump-only) the median is
    **1.691**. Charging cannot be the cause there. The ribosome constant is
    the published 12. Missing replication is the untested leading candidate.
  - In the full seven-module hybrid the median is **1.497** (before phase
    13a credited GDP and Pi, so 13b remeasures it). There the charged
    pool empties at rebuild instants and the constants fall about 20×, so
    charging plausibly contributes. How much is not measured.
  - **Remeasured 2026-09-25 by task 13.7** on the assembled model, every
    counter and credit wired: median **1.627** (min 1.215, max 2.193, one
    seed, job 17298242). The charged pool is under one particle at none of
    the 105 rebuild instants, so phase 11's clip-meets-rebuild mechanism is
    absent from this run.
  - Open for phase 14's F5, with §3's [1.7, 2.3] unchanged (amendment
    2026-09-24 F).
- **Does the lumped charging step's stoichiometry change the answer?** The
  published AMP-and-pyrophosphate form puts 553 events per second through the
  adenylate kinase; a two-ATP-to-two-ADP form puts zero. That changes how much
  the kinase's equilibrium constant can move the ATP/ADP ratio the rebuild reads,
  and so whether the kinase belongs in the target set. One cheap comparison run,
  in the charging phase. **After D14 it must be run at matched steady flux**
  rather than matched event count, because both forms are now mass-action rate
  laws with different self-limiting behaviour. **Partly answered 2026-09-23 by
  task 9.9; the ATP/ADP half is carried to phase 14 (task 14.9).** See the
  entry below.
- **Does the 60 s rebuild stay piecewise-constant?** Declarable, defaults to
  published. ~~The deviation is measured rather than argued (R11).~~
  **Amended 2026-09-05:** what is measured is the *cadence*, at a shorter
  piecewise-constant interval; a continuous cadence is refused at build time by
  phase 4 and so is labelled but not measured. See §12.
- **Is the promoter proxy inherited verbatim at Step 2?** Barely matters for
  synthetic recovery, where it only sets the truth we recover. Leaning toward
  inheriting it as the baseline the surrogate is measured against.
- ~~**[NEEDS CLARIFICATION: how large is the tRNA pool, and can the charged
  counter be made unable to clip?]** New with D14. The registry records no value
  for either tRNA species, and the pool size sets three separate things: the
  charging flux at fixed `k_chg`, the buffer that decides whether check 7's
  counter clips, and the lag on the dominant forward channel. One number, three
  consequences, and all three are ours. Check 1b bounds it from below and check 7
  from above; whether a value satisfies both is not yet known.~~ — **partly
  resolved 2026-09-23 by task 9.8** (`dev/scripts/trna_charging_diagnostics_result.md`,
  job 17044321), **and its check 7 half 2026-09-24 by task 11.7**: with
  translation consuming, 0.25 mM clips on none of a cycle's 6,300 drains and
  0.125 mM on 7, so **0.25 mM satisfies both checks** in phase 9's composition
  (`dev/scripts/translation_diagnostics_result.md`, job 17131232). The default is **0.25 mM at 0.8 charged**. At that fraction,
  check 1b's floor is ≈ 0.124 mM analytically (500 / (0.2 · 20,180.39)); the
  scan's uncharged minimum is 404 particles at 0.1 mM and 505 at 0.125 mM.
  0.25 mM clears it with 1,009 uncharged and 4,000 charged particles (cycle-end
  minima). **Check 7 cannot be decided before phase 11** (task 9.8's check 7
  half is deferred to 11.7). At cycle end the charged pool buffers **7.23 s** of
  demand at 0.25 mM, or 5.4 of the longest protein's worth of residues. Neither
  check gives an upper bound. What grows with the pool is the lag on the
  adenylate forward channel. Linearised, it is ≈ 6.0 s per mM (1.51 s at 0.25 mM);
  a 10% demand step measures **1.43 s**, against the 1 s handshake. At fixed
  `k_chg` the pool sets the flux sublinearly, from 56 /s at 0.025 mM to 1,873 /s
  at 1.0 mM, 33× for a 40× pool, because ATP falls as the flux rises.
- **Stoichiometry, partly answered by task 9.9 (2026-09-23, job 17044321).**
  At steady fluxes matched to 1.5e-9 (548.17 /s each; the two-ATP `k2`
  rescaled by 1.044), the published AMP + PPi form sends **548 /s** through the
  adenylate kinase and the two-ATP form sends **0**, and pyrophosphate settles at
  0.367 against 0.024 mM. Those are answered. **The ATP/ADP half is not.** Phase
  8's glycolytic double rephosphorylates ADP at a rate effectively first order
  in ADP alone (its phosphate factor is saturated), and both forms return two
  ADP-equivalents per event. So steady ADP is pinned by the double: 0.4317 mM
  in both runs. The measured 0.035% ATP/ADP difference is AMP moving within a
  fixed adenylate total, not the ratio responding. Nothing here says how far
  the kinase's equilibrium constant can move ATP/ADP under live glycolysis,
  whose PGK and PYK laws are reversible and read ATP. **Carried to phase 14
  (task 14.9):** rerun the comparison on the assembled model. (A first version
  also compared the forms at fluxes 0.86% apart and reported 1.02%, which was
  flux mismatch.)
- ~~**Whether the full-cycle checks belong in the default test suite.** They are
  the strongest evidence and the slowest thing added. Deciding needs the measured
  wall-clock.~~ — **partly resolved 2026-09-10 by phase 6's measurement, and left
  open for the assembled model.** One module's full-cycle checks cost six
  6,300 s stiff solves and about 1m40 to 2m10 of a 4m56 suite, which is
  affordable, so phase 6's stay in the default suite. That is a 13-state ODE
  block with the currencies held; the assembled model is a hybrid running ~6,300
  handshakes, and task 13.7's wall-clock is what decides for it. **Shortening the
  interval is not an option** — it is precisely the error of record.
  **Phase 7 adds a second measurement that splits the cost in two** (job
  16364852): its 6,300 s integration is **4.56 ms**, and the **19.67 s** around
  it is compiling the stiff-solver path, which a suite pays once for its first
  stiff solve whatever the horizon. So on a small ODE module the *horizon* is
  nearly free and what costs is compilation — which is why phase 6 can afford
  six such solves for barely more than one. Neither measurement touches the
  assembled model, where the cost is ~6,300 hook iterations and an SSA rather
  than one `solve`.

---

## 10. Risks

Each risk is paired with the signal that arrives *before* it costs weeks. Several
arrive in phase 3 or in the module phases rather than at the end, which is the
point of D0's reordering.

| # | Risk | Early warning | Pre-committed response |
|---|---|---|---|
| R1 | The reverse channel is too weak for the bidirectional claim | **Measure the charged-tRNA elasticity as soon as the translation rate law exists.** Below 0.02 and the reverse direction is bounded at a few percent by the structure of both published rate laws. **D13 delayed this warning and that is a real cost**: the charged pool is now an ODE state, so a jump module cannot produce it standalone. Either the translation phase carries a charged-pool double or the measurement waits for assembly. Carry the double | K2's response: reframe as a measured bound plus the replicate count needed, and promote Step 1b |
| R2 | The non-smooth drain operates in-regime, so the ODE block is not differentiable where it matters | **Any non-zero deficit carried on the first full-cycle run at published parameters** (check 7, ~~phase 3~~ **phase 14 — the assembled model; phase 3's toy exercises the mechanism only, see §12**). Do not wait to meet it as sampler divergences, which look like a step-size problem | Smooth the drain; validate the smoothed model against the clipped one. K5 |
| R3 | Integer rounding drift swamps the integrator and is read as a conservation leak | ~~The adenylate residual failing to shrink when tolerances tighten tenfold.~~ **Amended 2026-09-11 (see §12): a flat adenylate residual is the *healthy* signature for a moiety passing §3's exception, so it no longer discriminates.** The warning is now **the per-evaluation residual `|Σ nᵢ·duᵢ|` exceeding one ulp of the conserved sum**, which rounding drift does and exact conservation does not — and, at composition scope, the residual growing with *handshake count* rather than with tolerance (check 0 / output F4). **Amended 2026-09-25 (phase 14a):** read "one ulp" as `n` ulps for an `n`-term moiety, as §3's second gate now does. And at composition scope, growth *at roundoff* is expected: fractional carry's phosphate residual grows 9.4× over the cycle while staying below 10⁻⁸ of `tol_C`. The warning is growth to whole particles, which is what a rejected rounding policy does | Fractional-carry rounding (check 0). The magnitudes differ by ~140×, so this is not hypothetical |
| R4 | The reference system is skipped because production looks fine | **A calibration pass reported with no reference comparison beside it.** This is a process risk, so the warning is procedural | F10 is a phase deliverable, not an optional extra. A calibration claim without it is unattributable |
| R5 | ~~The message family's finite support truncates every hand-off~~ | **Retired by D13.** This risk was entirely about the density-estimate message family in the cut and iterative protocols, which pass nothing here and are §7 non-goals. It becomes live again only if the shared-parameter extension is taken, and the survey note carries it | None needed. Recorded rather than deleted so the extension inherits it |
| R6 | The six targets are not identifiable as a set | **A rank below six or a condition number above 1e6, computed before any sampling.** Costs one Jacobian | Reparameterise to the product plus an anchor. K6 |
| R7 | Calibration is unaffordable, so the claim cannot be made | ~~The phase 13 wall-clock exceeding 10 s per trajectory~~ **Amended 2026-09-24: the projected cost of one posterior** — task 13.7's per-trajectory wall-clock times the sequential evaluation count of the mechanism task 16.2 chooses — too long to run 50 replications side by side on the cluster; the bound is deferred until the first minimal-model inference lands (§8 K1) | K1's ladder, pre-committed so it is not negotiated under pressure |
| R8 | The lumped charging step — ours, not the model's — determines the answer. **Upgraded by D13** | Two signals. The stoichiometry comparison shifting any target's posterior mean by more than half a posterior standard deviation. And, near-certainly, **the charged-tRNA pool size appearing in the top three by sensitivity in F2**. This is no longer only about the reverse channel: after D13 the pool also gates the *dominant forward* channel, since ~81% of ATP turnover is charging and its flux now contains no stochastic-block quantity. A quantity we assert sits between the stochastic block and the adenylate pools | Report both channel gains as *functions* of the pool size rather than point values; make the pool size a first-class asserted quantity in T1 and T2, not a detail of `k_chg`; and treat check 1b as its lower bound ~~and check 7 as its upper~~ — **amended 2026-09-23:** check 7 bounds it from below too, since a larger pool buffers the counter; no check bounds it from above, and what grows with it is the forward-channel lag (task 9.8) |
| R9 | Core A′ is better-conditioned than the same subnetwork inside the full model, so success over-claims | **Measurable now, without waiting:** the posterior width on the pyruvate-kinase constant with the seven dropped sibling reactions stubbed in as a competing sink | Quote the factor alongside every recovery claim. K7 |
| R10 | Everything rests on synthetic data, so misspecification is untested by construction | **No internal early warning exists, and that is why this is a stated scope limit rather than a mitigated risk.** The nearest signal is F5: if either external check fails, the model is already misspecified against the data that does exist | State the limit in the abstract. Do not claim robustness that was not tested |
| R11 | The 60 s rebuild is an implementation artefact the posterior depends on | **Nominal trajectory at 60 s against 6 s rebuild cadence; any observed pool differing by more than one percent.** With the GTP pool turning over in half the interval, expect this to fire. **Amended 2026-09-05:** read "nominal trajectory" as D10 now reads it — a paired ensemble at aligned instants — because a rebuild also triggers the aggregation rebuild and so changes randomness consumption | A labelled decision in T2 with the posterior sensitivity reported. Turns an open question into a measurement |
| R12 | Protein fold change misses two, so the reduction throttles translation | **The arithmetic, run as soon as the translation rate law exists.** A first pass puts one gene near 1.3 against a published median near 2. Inherits R1's sequencing cost: the charged pool is an ODE state, so this needs the charging module or a double | Diagnose in the translation phase — the charging step, `k_chg`, the pool size or the ribosome constant — rather than explain it at recovery |
| R13 | The observation model's scale is wrong and every credible interval is misreported | **No noise-scale recovery test exists anywhere in `test/`.** A variance-versus-standard-deviation misreading of the multivariate normal constructor rescales every interval and passes every other check | Check 9, written before the first posterior |
| R14 | Summary statistics discard the low-copy information the reduction exists to exercise | **The ensemble observer averages replicates and the summaries keep per-time means only** — both true in the code today, both contradicting the scoping note's own recommendation | Carry transcripts as counts per replicate; add distributional summaries and a lag-one autocovariance, or supersede summaries with the exact likelihood |
| R15 | The seven module phases are written against a protocol that then changes again | **A framework need surfacing during a module phase.** That is an interface bug, and the wave plan's rule stands: it is a conversation on the main branch, not an edit on a module branch | Amend §12 with the change and its date, then land it as its own phase before the affected modules continue |

---

## 11. Task list

Seventeen phases plus a phase 0, a phase 5b, a phase 10b, a phase 11a, a phase 13a and a phase 14b, each one reviewable pull request. Ordering is
D0's: the two protocol changes, then the kill phase on a toy, then the drivers,
then the modules, then assembly, validation and inference.

**Parallelism.** Phase 1 runs alone. Then a four-way fan-out: the ODE track of
phases 6 to 9 runs concurrently with the framework track of phases 2 to 5. The
ODE track is **not** itself a clean fan-out after D13: phases 6, 7 and 8 are
mutually independent, but phase 9 depends on phase 8, because charging
contributes to the pools recycling owns. The stochastic track needs only phase
2, so phase 10 runs alongside the ODE track rather than after it; within that
track, phase 10 comes before phases 11 and 12, because both read the transcripts
phase 10 owns. Phase
11 also wants phase 9 landed, since its rate constant reads an ODE state; it can
proceed on a double, and R1 records that cost. Phases 13 to 17 are strictly
serial, because each validation check catches errors the next would mask.

**Phase 5b gates the fan-out, and that is the point of it.** Phase 7 needs the
inbound volume channel and R15 forbids building it on a module branch, so a
concurrent phase 7 would stall on framework code its author may not touch. 5b is
small and independent, so it costs one pull request to remove the one known
framework dependency from the module phases — after which 6, 7, 8 and 10 are
four genuinely independent branches. Three of their done-when clauses still are
not: task 8.4 needs phase 6 for the shared-enzyme assertion, task 8.5 composes
all four ODE modules, and task 10.7 checks copy numbers against the metabolic
modules'. Those belong to whichever pull request lands last, or to phase 13, and
carrying them as skipped tests naming the blocking phase is what keeps four
branches from stubbing each other.

### Phase 0 — Retire OpenSpec, keep its findings

**Goal:** make this spec the only live spec system without losing a derived
number.
**Done when:** every finding in the four drafted changes appears in this spec or
in a note, and `openspec/` is read-only history under `dev/archive/`.
**PR:** #41 (merged 2026-09-03)

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

**Goal:** make an ~~outbound~~ mass or currency edge execute, so a shared pool can
receive terms from every module that produces **into or draws from** it
(amended 2026-09-03; see §12).
**Done when:** a two-module ODE composition in which module A produces a species
module B owns shows B's state rising at A's declared rate, and every ~~outbound~~
mass or currency edge in the composition **on a species the declaring module
does not own, in either direction,** maps one-to-one onto an executed
contribution (amended 2026-09-03; see §12).
**PR:** #42 (merged 2026-09-03 UTC; 2026-09-04 AEST)

- [x] 1.1 Choose and record the mechanism — a second return value from
  `dynamics`, a separate protocol function, or a global-index write buffer —
  verify by the pull request recording what each costs at 32 states and by the
  chosen one adding no dependency that is absent from the Julia depot.
- [x] 1.2 Extend the protocol with the chosen channel, defaulting to empty —
  verify by `build_problem(TranscriptionTranslation())` producing byte-identical
  initial conditions, parameters and trajectory to the pre-change run, and by all
  829 existing tests passing unchanged.
- [x] 1.3 Accumulate contributions into the global derivative by
  registry-resolved index in `_build_rhs` — verify by a two-module test asserting
  the owner's derivative equals its own term plus the contributor's exactly, and
  that the contributor's own slice is unchanged.
- [x] 1.4 Reject a contribution to a chemostatted species, to a species no module
  owns, and to a name outside the registry — verify by three tests each throwing
  and naming both the species and the contributing module.
- [x] 1.5 Cross-check contributions against declared edges in both directions —
  verify by a test that deleting a contribution while keeping its outbound edge
  throws naming both, and that adding a contribution with no edge throws
  likewise; the drift check must be live, as the untyped-inputs check already is.
  **Amended 2026-09-03:** "outbound edge" reads "mass or currency edge in either
  direction on a species the module does not own" — a consumer of a foreign
  pool owes a negative term exactly as a producer owes a positive one (§12).
- [x] 1.6 Re-take the state-container decision at 32 states — the problem is
  currently built out-of-place over a static vector — verify by benchmarking
  static-out-of-place against mutable-in-place at 32 states and asserting the
  chosen path allocates nothing per right-hand-side call.
- [x] 1.7 Run the suite and hand off — verify by `sbatch test/run_tests.slurm`
  passing above 829 by the tests added, and by `docs/handoff.md` recording that
  outbound edges now execute rather than merely resolve. **Done:** job 16171912,
  932/932 (829 before the phase); the handoff records that mass and currency
  edges execute in both directions.

### Phase 2 — Several jump modules compose

**Goal:** make `_build_jump_problem` honour the index contexts it already
computes, so a jump module writes only its own global slice and can read and
write a declared peer's state.
**Done when:** two distinct jump modules — one owning a species, the other
incrementing and decrementing it — compose and produce a trajectory whose
statistics match a hand-written single-module equivalent within Monte Carlo
error.
**PR:** #43 (merged 2026-09-03 UTC; 2026-09-04 AEST)

- [x] 2.1 Pin the bug before fixing it — verify by a test composing two doubles
  with disjoint state names and asserting that *today* the second module's effect
  mutates the first module's state, so the fix has a witness rather than a claim.
- [x] 2.2 Give `reactions` the module's index context — verify by instrumenting
  each double's rate and effect functions and asserting every read and write
  falls inside that module's own state range. **Annotated 2026-09-04:** the
  context reaches the module as *views* rather than as indices — `reactions`
  returns `Reaction(rate, affect!)` in local coordinates and the orchestrator
  wraps them, mirroring `dynamics` — and the instrumented assertion is "own
  slice or a declared input", since a declared peer read or write is the point
  of tasks 2.4 and 2.5. Not an amendment; the mechanism choice and its rejected
  alternative are recorded in PR #43.
- [x] 2.3 Use the contexts in `_build_jump_problem`, which currently computes and
  discards them — verify by asserting the global initial-condition vector's
  length equals the sum of per-module state counts and that each module's counts
  sit at its own offsets.
- [x] 2.4 Let a jump module read a declared peer's state in a propensity — verify
  by a propensity proportional to the peer's count and a test where doubling the
  peer's initial count doubles the measured firing rate over many replicates.
- [x] 2.5 Let a jump module write a declared peer's state, gated on ~~an outbound
  edge~~ **a `written_states` declaration** — verify by a decay-shaped double
  decrementing a transcript it does not own, and by the same write throwing once
  the ~~edge~~ **declaration** is removed. **Amended 2026-09-04:** no edge kind
  can name a non-registry state (**since phase 11a, a `CatalyticEdge` can**,
  which writes nothing and so does not bear on this gate) and the transcripts
  are non-registry by task 10.2, so the gate is the declaration; where the written state *is* a registry
  species the resolver additionally holds the declaration to a mass or currency
  edge (§12).
- [x] 2.6 Preserve single-module behaviour exactly — verify by the existing
  stochastic and bursty gene-expression tests passing unchanged including the
  Fano-factor assertion, and by rewriting the "composition fails" test to assert
  what now actually fails (duplicate ownership) rather than what no longer does.
- [x] 2.7 Suite and handoff — verify by `sbatch test/run_tests.slurm` passing and
  by the handoff recording that the three gene-expression modules may now assume
  composability, correcting the drafted design that assumed it already held. **Done:** job 16185576, 1006/1006 (932 before the phase; 1002 before the review fixes); the handoff records the `Reaction` contract,
  `written_states`, and what phases 10–12 may now assume.

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
**PR:** #44 (merged 2026-09-05)

- [x] 3.1 Choose and record the composition mechanism — an outer split-operator
  loop over stepped integrators, a discrete callback on a single problem, or a
  jump process over an ODE problem — verify by the pull request recording, for
  each, whether the published piecewise-constant 1 s handshake is faithfully
  representable, and by the chosen one adding no dependency absent from the
  depot. **Done:** the outer loop, mirroring the published `hookSimulation` at
  `delt = 1.0` s. A `JumpProblem` over an `ODEProblem` makes propensities track
  the pools continuously, so phase 4 task 4.3 would have nothing to be true
  about; a discrete callback cannot advance the SSA half at all. `init`,
  `step!` and `reset_aggregated_jumps!` come from SciMLBase and JumpProcesses,
  both already direct dependencies, so `Project.toml` is unchanged.
- [x] 3.2 Replace the mixed-formalism refusal at `src/orchestrator.jl:27` with a
  hybrid build path — verify by a mixed composition returning a driver rather
  than throwing, and by the same composition *still* throwing with a named error
  when a required exchange declaration is absent. **Done:** `build_problem`
  returns a `HandshakeDriver`; the named refusals are asserted (a `param_slot`
  that is not the declaring module's own free parameter, a catalytic species no
  jump module owns, a debited pool no ODE module integrates, a counter no jump
  module owns, and a policy keyword on a homogeneous composition). Two refusals
  that did not exist before: a homogeneous `:sde` composition, which used to be
  reported as mixed (§2 G1); a deferred channel declared only by the pool's
  owner, or a catalytic edge declared only by the module owning the count,
  either of which would be declared and never executed; a counter that is a
  registry species, or an ODE module's own state; more than one pool credited
  from one counter, which would create matter; and **a state name owned in both
  blocks**, which
  would make every crossing that mentions it resolve to whichever block was
  asked first — phase 2's aliasing bug arriving *between* blocks rather than
  inside one. `_check_state_ownership` does span both blocks — it iterates every
  model regardless of formalism — but skips any name the registry does not know,
  and the transcripts phase 10 owns are non-registry by task 10.2, so that is
  exactly the case it misses. **Annotated 2026-09-05:** this is a new refusal,
  not a changed one, so it is recorded here and in the pull request rather than
  as a §12 amendment; nothing the spec said became false. Phases 6 to 12 inherit
  it, along with two more the same phase added: a module may not name a state
  another block owns, may not reach across the boundary through `inputs()` or
  `written_states()`, and may not fill a catalytic `param_slot` whose name a
  second ODE module also declares.
- [x] 3.3 Implement the count-to-concentration conversion through the registry at
  fixed volume with an explicit rounding policy — verify by check 0: exact-zero
  round-trip residual under fractional carry, and by the policy being a stated
  field on the driver rather than implicit in the arithmetic. **Done:** the
  factor is derived from Avogadro and the volume of a sphere and returns 20,180
  particles per mM at 200 nm, reproducing the scoping note rather than
  transcribing it. Check 0's round trip is exactly zero over 6,300 handshakes
  under fractional carry (`==`, not `≈`). All three signatures measured: carry
  drifts **0.0** particles at both 630 and 6,300 hooks; deterministic drifts
  **252.0 → 2,520.0**, exactly linear; stochastic RMS over 60 seeds goes
  **11.63 → 41.12**, a ratio of 3.53 against √10 = 3.16. Deterministic rounding
  is now rejected on measured evidence rather than on assertion.
- [x] 3.4 Isolate exchange error from integration error — verify by a toy whose
  ODE has zero derivative: after 600 handshakes every count is unchanged
  *exactly*, so any drift is attributable to the exchange alone. **Done:** every
  count and every pool is bitwise unchanged across all 600 handshakes, not only
  at the end, and the carried remainders stay below 1e-6. A second test closes
  the gap this one leaves open — with a zero derivative nothing would catch a
  `step!` that moved the clock without solving — by comparing a handshake's ODE
  step against an independent solve of the same problem over the same second.
- [x] 3.5 Debit the deferred counters under the published clamped policy — verify
  by a test where the cost exceeds the pool: the pool floors at zero, the carried
  deficit equals the shortfall exactly, and the next step debits it; and by the
  existing gradient report naming the edge when a differentiable sub-model is
  present. **Done:** a 250-particle cost against a 100-particle pool floors the
  pool at exactly 0.0 and carries a deficit of exactly 150.0, which the next
  handshake repays in full against a refilled pool. `check_gradient_safety`
  returns the one obstruction and names `M_atp_c`.
- [x] 3.6 Execute the enzyme-concentration channel — protein counts overwriting
  the ODE's enzyme concentration through a catalytic edge's parameter slot —
  verify by a test that a step change in the toy's protein count changes the ODE
  flux by the ratio of counts and changes nothing else. **Done:** 400 and 800
  copies fill the slot with exactly `n / 20180.39` mM, the flux ratio is exactly
  2.0 at a common state, and every other parameter is untouched. `param_slot`
  had been a required field no code read; it now has semantics and is resolved
  at build time.
- [x] 3.7 Run check 7's clipping census on the toy — verify by zero carried
  deficits at nominal parameters and by the fraction of 200 prior draws that clip
  being recorded, ~~since K5 fires here or nowhere~~ **(amended 2026-09-05: the
  census on a toy does not decide K5; see §12)**. **Done: K5 does not fire on the toy.**
  Zero carried deficits over 600 handshakes at nominal parameters, and **0 of
  200** prior draws clip at any handshake (0.0% against the 5% threshold). Read
  narrowly: the toy's pool is ~74,000 particles against a ~40 particle/s drain,
  so this bounds the *mechanism*, not Core A′'s charged-tRNA counter, which
  D13/D14 put at ~553 residues/s against a pool of order 10³ and which is
  re-scored in phase 14.
- [x] 3.8 Measure and record wall-clock — verify by a Slurm run reporting seconds
  per simulated second and the 6,300 s extrapolation, an explicit pass or fail
  against the phase's stated budget, and the number written into `docs/handoff.md`
  and into the scoping note's open-questions list, which carries it as unmeasured
  today. **Done:** job 16212707 on gina1, in
  `dev/scripts/bench_handshake_result.md`. 1.50e-06 s per simulated second with
  the pool frozen and 1.57e-06 s with the rate law live at a 600 s horizon. The
  frozen rows are steady across all three horizons (1.50–1.51e-06); the live
  rows run 2.16e-06 at 60 s and settle to 1.57e-06 by 600 s, so the **worst**
  row extrapolates to **0.014 s per 6,300 s trajectory against a 10 s budget —
  PASS** (0.010 s at the 600 s horizon). Recorded as a *floor*: the toy is one gene and two
  metabolite states against seventeen and thirty-two, so K1 is decided by task
  13.7, not here.
- [x] 3.9 Suite, count and verdict — verify by `sbatch test/run_tests.slurm`
  passing and by the kill-criterion verdict written down as a sentence rather
  than implied by the pull request merging. **The verdict: the project
  continues.** The 1 s handshake is expressible inside the sub-model protocol
  with no per-module special case — the driver reads only `formalism`,
  `coupling`, `states` and `parameters`, and neither toy module knows it is
  being handshaked — and the measured cost is three orders of magnitude inside
  the budget. Neither of the phase's two stated kill conditions fired, and K5
  did not fire on the toy. **Suite:** job 16213224, 1126/1126 (1006 before the phase); the handoff records the mechanism, the separate state vectors, the measured wall-clock and the cross-block ownership rule.

### Phase 4 — The 60 s rebuild

**Goal:** execute the reverse direction — live pools to recomputed jump rate
constants, piecewise-constant at the interval the edge declares. There are two
such channels in the assembled model, the nucleotide pools into transcription and
the charged-tRNA pool into translation; the toy exercises the mechanism, which is
shared.
**Done when:** the toy's propensities change only at interval boundaries, the
interval is read from the declared edge rather than hard-coded, and the measured
elasticity of a rate constant to its upstream pool is reported.
**PR:** #45 (merged 2026-09-05)

- [x] 4.1 Decide where mutable rate constants live — in the flat parameter vector
  the builder already returns as mutable, or as mutable state on the sub-model as
  the drafted transcription design chose — verify by the pull request recording
  the consequence for the parameter-substitution call every inference path uses,
  since a constant held on the struct is invisible to it. **Done:** the parameter
  vector, written by the hook exactly as phase 3's catalytic channel writes the
  ODE block's. A constant held on the struct is invisible to `remake(prob; p = θ)`
  *and* to `model_free_params`, so it could never appear in a posterior, in T1 or
  in K6's Jacobian; that is recorded in the pull request. Which slots and what
  fills them is a module-level declaration — `rebuilt_params(m)` and
  `rate_constants(p, t, m, pools)`, the rate-constant twins of phase 1's
  `contributed_states`/`contributions` — because `RateConstantEdge` carries no
  `param_slot` and one pool feeds many constants anyway. The pools arrive as an
  argument rather than through `inputs()`, which phase 3 forbids across the
  boundary. Mechanism recorded in the pull request, as tasks 2.2 and 3.1 recorded
  theirs; no §12 amendment, because nothing the spec said became false.
  **One consequence recorded rather than fixed here:** a rebuilt slot is still
  an entry of `model_free_params`, so `build_turing_model` and the ABC path
  sample it and the hook discards the draw at the first refresh. At seventeen
  transcription constants that is seventeen posterior dimensions coming back
  shaped like their priors, which must not be read as an identifiability
  result. Excluding derived names from the sampled set belongs to phases 15 to
  17, which are the first to compose a driver with `infer`; the consequence is
  documented on `rebuilt_params` where it will be met.
- [x] 4.2 Schedule the rebuild from the edge's declared interval — verify by a
  test that a 60 s edge refreshes ten times in 600 s and a 30 s edge twenty
  times, with the count read from the trajectory rather than asserted. **Done:**
  exactly **10** and **20**, counted by the handshakes at which the recorded jump
  parameter vector's rebuilt slot changed — `run_handshake!` now records that
  vector for the purpose. The value written at a refresh equals the module's own
  law at that handshake's recorded pool, bitwise. **Ten refusals**, each naming
  what is wrong: a rebuilt name that is not one of the module's own free
  parameters; a name another module in **either** block also declares, since
  across the boundary the two parameter vectors are separate and the hook would
  rewrite one copy and leave the other — which `_validate_shared_params` cannot
  catch, because the declared values agree; a pool no ODE module integrates;
  rebuilt parameters with no edge; an inbound edge with nothing to rebuild; an
  ODE module's outbound edge no jump module consumes; a *jump* module's
  **outbound** edge, which is the direction convention read backwards and was
  the one quiet mistake left, since the inbound form throws for missing
  `rebuilt_params` while the outbound form fell through the consumer filter and
  composed with the stochastic block silently keeping its nominal constants; an
  `:ode` module declaring `rebuilt_params`; two intervals on one module; and an
  interval that is not a whole number of handshakes, since the rebuild can only
  fire at one. The continuous cadence of 4.4 is counted separately. The last two
  of the ten — the outbound jump edge and the both-blocks name scan — came from
  the pre-merge review; phase 3's catalytic `param_slot` guard had the same
  one-block scan and is widened in the same pass.
- [x] 4.3 Assert piecewise-constancy — verify by sampling propensities between
  refreshes and asserting they are bitwise unchanged, so the coupling is the
  published piecewise-constant one and not an accidental continuous one.
  **Done:** over 180 s at a 60 s interval the rebuilt slot is bitwise unchanged
  across each inter-refresh window and changes at each boundary, while the
  upstream pool takes a different value at more than 150 of the 180 handshakes —
  so the constancy is not vacuously true of a channel that never fired. The
  propensity statement is made at a **fixed reference state**: a raw propensity
  moves with the jump state at every event, so only at a held state does
  "unchanged between refreshes" say anything about the coupling. **The toy
  rebuilds two constants for this test**, because transcription is zeroth order
  — in Core A′ as in the toy — so its propensity simply *is* its rate constant
  and asserting on that alone would restate the parameter-slot check. Rebuilding
  translation's constant as well gives a propensity first order in the
  transcript, which is the shape Core A′'s charged-tRNA channel has, and it is
  that one the reference-state assertion is made on.
- [x] 4.4 Handle a continuous cadence explicitly — verify by either implementing
  it and asserting `reduction_declarations` labels it a deviation, or rejecting it
  with a named error; silently discarding the cadence field is the failure to
  avoid. **Done: rejected, by name.** The outer loop of task 3.1 holds a constant
  between refreshes by construction — which is what makes 4.3 assertable at all —
  and propensities tracking a pool continuously would need the
  `JumpProblem`-over-`ODEProblem` shape 3.1 rejected. The error names the module,
  the species and that reason, and points at a shorter piecewise-constant
  interval, which is what §10 R11's cadence comparison actually varies.
  `reduction_declarations` still labels a declared-but-unexecuted continuous edge
  a deviation, unchanged.
- [x] 4.5 Measure the channel's gain on the toy — verify by an elasticity
  diagnostic and by the number recorded where phase 10 compares it against the
  0.044–0.051 the real transcription module predicts. **Done:**
  `rate_constant_elasticity` computes `d ln k / d ln pool` by central difference
  in log space. The toy's rebuild law is Michaelis-shaped, so its elasticity has
  the closed form `km/(km + pool)` and the diagnostic is **checked** against it
  rather than merely reported: fifteen (km, pool) settings in
  `dev/scripts/rebuild_channel_result.md`, every one
  agreeing to better than **7e-9** relative. The gain spans **0.0025 to 0.667**
  across that sweep, which is the point R8 makes — a channel gain is a function
  of the quantities we assert, not a constant of the model — and the row at
  `km = 0.175`, pool 3.6529 mM sits at **0.0457**, inside the 0.044–0.051 band
  phase 10 will compare the real transcription channel against.
- [x] 4.6 Run D10's granularity comparison in miniature — verify by nominal
  trajectories at 1 s, 5 s and 60 s drain granularity with the largest relative
  difference in any pool reported, so phase 13 inherits a measurement rather than
  an assumption. **Done:** a `drain_interval` on the driver, defaulting to the
  handshake interval, required to be a whole number of them, and labelled by
  `driver_declarations` when coarser. Measured over 300 replicates at a 600 s
  horizon, **paired per seed**, sampled at **drain-aligned** handshakes and
  reported one row per pool with nothing selected: the difference in the
  dynamics is **not resolved above Monte Carlo noise at either granularity, on
  either pool**. On the debited pool it is +0.02% ± 0.02% at 5 s and
  +0.01% ± 0.02% at 60 s; on the product pool, −0.17% ± 0.16% and
  −0.04% ± 0.17%. All four are far below D10's one percent.

  Reported separately, because it is not a granularity cost: the **sawtooth of
  the outstanding debit** on the debited pool, +0.65% ± 0.02% at 60 s against a
  closed form of +0.64%, and +0.06% ± 0.02% at 5 s against +0.04%. The closed
  form is `(drain − interval)·k_tl·(k_tx/γ)·cost / (factor·pool)`, computed from
  the toy's own declared parameters with nothing simulated. Between drains a coarse
  configuration holds up to `drain − interval` seconds of unpaid cost, so its
  pools sit high by exactly that; the quantity is deterministic and needs no
  simulation. **The first version of this measurement reported that sawtooth as
  the granularity cost**, sampling at every handshake so that the maximum
  landed at t = 599 — the instant furthest from the last drain. The pre-merge
  review of PR #45 caught it; §12's amendment B records it, because phase 13
  inherits the method and not only the number.

  **Phase 13 inherits three method points.** Sample at instants that are
  multiples of every drain compared, where each configuration has just debited
  the same total accrual. Pair the difference per seed and quote the standard
  error of *those*, since coarsening the drain changes how often the propensity
  aggregation is rebuilt and so consumes a different amount of randomness —
  three configurations at one seed give three different paths. And state the
  multiplicity behind any maximum, since a maximum over a noisy field is a
  selection statistic. Read narrowly, as the result file says: the toy's ATP
  pool is ~404,000 particles against a ~40 particle/s drain, a buffer of hours,
  where the pools D10 is about turn over in 109 s and 30 s — **so this bounds
  the mechanism and nothing else.**
- [x] 4.7 Suite and handoff — verify by `sbatch test/run_tests.slurm` passing and
  by the handoff recording that bidirectional coupling is now executed rather
  than declared. **Done:** job 16222708, **1230/1230** (1126 before the phase);
  the handoff records that five of the seven edge kinds now execute and that
  volume and clamped remain. Phase 3's wall-clock artefact was re-run in the
  same pass, since `run_handshake!` now copies the jump parameter vector per
  handshake: job 16221975, worst extrapolation **0.019 s** per 6,300 s
  trajectory against K1's 10 s budget, up from 0.014 s and still ~500× inside
  it.

### Phase 5 — Growth and volume

**Goal:** make the count-to-concentration conversion volume-aware, so growth
dilutes the ODE block as the published model's does.
**Done when:** a volume edge on the toy's membrane-protein state drives radius
and volume, both conversion directions read the live volume, and the published
initial surface area returns a 200 nm radius ~~exactly~~ *to the four
significant figures the source states it in — it is 200.03505 nm, and see task
5.2*.
**PR:** #46 (merged 2026-09-05)

- [x] 5.1 Promote a membrane-protein declaration to the protocol with an empty
  default, settling the open question the drafted transport design leaves for
  this phase — verify by a composition sweep finding the flagged states without
  naming any module type, and by a module declaring none composing unchanged.
  **Done:** `membrane_protein_states(m)`, empty default, promoted exactly as
  `add-pts-transport/design.md`'s D7 wanted and its Open Questions deferred to
  "wave 2's call". The flag is data on the module that *owns* the state, and
  `membrane_protein_states(models)` is the sweep — asserted at three
  composition orders, with no `isa` on any module type. **Either block may
  declare it**, because the toy and the real model differ: the phase-3 toy owns
  `M_ptsg_c` as a jump state and phase 11 task 11.4 makes translocation what
  increments it, while task 7.6 has `PtsTransport` own both phospho-forms in the
  ODE block. A jump state's count is read directly; an ODE state's is
  `u × factor` at the factor in force at that hook, continuous rather than
  rounded so the geometry adds no quantisation to check 0's round trip. A
  composition flagging none has an empty chain, a bitwise constant factor over
  100 handshakes and no growth label. **Fifteen refusals**, each naming what is
  wrong, and each asserted on its own message. From the sweep (2): a flag on a
  state its declarer does not own; the same species flagged twice, whose count
  would enter the area twice. From the lowering (8): an *inbound* edge, which is
  the half of this channel phase 5 does not build (see task 7.2) —
  **superseded by phase 5b**, which built that half; the refusal is now the
  constructor's, for an inbound edge carrying no `param_slot` and no
  `quantity`, and the count here stands as the record of what phase 5 did; a flag with no
  **outbound** `VolumeEdge`, which would execute undeclared; an **outbound** edge
  with nothing flagged, which would be declared and never executed; an
  **outbound** edge on a species the module does not flag — all three read
  *outbound* since phase 5b, which filtered this pairing to producers, because
  otherwise a module that both flags a protein and reads the geometry is refused
  for the wrong reason; a flagged state with no edge of its own — the trap a module
  flagging both ptsG phospho-forms and edging one would hit; and the three
  `extracellular_states` refusals of task 5.3. From the build (5): a growth
  keyword passed to a composition whose cell does not grow, and `radius_nm`
  passed to one whose cell does, either of which would state the geometry twice
  or not at all; a non-positive initial area; a non-positive footprint, which
  would leave the chain declared, labelled and frozen; and a derived baseline
  below zero, which would make area a super-linear function of count instead of
  the published law's affine one.
- [x] 5.2 Implement the surface-area, radius and volume chain — verify by
  asserting the published initial 502,831 nm² returns ~~exactly~~ 200.0 nm, and by
  the 28.0 nm² footprint carrying its calibrated-not-measured label, with the
  upstream docstring's contradictory 35 nm² recorded and the code's value taken.
  **Done, and "exactly" is the one word this task got wrong.**
  `sqrt(502831/4π) = 200.03505 nm`, not 200.000 — `4π(200 nm)²` is 502,654.8 nm²,
  so the published area is **176 nm² larger than an exact 200 nm sphere**, and
  both `dev/notes/well-stirred-minimal-cell.md:194` and the scoping note's "r =
  200.0 nm exactly" are four-significant-figure statements. The area stays the
  primitive, because task 5.5's arithmetic closes only with 502,831, and the
  published figure is asserted at the precision the source states it
  (`round(r; digits = 1) == 200.0`) with the full-precision value pinned beside
  it — phase 3's own precedent, where 20,180 particles/mM is asserted as
  `round(Int, ·) == 20180` against a true 20180.39. Not logged as a §12
  amendment: the number the spec asks for is reproduced, and it is the word
  "exactly" that is four significant figures. **The frozen baseline is derived,
  never typed:** 502,831 − 831 × 28.0 = **479,563 nm²**, asserted directly, because
  getting that offset wrong is what would make doubling ptsG double the whole
  cell. `driver_declarations` gains two labels whenever the chain is live —
  `:calibrated_constant` for the footprint, carrying `getProtSA`'s 35 nm²
  docstring against the 28.0 its code runs, and `:exogenous_growth` for the
  frozen non-ptsG baseline, which §6 T2 needs as a row.
- [x] 5.3 Route both conversion directions through live volume — verify by a test
  that a volume increase dilutes every ODE concentration by exactly the volume
  ratio with no change in any count. **Done:** `factor` was "particles per mM,
  at fixed volume" and is now recomputed at **step 0 of every handshake**, before
  the catalytic write, so every conversion in an exchange uses one volume. A
  bigger cell is *more* particles per millimolar, so the factor rises and every
  concentration at fixed count falls: at 831 → 1662 copies the ratio is
  **0.93440**, against `(200.035046/204.610919)³ = 0.93440` from the geometry.
  **Dilution is a fourth operation** — none of the three existing conversion
  sites implies it — and each state is rewritten as the count it held divided by
  the new factor, in that order, because going through the volume ratio instead
  would conserve the ratio and let the counts drift. Exactness is **one ulp, not
  bitwise**: `(u·f_old)/f_new · f_new` need not return `u·f_old` in binary
  floating point, so counts are asserted at `rtol = 1e-14` and the statement is
  made rather than hidden in a tolerance. The carried remainders are in particles
  and survive untouched, asserted. **Annotated 2026-09-05: "every ODE
  concentration" gains an exemption.** It is true of the toy and not of the
  assembled model — task 7.4 integrates external lactate as a dynamic state whose
  millimolar is a *medium* concentration through D6's 1e5 ratio, and diluting it
  by the cell's volume ratio would destroy lactate that has already left and open
  check 2. The registry has no extracellular marker (`M_lac__L_e` is an ordinary
  dynamic species there), so a second protocol declaration `extracellular_states`
  carries it, empty by default, with three refusals: a `:jump` module declaring
  it at all, since its states are counts and a count is referred to no volume; a
  state the declarer does not own; and a state that is both exempt and
  membrane-flagged, which would make the area depend on the volume it sets.
  Recorded here and
  in the pull request rather than as a §12 amendment, as tasks 2.2 and 3.2
  recorded theirs: nothing the spec said became false.
- [x] 5.4 Cap growth at exactly twice initial volume as the published model does —
  verify by a test driving counts past the cap and asserting volume stops rather
  than growing. **Done:** applied to the **volume**, as `in_out.py:107` applies
  it, leaving the radius uncapped exactly as upstream leaves it — an equivalent
  cap on area would carry a factor of 2^(2/3). Driven there artificially, since
  the cap needs ~11,400 membrane proteins against the 831 Core A′ starts with:
  at 50,000 copies the volume sits at exactly `2 × V(0)` and stays there at
  500,000 while the radius keeps rising. The published constant is a hard-coded
  6.70e-17 L, twice its own *rounded* 3.35e-17; ours is twice the computed
  3.35279e-17, and the difference is that rounding rather than a different rule.
- [x] 5.5 Reproduce the scoping note's ptsG arithmetic — verify by a test
  asserting that doubling 831 copies moves surface area from 502,831 to 526,099
  nm², radius from 200.0 to 204.6 nm, and volume to about 1.07×. **Done, all
  three, read off the growth trace rather than asserted from the counts that
  produced it:** area **502,831.0 → 526,099.0 nm² exactly**, radius
  **200.03505 → 204.61092 nm**, volume **1.07021×**. Checked first against a
  closed form written out longhand in the doubles file — as phase 4 checked its
  elasticity against `km/(km+pool)` — so the assertion is the published law and
  not the implementation restated, and only then against the note's three
  numbers.
- [x] 5.6 Refuse the doubling-time comparison in code — verify by the composed
  model exposing the reporting constraint that this is fractional growth or
  time-to-threshold, and by a test asserting the constraint is retrievable rather
  than living only in prose. **Done, in two halves because the task asks for
  two.** `reporting_constraints(driver)` returns the constraint as structured
  data — quantity, verdict, reason, and what to report instead — so a composed
  model is *asked* rather than read; and `doubling_time(driver)` exists and
  throws, naming the ~92% reduction. `growth_report` and `time_to_threshold` are
  the two admissible reportings, the latter read off the recorded trajectory
  rather than extrapolated from a rate, a rate being one algebraic step from the
  quantity being refused. The constraint does not depend on whether a given
  composition grows: it is a statement about the reduction. Phase 14 task 14.9
  consumes it.
- [x] 5.7 Suite and handoff — verify by `sbatch test/run_tests.slurm` passing.
  **Done:** job 16226536, **1381/1381** (1230 before the phase). Task 3.8's
  wall-clock artefact was re-run in the same pass, as phase 4 did, because step
  0 adds work to every handshake. `dev/scripts/bench_handshake.jl` gains a
  **matched pair** — the same gene at the same 831 initial copies, with and
  without the flag and the edge — because the pre-existing `live rate law` row
  starts at zero protein and that count fills the ODE rate law through a
  catalytic edge, so a gap against *it* would be a different trajectory as much
  as a volume chain. Each row is now the minimum of five repetitions with the
  spread reported beside it, because at this scale the scatter and the effect
  are the same size. Verdict unchanged: worst extrapolation well inside K1's
  10 s budget. **The chain's own cost is reported against its control and beside
  that spread, and is not claimed to be resolved above it** — and it is measured
  on two diluted ODE states against Core A′'s thirty-two, so it bounds the
  mechanism, as every toy number in phases 3 to 5 does.

### Phase 5b — The inbound volume channel

**Goal:** the other half of the volume channel — the cell's geometry re-entering
an ODE rate law as a parameter — so that no module phase carries a framework
dependency into the fan-out.
**Done when:** an inbound `VolumeEdge` resolves and executes; the slot it names
holds this handshake's geometry at every handshake, in a growing cell and in a
fixed one; the three quantities describe one sphere above the growth cap as
below it; and the four-way fan-out of phases 6, 7, 8 and 10 can start with no
module phase needing a framework change.
**PR:** #48 (merged 2026-09-09)

**Why this is its own phase.** R15's response is "amend §12 with the change and
its date, then land it as its own phase before the affected modules continue",
and that is exactly what this is. What did not happen is R15's *trigger*: the
need never surfaced during a module phase, because phase 5 diagnosed it, refused
the edge by name and wrote the diagnosis into this spec. So the only open
question was when, and the answer is before the fan-out rather than during it —
otherwise phase 7's author meets a live `ArgumentError` in framework code they
may not edit, and has to stop and wait for this phase anyway.
- [x] 5b.1 Give `VolumeEdge` a `param_slot` and a unit-bearing `quantity`,
  required inbound and refused outbound — verify by each of the five refusals on
  its own message; by `:radius`, `:volume` and `:radius_um` being rejected as
  hard as a nonsense name, since nothing else in the codebase records a unit and
  the symbol is the only place a nanometre is told from a centimetre; and by
  every existing outbound call site constructing unchanged.
  **Done:** five refusals, each on its own message. `VOLUME_QUANTITIES` sits
  beside `CLIP_POLICIES` and `RATE_CADENCES`; validation is in the inner
  constructor so positional construction cannot bypass it. Every existing
  outbound call site compiles unchanged, since neither field takes a
  direction-dependent default.
- [x] 5b.2 Require an inbound edge's `species` to be the declaring module's own
  state — verify by the refusal firing on `resolve_coupling` of the single
  module, without a composition. The geometry is a sum over every flagged state
  and so belongs to no module, but every edge must name a registry species
  (`src/resolver.jl`), so the species names *which of this module's rate laws is
  geometry-dependent*. Naming the membrane protein instead would be a claim that
  silently changes meaning once a second membrane protein is flagged — which
  task 7.6 plans, with `PtsTransport` owning and flagging both ptsG
  phospho-forms.
  **Done:** in `_resolve_edges`, beside the currency-pool check, so it fires
  on a single module. `src/resolver.jl` requires every edge to name a
  registry species, which is what rules out a pseudo-species and makes this
  the only honest reading of the field.
- [x] 5b.3 Lower it in `_lower_exchanges` rather than `_lower_growth` — verify
  by the slot resolving to the same index map the catalytic channel uses, since
  that is the only place with `ode_contexts`; by a slot that is not one of the
  module's own free parameters being refused; and by the cross-block
  slot-name-uniqueness rule applying unchanged.
  **Done:** `GeometryExchange`, which is `CatalyticExchange` minus `count_idx`
  and plus `quantity` — the absence is the channel, since the value comes from a
  driver field rather than a state vector. All three clauses check out: the slot
  resolves through the same `ode_contexts[i].param_idxs[j]` map the catalytic
  channel uses, asserted at `i > 1` by 5b.7's two-ODE-module composition; a slot
  that is not a free parameter is refused; and the cross-block
  slot-name-uniqueness rule is the catalytic one, reused unchanged.
- [x] 5b.4 Filter `_lower_growth`'s flag-to-edge checks to producers — verify by
  a module that both flags a membrane protein and reads the geometry building
  cleanly, which is the case an undirected check refuses with a message about
  the wrong problem; and by an inbound edge still not discharging a flag's
  obligation.
  **Done:** one `is_producer` in the comprehension. Without it the
  `e.species in names` check tells `ToyExportingPool` that its lactate edge
  names an unflagged species, which is why every test in
  `test_volume_inbound.jl` builds at all: testset 5b.4 composes the flagging
  form directly, and the rest use a non-flagging variant of the same double, so
  the filter is exercised where it bites while the growth arithmetic stays a
  function of one count.
  The flag-with-no-edge message gains the word *outbound*.
- [x] 5b.5 Refuse one slot with two writers, across the catalytic and geometry
  channels both — verify by the refusal firing, and note that this closes the
  same hole for two catalytic edges into one slot, which the catalytic channel
  has today.
  **Done:** `_claim_slot!`, shared by both channels. The catalytic hole it
  closes as a side effect was real: `slot_owners` guards two *modules*
  sharing a name and cannot see one module declaring two edges into one
  slot.
- [x] 5b.6 Refuse an inbound volume edge on a jump module — verify by the
  message naming the remedy that does exist, `RateConstantEdge` with
  `rebuilt_params`. Its parameters live in the propensity vector; resolved
  against the ODE block's layout instead, the write would be in range,
  type-correct and land in an unrelated rate law.
  **Done:** scanned after the ODE loop, in the shape of the catalytic mirror
  trap, with `RateConstantEdge` named as the remedy that exists.
- [x] 5b.7 Write at step 0 of every handshake, outside `_update_volume!` —
  verify by a fixed cell's slot holding the driver's geometry rather than the
  module's declared value, which is the case `_update_volume!`'s early return
  would silently skip; and by the value being this handshake's rather than the
  last one's, asserted against a count driven hard enough that a one-handshake
  lag is a whole percent rather than the 1e-5 it is at published parameters.
  **Done:** step 0b of `handshake_step!`. `ToyFixedExporter` declares
  `r_cell = 1.0` so a skipped write is a 200× error rather than a subtle one,
  and the lag test drives the count to 20,000 so a one-handshake lag is a whole
  percent against the 1e-5 it is at published parameters. That composition also
  puts the geometry slot **behind a second ODE module**, which is the only place
  `ode_contexts[i].param_idxs[j]` is exercised at `i > 1`: local index 2 against
  global 5, so a mapping that confused the two would write into the wrong
  module's rate law, in range and silent.
- [x] 5b.8 **Derive the rate law's geometry from the capped volume**, and label
  it — verify by all three quantities describing one sphere above the cap, by
  the reported `radius_nm` still being the published uncapped value, and by
  `driver_declarations` carrying the departure. The cap is on the volume and the
  radius is left uncapped, as `in_out.py:107` leaves it; that is harmless while
  the radius is only reported and wrong once it drives `3P/r`, whose whole
  content is that `3/r` is a sphere's surface-to-volume ratio. §2 records
  that the
  cap is never approached at published parameters — it needs ~11,400 membrane
  proteins against Core A′'s 831 — so this is not a regime the nominal model
  visits. It is fixed rather than argued away because a rate law reading a
  geometry that is not a sphere is wrong wherever it happens, and because
  nothing bounds a prior draw away from it.
  **Done, and the first Slurm run corrected it.** Deriving unconditionally
  round-trips radius → volume → radius and loses an ulp, so a fixed cell at
  exactly 200 nm handed its rate law 200.00000000000003. Below the cap the
  reported geometry is already one sphere — volume was computed from radius,
  radius from area — so it is passed through and only the capped branch
  reconstructs. The below-cap assertion is `==`, which is what says so.
- [x] 5b.9 Make the channel observable from outside — verify by
  `growth_census(driver).consumers` naming the species, quantity, slot, declarer
  and current value; by `run_handshake!` recording `ode_p` per handshake, so a
  write can be read off the trajectory rather than inferred from the schedule,
  which task 13.2's "every declared edge is executed" needs; and by
  `driver_written_params(driver)` enumerating all three channels' slots.
  **Done:** `growth_census(driver).consumers`, `run_handshake!`'s `ode_p`, and
  `driver_written_params(driver)` across all three channels, the rate-constant
  one asserted on phase 4's own composition. The census row is **read back out
  of the parameter vector** rather than recomputed: recomputing would report the
  right geometry even if the write never landed, which is the one thing the row
  exists to witness. A test pokes the slot and watches the census follow.
- [x] 5b.10 **Record what this does not fix.** A `param_slot` must be a *free*
  parameter, and free non-observation parameters are exactly what inference
  samples, so a driver-written slot is sampled and then overwritten. This is
  pre-existing — `rebuilt_params` documents it at `src/interface.jl`, and the
  catalytic channel has it undocumented — but it is sharper here: a geometry
  slot is not an unknown at all, and in `3P/r` its Jacobian column is exactly
  proportional to the permeability's, so a rank deficiency of one is an artefact
  of the declaration rather than a finding. Verify by
  `driver_written_params(driver)` enumerating the slots of all three channels,
  which is the checkable half; the warning itself is `rebuilt_params`' restated
  on `VolumeEdge`, with the two reasons it is sharper here. **Excluding them
  from the
  sampled set belongs to the inference phases (15 to 17)**, which is where
  `rebuilt_params` already assigns it and where this spec keeps it. It is
  explicitly **not** task 13.3: 13.3 is hybrid *dispatch* — which code path runs
  — and its test, a mixed composition dispatching or throwing, cannot catch
  this, because the composition that silently samples a geometry slot is the
  *homogeneous* ODE block.
  **Done:** the `rebuilt_params` warning carried verbatim onto `VolumeEdge`,
  with the two reasons it is sharper here — a geometry slot is not an
  unknown, and in `3P/r` its Jacobian column is exactly proportional to the
  permeability's. Recorded in §12 that task 13.3's own test does not cover
  this, since the composition that silently samples one is the homogeneous
  ODE block.
- [x] 5b.11 Suite, figures and handoff — verify by `sbatch test/run_tests.slurm`
  passing, and by the progress figures rebuilding with phase 5b parsed rather
  than silently skipped.
  **Done:** job 16341043, **1468/1468** (1381 before the phase). Both
  figures rebuilt; the parser now reads a lettered phase and refuses a
  heading it cannot parse rather than skipping it in silence — a guard worth
  having independently. The old regex did not skip this heading quietly — it
  attributed 5b's `**PR:**` line to phase 5 and tripped the pre-existing
  two-PR-lines assertion — but it would have gone quiet for a phase whose
  heading carried no PR line, and the guard closes that. **A trap worth
  recording:** the figures are generated from this section, so ticking a box
  without re-running `make_progress.py` leaves them understating the work by
  exactly that many tasks, silently. It happened once in this phase and the
  pre-merge review caught it; the figure README now says so.

### Phase 6 — Central glycolysis

**Goal:** the ten reactions from glucose-6-phosphate through lactate, owning the
eleven glycolytic intermediates and the two redox species.
**Done when:** the thirteen owned states integrate non-negatively over a full
cycle from the registry's initial conditions at published parameters, the redox
pair's residual ~~is conserved to a bound that shrinks with the solver
tolerance~~ **sits at the floating-point floor and does not move with the solver
tolerance, which is what an exactly conserved moiety looks like and what check 3
predicted** (amended 2026-09-10; see §12), and a
mutation to the dehydrogenase's stoichiometry makes that check fail.
**PR:** #50 (merged 2026-09-10)

The drafted task list at
`dev/archive/openspec/changes/add-central-glycolysis/tasks.md` holds the
fine-grained sub-checks and the spot values to assert; it is the reference for
this phase's detail. Three amendments, from D0: the outbound currency edges now
**execute**, so add an assertion that each moves the mass its rate law says; drop
the "declared but unexecuted" clause from its handoff note; and use the real
recycling module in place of its held-energy double wherever phase 8 has landed,
keeping the double only for standalone runs.

- [x] 6.1 Vendor the derived extract of ~65 rows and its regeneration script —
  verify by the row classes and counts matching (10 forward and 10 reverse
  catalytic constants, 32 Michaelis constants, 13 concentrations), by four spot
  values matching upstream byte for byte, and by the README's command
  round-tripping to an unchanged file.
- [x] 6.2 Implement the sub-model, its thirteen states and the modular rate law
  generic over substrate and product counts — verify by the state set equalling
  the registry's glycolytic and redox groups with strictly increasing indices, and
  by unit tests that equal forward and reverse terms give zero net rate and that
  the per-reaction term counts sum to 32.
- [x] 6.3 Import every value through the loader with log-normal priors from mode
  and geometric standard deviation, fixed by default — verify by every parameter
  reporting the central file as its source, by the two prior-default species
  keeping their width, and by the registry-agreement check firing when a
  concentration row is mutated in a temporary copy.
- [x] 6.4 Set the ten enzyme concentrations from copy number at the registry's
  volume, marked nominal and overridable, and declare the ten protein counts as
  inputs so translation later supersedes them — verify by each value equalling
  copies over ~~20,180 to six decimals~~ **`corea_particles_per_mM()`, which is
  20,180.39 and not the scoping note's rounded 20,180 — the two agree to five
  decimals and differ in the sixth for PFK, GAPD, PGK, ENO and LDH_L** (amended
  2026-09-10; see §12), by none being the published no-rule default
  of 0.001 mM, and by overriding one scaling exactly the rates that enzyme
  catalyses.
- [x] 6.5 Declare the boundary: currency edges on the energy species and mass
  edges on the shared intermediates, every peer unnamed — verify by the edge
  count, kinds and directions matching, by no edge naming a redox species, and by
  standalone resolution succeeding while listing the energy species as unowned.
- [x] 6.6 Register the reduction notes: the dropped oxidase with its three
  reasons, the Michaelis-constant column choice naming all eight differing
  constants with both values, the held currencies, and the nominal enzyme
  concentrations — verify by `reduction_declarations` returning all four as
  sentences naming the affected reactions.
- [x] 6.7 Integrate and check redox balance — verify by non-negativity at every
  save point, by the redox residual ~~satisfying the tolerance principle across two
  tolerance settings~~ **measured across a ladder of tolerances spanning eight
  decades and asserted flat at the floating-point floor** (amended 2026-09-10;
  see §12), and by the mutation test failing and naming the conserved
  pool and the drift.

### Phase 7 — Phosphotransferase transport and lactate export

**Goal:** the five-step cascade that imports glucose against phosphoenolpyruvate,
plus passive lactate export, owning the eight carrier phospho-states and external
lactate.
**Done when:** the four carrier sums are independently conserved to
tolerance-derived bounds, each carrier's two forms sum to its published copy
number exactly, and the export rate matches the permeability law at the
registry's radius.
**PR:** #49 (merged 2026-09-10)

Reference detail at
`dev/archive/openspec/changes/add-pts-transport/tasks.md`. One amendment: the
membrane-protein declaration is the protocol function phase 5 promoted, not a
plain function on this module, which settles that design's own open question.

Three more things moved under that design since it was written, all of them the
protocol rather than the chemistry: phase 1 made the contribution channel
execute, so this module implements `contributed_states` and `contributions`
rather than declaring production and leaving it unexecuted; phase 5b delivered
the inbound volume channel task 7.2 needs; and phase 5 added
`extracellular_states`, which task 7.4 uses. The phase also amended §3's
tolerance principle — see task 7.7 and §12.

- [x] 7.1 Vendor two extracts — eleven rate constants with no uncertainty column
  and nine initial conditions with their derivation — verify by every constant
  loading with informedness `asserted` and appearing in the composed model's
  asserted-prior enumeration, and by the balanced glycolytic parameters not
  appearing there.
  **Done:** `pts_transport.tsv` (11 rows, no `GeometricStd` column) and
  `pts_initial_conditions.tsv` (9 rows), regenerated by
  `dev/scripts/extract_pts_transport.py` and byte-identical on a re-run. All 11
  load `:asserted` from the table's *shape*, not a declaration. The
  asserted-prior enumeration returns **23** for this module — the 11 constants,
  the radius, the volume ratio, the external-glucose clamp and the 9 initial
  conditions — and a balanced glycolytic parameter composed beside it does not
  appear. The clamp is in that list and belongs there: 40 mM is a point value
  with no quantified uncertainty, which is what the registry records for it too,
  and that is a different claim from the *clamp* being ours, which task 7.5
  checks separately. `ambiguity_report`
  over the two tables is empty.
  **The copy numbers are read from `proteomics.xlsx` with the Python standard
  library**, `zipfile` plus `ElementTree`: `openpyxl` is in no environment on
  this machine and compute nodes have no network. The four carriers are keyed by
  AOE id — `AOE93321.1`, `AOE93571.1`, `AOE93322.1`, `AOE93596.1` — which skips
  the published loader's second workbook, and the script asserts each row's NCBI
  description as well as its rounded count. Phase 10's promoter proxy reuses the
  reader.
  **One divergence from the archived design.** Its D3 table used a rounded
  20180 particles per mM; this uses the 20180.3873819 the code derives from
  Avogadro, so **five** of the eight phospho values differ in the seventh
  decimal, by at most 7e-7.
  The derived factor is the one the values are read back with, so the carrier
  totals land on 353, 314, 290 and 831 to floating-point rather than on
  353.007 — ptsI sums to 352.99999999999994, which is why the test asserts
  `rtol = 1e-12` and not equality.
- [x] 7.2 Implement the cascade as mass-action forward and reverse pairs and the
  export law with permeability and radius carried separately — verify by a
  hand-checked derivative at one state vector, and by the export rate equalling
  the permeability law rather than a folded constant, so phase 5's growing radius
  changes it. **Phase 5b delivers the channel this needs**, so this is ordinary
  module work: declare an inbound `VolumeEdge` on `M_lac__L_e` with
  `param_slot` naming the radius in this module's own rate law and
  `quantity = :radius_nm`, and the driver fills it at every handshake. Two
  things to know. The geometry a rate law receives is derived from the *capped*
  volume rather than from the reported uncapped radius, so above the growth cap
  it is the smaller of the two — a labelled departure, see §12. And a standalone
  homogeneous build executes no handshake, so the slot keeps this module's
  declared value; declare it as the registry's initial radius, and note that the
  driver writes 200.03505 nm rather than a round 200.0.
  **Done:** all nine derivatives hand-checked against the signed sum of the
  rates touching them at one state vector. The export rate is
  `P·(lac_c − lac_e)·3/(r_nm·1e-9)`: zero at equal pools, negative when external
  exceeds cytosolic, `0.075·(lac_c − lac_e)` at 200 nm, and **exactly doubled at
  100 nm**, which is what says the radius is not folded into a constant. The
  slot is in nanometres and the permeability in metres per second, so the
  conversion is in the rate law and the unit-bearing `quantity` is what makes it
  checkable. Declared at `COREA_INITIAL_RADIUS_NM`. No saturation: with
  pyruvate absent, doubling both of GLCpts0's substrates quadruples it.
- [x] 7.3 Set the eight initial conditions from copy number times proteomics
  fraction (D7) — verify by each carrier's forms summing to 353, 314, 290 and 831
  copies exactly, and by both the totals and the split recorded as the published
  model's rather than asserted.
  **Done:** the four sums round-trip to 353, 314, 290 and 831, and to a
  relative 1e-12 rather than merely to the nearest particle, because the extract
  is built with the same conversion factor it is read back with. Both published
  inputs travel with each row — the extract's `Copies` and `ProteomicsFraction`
  columns and an `UpstreamRow` naming `protein_metabolites_frac.csv` or
  `membrane_protein_metabolites.csv` — so the split is recorded as the published
  model's. The *prior* on each is `:asserted`, which is a different claim and
  the true one: no upstream file quantifies uncertainty on a copy number.
- [x] 7.4 Integrate external lactate with the volume ratio (D6) — verify by
  external lactate after a full cycle staying below one percent of steady
  cytosolic lactate at the default ratio, by a ratio of one visibly saturating
  export, and by the ratio registered as ours. **Declare `M_lac__L_e` in
  `extracellular_states`** — phase 5's exemption from the dilution growth applies
  to every other ODE concentration. Its millimolar is a medium concentration, so
  diluting it by the *cell's* volume ratio would destroy lactate that has already
  left and open check 2; verify by a growing composition leaving it untouched
  while its cytosolic neighbours dilute.
  **Done, and it needed a lactate source to be posable at all.** Steady
  cytosolic lactate is glycolysis's, and glycolysis is phase 6, so the
  `HeldMetabolites` double supplies it at the published two-lactate-per-glucose
  rate. D6's own figures come back: cytosolic lactate settles at
  **1.46836 mM** and external reaches **0.0068918 mM** after 6,300 s, which is
  **0.469%** against D6's stated 0.47%. A ratio of 1e4 still misses the
  criterion, at **4.50%** — D6's 4.7% is a naive tenfold of the 1e5 figure and
  ignores that cytosolic lactate rises with the external pool, so the decision
  D6 records stands while its second percentage does not reproduce. The settling value is `production/0.075` **plus** the
  external pool, because export is driven by the difference between them; the
  0.47% by which it exceeds 1.4615 mM is exactly the lactate that has left,
  which is the coupling working rather than an error in it. A first draft
  asserted that difference away and failed, correctly.
  **The 6,300 s integration costs 4.56 ms** (job 16364852: best of five, 286
  accepted steps and 3 rejected over 64 save points). It runs in the default
  suite rather than behind a flag. What the suite actually pays is **19.67 s of
  compilation** on the first call, and that is the stiff-solver path being
  compiled once, not the horizon — a 600 s integration pays the same. An earlier
  draft of this note quoted ~19.5 s as the cost *of the full cycle*, which had
  the magnitude right and the meaning wrong; the pre-merge review caught that
  the figure came from a failing run of a superseded test, and separating the
  two is what answers §9 rather than restating it.
  At a ratio of one, total lactate is conserved across the two pools and export
  collapses to zero once they equalise. `M_lac__L_e` is declared in
  `extracellular_states`, and the ratio is registered in `reduction_notes`.
- [x] 7.5 Declare the boundary, including external glucose clamped at 40 mM with
  published origin — verify by that clamp *not* appearing in the
  what-is-ours enumeration while the chemostatted nucleotide and amino-acid pools
  do, since the published model clamps external species too.
  **Done:** ten edges, every `peer` unnamed. The clamp registers no `:clamp`
  deviation, because its origin is `:published`; the volume ratio does register.
  A rival clamp at the transport table's 42.77 mM throws against the registry's
  40. Standalone resolution succeeds and reports `M_pep_c` and `M_lac__L_c`
  among the unowned states rather than failing.
  **Ten edges, not the archived design's five.** Two of the five extra are the
  reverse readings: both reversible steps read the pool their forward term
  fills, and `inputs` — the only channel that wires a foreign state into
  `dynamics` — is held to the *inbound* edges, so `M_pyr_c` and `M_g6p_c` each
  carry an edge in both directions and `contributed_states` returns one net
  signed term for each. Phase 1 is why that matters: before it, the archived
  design could declare production and leave it unexecuted. The other three are
  the volume edges the archived design predates — two outbound on ptsG's
  phospho-forms (task 7.6) and one inbound on `M_lac__L_e` (task 7.2).
- [x] 7.6 Flag ptsG through phase 5's protocol function with its 28.0 nm²
  footprint — verify by a composition sweep finding both phospho-forms without
  naming this module's type.
  **Done:** `membrane_protein_states(models)` returns exactly
  `[:M_ptsg_c, :M_ptsg_P_c]`, found by sweeping the composition. Each carries
  the outbound `VolumeEdge` the flag owes. The 28.0 nm² footprint is the
  driver's `:calibrated_constant` declaration from phase 5, not restated here.
- [x] 7.7 Check carrier conservation as four separate assertions — verify by four
  independent bounds each naming its own carrier, by a mutation to one cascade
  step failing that carrier's check while the other three still pass, and by the
  bounds satisfying the tolerance principle.
  **Done, under §3's exact-invariant exception rather than the fivefold fall.**
  Each cascade step transfers one phosphate between adjacent carriers, so the
  composed right-hand side returns each pair's two derivative terms as
  bit-for-bit negatives. The test asserts §3's mechanical gate —
  `du[unphos] + du[phos] === 0.0`, bitwise, at every saved state of the
  trajectory rather than at one point — and then the ladder the exception
  requires: six rungs spanning ten decades of `(abstol, reltol)`, every rung
  within 100 ulps of the conserved sum and the largest within 100× of the
  smallest. Four bounds are still reported, each
  `2·max(abstol, reltol·maxₜ|x|)` with `N_restarts = 1`, since a standalone
  solve has no handshakes.
  **Phase 7 met the exception independently of phase 6, on different
  chemistry**, which is the useful part: measured 3.47e-17, 2.78e-17, 8.67e-18
  and 1.32e-16 at `(1e-10, 1e-8)` for ptsI, ptsH, Crr and ptsG against bounds
  seven orders larger, and at `(1e-11, 1e-9)` two rose and two fell with none
  falling by the required fivefold. A first draft of this phase proposed its own
  two-class amendment to §3; phase 6 had already landed the same finding with a
  stricter, mechanically checkable rule, so that draft was dropped on rebase and
  this check now asserts phase 6's criterion. See §12, 2026-09-10.
  A structural leak still fails, which is the part that matters:
  `PtsCarrierLeak` writes GLCpts1 as creating phospho-HPr rather than
  transferring it, and the ptsH residual then exceeds its bound while the other
  three stay under theirs.
  **Suite:** `sbatch test/run_tests.slurm` job **16374055**, **2023/2023** in
  4m58.2s on the tree rebased onto phase 6 (#50). The full-cycle timing is job
  16364852. Before the rebase this phase stood at 1654/1654 on the 1468-test
  base (job 16364952); earlier counts of 1619 and 1651 were measured against
  dirty or superseded trees and are not the ones it stands on.

### Phase 8 — Nucleotide recycling

**Goal:** the five reactions that make GTP and close the adenylate, guanylate and
phosphate moieties, owning the four adenylate, three guanylate and one
pyrophosphate species — the pools every other module routes energy through.
**Done when:** over a full 6,300 s cycle against a charging drain, adenylate and
guanylate are each conserved; removing the adenylate kinase drives ATP below 1%
of its initial value on the timescale the scoping note computes, and removing the
pyrophosphatase ~~leaves pyrophosphate unbounded~~ **strands the phosphate
moiety and stalls the pathway** (amended 2026-09-10; see §12 and task 8.6).
**PR:** #52 (merged 2026-09-17)

Reference detail at
`dev/archive/openspec/changes/add-nucleotide-recycling/tasks.md`. Its charging
drain double is the acceptance criterion and not scaffolding; it gains one task,
recording its stoichiometry so phase 9 can assert the real module reproduces it.

- [x] 8.1 Vendor two extracts — the governing nucleotide file and the central
  file's rival values for the same identifiers — verify by every one of the 35
  identifiers being held by both, so a truncated rival file fails loudly rather
  than quietly disarming the ambiguity check.
- [x] 8.2 Declare a governing file for all five reactions and all eight initial
  conditions (D2) — verify by loading failing when a governing declaration is
  removed, by `governing_choices` returning each value with the file chosen and
  the file rejected, and by the pyrophosphatase constant being the nucleotide
  file's 646.73 and not the central file's 583,611.
- [x] 8.3 Implement the five reactions with repeated stoichiometry expanded as the
  published builder does — verify by the two reactions with a doubled product
  each naming that species twice in the rate law, giving 17 Michaelis constants
  across five reactions rather than 19, and by the pyrophosphatase carried
  reversibly.
- [x] 8.4 Set the three new enzyme concentrations, and assert the two shared with
  glycolysis agree — verify by a cross-module test asserting one enzyme is not
  run at two concentrations, since two modules carry a nominal for each of the
  two shared genes. **Not independent of phase 6.** In a fan-out this assertion
  belongs to whichever pull request lands second; carry it here as a skipped
  test naming phase 6 rather than stubbing glycolysis to make it pass.
  **Annotated 2026-09-11: phase 6 landed first as PR #50, so phase 8 is the
  second and the assertion is written, not skipped.** It was worth writing —
  it failed on the first run, because this module divided copy numbers by a
  transcribed 20,180 where `central_glycolysis.jl` uses
  `corea_particles_per_mM()` = 20,180.39, putting `JCVISYN3A_0606` at 0.020367
  mM in one module and 0.020366 in the other. It is asserted as an equality
  rather than a tolerance for that reason.
- [x] 8.5 Declare the boundary: mass edges for species not owned, currency edges
  where this module is the **principal** producer or consumer of a pool it owns.
  The original rule said *sole*, and D13 broke it: charging also produces AMP and
  pyrophosphate, so recycling is no longer the only holder of those crossings —
  verify by the edge count and kinds matching the restated rule, by a test
  composing **all four** ODE modules' declarations asserting no species is
  described as both mass and currency in one direction, and by the restatement
  recorded so a reader does not apply the stale rule. **Annotated 2026-09-10:
  eleven edges, not the nine of the archived design.** The published PGK3 and
  PYK3 rate laws are reversible and read `M_3pg_c` and `M_pyr_c` twice over —
  once in the reverse numerator, once in the denominator's product bracket — so
  a declaration carrying them outbound only cannot evaluate its own rate law,
  and an ODE module reaches a foreign state through `inputs()` alone, which
  `src/resolver.jl` holds to an *inbound* edge. Each product therefore carries a
  mass edge in both directions, which is phase 1's rule of one edge per
  (species, direction) that actually occurs and is a true statement about a
  reversible reaction; `contributed_states` still names each species once and
  carries the net signed rate. Recorded here rather than in §12, since nothing
  this spec said became false — the rule above fixes no direction. **The four-module
  composition is not independent of phases 6, 7 and 9.** In a fan-out it belongs
  to the last pull request to land, or to phase 13; carry it here as a skipped
  test naming what it waits on. **Annotated 2026-09-11:** phases 6 and 7 have
  landed, so only phase 9 is outstanding and the skip names it alone; the
  three-module version of the assertion passes today.
- [x] 8.6 Build the charging drain double and run three full-cycle configurations
  — verify by the all-reactions run conserving both moieties with ATP positive and
  pyrophosphate settling bounded; by the kinase-removed run driving ATP below 1%
  of initial within an order of magnitude of the scoping note's 144 s; and by the
  pyrophosphatase-removed run ~~rising without bound~~ **stranding the phosphate
  moiety**. **Assert a threshold
  crossing, not exhaustion**: the note's 144 s assumes a constant drain, and the
  real charging module of phase 9 is mass-action, so ATP decays exponentially and
  never reaches zero. Record both numbers and which model each belongs to.
  **Annotated 2026-09-10:** pyrophosphate *cannot* rise without bound in a model
  whose phosphate is closed, which check 4b says Core A′'s is. Removing the
  enzyme instead sequesters the moiety — measured, a 128-fold rise holding 72.9%
  of the phosphate budget — and ATP regeneration then collapses for want of free
  phosphate, so charging stops and the rise flattens. That is a stalled pathway
  and it is the stronger demonstration that the reaction is required, as well as
  the one that survives assembly. The note's 173 mM is open-pool arithmetic:
  one pyrophosphate per charging event for a whole cycle, which assumes a
  phosphate supply the closed model does not have. Both are recorded with the
  model each belongs to. See §12.
- [x] 8.7 Check phosphate closure in both forms — verify by exact invariance with
  the GTP-branch reactions inactive, and by the flux-corrected form with them
  active *subtracting* the inbound flux rather than relaxing the bound.
  **Annotated 2026-09-11:** this task, and task 8.6's two moiety checks, assert
  the **second gate** of §3's exact-conservation exception — `|Σ nᵢ·duᵢ| ≤ 1
  ulp of the conserved sum` at sampled states, then the ≥6-decade ladder — and
  not the fivefold fall, which none of the three meets and none should. Phase
  6's first gate is unavailable here: the terms arrive from three modules and
  do not cancel bitwise. Measured 1.000, 0.125 and 0.115 ulps per evaluation
  against a mutated 1.7e16. The flux-corrected form is inside the gate rather
  than outside it: the correction is `u[GTP] − u0[GTP]`, a state difference and
  not a time-integrated quadrature, so the corrected closure is still a plain
  linear functional of the state with GTP weighted 2. **That stops being true
  at assembly**, where translation consumes GTP and the proxy fails, so check
  4b's assembled restatement needs a real flux accumulator — task 14.4's, not
  this one's. See §12, 2026-09-11.
- [x] 8.8 Record the drain's stoichiometry and rate for phase 9 — verify by the
  recorded rate being 553.1 per second, derived from 3,484,518 residues over
  6,300 s, and by it being recorded as the *demand* the real module must meet
  rather than as a value phase 9 may calibrate against and then re-check, which
  would be circular (D14).

### Phase 9 — Lumped tRNA charging, in the ODE block

**Goal:** the one reaction that is ours rather than the published model's,
integrated deterministically, owning the charged and uncharged tRNA pools and
contributing to the energy pools nucleotide recycling owns.
**Done when:** the tRNA pair is conserved over a full cycle against a
translation-demand double, ~~the charging flux reproduces the published residue
demand from independently derived inputs~~ **the charging flux's drift from the
independently derived residue demand is measured and reported as
self-consistency, not validation** (amended 2026-09-23, with task 9.7), and
`reduction_declarations` reports both the lumping and the formalism.
**PR:** #59 (merged 2026-09-23)

This is an ODE module, per D13. It was specced as a stochastic one, which
contradicted both the frozen registry and the scoping note; §12 amendment 1
records the correction. It is the fourth ODE module and it runs in the ODE track,
but it depends on phase 8, so the track is not a clean fan-out.

- [x] 9.1 Implement the single reaction as an ODE sub-model owning both tRNA
  states, with the mass-action rate law of §3 — verify by the state set equalling
  the registry's tRNA group with strictly increasing indices, by `formalism`
  reporting `:ode`, and by a hand-checked derivative at one state vector.
- [x] 9.2 Assert the total tRNA pool size and derive `k_chg` from it (D14) —
  verify by both loading with informedness `asserted` and appearing in the
  composed model's asserted-prior enumeration; by `k_chg` being computed from the
  published residue demand and the pool size rather than typed in; and by the
  derivation recorded so the pool size can be changed and `k_chg` follow.
  **Annotated 2026-09-23:** the derivation also needs the nominal charged
  fraction, a third asserted quantity; defaults 0.25 mM and 0.8. See §12.
- [x] 9.3 Declare the three currency edges on ATP, AMP and pyrophosphate, which
  recycling owns — verify by the kinds and directions matching, by composing with
  phase 8 asserting no species is described as both mass and currency in one
  direction, and by the contributions actually **executing** through phase 1's
  mechanism rather than resolving and doing nothing, which is what they would
  have done in the stochastic block.
- [x] 9.4 Register two declarations, not one — verify by `reduction_declarations`
  returning the lumping under `lumping`, naming what it replaces (twenty chains
  of five reactions) and why the published product stoichiometry was preferred
  over a two-ATP form that closes the same moieties by fiat; **and** returning
  the deterministic formalism separately, since the placement is the scoping
  note's but integrating it deterministically is ours.
- [x] 9.5 Add a translation-demand double to this phase's tests — verify by it
  consuming charged tRNA at the published residue rate and returning the
  uncharged form, since without a consumer the pool saturates, the flux goes to
  zero and every check below passes trivially.
- [x] 9.6 **This module's own check: tRNA conservation.** Assert the pair's sum
  equals its initial value at every save point over a full cycle against the
  double — verify by the check satisfying the tolerance principle, and by a
  mutation that creates tRNA rather than transferring it making it fail.
  **Annotated 2026-09-23:** asserts §3's first gate (bitwise `0.0` at all 106
  save points) with the `tol_C` ladder the same day's amendment allows, over
  six decades from (1e-4, 1e-2) to the pinned (1e-10, 1e-8): 12.2 to 4.75 orders
  below `tol_C`, 70× spread (job 17044321, all nine rungs reported). See §12.
- [x] 9.7 Show the flux is right without assuming it — verify by the steady
  charging flux reproducing 553.1 per second where that figure is computed
  independently, from 3,484,518 residues over 6,300 s, with `k_chg` reported as
  the derived quantity; and by a test asserting the check is not circular, since
  calibrating `k_chg` to 553.1 and then checking 553.1 verifies arithmetic (D14).
  **Amended 2026-09-23:** the first clause is circular under any demand double,
  so it is restated. Report the composed steady flux and the ATP it settles at
  against the independently derived 553.1 as a **self-consistency drift**, not
  as validation. The non-circularity test stays. See §12.
- [x] 9.8 Report the pool size's three consequences — verify by a diagnostic
  giving, for a range of pool sizes, the charging flux, the buffer against
  check 7's charged-tRNA counter, and the lag on the adenylate forward channel;
  and by the §9 open question on the pool size being replaced with a value that
  satisfies both check 1b and check 7, or by a statement that none does.
  **Annotated 2026-09-24: answered by task 11.7** — 0 clips per cycle at
  0.25 mM and above, 7 at 0.125 mM (job 17131232).
  **Annotated 2026-09-23: the check 7 half is deferred to task 11.7**, which
  already feeds its census back here. Check 7 counts clips of translation's
  charged-tRNA debit, and there is no translation until phase 11. This phase
  delivers the scan, check 1b's floor and the analytic buffer. See §12.
- [x] 9.9 Answer the stoichiometry question at matched steady flux — verify by a
  comparison run under both lumpings, matched on flux rather than event count
  (D14), recording the kinase traffic and the resulting ATP-to-ADP ratio the
  rebuild reads, and by the scoping note's open-questions entry being replaced
  with the measured answer.
- [x] 9.10 Suite and handoff — verify by `sbatch test/run_tests.slurm` passing,
  and by the handoff recording that phase 8's drain double is now superseded for
  composed runs and kept only for standalone ones.

### Phase 10 — Transcription

**Goal:** seventeen genes, one constitutive transcription reaction each, with the
corrected base mapping, the five cost counters, and rate constants that are
recomputable but never self-recomputed.
**Done when:** seventeen transcripts simulate over a full cycle as non-negative
integers, the time-averaged counts reproduce the measured means to §3's external
bound — Spearman at least 0.7 across the seventeen and within a factor of two for
at least fifteen of them — the seventeen rate constants fall in 1.26e-3 to
8.29e-3 per second, and the GTP elasticity reproduces 0.0079 to 0.0117 under the
corrected mapping.
**PR:** #51 (merged 2026-09-17)

Reference detail at
`dev/archive/openspec/changes/add-corea-transcription/tasks.md`. One amendment:
its assertion that simulating past the declared cadence leaves the constants
unchanged becomes "phase 4's driver refreshes them, and the module still never
refreshes itself".

- [x] 10.1 Vendor the per-gene extract from ~~three~~ **five** upstream sources
  (see amendment 2026-09-10) — verify by
  seventeen rows whose four base counts sum to the transcript length, by the
  totals matching A 7236, C 2078, G 3094, U 5868, and by the generator failing
  loudly on a missing locus rather than defaulting, since a silently absent gene
  shows up only as a model with sixteen transcripts. The extract also carries
  the transcript's **first two bases**, which the rate law reads as `C₁` and
  `C₂`, and which the archived `add-corea-transcription` design's D8
  header omits (that archive's D8, not §4's).
  **Done:** 17 rows; every gene's four base counts sum to its length; totals
  A 7236, C 2078, G 3094, U 5868. Four spot rows exact: `JCVISYN3A_0445`
  (1284; A 521, C 125, G 197, U 441; 266; 0.4403), `JCVISYN3A_0607` (1017;
  1355; 2.1781), `JCVISYN3A_0779` (2238; 831), `JCVISYN3A_0694` (270; 290).
  All seventeen copy numbers equal the scoping note's table, derived
  independently through the five-file chain. Two runs produced identical bit
  patterns. Removing `JCVISYN3A_0694` from `mRNA_counts.csv` aborts with
  "locus JCVISYN3A_0694 (ptsH) is absent from mRNA_counts.csv" rather than
  writing sixteen rows.
- [x] 10.2 Implement the seventeen reactions, genes carried as fixed quantities
  rather than states — verify by the state count being 17 transcripts plus 5
  counters, by none being a registry species, and by a recorded note that adding
  replication later must promote the genes to states.
  **Done:** `length(states(m)) == 22`, `!any(is_registered, states(m))`, 17
  reactions. Firing GAPD's raises its transcript by one and no other, `ATP_trsc`
  by 1017, and each monomer counter by that gene's base count; the propensity is
  the rate constant at any state and any time, which is zeroth order made
  assertable.
- [x] 10.3 Implement the rate constant with the corrected base mapping as default
  and the published permutation behind a keyword (D4) — verify by the correction
  appearing in `reduction_notes` while the published mapping does not, by the
  turnover cap never binding (the largest is 8.85 against a ceiling of 180), and
  by both mappings computable so the difference is one argument away.
  **Done:** the seventeen constants span **1.2579e-3 to 8.2909e-3** per second.
  Largest promoter-scaled turnover is GAPD's **8.8516** nt/s against the ceiling
  of 180, reported by `turnover_headroom` rather than assumed. Under
  `:corrected` one more cytosine costs more than one more guanine, and under
  `:published` the reverse — the permutation made observable rather than
  asserted. `reduction_notes` carries the correction; the published mapping
  does not appear there.
- [x] 10.4 Take all four nucleotide concentrations from balanced tables, never the
  first-minute setup constants — verify by the four values and their widths
  matching the balanced files, and by a recorded note that the nucleotide file
  repeats two of the setup constants, which makes the wrong choice look
  corroborated.
  **Done:** ATP 3.6529/1.2825 central, GTP 1.6627/1.5684, CTP 0.6874/2.0574,
  UTP 2.7681/1.3664 nucleotide, each reporting its file and `:balanced`. None
  equals its setup constant (1.04, 0.34, 0.68, 0.68). ATP and GTP agree with the
  registry; CTP and UTP have no registry value, which is why this module carries
  theirs.
- [x] 10.5 Expose one recomputation entry point and confirm the module never calls
  it — verify by simulating past the declared cadence with the constants
  unchanged, and by phase 4's driver changing them when it is present.
  **Done:** `rate_constants` reproduces the module's own law to 1e-12 at the
  declared pools and returns 17 values. Scaling all seventeen `k_tx` slots by
  1000 changes nothing it returns, so the law cannot compound its own previous
  value; doubling one promoter strength doubles exactly that gene's constant and
  leaves the other sixteen bitwise unchanged.
- [x] 10.6 Declare the five deferred counters and the two clamped nucleotide pools
  with origin ours — verify by the counter table matching what each drains into,
  by all five taking the published clamped policy so none is a labelled
  deviation, and by three kinds coexisting on ~~ATP~~ **CTP and UTP** inbound as
  the contract requires (ATP carries three edges of two kinds; see the
  2026-09-10 amendment).
  **Done:** exactly 11 edges — 4 rate-constant at 60 s, 5 deferred counters, 2
  clamped at `:ours` holding 0.6874 and 2.7681. `inputs(m)` is empty.
  `resolve_coupling` accepts two counters plus a rate-constant edge on
  `M_atp_c` inbound, and counter-plus-rate-constant-plus-clamp on CTP and UTP.
  None of the five counters deviates from published; both clamps do.
  **Standalone only.** That resolution is `resolve_coupling([m])`; the CTP and
  UTP declarations cannot be satisfied by any *composition*, which amendment E
  records and a test pins.
- [x] 10.7 Register both promoter-proxy declarations (D5) — verify by
  `reduction_declarations` returning one entry saying it is a proxy and a second
  naming the circularity and the seventeen parameters affected, and by a test
  asserting this module's copy numbers agree with the metabolic modules', since
  one number has three consumers. **That cross-check is not independent of
  phases 6 and 7.** In a fan-out it belongs to whichever pull request lands
  later; carry it here as a skipped test naming them.
  **Done:** two declarations, the proxy and the circularity, the second naming
  all seventeen `k_tx_*`. Values match to three decimals — PGI 1.478, GAPD
  7.528, ptsG 4.617, GK1 1.033 — and GAPD/PGI equals 1355/266 exactly. The
  cross-check against the metabolic modules is carried as a skipped test naming
  phases 6 and 7, not stubbed.
- [x] 10.8 Simulate a full cycle **against a transcript-decay double** and report
  the elasticities — verify by
  non-negative integer counts throughout, by the time-averaged counts meeting
  §3's external bound (Spearman ≥ 0.7, and within a factor of two for at least
  fifteen of seventeen), by the elasticity to all four pools falling in 0.044 to
  0.051, and by the GTP-alone elasticity being reported under both mappings so
  the 1.9× effect is measured rather than argued. **The double is the acceptance
  criterion, not scaffolding** — without a consumer the transcripts only
  accumulate and every average is an order of magnitude high, so the check would
  be vacuous. It uses the published per-gene law, derivable from the extract's
  own length column; phase 12 supersedes it for composed runs and keeps it for
  standalone ones. See amendment 2026-09-10.
  **Done:** eight replicates over a 6,300 s cycle against the decay double.
  Counts non-negative integers throughout and the five counters monotone, since
  nothing debits them without a driver. Spearman **0.8701** against the measured
  means and **16 of 17** within a factor of two — FBA the outlier at **3.67×** on
  the simulated time averages (the deterministic `k_g / k_deg` band below puts
  it at 3.96×; both are reported, neither is the other), for the reason
  amendment 2026-09-10 B records. Elasticity to all four pools
  **0.043742 to 0.050579** against 0.044 to 0.051; GTP alone **0.0078915 to
  0.0117062** corrected against 0.0079 to 0.0117, and **0.015617 to 0.019536**
  published against 0.0156 to 0.0196. Per-gene ratio of the two spans 1.53 to
  2.45, median **1.85** — D4 states the effect as a flat 1.9×, where it is in
  fact the per-gene U-to-G count ratio and 1.85 is its median. Suite: job
  16614893 on the rebased tree, **2513 passed, 0 failed, 2 broken**, 4m40.4s
  (2265 before the phase, which is phase 8's recorded count on this tree; the two
  broken are the repo's only two `@test_skip`s — phase 8's and task 10.7's).
  The pre-rebase figure of 1708/1708 against a 1468 base belongs to job
  16373862 and to a tree without phases 6, 7 and 8.

### Phase 10b — Framework and phase 10 fixes from the 2026-09-22 review

**Goal:** fix the review findings that are small, independent and outside any
module phase's remit (§12, 2026-09-23 C to F and H), so that phase 13 meets
none of them.
**Done when:** every returned ABC particle lies within the reported tolerance,
shared priors are compared by value, a named peer must supply the quantity, an
empty boundary is refused, and transcription's states and provenance follow its
own configuration.
**PR:** #57 (merged 2026-09-23)

Independent of phases 11 and 12, so it can land at any point before phase 13.

- [x] 10b.1 Report the tolerance the returned population was accepted under —
  verify by a seeded run in which every returned particle's distance is at or
  below the reported tolerance, including at `n_populations = 1`, where
  population 1 must itself be filtered or the call refused.
  **Done:** `abc_smc` returns the threshold its final population was accepted
  under. The test uses a deterministic simulator, so every distance can be
  recomputed, and runs at `n_populations` 2 and 5. Every returned particle is
  at or below `tolerance`, and the largest is above `tolerance / 2`, so the
  label is tight rather than inflated. `n_populations = 1` is **refused**
  with an `ArgumentError`, not filtered.
- [x] 10b.2 Compare shared priors by value, not by type — verify by
  `LogNormal(0, 0.1)` against `LogNormal(5, 2)` under one name throwing and
  naming both modules, and by identical priors still composing.
  **Done:** Distributions' `isequal`, which compares family and fields. Under
  `:k_tx`, the pair throws, naming `:txl`, `:metab` and both priors.
  `LogNormal(0, 0.1)` against `LogNormal(0.0, 0.1)` composes.
  `truncated(Normal)` against `truncated(LogNormal)` throws. An earlier
  `nameof`/`params` fallback let that pair through, and `/check-PR` caught it.
- [x] 10b.3 Require a named peer to own or supply the edge's quantity — verify by
  an owner, a consumer and a bystander, with the edge naming the bystander,
  throwing and naming the actual owner; and by every existing module's edges
  still resolving.
  **Done:** a named peer must integrate the species or declare the opposite
  direction of the same crossing. The bystander case throws with "integrated
  by Central". Naming the owner resolves, and so does a producer naming its
  consumer. No Core A′ module names a peer, and all five named peers in the
  tests pass. A clamp or an inbound volume edge has no counterpart, so its peer
  stays unnamed, which the `CouplingEdge` docstring now says.
- [x] 10b.4 Refuse an empty boundary in `iterative_infer` — verify by two blocks
  with no shared free parameter throwing before any sampling, instead of
  returning `n_iters = 2` with a zero KL.
  **Done:** `TranscriptionTranslation([:g1])` against
  `StochasticGeneExpression()` throws an `ArgumentError` ("share none") at a
  guard placed before the iteration loop. The message points a state-coupled
  composition at D10.
- [x] 10b.5 Make transcription's states and reduction notes follow its
  configuration — verify by a three-counter construction exposing exactly three
  counter states, and by `reduction_notes` under `base_mapping = :published`
  describing the published permutation, not the corrected mapping.
  **Done:** a three-counter build (`ATP_trsc`, `GTP_mRNA`, `UTP_mRNA`) has 20
  states, with those three at 18:20. **Wider than the review said:**
  `reactions` wrote five hard-coded positions and `counter_drains` returned
  the default five, so both now follow the configured counters too. Firing
  GAPD writes its length, G count and U count to 18:20. An unknown counter
  is refused, and so is a counter charged against the wrong pool. Under
  `:published` there is one mapping note, it says "is the published
  permutation", and no note claims the correction. The ATP double-debit note
  appears only when both ATP counters are configured.
- [x] 10b.6 Replace the 10.7 `@test_skip` with the executed check — verify by the
  seventeen transcription loci's copy numbers agreeing with those phases 6 and 7
  declare wherever a locus appears in both, and by the suite's broken count
  falling from 2 to 1.
  **Done:** equal to phase 6's per-locus `copies`. Phase 7's phospho-state
  pairs times particles per mM agree to rtol 1e-9, and phase 8's per-locus
  `copies` match too, which is what covers ADK1, PPA and GK1. The set of
  compared loci equals all seventeen. The broken count went **2 → 1**; the one
  left is phase 8's skip, waiting on phase 9.
- [x] 10b.7 Suite and handoff — verify by `sbatch test/run_tests.slurm` passing.
  **Done:** job 17021441 at `b04d467` gave **2580 passed, 0 failed, 1
  broken** in 4m56.0s, against 2513 and 2 broken before the phase (job
  16614893). Job 17019818 at `f611345`, before the review fixes, gave 2577,
  0 and 1. `docs/handoff.md` is updated. ~~CI didn't run, because the Actions
  quota is exhausted until 2026-10-01.~~ **Corrected 2026-09-23:** CI did run.
  The Test workflow passed on this branch (runs 35793639222 and 35795062475),
  because the repo had gone public and public repos run Actions free.

### Phase 11a — Catalytic edges on jump-owned protein counts

**Goal:** give the translation phase a channel to write into. A catalytic edge
may name a protein count the jump block owns, which is not a registry species,
and the two metabolic modules whose enzymes translation makes expose those
enzymes as slots the edge can fill.
**Done when:** a catalytic edge on a non-registry protein count resolves and
executes, every other edge kind still refuses a non-registry species,
`CentralGlycolysis` and `NucleotideRecycling` each have an opt-in mode whose
enzyme concentrations are filled by catalytic edges from `P_<locus>` counts,
their default mode is unchanged, and a jump module's producer counter is shown
to raise a PTS carrier and, through it, the cell's surface area.
**PR:** #62 (merged 2026-09-24)

Added 2026-09-23 because phase 11 could not be built as written (§12, same
date). Framework and two module files, so §10 R15 puts it in its own pull
request, as 5b and 10b were.

- [x] 11a.1 Let the resolver accept a `CatalyticEdge` on a non-registry species —
  verify by a toy hybrid in which a catalytic edge on a jump-owned `:P_toy`
  resolves, lowers and fills its slot at every handshake; by a mass, currency,
  deferred-counter, rate-constant, volume or clamped edge on the same name still
  throwing; and by a catalytic edge on a count no jump module owns throwing at
  build time and naming it.
  **Done:** `ToyMembranePool`'s edge on `:P_toy` resolves with `state_index`
  0 and no dead end. Mass, currency, deferred-counter, rate-constant, volume
  and clamped edges on `:P_toy` each throw, naming it and "Only a
  CatalyticEdge". `:P_missing` with no owner throws at `build_problem`, naming
  it. Composed with a `ToyProteome` making `:P_toy` at 2/s, the slot equals
  `counts_to_mM` of the count read at that handshake.
  **Added on review (/check-PR of #62):** before this, the two translated
  modules composed with no jump block built as a plain `ODEProblem` with 15
  free enzyme slots that nothing ever filled. A build with no jump block now
  refuses any `CatalyticEdge`, naming the slot and the count. The nominal pair
  is still refused through its inputs, and `ToyPool` without the edge builds.
- [x] 11a.2 Give `CentralGlycolysis` a translated-enzyme mode — verify by its ten
  enzyme concentrations becoming free slots `enz_R_<reaction>` at the nominal
  values, by ten catalytic edges from `P_<locus>` to those slots, by the protein
  counts leaving `inputs()` in that mode only, by the derivative at the registry
  state being bitwise equal to the default mode's, and by the default mode's
  parameters, edges and inputs being unchanged.
  **Done:** `enzymes = :translated` gives exactly ten free parameters,
  `enz_R_PGI` to `enz_R_LDH_L`, equal to `enzyme_conc` and `:asserted`. Ten
  `CatalyticEdge`s run from `default_protein_sources()`, the nine boundary edges
  are unchanged, and `inputs` is the three currencies. `dynamics` and
  `contributions` at the registry state are `===` the default mode's. Doubling
  `enz_R_GAPD` doubles R_GAPD and leaves the other nine bitwise unchanged. The
  default mode has no free parameter and no catalytic edge, and keeps the ten
  counts in `inputs`. An unknown mode is refused.
- [x] 11a.3 Give `NucleotideRecycling` the same mode — verify by its five
  `enz_R_*` slots freed and filled by catalytic edges, by PGK3 and PYK3 reading
  the same `P_<locus>` counts as glycolysis's PGK and PYK, and by one count
  filling both modules' slots in a hybrid build.
  **Done:** five free `enz_R_*` slots, each with a catalytic edge from
  `P_<locus>`. The PGK3 and PYK3 counts are the ones glycolysis's PGK and PYK
  name. Built with both translated modules and a 13-count `ToyProteome`, the
  driver lowers **15** catalytic exchanges. After one handshake every slot
  equals its count's `counts_to_mM`, PGK3 equals PGK, PYK3 equals PYK, and
  ADK1 equals its nominal concentration to rtol 1e-12. **Added on review:**
  a reaction disabled through `reactions` is neither freed nor wired. Without
  ADK1 the mode frees and wires exactly PGK3, PYK3, GK1 and PPA.
- [x] 11a.4 Pin the carrier path translation will use — verify by a jump module's
  producer `DeferredCounterEdge` raising a PTS carrier state by exactly its
  accrual per drain, and by the volume chain reading the larger surface area
  through the carrier owner's existing membrane flag.
  **Done:** a `ToyProteome` making `:P_toy` at 5/s credits `M_ptsi_c` in
  `ToyMembranePool` through a producer counter. Over 100 handshakes, each drain
  raises the carrier by exactly the pending accrual, in particles, to rtol
  1e-12, and carriers plus pending equal the initial carriers plus proteins
  made. `growth_census` counts the credited carriers through the owner's flag,
  and the area grows. No framework change was needed for this path.
- [x] 11a.5 Suite and handoff — verify by `sbatch test/run_tests.slurm` passing.
  **Done:** job 17057406 at `336e504`, **2976 passed, 0 failed, 0 broken**,
  6m37.6s. CI on PR #60's branch gave 2803 for phase 12 (run 35814484827). The
  first run, job 17057084 at `0121176`, errored once in 11a.1. The toy kept its
  membrane flag after its volume edge was overridden away, and the driver
  correctly refused that. The fix was in the test (`336e504`). After the
  review fixes: job **17127107** at `5e25125`, **2985 passed, 0 failed, 0
  broken**, 7m55.2s.

### Phase 11 — Translation

**Goal:** one translation reaction per transcript plus ptsG translocation, with
the rate constant built from the live lumped charged-tRNA pool, publishing the
seventeen protein counts the ODE block and the growth channel read.
**Done when:** seventeen protein counts simulate over a full cycle as
non-negative integers, the energy counter equals exactly twice the residues
translated, and the counts feed the metabolic modules' enzyme concentrations
through executed catalytic edges rather than nominal stand-ins.
**PR:** #64 (merged 2026-09-24)

- [x] 11.1 Extract per-gene amino-acid counts and residue totals for the
  seventeen loci from the genome record — verify by the residue counts summing to
  3,484,518 for a full proteome doubling, matching the scoping note's own figure,
  and by four spot lengths matching the recorded ~~746, 574, 155 and 90~~
  **745, 573, 154 and 89** residues. **Amended 2026-09-23:** 3,484,518 is
  Σ copies × (length/3 − 1), so it excludes the stop codon, and the recorded
  spot lengths include it. Residues exclude the stop throughout, which is the
  convention `k_chg` was calibrated on in phase 9. Upstream charges GTP on
  `len(aasequence)`, which counts the `*`, so it charges two more GTP per
  protein; that difference is recorded as ours (§12).
  **Done:** the gene extract gains `!Residues`, each transcript translated
  under NCBI table 4 as `MinCell_CMEODE.py:121` does, less its one stop; the
  generator asserts a single terminal stop, so it equals `length/3 − 1` by
  translation. Σ copies × residues = **3,484,518**; ptsG, ptsI, Crr, ptsH have
  **745, 573, 154, 89**. The first nine columns are byte-identical and two runs
  give identical files. **Per-amino-acid counts are not extracted** (amendment
  2026-09-24 C): under the lumped law they drop out, since Σ n_aa = residues.
- [x] 11.2 Implement seventeen jumps catalytic in the transcript count, reading
  transcripts as a phase-2 peer state — verify by firing one gene's reaction
  raising only that gene's protein by one, leaving the transcript unchanged, and
  raising the energy counter by exactly twice that gene's residue count.
  **Done:** `CoreATranslation`, 18 jumps (17 plus translocation), each
  `k_tl_g · mRNA_g` read through `inputs`. For all seventeen genes one firing
  raises that gene's count by one (ptsG's cytosolic count for ptsG), leaves
  the transcripts unchanged, and adds exactly 2r_g to `GTP_translat` and r_g
  to `tRNA_translat`. It composes jump-only with `CoreATranscription`.
- [x] 11.3 Implement the rate constant with the lumped pool substituted for the
  twenty per-amino-acid pools — verify by the three polymerase-capacity constants
  carried with the note that all are computed at initial volume and are among the
  genuinely frozen quantities, and by the lumping **not** being registered here:
  phase 9 owns that declaration after D13, and registering it twice would
  double-count it in phase 13's enumeration.
  **Done, with two choices approved 2026-09-24 (amendment A, B):** the restart
  law (`translation_rate_restart.py`: riboKcat 12, riboKd 1e-3, kcat_mod
  (0.25n + 0.2)·riboKcat), and each of the twenty-one tRNA concentrations read
  as the lumped pool's per-amino-acid share, [M_trna_chg_c]/20. Hand-checked at
  GAPD to 1e-14, and equal to the upstream law at upstream's own 150-copy
  pools. **Of the three constants only the ribosome concentration is computed
  at initial volume** (503 copies, as upstream computes `ribosomeConc` once);
  K₀ and K_d are dissociation constants in mM. All three, and riboKcat, are
  fixed `:asserted` parameters citing the restart file. `reduction_notes`
  carries the share and the restart law, and no lumping declaration.
- [x] 11.4 Add the ptsG translocation reaction — verify by only that locus
  carrying it, by translocation being what increments the state phase 5's volume
  edge reads, and by no cytosolic protein having it.
  **Annotated 2026-09-23:** `PtsTransport` owns that state and its volume edge,
  so translocation increments it through a producer deferred counter crediting
  `M_ptsg_c` (task 11a.4), not by owning it. The other three carriers are
  credited the same way. Upstream translocation also charges `int(len/10)` ATP
  (`ATP_transloc`, `MinCell_CMEODE.py:1012`), which this task did not name:
  declare it as a counter or record its omission as ours (§12).
  **Done:** one translocation jump, `50/746 · Pcyto_ptsG`, and only ptsG has
  it. A firing moves one ptsG into `P_JCVISYN3A_0779`, accrues **74**
  `ATP_transloc` (declared, debiting ATP, products pending 13.10) and one
  `ptsG_transloc`, which credits `M_ptsg_c`. No other reaction lowers a
  protein count. In the seven-module composition, over 120 handshakes, PTS's
  ptsG pair plus the hook's carried remainders rises by exactly the pending
  credit at every handshake (rtol 1e-9), and the area grows through PTS's flag.
- [x] 11.5 Declare the boundary, **rewritten by D13 and by the defect below**: a
  60 s inbound rate-constant edge on the charged pool, a deferred counter on GTP
  under the published clamped policy, **a paired deferred counter debiting the
  charged pool and crediting the uncharged one**, seventeen outbound catalytic
  edges naming each rate law's parameter slot, and a volume edge on ptsG
  (**amended 2026-09-23:** the catalytic edges are declared by the consuming ODE
  modules through phase 11a's mode, since only the consumer has a slot, and
  there are thirteen distinct enzyme counts behind fifteen slots, not seventeen;
  the four PTS carriers have no enzyme slot and take producer counters instead;
  translation declares no volume edge, because `PtsTransport` owns that one; and
  the GTP counter debits GTP alone until task 13.10 gives it GDP and Pi —
  see §12) —
  verify by the exact edge count, kinds and directions; by the catalytic edges
  carrying no mass, with an attempt to include one in a conservation check
  rejected rather than counted as zero; and by no tRNA edge being declared as
  *mass*, since a jump module cannot continuously write an ODE state, which is
  the same constraint D13 applies to charging's currency edges.
  **Done:** nine edges, exactly. One 60 s inbound rate-constant edge on
  `M_trna_chg_c`, and eight counter channels: `GTP_translat` debiting GTP;
  `tRNA_translat` debiting `M_trna_chg_c` and crediting `M_trna_c`;
  `ATP_transloc` debiting ATP; four carrier credits. No catalytic, volume, mass,
  currency or clamped edge. Composed with all seven modules the driver lowers
  15 catalytic exchanges from 13 of translation's counts, the eight channels
  and one rebuild. All 15 catalytic edges refuse `mass_contribution`.
- [x] 11.6 **Fix the missing debit, which is a defect in the committed spec.**
  As originally written this module credited uncharged tRNA and never debited the
  charged pool: it read the charged pool through a 60 s rate-constant edge rather
  than consuming it. That is unbalanced — with charging producing and nothing
  consuming, all tRNA charges, the uncharged pool empties, the charging flux
  collapses and the composed model has no steady state. Verify by a composed run
  with phase 9 reaching a steady tRNA split rather than a fully charged pool, and
  by a test that removing the debit reproduces the collapse, so the fix has a
  witness rather than a claim.
  **Done, with one correction (amendment 2026-09-24 D):** composed with
  phase 9's modules and live transcripts, the split is a stochastic steady
  state. The test asserts each 300 s mean over 1,200 s within (0.6, 0.95).
  An ad hoc full-cycle run on the login node, recorded in no artefact, put
  the 300 s means at 0.71 to 0.90 with no trend. tRNA plus carries equals 0.25 mM to
  1e-9. **With no transfer at all the pool collapses**, fully charged with
  uncharged tRNA at 2e-156 mM. **Removing only the debit does not collapse**:
  the credit then creates tRNA, 276,128 particles in 600 s. Both are pinned.
- [x] 11.7 Report whether the charged-tRNA counter can clip — verify by the
  census of check 7 run on this counter specifically, since it debits ~553
  residues per second against a pool of order 10³ particles and is now the most
  likely clip in the model (K5); and by the result feeding phase 9's pool-size
  diagnostic rather than being reported in isolation.
  **Done** (`dev/scripts/translation_diagnostics_result.md`, job 17131232).
  With translation as the consumer in phase 9's composition, over one full
  cycle per pool, the counter clips on **0** of 6,300 drains at 0.25, 0.5 and
  1.0 mM, **7** at 0.125, 17 at 0.1, 253 at 0.05 and 1,009 at 0.025. At
  0.25 mM the charged pool's minimum is 1,144 particles. The suite pins zero
  clips over 1,200 drains at 0.25 mM, and 77 of 120 at 0.01 mM so the census
  can fire. **The first full seven-module cycle clips it on 911 drains**,
  with ATP falling to about 0.44 mM: an early warning for K5, carried to 14.8
  (amendment E).
- [x] 11.8 **This module's own check: residue-to-energy closure.** Assert the
  accumulated counter equals exactly twice the residues translated at every write
  point — verify by the assertion passing and by a mutation halving one gene's
  residue count failing and naming that gene.
  **Done:** `assert_residue_energy_closure` takes its residues from the
  extract, not the module. Over a full jump-only cycle with transcription and
  decay it is exactly zero at all 106 save points, and every protein count is
  a non-negative integer. Halving GAPD's residues gives −2 × made × 169 and
  an error naming `JCVISYN3A_0607 (GAPD)` and no other gene.
- [x] 11.9 **Resolve whether protein degradation is a reaction at all, then
  implement or drop it.** The published model builds three reactions per locus
  plus one translocation, which is the 52 the scoping note counts, and carries
  `ptnDegRate` as a constant with no matching per-gene reaction in that list. So
  adding seventeen degradation jumps is **ours**, and the committed spec added
  them without saying so — verify by a written determination citing where
  `ptnDegRate` is consumed upstream; if it is a reaction, implement it with the
  count in §2 corrected and the departure labelled if it falls outside the
  published 52; if it is not, drop it and record why. Either way §2's count and
  the count the tests assert must agree. If kept, degradation returns nothing to
  the chemostatted amino-acid pool, and that exemption is recorded so a later
  change making the pool live reinstates the check.
  **Determined 2026-09-23 (§12):** not a reaction. `ptnDegRate` is assigned at
  `MinCell_CMEODE.py:273` and `:343`, `translation_rate_start.py:29`,
  `translation_rate_restart.py:30` and `MinCell_restart.py:440` and `:513`, and
  no reaction reads it; `ptnDegProd = []` at `MinCell_CMEODE.py:1006` is never
  used. Dropped, and §2 corrected. What remains for this task is the test
  asserting 52.
  **Done:** 17 + 17 + 18 = 52 reactions across transcription, decay and
  translation, asserted; no translation reaction lowers a protein count.
- [x] 11.10 **Settle two open questions from §9 here.** Measure the elasticity of
  the translation rate constant to the charged pool, and the protein fold change
  over a full cycle — verify by both numbers recorded against the predictions
  (a reverse channel bounded at a few percent, and a fold change that may fall
  short of two), and by either shortfall diagnosed to the charging step or the
  ribosome constant *in this phase*, since R1 and R12 both say diagnosing here is
  far cheaper than explaining at the end.
  **Annotated 2026-09-23:** `riboKcat` is 12 in the two translation-rate files
  that run and 10 in `MinCell_CMEODE.py:332`; `riboKd` is 1e-4 in the start
  file and 1e-3 in the restart file. Record which value the module uses and
  why before either number here is quoted.
  **Done; the diagnosis is "neither" in the jump-only runs and unmeasured in
  the hybrid (amendment 2026-09-24 F):**
  riboKcat 12 and riboKd 1e-3, the restart law (11.3). **Elasticity
  0.0909 to 0.0911**, about 1.9× transcription's 0.044 to 0.051. So the
  reverse channel is **stronger** than predicted, and R1 does not fire.
  **Fold change, two numbers, and they differ.**
  - **Jump-only**, eight seeds, constants held at the nominal pool: median
    **1.691**, 1.604 to 1.822, none below 1.5, log-log slope against length
    −0.012. Here it is neither candidate. The pool is nominal by construction,
    and the ribosome constant is the published 12. Reaching 2 needs 1.45× more
    translation output (riboKcat ≈ 17.4). The **untested** leading candidate
    is that Core A′ cuts replication, so gene dosage never rises; a 1.5× mean
    dosage gives 2.04.
  - **Full seven-module hybrid**, one seed: median **1.497**, minimum
    **1.232**. That fails §3's median and its no-gene-below-1.5 rule. Here
    charging plausibly does throttle translation. The charged pool reads
    0.0000 mM at some rebuild instants: alternate ones among the last five
    printed samples, with how many over the whole cycle not recorded. After
    such an instant the floored law puts every constant about 20× below
    nominal for the next 60 s.
  - **How much of the drop from 1.691 to 1.497 is due to that is not
    measured.** It needs a hybrid rerun with the constants frozen at 0.2 mM,
    or with a 6 s rebuild.
  - The result file's sentence "The charged pool is not what holds it back",
    and its 1.53× and riboKcat ≈ 18.4, are **superseded** by this entry. They
    come from the analytic baseline, and the file is left as job 17131232
    wrote it.
  - Carried to phase 14's F5, with §3's bound unchanged.
- [x] 11.11 Suite and handoff — verify by `sbatch test/run_tests.slurm` passing and
  by the handoff recording that the metabolic modules' nominal enzyme
  concentrations are now superseded by live counts.
  **Done:** job **17131233** at `9cd7b66`, **27,305 passed, 0 failed, 0
  broken**, 9m40.8s, against 2,987 on the merged phase 11a tree (job
  17127716). `docs/handoff.md`
  records the supersession.

### Phase 12 — Transcript decay

**Goal:** one decay reaction per transcript at the published length-dependent
global rate, returning the four monomers to the pools that can accept them.
**Done when:** each decay returns exactly its own base counts, the two
recyclable monomers reach the pools the recycling module closes, the two
chemostatted monomers are exempt with the exemption recorded, and the guanylate
return over a cycle ~~matches the scoping note's ~24,000 particles~~ **matches
the analytic `Σ_g k_g·T·G_g` — 61,531 particles, about 1.5× the guanylate pool —
within Monte Carlo error** (amended 2026-09-23; see §12).
**PR:** #60 (merged 2026-09-23)

- [x] 12.1 Implement seventeen decay jumps with the published single global
  constant over transcript length — verify by exactly one rate parameter across
  all seventeen genes, by each half-life being a function of length alone, and by
  the upstream hand-tuning comment recorded rather than silently inherited.
  **Done:** `CoreATranscriptDecay` has one free parameter, `:krnadeg` =
  `(18/452)·88` = **3.5044** nt/s, across 17 reactions; `k_g·n_g` equals it for
  every gene, and the half-lives sort exactly by length. `RNADEG_KCAT`'s
  docstring carries the upstream `# INSTEAD OF 18 or 20` and records that
  upstream's variable *named* `krnadeg` is the dead `0.00578/2` at line 272,
  beside the equally dead `rnaDegRate` at line 295. Phase 10's
  `transcript_decay_constant` is pinned equal to it.
- [x] 12.2 Decrement the transcript this module does not own, through phase 2's
  declared peer write — verify by firing decay lowering only that transcript, by
  it being unable to fire at zero copies, and by the write throwing once the
  ~~outbound edge~~ **`written_states` declaration** is removed (amended
  2026-09-04; see §12 — a transcript is not a registry species, so no edge can
  name it).
  **Done:** firing genes 1, 5 and 17 lowers that transcript alone; every
  propensity is exactly `0.0` at zero copies, and a 20,000 s run from one copy
  each ends at zero with no count ever negative. With `written = Symbol[]` the
  first firing throws, naming `CoreATranscriptDecay`, the transcript and
  `written_states`.
- [x] 12.3 Return the four monomers per firing — verify by the returned counts
  equalling that gene's four base counts exactly, cross-checked against phase 10's
  extract so the two modules cannot disagree on a base count.
  **Done:** for all seventeen genes one firing writes exactly `[n_g, A, G, C,
  U]` to the five counters, read from `transcription_genes(CoreATranscription())`
  — the module's genes default to that extract, so there is one source.
- [x] 12.4 Declare outbound counters on the two recyclable monomers and take the
  chemostat exemption on the other two — verify by `resolve_coupling` reporting
  the exempt pair among its chemostat exemptions with the exemption recorded.
  **Annotated 2026-09-23:** Core A′ has no CMP or UMP species, so "the other
  two" credit the **CTP and UTP chemostats**, which absorb them; that is ours and
  `reduction_notes` says so. No hybrid can build with those two counters until
  task 13.9 (§12).
  **Done:** five `DeferredCounterEdge`s — `ATP_mRNAdeg` `:in` on ATP; AMP, GMP,
  CTP, UTP `:out`. Resolved beside a transcript-source double — not beside
  transcription, whose own CTP/UTP counters would supply the exemptions
  whatever decay declared — decay alone puts `M_ctp_c` and `M_utp_c` among the
  chemostat exemptions, and its four credits resolve on AMP, GMP, CTP and UTP;
  the hybrid throw is pinned. A counter aimed at the wrong pool, a duplicate
  or an unknown name is refused at construction.
  **Executed, not only declared:** with the ATP, AMP and GMP counters, decay
  composes with `NucleotideRecycling`, the phase 8 glycolytic double and a
  transcript-source double; over 600 handshakes 75 decays returned 13,162 G and
  the guanylate pool rose by **13,162.08** particles, while adenylate moved by
  **−46,929.00** against the −46,930 its AMP credit and uncredited ATP debit
  predict. **Superseded 2026-09-24 by phase 13a:** the pinned throw is replaced by a build with all
  five counters, the chemostats' census matches the C and U returned, and
  adenylate now moves by the AMP credit alone, since `ATP_mRNAdeg` credits ADP.
- [x] 12.5 Add the decay energy cost counter the published hook drains — verify
  by it accruing per firing and by the composed model reporting which registry
  species it debits.
  **Done:** `ATP_mRNAdeg` accrues `n_g` per firing (one ATP per nucleotide,
  `in_out.py:178`); `counter_drains` and the resolved edge report `:M_atp_c`,
  declared by `CoreATranscriptDecay`, with products ADP and Pi recorded and not
  yet credited — task 13.10 now covers it (§12). **Superseded 2026-09-24 by phase 13a:** the products
  are credited, one each per ATP paid.
- [x] 12.6 **This module's own check: monomer closure.** Assert that over a
  composed transcription-and-decay run the total monomers returned equals the
  total bases polymerised, per moiety — verify by the assertion passing and by a
  mutation to one gene's returned counts failing and naming the moiety.
  **Done:** read as *polymerised + initial stock = returned + final stock*,
  since transcripts alive at either end sit on one side only. Over eight 6,300 s
  jump-only runs the residual is **exactly 0** for A, C, G and U, and
  `ATP_mRNAdeg` equals the monomers returned. One extra guanine in GAPD's
  returned counts gives a negative G residual, zero elsewhere, and
  `assert_monomer_closure` throws naming G alone.
- [x] 12.7 Record the leak this module closes — verify by the composed model
  reporting the guanylate return over a cycle against the guanylate pool, about
  ~~60%~~ **150% (amended 2026-09-23; see §12)**, as the recorded reason the
  guanylate kinase is in the core.
  **Done:** **61,219** GMP returned per cycle over eight replicates, against the
  analytic **61,531** and a guanylate pool of **39,806** particles — a ratio of
  **1.54**. Without GK1 one cycle strands more guanylate than the whole pool.
- [x] 12.8 Suite and handoff — verify by `sbatch test/run_tests.slurm` passing.
  **Done:** job 17023839 at `a6b7abd`, **2703 passed, 0 failed, 1 broken**,
  5m20 (2580 before the phase, 10b's recorded count; the 123 added are this
  phase's; the one broken is the repo's remaining `@test_skip`). The
  `/check-PR` fixes added five assertions (128 in the phase file), run on the
  login node and in the PR's CI, not in a further Slurm job. CI on the tree merged
  with phase 9 (run 35814484827): **2803 passed, 0 broken**.

### Phase 13a — The framework fixes assembly needs

**Goal:** remove the three framework gaps that stop the seven modules composing
into a model whose moieties can balance, so that 13b assembles and measures the
real model rather than one with GTP starved from 44 s.
**Done when:** transcription and recycling build together with the chemostatted
CTP and UTP taking ownerless debits and credits; every deferred counter credits
its products per unit actually paid; and a mixed composition either dispatches
on the composition or throws a named error.
**PR:** #66 (merged 2026-09-24)

Split out of phase 13 on 2026-09-24 (§12, same date). Framework files, so §10
R15 puts them in their own pull request, as 5b, 10b and 11a were. The tasks are
13.9, 13.10 and 13.3 as written in 13b's list, which keeps their ids and their
annotations; they are ticked there.

- [x] 13a.1 = task 13.9, the ownerless path for chemostatted pools. Its
  two-module build carries phase 8's `HeldGlycolytic` double as a third module
  (§12, 2026-09-24 C).
- [x] 13a.2 = task 13.10, per-product stoichiometry on deferred counters. Its
  last clause, task 13.2 seeing every product edge executed, is verified in 13b,
  where 13.2 is built.
- [x] 13a.3 = task 13.3, hybrid inference dispatch.

**Done (#66):** Slurm job 17163832 at `20dce20`, **27,382/27,382**. The later
commit `54f4b23` touched only prose; CI's unit tests passed on it.
- 13.9: `NucleotideRecycling + HeldGlycolytic + CoreATranscription` builds. The
  rebuild holds CTP at 0.6874 and UTP at 2.7681, and the 60 s refresh
  reproduces `rate_constants` at those values to rtol 1e-14. Decay with all
  five counters composes, and its chemostat census equals Σ fired × C and
  Σ fired × U. Two pinned throws are replaced.
- 13.10: on a toy with ATP clipping (100 paid of 250), ADP and Pi each rise by
  100 and both moieties balance to 1e-6 particles. With transcription and
  `ATP_trsc` forced to clip, adenylate stays within 2 particles. In decay's
  executed-credit test, the adenylate change now equals the AMP returned.
- 13.3: a hybrid is refused by name in either module order, and mixed
  inference modes are refused.
- /check-PR added two things (§12, 2026-09-24 E and F): products are refused
  on a counter with two consumers, and inert clip fields are no longer
  reported. A pool owner that mirrors a product credit as `:in` would still
  be reported as clipping. No module does that; it is latent.

### Phase 13b — Assemble Core A′ and assert structural completeness

Was phase 13 until the 2026-09-24 split. Task ids are unchanged; 13.3, 13.9 and
13.10 land in 13a.

**Goal:** compose all seven modules into one hybrid model and make the assertion
the wave-0 contract deliberately withholds.
**Done when:** every registry dynamic state is owned, no dead ends remain, every
declared edge is executed, the model runs a full 6,300 s cycle, and the
per-trajectory wall-clock is recorded.
**PR:** #68 (merged 2026-09-25)

**Done (#68):** Slurm suite **27,453/27,453** at `a777a03` (job 17315831);
CI unit tests passed.
- 13.1: `build_corea()` resolves in completeness mode with nothing unowned and
  no dead ends. Four PTS carriers are reported as `accumulating`. Dropping
  recycling or charging fails, naming the states and the stranded moieties.
- 13.2: 0 of 89 declared edges unexecuted. An inert outbound currency edge
  fails the build, naming it. Every product edge of the four product-bearing
  counters is credited.
- 13.4: the drain stays at the published 1 s. The worst final paired
  difference is −52.9% ± 9.6% at 5 s and +288% ± 29% at 60 s (job 17298253).
- 13.5: a full 6,300 s cycle at Rodas5P, 1e-10 and 1e-8 with no solver
  failure, read back by `solver_settings`. It runs in the default suite.
- 13.6: the lumping is reported once and its formalism apart. The fourteen
  are among 52 asserted priors, with 16 discarded driver-written priors.
- 13.7: 32.2 s warm per cycle, 5.11 ms per handshake, plus a 28.1 s build and
  an 81.8 s one-off compile (job 17298242).
- /check-PR verdict MERGE, with nine MINOR should-fix items left open in
  `PR_REVIEW_68_2026-09-25.md`, which is not tracked. The most
  consequential:
  - `accumulating` would also file a stranded phospho-carrier as biomass;
  - the `:driver_policy` drain row calls any handshake interval the
    published one;
  - the full-cycle driver seeds after the build, so a `DRAIN_S≠1` rerun is
    not deterministic.

- [x] 13.1 Add a completeness mode — every dynamic state owned and no dead ends,
  as a *failure* rather than a report — verify by the assembled composition
  passing it and by a composition with one module removed failing it and naming
  the unowned states and stranded moieties. The contract is explicit that a
  successful resolve is not evidence of closure, so **this assertion must be
  shown to fail** or it is not an assertion.
  **Annotated 2026-09-25:** the real assembly had eight dead ends. Four were
  missing owner-mirror edges, which are now added. Four are the PTS carriers,
  reported as `accumulating` biomass by a `:protein`-regime rule (§12,
  2026-09-25 A and B).
- [x] 13.2 Assert every declared edge is executed — verify by a test
  cross-checking the resolved graph against the driver's registered
  contributions, debits, rebuilds and volume reads, so a declared-but-inert edge
  fails the build rather than passing silently.
  **Annotated 2026-09-25:** a clamp is executed when it is held and its value
  is read. An owner's mirror edge is executed by the owner's own term, and a
  missing peer is 13.1's to catch (§12, 2026-09-25 C).
- [x] 13.3 (in 13a, as 13a.3; #66) Fix inference dispatch for a hybrid composition — verify by a test that
  a mixed composition either dispatches on the composition or throws a named
  error, since dispatch currently reads only the first module in the vector and
  would silently take whichever path that module declares.
- [x] 13.4 Settle the drain granularity with D10's measurement at full scale —
  verify by nominal trajectories at 1 s, 5 s and 60 s granularity, by the largest
  relative difference in any observed pool reported, and by the chosen
  granularity registered as a labelled reduction with its *measured* cost.
  **Annotated 2026-09-25:** the chosen granularity is the published 1 s, so
  there is no reduction to register. It is recorded as a `:driver_policy` row
  quoting what the coarser drains cost (§12, 2026-09-25 G).
- [x] 13.5 Run a full cycle at published parameters with a stiff solver and pinned
  tolerances — verify by a trajectory over the full interval with no solver
  failure, and by the solver and tolerances recorded on the model rather than in
  a test file.
- [x] 13.6 Enumerate what is ours across the whole composition — verify by
  `reduction_report` returning the lumped charging step **and its deterministic
  formalism as two separate declarations**, `k_chg`, the tRNA pool size and its
  charged fraction, the
  chemostatted nucleotide and amino-acid pools, the corrected base mapping, the
  ~~thirteen~~ fourteen asserted priors (amended 2026-09-23), the Michaelis-constant column choice, the lactate
  volume ratio, the rounding policy, the drain granularity, the capped rate-law
  geometry and exogenous membrane
  growth; by the count matching the phases that registered them; and by the
  lumping appearing exactly once, since phase 11 no longer registers it.
  **Annotated 2026-09-23 (phase 11a):** the translated-enzyme modes add fifteen
  enzyme slots with `:asserted` priors, beside the five nominal recycling slots
  phase 8 already registered as asserted and not in the fourteen. Decide here
  whether driver-written slots count, and record the decision.
  **Decided 2026-09-25:** they are reported apart as `:discarded_prior`.
  There are 16, not 15: the inbound volume channel also overwrites
  `r_cell_nm`. The fourteen stay the asserted rate priors an inference might
  free. The report gains `:model_note`, `:formalism` and `:driver_policy`
  (§12, 2026-09-25 D).
- [x] 13.7 Measure and record per-trajectory wall-clock at full scale — verify by
  a Slurm run reporting seconds per cycle, comparison against phase 3's
  extrapolation, and an explicit statement of how many trajectories the inference
  budget affords. ~~**K1 and K7 are decided by this number.**~~ **Amended
  2026-09-24:** this number feeds K1, which is scored per posterior in phase 16
  (§8); it does not decide K1 alone, and it never bore on K7, which is about
  posterior width.
  **Annotated 2026-09-24 (phase 11):** an early reading exists. The seven
  modules composed as phase 11 composes them, without transcription's
  counters, take **17.9 s** for 6,299 handshakes after an 80 s first-handshake
  compile (2.84 ms each, job 17131232). That is above K1's roughly 10 s, and
  it is R7's early warning. It is not the phase 13 number. *(The 10 s bound
  was retired the same day, amendment 2026-09-24 A; the reading stands as a
  measurement.)*
- [x] 13.8 Suite and handoff — verify by `sbatch test/run_tests.slurm` passing,
  with the full-cycle run marked and Slurm-gated if it is too slow for the default
  suite. Shortening the interval is forbidden; it is precisely the error of
  record. **Annotated 2026-09-25:** no gate. A warm cycle is 32.2 s, so the
  full cycle runs in the default suite (§12, 2026-09-25 E).
- [x] 13.9 (in 13a, as 13a.1; #66) Give a chemostatted pool both an ownerless debit path and a rebuild
  that reads its held value (§12, 2026-09-10 E, widened by 2026-09-23 A) —
  verify by `build_problem([CoreATranscription(), NucleotideRecycling()])`
  building, by the transcription rebuild reading 0.6874 and 2.7681 mM for CTP
  and UTP, and by the two pinned throws in `test/test_corea_transcription.jl`
  being replaced by these assertions. **Widened 2026-09-23 by phase 12:** decay's
  `CMP_mRNAdeg` and `UMP_mRNAdeg` credit the same two chemostats, so the
  ownerless path must take a credit as well as a debit, and the pinned throw in
  `test/test_corea_transcript_decay.jl` is replaced too. *(13a: transcription
  had one pinned build throw, not two; the `CtpOwner` refusal beside it is
  kept, as the reason the path must be ownerless. Two throws were replaced in
  all.)*
- [x] 13.10 (in 13a, as 13a.2; #66; its last clause, 13.2 seeing every product edge executed, is verified in 13b — done in #68) Give a deferred counter per-product stoichiometry and wire the five
  transcription counters' products (§12, 2026-09-23 B) — verify by `ATP_trsc`
  crediting ADP and phosphate one each per ATP actually paid, by the adenylate
  and phosphate moieties balancing across a handshake on which ATP clips, and by
  task 13.2 seeing every product edge executed. **Widened 2026-09-23 by phase
  12:** decay's `ATP_mRNAdeg` has the same two products and the same gap, and
  its credits are wired here too. **Widened 2026-09-23 by phase 11a:**
  translation's GTP counter is GTP → GDP + Pi upstream, two per residue, the
  same shape again, and its credits are wired here too.
  **Annotated 2026-09-24 (phase 11):** without this, no assembled run means
  anything for guanylate. In the first full seven-module cycle, GTP sits at
  zero from **44 s** and `GTP_translat` clips on 6,257 of 6,300 drains (job
  17131232). `ATP_transloc`'s ADP and Pi join the list here too. *(Those
  numbers, and 11.10's 1.497, predate these credits; 13b remeasures them.)*

### Phase 14a — The balance checks, in the scoping note's order

Was phase 14 until the 2026-09-25 split (§12, same date). Task ids are
unchanged. 14a takes 14.1, check 1's half of 14.2, 14.3 to 14.7, and 14.10 for
those checks. 14b takes the rest.

**Goal:** run the balance checks on the assembled model in the order in which
each catches errors the next would mask.
~~**Done when:** checks 0 through 8 all pass over a full cycle, and each has a
mutation test showing it can fail.~~
**Done when (amended 2026-09-25):** checks 0, 1, 2, 3, 4, 4b and 5 pass over a
full 6,300 s cycle of the assembled model, and each has a mutation test that
fails it naming the intended quantity.
**PR:** #70 (merged 2026-09-26)

**Done (#70):** the suite passes 27,556/27,556 at `3bee403` (Slurm job
17539617), and CI's unit tests pass at `289c491`, which touched only prose and
one test comment. The run of record is job 17539616,
`dev/scripts/corea_validation_result.md`.
- Every whole-cell closure is at 3e-10 to 9e-9 of its `tol_C`, over the full
  cycle and a ladder of seven decades. Nothing falls: the residuals are
  roundoff.
- Carbon and phosphate pass the exception at `n` ulps (3.94 of 13, 4.50 of
  21), and the other seven bitwise or at one ulp.
- Each stoichiometry mutant fails its moiety at 1.6e4 to 1.7e5 × `tol_C`
  and misses the gate by 11 to 18 orders. The boundary mutant fails adenylate
  alone, at 89.6×.
- With the kinase removed, ATP crosses 1% at 246 s. The homolactic ratio is
  2 − 6.0e-13, and 99.86% of the lactate formed is exported.
- /check-PR: MERGE AFTER FIXES (5 major, 9 minor). Every item was fixed or
  stated, then re-reviewed. The fixes also raised CI's timeout from 40 to 90
  minutes, since its first run took 37m53.
- Split tasks: 14.2's check 1 half and 14.10's part for checks 0 to 5 are
  done here, and their boxes stay open for 14b.

- [x] 14.1 Re-assert check 0 on the assembled model — verify by the round-trip
  residual showing the signature its policy predicts across 6,300 handshakes, so
  no downstream residual can be a rounding artefact.
  **Annotated 2026-09-25:** on the assembly the rejected policies' growth laws
  do not separate. So 14.1 asserts the carry within 10⁻⁶ of `tol_C`, and
  deterministic rounding at least 10⁶ times it. The driver records stochastic
  rounding (§3, §12 2026-09-25 B).
- [ ] 14.2 (check 1's half done in 14a, #70; 1b in 14b) Check 1 and 1b — verify by non-negativity naming the first state and
  time of any violation rather than reporting a global failure, and by the
  particle-floor report flagging every state below 500 particles and
  cross-checking the three smallest against a chemical-Langevin ensemble, with
  any observable outside the band excluded from the likelihood.
  **Split 2026-09-25:** check 1 is 14a's; check 1b is 14b's. 1b's
  chemical-Langevin ensemble is a **test-local double**: Euler–Maruyama over
  the ODE block's reactions with the jump path frozen, used by this check
  only. The `:sde` formalism stays refused in `src`, and §7's non-goal
  stands (§12, 2026-09-25). **Amended 2026-09-26:** not Euler–Maruyama but
  its local-linearization form, at h = 0.01 s, since the block is too stiff
  for an explicit step. The cross-check is on the three pools with the smallest
  cycle medians of at least 10 particles; a pool below that is excluded
  outright (§12).
  **Check 1's half done in 14a:** no state falls below its bound over the
  pinned cycle. With translation's debits unclamped, check 1 names `M_gtp_c`
  at 4,914 s (job 17539616). 1b's half stays open for 14b.
- [x] 14.3 Check 2, carbon balance — verify by closure satisfying the tolerance
  principle, by the homolactic ratio equalling 2.000 to integrator tolerance
  since it is analytically exact, and by a mutation removing lactate export
  failing with cytosolic lactate heading for the recorded ~691 mM.
  **Annotated 2026-09-25:** glucose enters through a clamp, so the closure
  needs `MeteredPtsTransport`'s accumulators (§12, 2026-09-25 D). The export
  mutation leaves carbon closed and heads for 643 mM, and the check that fails
  is the exported fraction. It is run by the driver only: a zero-permeability
  module is a distinct type, and compiling it would add about five minutes to
  the suite. The stoichiometry mutation (LDH making two lactate)
  fails carbon alone.
- [x] 14.4 Check 3, redox — verify by the module-local check from phase 6
  restated at composition scope rather than reimplemented, so there is one
  implementation and one bound.
- [x] 14.5 Check 4, adenylate and guanylate **over a full cycle** — verify by both
  moieties closing over 6,300 s, by each being asserted separately so a failure
  names the moiety, and by an explicit assertion that the interval is a full
  cycle, since 144 s of it looked fine before the dead end was found. After D13
  charging's transfer is internal to the ODE block and conserves adenylate, so
  there is no drain to correct for; the kinase-removed configuration asserts a
  threshold crossing rather than exhaustion, per phase 8.
- [x] 14.6 Check 4b, phosphate closure — verify by closure after accounting for
  transcription's pyrophosphate, which the scoping note's own accounting omits,
  and by the flux correction being a subtraction rather than a relaxed bound.
  Charging's pyrophosphate needs no correction after D13: ATP's three phosphates
  become AMP's one plus pyrophosphate's two, exactly closed inside the block. A
  test should assert that rather than leaving a reader to assume it.
  **Annotated 2026-09-24 (phase 13a):** phosphate now also crosses the CTP and
  UTP chemostats. `CTP_mRNA` and `UTP_mRNA` are paid by the chemostat and credit
  PPi to a live pool, and `CMP_mRNAdeg` and `UMP_mRNAdeg` return one phosphate
  each into the chemostat, not into any integrated pool. The closure needs `chemostat_census`, each row weighted by what
  its counter carries — three phosphates for a supplied NTP, one for a returned
  NMP — not by the species.
- [x] 14.7 Check 5, carrier conservation — verify by four independent bounds
  naming the carrier that drifts, restated as conserved up to what translation
  adds, again by subtraction.
- [ ] 14.8 Check 7 and check 8 — verify by the clipping census reporting zero
  carried deficits at published parameters and the clipping fraction across 200
  prior draws (K5), and by the metabolic control analysis identities holding
  within 1e-6, which is the one analytic result the implementation must reproduce.
  **Annotated 2026-09-24 (phase 11), R2's early warning:** the first full
  seven-module cycle clips the charged-tRNA counter on **911** of 6,300 drains,
  first at 199 s, as ATP falls to about 0.44 mM (job 17131232). That run lacks
  13.10's credits and transcription's costs, so it is not the K5 score. In
  phase 9's composition the same counter never clips at 0.25 mM (task 11.7).
  **Annotated 2026-09-25 (phase 13b):** on the assembled model with every
  credit, one seed, `GTP_translat` clips on 253 of 6,300 drains, first at
  953 s, and `tRNA_translat` on 20 (job 17298242). This is still non-zero, so
  K5 is expected to fire at published parameters unless phase 14 finds
  otherwise.
  **Amended 2026-09-25 (14b):** check 7 is **run and scored, not gated**. The
  carried-deficit count at published parameters and the 200-draw fraction go
  on T3 as K5's measured values whatever they show, and phase 14b is done
  when they are recorded, not when they are zero. What to do if K5 fires is a
  separate decision.
- [ ] 14.9 Check 6, the nominal trajectory — verify by the trajectory produced and
  reported as fractional growth or time-to-threshold, with the doubling-time
  comparison refused in code as phase 5 established, not merely in prose.
  **Annotated 2026-09-23 (from phase 9):** also rerun task 9.9's stoichiometry
  comparison on the assembled model, matched on steady flux, and record the
  ATP/ADP ratio under both lumpings. Phase 9 could not answer this, because its
  glycolytic double pins ADP (§9, §12 2026-09-23).
- [ ] 14.10 (checks 0 to 5 done in 14a, #70; 1b, 7 and 8 in 14b) Show each check can fail — verify by one mutation test per check, each
  failing the intended check and naming the intended quantity, collected into the
  mutation table of output F3.
  **Split 2026-09-25:** 14a lands the mutations for checks 0 to 5 and the
  F3 table. 14b adds checks 1b, 7 and 8 and completes F3's table.
  **Annotated 2026-09-25:** the table is `dev/scripts/corea_validation_result.md`
  and not a function. It gives every mutant's residual over `tol_C` for every
  moiety, so a mutation shows which checks it fails as well as the one it
  targets.

### Phase 14b — The derivative, ensemble and census checks

Split from phase 14 on 2026-09-25 (§12, same date). Its tasks are 14.2's check
1b half, 14.8, 14.9 and the rest of 14.10, as written in 14a's list, which
keeps their ids and their annotations. They are ticked there.

**Goal:** finish the checks that each need new machinery: an ensemble for 1b,
a steady-state solve with derivatives for check 8, and a prior-draw census on
Slurm for check 7.
**Done when:** check 1b's floor report and its three smallest pools are
cross-checked against the test-local chemical-Langevin double; check 7 is run
at published parameters and across 200 prior draws and scored on T3 as K5,
pass or fail; check 8's summation identities hold within 1e-6; check 6's
nominal trajectory is reported with task 9.9's comparison rerun; and F3's
mutation table covers every check.
**PR:** _not started_

- [ ] 14b.1 = task 14.2's check 1b half.
- [ ] 14b.2 = task 14.8.
- [ ] 14b.3 = task 14.9.
- [ ] 14b.4 = task 14.10's remainder, for checks 1b, 7 and 8.
- [ ] 14b.5 Output F5, the two external comparisons of §3 (added 2026-09-26)
  — verify by the assembled model's time-averaged transcripts against the
  measured means, over several seeds, with Spearman at least 0.7 and at least
  15 of 17 genes within a factor of two; and by protein fold change over one
  cycle with a median in [1.7, 2.3], no gene below 1.5 or above 3.0, and a
  negative fold-against-length slope. A miss is reported as a divergence, not
  absorbed into the bounds.

**Placed 2026-09-26:** output F5 was unassigned in either half. §9 and tasks
11.10 and 13.7 defer the fold-change comparison to "phase 14's F5", and the
original phase 14 had no F5 task either. It is 14b.5 (§12, same date).

**Decided at planning, 2026-09-26 (§12):** check 7's prior draws vary every
ODE-block kinetic constant with an informed prior. Asserted and uninformed
priors stay at their values, and a slot the driver writes is never drawn.
Task 9.9's rerun matches the two lumpings on cycle-mean charging flux, since
a growing cell has no steady state to match at.

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
- [ ] 15.9 Map observations to states by species name, and keep ensemble spread
  in the ABC summary (§12, 2026-09-23 G) — verify by permuting the observed
  species leaving the log-density unchanged, by a subset of states building and
  scoring, by an unknown species name throwing, and by the time-series summary
  carrying a per-time spread as well as the mean.

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
- [ ] 16.2 **Choose the path-update mechanism and record why.** D10 commits to the
  three-block conditional scheme and leaves this one decision open. Choose among
  particle Gibbs with ancestor sampling, exact forward filtering over a truncated
  state space, and no augmentation at all — verify by the choice justified against
  the granularity measurement of phase 13, against the first-order structure D13
  established, and against what the gradient report says about the assembled
  model's non-smooth boundary edges; and by the two installed constraints of §5
  being addressed explicitly if a particle route is chosen, since one of them
  silently shares a solver object across forked particles rather than erroring.
- [ ] 16.3 Confirm the first block is not smooth just because the path is fixed —
  verify by a test that a clipped drain under a fixed path still has a
  discontinuous derivative in an ODE parameter, so the gradient-based step meets
  K5 inside the conditional scheme exactly as it does outside it, and by check 7's
  census being the thing consulted rather than an assumption.
- [ ] 16.4 Recover the tight-prior control alone — verify by the posterior
  concentrating on truth well inside its prior, by a deliberate perturbation of
  the truth moving the posterior with it, and by shrinkage below 0.9. **If the
  control does not recover, the machinery is wrong and this phase stops** (K4).
- [ ] 16.5 Recover the full six-parameter set — verify by every marginal
  containing truth and by the run's cost reported against phase 13's measured
  per-trajectory budget. **Annotated 2026-09-24:** the per-trajectory budget is
  retired (§8 K1); report the cost per posterior, which is what K1 is judged
  on.
- [ ] 16.6 Compute coverage over repeated datasets — verify by nominal 50, 80, 90,
  95 and 99 percent intervals covering truth at those rates within binomial
  error, per parameter, with the replicate count stated, and by the tight-prior
  control falling inside K4's band.
- [ ] 16.7 Produce the shrinkage table, output F6 — verify by all six targets
  reported under transcripts-only, metabolites-only and joint data, and by the two
  load-bearing cells present: the ptsG promoter under metabolites only and the
  tight-prior control under transcripts only. **These two cells are what "crosses
  the boundary" means operationally.**
- [ ] 16.8 Compare against the reference — verify by output F10 reporting the
  divergence between the production and reference posteriors on the two-gene
  system. There is no cut arm to compare against: D13 shows the existing cut
  passes nothing here, so F11 is instead the composition-and-audit figure §6
  describes, and this task does not produce it.
- [ ] 16.9 Decide K2 — verify by the prior-to-posterior divergence under
  transcript-only data reported for every metabolic parameter against the 0.05-nat
  threshold, and by the verdict written as a sentence with its consequence for how
  the result is framed.
- [ ] 16.10 Label what the result does and does not license — verify by the report
  carrying the scoping note's own caveats: that Core A′ is the best-measured
  region of the network, that enzyme competition is partly removed, that success
  is evidence the architecture works rather than that inference on the full model
  is well-posed, and K7's measured factor if it fired.
- [ ] 16.11 Fill the kill-criteria scoreboard, output T3 — verify by all seven
  criteria carrying a threshold, a measured value and a verdict, **published
  whether or not everything passed.** **Annotated 2026-09-24:** K1's threshold
  is set from the first minimal-model inference, before its verdict is
  written, and T3 records that it was set after the measurement.

### Phase 17 — Per-module calibration

**Goal:** check that each module's posterior is calibrated rather than merely
centred, and that a failure can be attributed to the composition rather than the
sampler.
**Done when:** rank statistics are uniform within test for the ODE block given
the stochastic path and for the stochastic block given the rate-constant
sequence, the joint is checked over the whole target set, and coverage on the
quantities whose information crosses the boundary is reported separately. The
joint check is not about a shared parameter, because D13 establishes there is
none; it is about the composition and the path update.
**PR:** _not started_

**Tasks are written when phase 16 has run.** This is deliberate rather than
lazy. Calibration for a modular system has no established protocol — the survey
states that one will likely have to be built — and its shape depends on the
coupling structure phases 3 to 5 actually produce and on which path update phase
16 settles on. What is committed now is the shape in D12: the three checks, the
attribution rule, the weighted-rank mechanic, and the rank-ECDF test with
simultaneous bands rather than a bare histogram. **The survey's natural
formulation, calibrating each site's outgoing message, is not available here**:
with no shared parameter there is no outgoing parameter message to rank, which
D12 records. Anything more detailed written now would be invented, and an
executor would believe it.

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

### 2026-09-26 — phase 14b planning: check 8's perturbation set, what check 1b can exclude before phase 15, and F5's home

**Trigger:** planning phase 14b against the code. Check 8's summation
theorems hold over a set of parameters that each multiply one reaction's whole
rate, and must cover every reaction. §3's "17 enzymes plus `k_chg`" is not
such a set on this model:
- only thirteen of the seventeen genes multiply a rate: the ten glycolytic
  enzymes, ADK1, GK1 and PPA;
- PGK and PYK each multiply two reactions (the glycolytic step and its GTP
  branch, `nucleotide_recycling.jl` slots `enz_R_PGK3` and `enz_R_PYK3`);
- the four PTS genes are carrier totals, and the cascade's steps are bilinear
  in two carriers, so scaling a total scales no rate proportionally;
- lactate export and the demand forcing that the frozen variant needs have no
  enzyme.

Separately, check 1b excludes an observable "by more than the assumed
observation noise", and that noise is task 15.4's. F5 had no task.

**Change:**

**A — check 8 perturbs one multiplier per reaction.** Approved 2026-09-26.
The set is the 22 ODE reactions plus the frozen variant's tRNA and GTP
demands. Gene-level coefficients are reported by summing a gene's
reactions, and carrier-total responses are reported apart and are not in the
sums. The coefficients come from the implicit function theorem on
automatic-differentiation Jacobians at the solved steady state. That is
exact to the solve and needs no new dependency. A reaction whose steady flux
is zero has no flux control coefficient, so it is named and left out of the
flux sum. It stays in the concentration sums.

**B — check 1b reports the excursion.** Approved 2026-09-26. Per pool, the
largest relative excursion of the ODE trajectory outside the Langevin
ensemble's 90% band. That is the noise level below which phase 15 must
exclude the pool.

**C — F5 is task 14b.5.** Approved 2026-09-26, with §3's thresholds
unchanged. 13b's median fold change of 1.627 predicts a miss. A miss is a
divergence to report, not a band to widen.

**D — check 7's draws and task 9.9's matching.** Approved 2026-09-26. The
draws vary the informed ODE-block kinetic priors only. Task 9.9 matches on
cycle-mean charging flux.

**E — check 1b's ensemble is a local-linearization CLE, at h = 0.01 s.**
Approved 2026-09-26. Check 8's frozen steady state has a fastest mode of
−2.09e4 /s (job 17545158), so explicit Euler–Maruyama needs steps under
5e-5 s, about 1.3e8 per cycle. A drift-implicit step damps the noise on the
fast near-equilibrium modes by about 1/(1 + λh), and those modes set the small
glycolytic pools. Each step is instead the exact solution of the SDE
linearised at its start. The noise is summed over each reaction's forward and
reverse channel separately, and its covariance integral is built by scaling
and squaring. Over one interval this integrator misses the pinned Rodas5P by
up to 24% on the lower-glycolysis pools at h = 0.05 s and 12% at 0.025 s,
first order in h (job 17545342). At 0.01 s that is about 5%. The frozen path
absorbs it at each handshake, and the driver reports the h-halving.

**F — check 1b's cross-check skips pools at about one particle.** Approved
2026-09-26. Over the pinned cycle at seed 1410 (job 17545535), 22 of the 32
ODE states fall below 500 particles and 10 fall below one. Several sit at
about one particle for most of the cycle:

| pool | start | median | time under one particle |
|---|---|---|---|
| PPi | 2,019 | 1.0 | 42% |
| 13DPG | 198 | 0.7 | 73% |
| phospho-EI | 336 | 0.4 | 99% |

NADH, phospho-HPr and phospho-Crr have medians of 2.5, 3.4 and 4.3.

A Langevin band around a pool of one particle compares two invalid
descriptions. So a pool whose cycle median is under 10 particles is excluded
outright, and the cross-check runs on the three with the smallest medians of
at least 10. Ranked by median, not minimum: GTP, GDP and GMP each reach zero
for one handshake when a clipped debit empties them, and hold thousands of
particles otherwise. For phase 15: **13DPG, one of D8's informative small pools, cannot
be an observable on the assembled model.**

**Sections touched:** header; §3 (checks 1b and 8); §11 (14b.5 added, 14b's
decisions recorded, 14.2 annotated); §12.

### 2026-09-25 — phase 14a: whole-cell closures, an n-ulp gate, and what check 0 can assert on the assembly

**Trigger:** implementing phase 14a. The assembled checks were run at
full scale, first in Slurm job 17506883 at `fa3bedc`. After /check-PR's
fixes they were run again in job 17539616 at `3bee403`, which is the run of
record; the suite at that commit passes 27,556/27,556 (job 17539617). Every
residual below is unchanged between the two runs. Two of §3's criteria did not fit
what the assembled model does.

**Change:**

**A — the second gate is `n` ulps.** Approved 2026-09-25. Carbon and phosphate
close to 2.3e-9 and 8.7e-9 of their `tol_C`. Their composed right-hand side
cancels to 3.94 and 4.50 ulps per evaluation. That misses §3's one ulp, but it
is within 0.043 and 0.018 of the summation bound `n·ε·Σ|wᵢ·duᵢ|`. The sums
have 13 and 21 terms from several modules, against phase 8's three. So the
second gate is now `n` ulps (§3). Redox and the four carriers still pass
bitwise, and adenylate and guanylate at one ulp. The ladder over (1e-4, 1e-2)
to (1e-10, 1e-8) keeps every rung at least seven orders below `tol_C`. Its
largest rung is within 79× of its smallest (phosphate), or 30× for carbon.
Nothing falls. The residual is roundoff throughout.

**B — check 0 on the assembly asserts the carry and the gap, not the growth
laws.** Approved 2026-09-25. Fractional carry holds every remainder within
half a particle. Over all nine moieties its whole-cell residuals are 2e-11
(Crr) to 9.5e-6 (carbon) particles, and 4e-11 to 1.2e-6 over the four check 0
tabulates. They stay at roundoff but are not flat. Phosphate's largest
residual grows 9.4× from the first 630 handshakes to the whole cycle, the ulp
per handshake that `N_restarts` puts in `tol_C`. The rejected policies inject
3 to 485 particles. Over the same span, the running maximum grows:
- deterministic: 2.4× to 4.8× by moiety, one seed;
- stochastic: 1.3× to 7.0× across five seeds, median 3.6 against √10 = 3.16.

So stochastic rounding keeps the toy's √N law, and deterministic rounding does
not show the toy's linear one. The pool fractions on the assembled model are
not systematically biased, so its bias does not accumulate. Both policies grow
roughly as √N here, and the growth does not separate them. That toy
remains the signature evidence. On the assembly, 14.1 asserts the carry
within 10⁻⁶ of `tol_C` and deterministic rounding at least 10⁶ times it. The
driver records stochastic rounding. Against the carry's whole-cycle maximum,
its smallest margin is 3.6e7 over the first 630 handshakes (seed 1410,
phosphate) and 1.3e8 over the whole cycle (seed 1413, phosphate). See §3,
check 0's row.

**What `tol_C` does and does not discriminate.** At the pinned pair `tol_C` is
143 particles for phosphate and 4,140 for carbon, because it carries
`N_restarts`, and for carbon the growing meters as well. So an
`assert_conserved` at `tol_C` would admit a leak of whole particles.
Stochastic rounding's phosphate residual is 44 to 107 particles over the first
630 handshakes, under the 143, and 165 to 485 over the cycle. The discriminating asserts are check 0's `10⁻⁶·tol_C` bound,
which every closure meets with at least 100× to spare, and the mutants, which
exceed `tol_C` by 90 to 1.7e5.

**C — closures are whole-cell, so 4b needs no quadrature.** Each moiety
counts four things:
- the ODE pools in particles, with the carried remainder;
- the transcripts, by base content from the extract;
- the CTP and UTP chemostat rows, weighted by what each counter carries
  (three phosphates for a supplied NTP, one for a returned NMP);
- the accrual not yet applied.

Every boundary term is an exact integer. So the prediction in 4b's row, that
an accumulated GTP flux needs a quadrature and so the fall, does not arise.
Because the transcript side reads the extract and not the module, a mutation
that makes GAPD's counters charge ten A for ten U fails adenylate at 89.6×
`tol_C`, and nothing else.

**D — carbon needs a meter, and external lactate is not an exact record.**
Glucose crosses in through a clamped pool, so nothing recorded it.
`MeteredPtsTransport` integrates GLCpts4 and the exporter into two diluted
accumulator states. External lactate is referred to a medium of `R` times the
*live* cell volume and is never diluted. So `lac_e·R·factor` differs from the
metered export by +1.8% at the final factor and −1.9% at the initial one over
a cycle. This changes no rate law materially, since lac_e is under 1% of
cytosolic lactate. It is recorded, not changed. The homolactic ratio is 2 −
6.0e-13. That is the carbon closure written as a ratio, not an independent
result. The independent quantity is that 99.86% of the lactate formed leaves
the cell.

**E — what the configurations and mutations showed.**
- Kinase removed: ATP falls below 1% at **246 s**, against the scoping note's
  144 s for a constant drain. A mass-action drain slows as ATP falls.
  Adenylate, guanylate and phosphate still close.
- With the exporter's permeability zeroed, cytosolic lactate rises at 0.102
  mM/s, heading for 643 mM over a cycle, against the 691 mM recorded.
- Each of the six stoichiometry mutations fails its own moiety at 1.6e4 to
  1.7e5 × `tol_C`. ADK1, GK1 and the carrier leak also break phosphate, which
  is what their chemistry does. Redox, carbon and PPA break nothing else.
- Unclamped translation debits first take GTP negative at **4,914 s** (seed
  1410), and check 1 names it.

**F — cost.** On the node of record, which the suite job was sharing at the
time, the first metered cycle takes 381 s including the build and compilation. A warm cycle takes 12 to 40 s from
(1e-4, 1e-2) to the pinned pair, and 55 s one decade tighter. Two mutant
module types cost about five minutes each to compile, about 10 of the 22
minutes the phase adds to the suite, which now takes 33 to 38 minutes (jobs
17506884, 17539617). `run_tests.slurm`'s limit rises
to 60 minutes. CI took 37m53 against a 40-minute timeout, so the workflow's
timeout rises to 90.

**Sections touched:** §3 (the exception's second gate; check 0's, 4b's and
7's rows); §6 F3; §10 R3; §11 (14a's tasks annotated; 14b notes F5); §12.

### 2026-09-25 — phase 14 splits into 14a and 14b; 1b's ensemble is a test double; check 7 is scored, not gated

**Trigger:** planning phase 14. Three of its ten tasks need machinery that
does not exist: a chemical-Langevin ensemble for check 1b, which §7 lists as a
non-goal and the orchestrator refuses as `:sde`; a steady-state solve with
derivatives for check 8, for which nothing and no dependency exists; and a
200-draw prior census for check 7, about 1.8 CPU-hours on Slurm. Separately,
13b's one assembled seed already clips `GTP_translat` on 253 drains and
`tRNA_translat` on 20 (job 17298242), so "checks 0 through 8 all pass" cannot
be the done-when for check 7, whose zero requirement is K5's threshold.

**Change:**

**A — phase 14 splits.** 14a is the balance checks: 0, 1, 2, 3, 4, 4b and 5,
with their mutation tests. 14b is checks 1b, 7, 8 and 6, and completes F3's
mutation table. Task ids are kept and annotated in place, as in the 13a/13b
split. Approved 2026-09-25.

**B — check 1b's ensemble is a test-local double.** Euler–Maruyama over the
ODE block's reactions, with the jump path frozen, used only by check 1b. The
`:sde` formalism stays refused in `src`, so §7's non-goal is unchanged.
Approved 2026-09-25.

**C — check 7 is scored, not gated.** The deficit count at published
parameters and the fraction over 200 prior draws are K5's measured values and
go on T3 whatever they show. 14b is done when they are recorded. Whether to
smooth the drain or resize a pool if K5 fires is a separate decision, not
taken here. Approved 2026-09-25.

**Sections touched:** header; §11 (preamble count, phase 14 retitled 14a with
its done-when amended, new phase 14b, tasks 14.2, 14.8 and 14.10 annotated);
§12.

### 2026-09-25 — phase 13b: what closure means on the real assembly, what "executed" means per edge, and how the report is typed

**Trigger:** implementing phase 13b. The assembled model had eight dead ends
where the spec expected none. 13.2 and 13.6 each needed a decision the spec
left open, and the first 13.2 rule failed eight tests of partial compositions
(job 17293734).

**Change:**

**A — the four owner-mirror edges that were missing.** Completeness mode on
`corea_models()` found nothing unowned, but it found eight dead ends. Four of
them were owners that never declared their side of a peer's flow:
- glycolysis on 13DPG, which recycling's PGK3 draws;
- recycling on GDP, which translation's GTP counter credits;
- charging on the tRNA pair, which translation's counters draw and return.

Glycolysis and recycling already declare such edges for g6p, pyr, pep, lac
and the nucleotide currencies. The four are added, and no flow changes.
Approved 2026-09-25.

**B — translated protein accumulates, and that is not a dead end.** The other
four dead ends are the PTS carriers. Translation produces them, nothing draws
them down, and the cascade only moves them between phospho-states. Core A′
has no protein degradation (task 11.9), so their only sink is growth
dilution, and that moves no mass between modules. The rule is mechanical: a
registry `:protein`-regime species produced with no consumer is reported in
`CouplingGraph.accumulating`, not in `dead_ends`. A protein consumed with no
producer is still a dead end. So 13b's "no dead ends remain" holds with four
species accumulating. Approved 2026-09-25.

**C — what "executed" means for each edge kind (task 13.2).**
`unexecuted_edges` maps each resolved edge to a record on the driver:
- contributions and peer writes: the driver now keeps these, as
  `ContributionRecord`;
- debits and credits, sign-matched on the accruing side;
- catalytic and geometry exchanges, rebuilds and growth chains.

Two rules were decisions:
- **A clamp is executed when it is held** (approved 2026-09-25). No module
  integrates or writes the species, and its value is read, either by a
  rebuild through the chemostat path or as a fixed parameter sourced from the
  species.
- **An owner's mirror edge is executed by the owner's own term.** The first
  rule asked for a peer's flow in the opposite direction. That refused
  legitimate partial compositions: 8 errors in job 17293734, for example
  recycling's GMP edge with no transcription composed. Whether the peer is
  there is closure, and A's completeness mode already fails it.

**D — the report is typed (task 13.6).** Approved 2026-09-25, with the code
fixed rather than the list amended. Four categories are added:
- `:model_note`, for a module's own simplifications. Every note used to be
  filed as `:lumping`, so "the lumping exactly once" could not be counted.
- `:formalism`, for the charging step's deterministic integration.
- `:discarded_prior`, for an asserted prior on a slot the driver overwrites.
- `:driver_policy`, which records the 1 s drain and fractional carry even at
  their defaults, so both reach every report. The report prints them apart
  and leaves them out of its departure count.

`reduction_notes` entries may now be `category => text` pairs. The driver-written
slots are reported apart from the fourteen: there are **16, not 15**, because
the inbound volume channel also overwrites `r_cell_nm`. The amino-acid
chemostat had no declaration anywhere. `TrnaCharging` now carries one as a
note. A `ClampedEdge` would have been inert, since no reaction reads the pool,
and C's rule refuses an inert clamp.

**E — 13.8 needs no gate.** A warm full cycle is 32.2 s of stepping (job
17298242), so the 6,300 s run sits in the default suite.

**F — what the assembled cycle shows, recorded for phase 14.** One seed, job
17298242:
- `GTP_translat` still clips, on 253 of 6,300 drains, first at 953 s. Phase
  11's partial run clipped it on 6,257, from 44 s.
- `tRNA_translat` clips on 20 drains, down from 911.
- Check 7 requires zero at published parameters (§8 K5). That is scored in
  task 14.8, and this is its early warning.
- The fold-change median is 1.627, and the charged pool never reads under one
  particle at a rebuild instant (§9).

**G — the drain stays at 1 s (task 13.4).** Over 50 paired seeds at aligned
instants (job 17298253):
- At 5 s, 26 of the 29 compared pools have a final point estimate beyond
  1%, though not all of them outside their standard error. The worst is
  **−52.9% ± 9.6%** on 3PG. The phospho-carriers, 2PG and NADH move by
  48–50%.
- 60 s reaches **+288% ± 29%** on phosphate, and several glycolytic pools
  fall about 100%.

D10 expected 60 s to fail and 5 s to pass. Only the first came true. So the
published 1 s drain is used and no reduction is registered. D10's
"no augmentation" route for task 16.2 needs a passing aggregate drain, so it
is closed for this model. A coarse drain cuts a warm trajectory from 29.7 s to
16.6 s (5 s) or 4.6 s (60 s), but that saving cannot be claimed.

**Sections touched:** header; §4 D10 (the path-update table); §9 (the
granularity and fold-change entries); §11 (tasks 13.1, 13.2, 13.4, 13.6 and
13.8 annotated, 14.8 annotated); §12.

### 2026-09-24 — K1's 10 s proxy is retired, and phase 13 splits into 13a and 13b

**Trigger:** planning phase 13. Phase 11's partial seven-module run took 17.9 s
per cycle (job 17131232), above K1's roughly 10 s, before transcription's costs
or any product credit exist. The user's view: a 60 s trajectory is affordable,
because calibration replications are independent and a cluster runs them in
parallel. Separately, phase 13's ten tasks mix three framework fixes with the
assembly and its Slurm measurements, and the measurements mean nothing until
13.10 stops GTP sitting at zero from 44 s.

**Change:**

**A — K1 keeps its criterion and drops its proxy.** Fewer than 50 joint
calibration replications still fires K1, for the same reason: below it the rank
test cannot tell 95% coverage from 85%. The "roughly 10 s per trajectory"
translation is struck. It was never derived from an evaluation count; it appears
only as a flat budget in phase 3, 4 and 5b's notes and in phase 11's readings
(13.7's annotation, and entry E of the phase 11 amendment), which stay as the record of
what they passed. What binds instead is the cost of one posterior: sequential
evaluations times per-trajectory wall-clock. Replications parallelise; D10's
alternating path and parameter updates do not, within a chain. So K1 is scored
in phase 16, once task 16.2 fixes the mechanism and its evaluation count. R7's
early warning becomes that projection. Task 13.7 still measures and records, and
no longer "decides K1 and K7"; the K7 half was wrong, since K7 is about posterior
width. Approved 2026-09-24.

**B — phase 13 splits.** 13a takes tasks 13.9, 13.10 and 13.3, all framework
code, so §10 R15 wants them in their own pull request. 13b is the old phase 13
with the rest. Task ids are kept and annotated in place rather than renumbered.
The parameter-accessor refactor that 2026-09-11 left to "phase 13 or 14" is not
in 13a; it stays unassigned. Approved 2026-09-24.

**C — task 13.9's build needs a third module.** Its clause says
`build_problem([CoreATranscription(), NucleotideRecycling()])` builds. With the
ownerless path in place that call gets past every CTP and UTP check and then
fails for an unrelated reason: recycling reads `M_13dpg_c` and three other
glycolytic species, and in a two-module composition nothing owns them. The test
adds phase 8's `HeldGlycolytic` double, which is what phase 12's tests already
compose recycling with. Nothing about the chemostat path changes. Recorded
2026-09-24, during 13a.

**D — where a chemostat's held value comes from.** 2026-09-23 A says the rebuild
"reads the registry's held value". The registry imports none for CTP or UTP
(`:not_imported`), so the value is the composition's `ClampedEdge` on the
species — transcription's 0.6874 and 2.7681 — with the registry as fallback and
a named error when neither exists. The resolver already holds clamps to the
registry and to each other.

**E — a counter with two consumers may not credit products.** Each consumer pays
the whole accrual and the credit is their sum, so relaxing the one-producer rule
let A + B → C credit C once per reactant (/check-PR of #66: 350 from an accrual
of 250). The driver refuses products on such a counter. No Core A′ counter has
two consumers.

**F — K1 has no numeric bound for now.** Approved 2026-09-24: build and run the
first inference on the minimal model, then set the bound from where it lands,
before phase 16 scores K1. Tasks 16.5 and 16.11 are annotated. This is recorded
so the missing number reads as a decision rather than an omission.

**Sections touched:** header; §5 (seeds and cost); §8 K1; §9 (the fold-change
entry); §10 R7; §11 (preamble count, new phase 13a, phase 13 retitled 13b, tasks
12.4, 12.5, 13.3, 13.7, 13.9, 13.10, 14.6, 16.5 and 16.11 annotated); §12
(this entry; 2026-09-23 C annotated).

### 2026-09-24 — phase 11: the lumped pool's per-amino-acid share, the restart law, and what the first full cycle showed

**Trigger:** implementing phase 11. Two things the spec left open had to be
chosen, two task clauses did not match what the code showed, and diagnostic
job 17131232 produced four results the spec did not predict.

**Change:**

**A — the lumped pool enters the rate law as its per-amino-acid share.**
Task 11.3 said to substitute "the lumped pool" for the twenty per-amino-acid
pools. Read literally, each of the twenty-one concentrations in the published
law would be the whole 0.2 mM pool, 27× upstream's 150 copies (0.0074 mM).
That would cut the charged-pool elasticity from about 0.09 to 0.005 for no
reason but the lumping. Each concentration instead reads [M_trna_chg_c]/20,
0.01 mM at the nominal pool. **Ours**, approved 2026-09-24, and in
`CoreATranslation`'s `reduction_notes`.

**B — the restart law throughout.** Upstream runs `translation_rate_start.py`
(riboKd 1e-4, kcat_mod +0.25) for the first minute and
`translation_rate_restart.py` (riboKd 1e-3, +0.2) at every 60 s rebuild after
it. The module's rebuild is that channel, so it takes the restart law from
t = 0. riboKcat is 12 in both files. The 10 at `MinCell_CMEODE.py:332` is a
variable no rate law reads. Recorded in `reduction_notes`. Approved 2026-09-24.

**C — task 11.1 extracts residues, not per-amino-acid counts.** Under A the
composition drops out of the law, since Σ n_aa = residues. So the extract
gains only `!Residues`, counted by translating each transcript under table 4.

**D — task 11.6's witness has two halves.** It said removing the debit
"reproduces the collapse". It does not. With the credit kept and the debit
removed, translation creates tRNA, 276,128 particles in 600 s. The collapse,
a fully charged pool with no charging flux, needs the whole transfer removed.
Both are pinned. The fix's criterion, a steady split, is a stochastic steady
state: in an ad hoc full-cycle run its 300 s means wandered 0.71 to 0.90 with
transcript noise, and the test asserts a band over 1,200 s.

**E — early warnings from the first full seven-module cycle**, which is not the
assembled model: it has none of 13.10's product credits and none of
transcription's costs.
- `GTP_translat` debits GTP with no GDP credit. GTP sits at zero from 44 s,
  and the counter clips on 6,257 of 6,300 drains. 13.10 is annotated.
- The charged-tRNA counter clips on 911 drains, first at 199 s, as ATP falls
  to about 0.44 mM. This is R2's early warning, and 14.8 is annotated. In
  phase 9's composition, with ATP held by the double, it never clips at
  0.25 mM (task 11.7).
- 17.9 s per cycle after an 80 s compile, 2.84 ms per handshake, against K1's
  roughly 10 s. This is R7's early warning, and 13.7 is annotated.
- **The clip meets the rebuild.** In that run the charged pool reads 0.0000 mM
  at some rebuild instants: alternate ones among the last five printed
  samples, with the count over the cycle not recorded. At those instants the
  60 s rebuild evaluates the translation law at its one-particle floor, and
  every constant sits about 20× low until the next rebuild. ATP alternates in step. That is R11's rebuild-cadence
  concern made concrete: a pool the counter clips is read at the instant it
  is pinned. The protein fold-change median falls to 1.497 (F).

**F — task 11.10's diagnosis: "neither" in the jump-only runs, unmeasured in the hybrid.**
- The elasticity is **0.091**, about 1.9× transcription's, so §9's first-pass
  prediction fails the other way. The reverse channel is stronger, and R1
  does not fire.
- The fold change has two answers, depending on the composition.
  - **Jump-only**, with the constants held at the nominal pool: median
    **1.691** against §3's [1.7, 2.3]. There it is neither of the task's two
    candidates. The pool is nominal by construction, and the ribosome
    constant is the published 12. Reaching 2 needs 1.45× more translation
    output (riboKcat ≈ 17.4). The untested leading candidate is that Core A′
    cuts replication, so gene dosage never rises; a 1.5× mean dosage gives
    2.04.
  - **Full seven-module hybrid**: median **1.497**, minimum 1.232. There
    charging plausibly contributes, through the clip at rebuild instants (E).
- **The split between those two causes is not measured.** Measuring it needs
  the hybrid rerun with the constants frozen at 0.2 mM, or a 6 s rebuild.
- §3's bound is **not relaxed**. The comparison belongs to phase 14's F5. It
  should expect a clear miss on the hybrid unless 13.10 and a fix for the
  clip change the picture.
- *Corrected before merge (/check-PR of #64):* an earlier version of this
  entry said "not the charging step" on the strength of the jump-only runs,
  which cannot test it. It also quoted 1.53× and riboKcat ≈ 18.4, which came
  from an analytic median (1.653) built on the measured rather than the
  simulated transcript means.

**Sections touched:** header; §9 (the channel-strength, fold-change and
pool-size entries); §11 (tasks 11.1 to 11.11 annotated in place and ticked;
9.8, 13.7, 13.10 and 14.8 annotated); §12.

### 2026-09-23 — phase 11 cannot be built as written: catalytic edges, enzyme slots, ptsG, residues and degradation

**Trigger:** confirming phase 11 against the source before starting it. Six
things the phase assumes are false, each checked in code or upstream.

**Change:**

**A — a catalytic edge cannot name a protein count.** The resolver refuses any
edge on a non-registry species (`src/resolver.jl:401`). The registry's only
proteins are the eight PTS carrier forms, so thirteen of the enzyme counts
translation publishes have no species to name. The count is not a registry
species and should not become one: it is jump-block state, like a transcript.
**Phase 11a relaxes the rule for `CatalyticEdge` alone**, since a catalytic edge
carries no mass and so touches no moiety the registry exists to track. Every
other kind still needs a registry species.

**B — the metabolic modules have no slot to fill.** A catalytic edge fills one of
the consuming ODE module's free parameters (`src/handshake.jl:697`).
`CentralGlycolysis` holds its enzyme concentrations in a struct field, and
`NucleotideRecycling` holds them as fixed parameters. **Phase 11a gives both an
opt-in translated-enzyme mode** whose enzyme concentrations are free slots filled
by catalytic edges. The default mode stays as it is, so phases 6 and 8 keep
their standalone results. Only the consumer knows its slot, so the edges are
declared by these two modules and not by translation, which is how the channel
already lowers. PGK and PGK3 share locus 0606, and PYK and PYK3 share 0221, so
thirteen counts fill fifteen slots.

**C — the PTS carriers take protein by credit, and translation declares no
volume edge.** The four carriers are reactants in `PtsTransport`, with no enzyme
factor, so no catalytic edge can reach them. `PtsTransport` already owns
`M_ptsg_c` and `M_ptsg_P_c`, flags them as membrane proteins and declares their
volume edge, and a second owner, flag or outbound volume edge is refused
(`handshake.jl:1142`, `:1271`; `interface.jl:343`). So translation credits each
carrier through a producer deferred counter, the channel transcript decay
already uses to credit AMP and GMP, and the volume chain sees new ptsG through
`PtsTransport`'s flag. Check 5's "conserved up to what translation adds" is
exactly this credit. Task 11a.4 pins the path, and tasks 11.4 and 11.5 are
annotated.

**D — residues exclude the stop codon.** 3,484,518 = Σ copies × (length/3 − 1),
exactly, from `transcription_genes.tsv`. The spot lengths 746, 574, 155 and 90
are length/3, so they include the stop. Phase 9 calibrated `k_chg` on the
former, so 11.1's residue counts and 11.8's closure use it too, and the spot
checks become 745, 573, 154 and 89. Upstream charges GTP on
`len(aasequence)` (`MinCell_CMEODE.py:1008`), which counts the `*`, so it
charges two more GTP per protein. That is recorded as ours.

**E — protein degradation is not a reaction (resolves task 11.9).**
`ptnDegRate` is assigned six times across four upstream files and read by no
reaction; `ptnDegProd = []` is never used. The seventeen degradation jumps are
dropped, and §2's count of 52 stands without the "ours if kept".

**F — two upstream facts the phase did not name.** Translocation charges
`int(len/10)` ATP (`ATP_transloc`, `MinCell_CMEODE.py:1012`), which task 11.4
must declare or record as omitted. And `riboKcat` is 12 in the translation-rate
files that run but 10 in `MinCell_CMEODE.py:332`, while `riboKd` differs between
the start and restart files. Task 11.10 records which is used before quoting a
number.

**Not changed here:** the GTP counter still cannot credit both GDP and Pi,
because of the one-producer rule (`handshake.jl:852`). Task 13.10 owns the
per-product mechanism, and **is widened** to wire translation's counter too, as
phase 12 widened it for decay's. Task 11.5 records it.

Approved 2026-09-23, before phase 11a was started.

**Sections touched:** header Status; §2 (the stochastic-block row); §9 (the
D13 entry's reaction list); §11 (new phase 11a; tasks 11.1, 11.4, 11.5, 11.9
and 11.10 annotated in place; task 13.10 widened); §12.

### 2026-09-23 — phase 12: the guanylate leak is 1.5× the pool, not 60%, and decay's CMP and UMP credit the chemostats

**Trigger:** implementing phase 12. One acceptance number was wrong, and one
task named a species the registry does not have.

**Change:**

**A — the guanylate return is about 61,500 particles per cycle, not ~24,000.**
The phase's done-when and task 12.7 took the scoping note's ~24,000 GMP per
cycle, about 60% of the ~39,800-particle guanylate pool. At the rate constants
phase 10 landed, the steady return is the transcription flux's guanine,
`Σ_g k_g · 6,300 · G_g` = **61,531**, whatever decay's constant is. Eight
simulated cycles give **61,219**. That is **1.54×** the pool. The note's
derivation is not recorded, and no route reproduces it: the measured transcript
means under the published per-gene law give 53,110, and under the dead
`rnaDegRate` of 2026-09-10 B they give 51,055. Every estimate lands at 51,000 to
62,000. The consequence **strengthens** the note's argument for GK1: without it,
one cycle strands more guanylate than the whole pool, which removes the
guanylate pool rather than biasing it. The done-when now asserts the analytic
figure within Monte Carlo error, and 12.7 reads "about 150%".
Because §0 lets the note win on model content, the note is corrected too, in
its own convention: a Was/Actually row of record beside its second-review
table, and an inline correction to the GK1 paragraph. Neither deletes the
original figure.

**B — CMP and UMP credit the CTP and UTP chemostats.** Task 12.4 said to "take
the chemostat exemption on the other two". Core A′ has no CMP or UMP species at
all, chemostatted or otherwise, and the resolver refuses an edge that names a
non-registry species. So decay's `CMP_mRNAdeg` and `UMP_mRNAdeg` are declared as
outbound counters on `M_ctp_c` and `M_utp_c`. The chemostat absorbs them, and
`resolve_coupling` records the exemption the task asks for. This is **ours**: a
monophosphate credited to a triphosphate pool is an exemption, not a reaction,
and `reduction_notes` says so. It inherits 2026-09-10 E's defect from the
producer's side, so no hybrid can build with those two counters until task 13.9.
**Task 13.9 is widened** so that the ownerless path takes credits as well as
debits. The rejected alternative was to accrue CMP and UMP into counters with no
edge. That verifies nothing, and 12.4 would have had to be rewritten rather
than annotated.

**C — decay's energy counter joins task 13.10.** `ATP_mRNAdeg` is ATP → ADP +
Pi upstream (`in_out.py:178`), the same two-product shape as `ATP_trsc` in
2026-09-23 B. It is declared as a debit on ATP alone with its products
recorded, and **task 13.10 is widened** to wire them. Until then the adenylate
moiety loses one ATP per decayed nucleotide, and the phase's executed-credit
test asserts exactly that loss. *(Superseded 2026-09-24: phase 13a wired the
credits, and the test now asserts no loss.)*

**Sections touched:** §11 (Status line; phase 12's done-when and tasks 12.4 and
12.7, annotated in place; tasks 13.9 and 13.10 widened), §12; and
`dev/notes/reduced-syn3a-scoping.md` (a correction-of-record row and the GK1
paragraph).

### 2026-09-23 — task 9.8's check 7 half waits for phase 11, and 9.9 is rematched

**Trigger:** the /check-PR review of #59.

**9.8.** Its verify clause asks for a pool size satisfying checks 1b and 7, or a
statement that none does. Check 7 counts clips of translation's charged-tRNA
debit, and there is no translation until phase 11. So this phase delivers the
scan, check 1b's floor (≈ 0.124 mM at 0.8 charged) and the analytic buffer.
The check 7 half is **deferred to task 11.7**, which already feeds its census
back to this diagnostic. R8's "check 7 as its upper" bound was wrong: a larger
pool buffers the counter, so check 7 bounds from below. R8 is corrected in
place.

**9.9.** The first comparison derived the two-ATP form's constant at nominal
ATP. ATP settles below nominal, and that law is second order in ATP, so the
two ran at fluxes 0.86% apart. The 1.02% ATP/ADP difference it reported was
almost all that mismatch. Rescaling `k2` by 1.044 matches the steady fluxes to
1.5e-9, and the difference is then 0.035% (job 17044321). **But that figure
does not answer the ATP/ADP question either,** as the second review round
found. The glycolytic double pins steady ADP (0.4317 mM in both runs), so the
ratio cannot respond. 9.9 answers the kinase traffic and the pyrophosphate
level. The ATP/ADP half moves to task 14.9, on the assembled model.

**The pool size and charged fraction are fixed at construction.** They set the
initial conditions and `k_chg` once, when `TrnaCharging` is built. Varying
either, in R8's sensitivity or in any scan, means rebuilding the module, not
freeing a parameter. D7's "fourteen an inference phase might free" includes
two that can only be varied that way. The charged-fraction prior is truncated
to [0, 1], so a draw from the asserted priors cannot reach the value where
the derivation throws.

**Sections touched:** §4 D7, §9, §10 R8, §11 tasks 9.8 and 14.9, §12.

### 2026-09-23 — a small moiety fails the first gate's flat ulp bound without leaking

**Trigger:** task 9.6. Measured on the four-module composition (charging,
recycling, the glycolytic double and the translation-demand double) over a full
cycle, the tRNA pair's summed derivative is bitwise `0.0` at **106 of 106** save
points, which is the first gate's criterion. Over nine rungs from (1e-4, 1e-2)
to (1e-12, 1e-10), Slurm job 17044321 (dave29) measures an integrated drift of
29, 613, 1,061, 2,043, 1,077, 383, 805, 391 and 217 ulps of the 0.25 mM sum
(`dev/scripts/trna_charging_diagnostics_result.md`). That is flat, so the
fivefold fall fails, as it should. It also fails the first gate's
every-rung-within-100-ulps bound at eight of nine rungs. In absolute terms the
drift is 1.6e-15 to 1.1e-13 mM, which is implicit-solver roundoff at the scale
of the composition's large states: a ~35 mM phosphate moiety, ~18 mM free
phosphate, ~3.6 mM ATP. It only looks large in ulps of a much smaller sum. §3
already made this argument for the second gate: a flat ulp count is fitted to
one composition.

**The residual depends on the machine, and a first draft of this entry quoted
the wrong run.** It gave a login-node ladder (3,097 down to 60 ulps, and two
tight rungs at 4.2 and 2.58 orders below `tol_C`) as the measurement. The
suite on a compute node logged different values, and the /check-PR review of
#59 caught the mismatch. The numbers above come from the driver under Slurm,
and the suite, run on the same node, logs the same seven asserted rungs. The
conclusion did not change, but the tight-rung reason the draft gave for
stopping the asserted ladder at 1e-10 did not survive: on dave29 all nine
rungs clear three orders.

**Change:** a composition passing the first gate's bitwise test may assert the
second gate's `tol_C` ladder in place of the flat 100-ulp bound. The asserted
ladder ends at the pinned tolerances, and every rung run is reported from a
Slurm run (§3). Task 9.6 asserts six decades, (1e-4, 1e-2) to (1e-10, 1e-8), at
**12.2 down to 4.75 orders** below each rung's `tol_C`, with a **70× spread**
(bound 100×). The margin on the spread is thin, and machine-dependent: the
login node gave 51×. If a future node breaks 100×, that is roundoff moving, not
a leak, and the per-evaluation bitwise check is the guard that tells the two
apart. Reported, not asserted: 4.06 orders at 1e-11 and 3.32 at 1e-12. The
mutation that creates tRNA moves the per-evaluation sum by 0.0274 mM (bitwise 0
unmutated) and the sum by 15.9 mM in 600 s.

*Who else this reaches:* nobody retroactively. Phase 6's redox pair and phase
7's carriers pass the flat bound and keep it.

Approved at implementation time, before 9.6 was written; the two conditions and
the corrected figures were added on review.

**Sections touched:** §3 (the exception's first-gate ladder bound), §11 task
9.6, §12.

### 2026-09-23 — the charging derivation needs a charged fraction, and its flux check is a drift

**Trigger:** writing task 9.2. D14 derives `k_chg` from the residue demand and
the asserted pool size, but the rate law is `k_chg·[M_trna_c]·[M_atp_c]`, and
`[M_trna_c]` is the *uncharged* pool. The total does not fix it. The spec named
neither a value for the pool nor a split, and the registry records no initial
value for either species.

Working that through also showed that task 9.7's first clause cannot be a check.
Against a demand double, the steady charging flux equals what the double
consumes, by mass balance. `k_chg` is calibrated so that happens at the nominal
pools. So "the steady flux reproduces 553.1" holds by construction wherever the
composition settles at nominal ATP, and away from nominal ATP the gap measures
the calibration drifting, not the module being right. D14 already made
conservation the acceptance test. This makes 9.7 say what it measures.

**Change:**

- **A third asserted quantity: the nominal charged fraction.** Defaults: a total
  pool of **0.25 mM** (5,045 particles at 20,180.39 per mM) at **0.8 charged**,
  i.e. 1,009 uncharged and 4,036 charged particles. Both are `:asserted`, and
  **both are recalled rather than cited**: the total as the order of bacterial
  total-tRNA concentrations (a few hundred µM, taken as a concentration and
  applied at Syn3A's volume only to count particles), the split as a typical
  bacterial charged fraction. Neither is from the published model, and **each
  still needs a source** before a result depends on it. At these nominal
  values the charged pool buffers **7.30 s** of the 553.10/s demand (at cycle
  end, where it has settled to 4,000 particles, 7.23 s), the
  uncharged pool clears check 1b's 500-particle floor, and `k_chg =
  0.027408 / (0.05 · 3.6529) = 0.15006 /(mM·s)`. Task 9.8 scans both.
- **Task 9.7's first clause is restated as a drift.** It now reports the composed
  steady flux, and the ATP it settles at, against 553.10/s derived in phase 9's
  own suite from 3,484,518 residues over 6,300 s. Where ATP settles also
  reflects phase 8's glycolytic double, whose rephosphorylation is calibrated to
  the same demand, so the drift is the composition's and not the module's
  alone. The gap is labelled
  self-consistency and not validation. The non-circularity test is unchanged.

Approved at implementation time, before task 9.2 was written.

**Sections touched:** §3 (assumptions table), §4 D14, §11 tasks 9.2 and 9.7,
§12.

### 2026-09-23 — what the 2026-09-22 correctness review found that §12 did not already hold

*Trigger:* a whole-project correctness review
(`docs/reviews/full-project-correctness-review-2026-09-22.md`) raised eight
findings and three smaller ones. Each was checked against the source on
`main` at `fb2d37b`. All are real, but about half were already recorded here
with an owner, and the review cites neither this log nor the phases that own
them. This entry separates the two groups and assigns each new item a task.
The review's links point to `src/framework/` and `src/inference/`, which do not
exist (the files are flat under `src/`), and its "through Phase 10" scope skips
the fact that phase 9 has not started.

**Already recorded, not re-assigned.** CTP and UTP counters that need an owner
the registry forbids (2026-09-10 E, phase 13). Transcription products declared
but never credited (the same entry). Inference dispatch on `models[1]` (task
13.3). A freed initial condition sampled and ignored (2026-09-10, phases 15 to
17). A driver-written slot that must be free and is therefore sampled, which is
`:r_cell_nm` in phase 7 (2026-09-09, phases 15 to 17). A single scalar noise
term and a replicate-mean transcript observer (tasks 15.3 and 15.4).

**A — entry E is wider than it says: the rebuild needs an owner too.** Entry E
names the two deferred counters on CTP and UTP. Transcription also declares an
inbound `RateConstantEdge` on each, and the rebuild lowering throws in exactly
the same way when no ODE module integrates the pool (`src/handshake.jl:1013`).
So E's module-local alternative, declaring no counter on CTP and UTP, would
not make the composition build. The framework has to give a chemostatted pool
two paths: a debit with no owner behind it, and a rebuild that reads the
registry's held value. **Task 13.9.**

**B — one accrual cannot credit two products, and splitting the counter does
not fix it.** `ATP_trsc` turns ATP into ADP and phosphate, one of each per
event. `src/handshake.jl:853-863` refuses a counter with more than one producer,
on the grounds that each would receive the whole accrual and create matter. For
this reaction, crediting both pools in full is the correct stoichiometry: the
check counts molecules as if they were mass. The workaround the error message
suggests, a second counter that only produces phosphate, has no consumer to
match. It is therefore credited `accrued_now` instead of what ATP actually paid,
so it creates phosphate on every handshake where ATP clips. The deferred-counter
contract needs a per-product stoichiometry, with credits scaled from what the
consumer paid. **Task 13.10**, which also wires the five counters' products,
the second half of entry E.

**C — the resolver accepts a named peer that does not supply the quantity.**
`src/resolver.jl:397` checks only that the named peer is present in the
composition. An edge naming a bystander resolves, keeps the wrong peer in its
metadata, and executes against the real owner. The declared topology then
disagrees with what runs, and task 13.2 cross-checks against that declared
topology. **Task 10b.3.**

**D — shared-parameter validation compares prior types, not priors.**
`src/orchestrator.jl:87` compares `typeof(p.prior)`, so `LogNormal(0, 0.1)` and
`LogNormal(5, 2)` pass as the same prior, and module order decides which one is
kept. D13 means Core A′ has no shared parameters today, so nothing currently
depends on this. It is still a silent wrong answer. **Task 10b.2.**

**E — ABC-SMC returns a tolerance it never applied.** `src/abc_smc.jl:94` sets
the next threshold from the accepted distances and returns it with the
population that was accepted under the previous threshold. At `alpha = 0.5`
about half the returned particles exceed the reported tolerance (the review
measured 50 of 100). The population itself is sound; the label is wrong, and
that label is the number a paper would quote. Population 1 is never filtered,
so `n_populations = 1` returns the prior. **Task 10b.1.**

**F — `iterative_infer` converges on nothing.** With no shared parameters,
`chain_kl` sums over an empty set and returns `0.0`, and the loop reports
convergence at iteration 2 (`src/boundary.jl:263`, `:219`). D13 makes that set
empty for Core A′ by design, and D10's conditional scheme, not this loop, is the
production path. But an empty boundary should be an error, not a success.
**Task 10b.4.**

**G — observations are matched to states by row position, and the ABC path
throws away variance.** The NUTS likelihood compares `data_obs[:, i]` with
`sol[:, i]` and never reads `data.species`. Permuting the observed species
silently compares data against the wrong states, and observing a subset fails
on dimension. Separately, the ABC path calls `compute_summary_stats(...; times =
data.times)` (`src/inference.jl:137`), whose time-series branch keeps means only.
Its final-state branch does keep variance and Fano factors, so the review's
"means only" is half right. Tasks 15.3 and 15.4 cover the noise model and the
observer, but neither covers species-to-state mapping. **Task 15.9.**

**H — three phase 10 defects.** `states(m::CoreATranscription)`
(`transcription.jl:464`) takes its counters from `TRANSCRIPTION_COUNTERS`
rather than from the constructor's `counters`, so a custom counter set gets the
default states. `reduction_notes` (`:541`) says the mapping is corrected even
under `base_mapping = :published`. And the 10.7 `@test_skip` still says it is
blocked on phases 6 and 7, both of which merged on 2026-09-10. §11's preamble
gives that check to whichever pull request lands last, and phase 10 did. The
review describes it as a cross-module integration test; it is a copy-number
consistency check. **Tasks 10b.5 and 10b.6.**

*Why a phase 10b:* C to F and H are small, independent of each other and of
phases 11 and 12, and C to F sit in framework files that §10 R15 freezes against
a module branch. Phase 5b is the precedent for landing a framework fix as its
own pull request rather than waiting for a module phase to trip over it. A
and B stay in phase 13, because only the assembled composition exercises them
and task 13.2 is the check that would catch them.

**Sections touched:** §11 (preamble phase count; new phase 10b; tasks 13.9,
13.10 and 15.9), §12.

### 2026-09-11 — the exact-conservation exception gains a second gate, and a closed moiety cannot diverge

Two findings from phase 8, both from the same full-cycle run.

**Amendment A — the exact-conservation exception gains a second gate, for a
moiety conserved by the right-hand side but not bitwise.**

*Trigger:* §3 requires every conservation residual to fall at least fivefold
when `(abstol, reltol)` are tightened tenfold. Measured on phase 8's
composition over a full 6,300 s cycle, across a ladder of nine rungs from
`(1e-4, 1e-2)` to `(1e-12, 1e-10)`: adenylate wanders between 8.3e-14 and
5.5e-13 mM with no trend, guanylate between 4.7e-15 and 1.4e-14, flux-corrected
phosphate between 2.6e-13 and 7.9e-13. Not one falls fivefold, and none should.

*What this is not:* it is **not** the exception phase 6 landed the day before
(§12, PR #50). That exception is gated on the composed right-hand side
returning `sum(nᵢ · duᵢ) === 0.0` bitwise, and phase 8's composition does not:
adenylate is bitwise zero at only 68–80 of 200 sampled states, guanylate at
117–126, flux-corrected phosphate at 0–32 of 106. Phase 6's entry names task
8.7 among those "not pre-granted the exception", and it was right to. **A first
draft of this amendment claimed the fivefold rule binds only where the
trajectory is restarted, and exempted these three moieties by argument. That
draft was wrong twice over** — it reinstated a carve-out phase 7 had already
withdrawn on rebase in favour of phase 6's stricter rule (task 7.7), and its
conclusion did not follow from its premise: if the solver preserves a linear
invariant exactly, it preserves it exactly after a restart too, so the residual
in the assembled model is set by rounding concentrations to integer particle
counts — half a particle is ~2.5e-5 mM at 20,180.39 particles/mM, eight orders
above the round-off floor — and does not scale with tolerance either. It would
have handed task 14.4 a rule that cannot pass.

*Why a second gate is nevertheless right:* phase 6's gate is written for a
composition whose two derivative terms are exact IEEE negatives, which happens
when both come from one reaction in one module. Phase 8's adenylate is exactly
conserved in real arithmetic — `nᵀf ≡ 0`, checked term by term across all three
composed modules — but is summed as `(−adk1 − gk1) + (2·adk1 + gk1) + (−adk1)`
with two further modules adding into the same slots, and that does not cancel
bitwise. It contains no local truncation error either, so it cannot fall when
the solver is tightened. The gate that separates it from a real leak is the
*per-evaluation* residual, and it is tight: measured over 106 trajectory
states, adenylate is at exactly **1.000 ulp** of the 3.9539 mM sum, guanylate
**0.125**, flux-corrected phosphate **0.115**, while a mutated ADK1
stoichiometry puts adenylate at 3.70 mM — 1.7e16 ulps.

*Change:* §3's exception now has two gates. The first is phase 6's, unchanged.
The second requires `|Σᵢ nᵢ · duᵢ| ≤ 1 ulp of the conserved sum at every
sampled state`, then the same ≥6-decade ladder, with the per-rung bound the
accumulated `n_evaluations · eps(sum)` rather than a flat 100 ulps. Both gates
are assertions about a composition, neither may be claimed by argument, and a
phase that passes neither keeps the fall. The mutation test is unchanged and is
still what shows the check can fail.

*Who else this reaches:* **nobody retroactively.** Phase 6's redox pair and
phase 7's four carrier sums already assert the first gate and are untouched —
task 7.7 in particular already carries a six-rung ten-decade ladder within 100
ulps, which is strictly stronger, and inherits nothing from this entry. Phase
14's checks 2 to 5 keep what phase 6 assigned them: restated in particles with
the `N_restarts` factor, which is the assembled-model statement and is not this
gate. Phase 9's task 9.6 must still assert a gate or keep the fall.

**Two things recorded rather than fixed here, both framework.**

*A `fixed` parameter has no channel into the composed vector, and three modules
now work around it separately.* `model_free_params` (`src/parameters.jl`)
filters `fixed` out and `src/orchestrator.jl` builds `p` from that alone, so a
module cannot read its own held constants from the parameter vector at all.
Phase 6 carries `values` + `free_slot`, phase 7 `q` + `free_slots`, phase 8
`held` + `pidx` — one workaround in three spellings, and freeing a parameter
shifts every later index. The root fix is in the framework: resolve all of a
module's parameters and hand `dynamics` an accessor that reads free ones from
`p` and fixed ones from a composition-owned constant vector, which deletes the
extra field, the slot map and the `_k`-style helper from every module. §10 R15
freezes framework files against a module branch, so it is not phase 8's to do —
but it is now duplicated three times rather than once, which is the threshold
at which it stops being a convention and starts being a defect. **Phase 13 or
14 owns it.**

*The anti-circularity guard enforces a weaker proposition than the one that
matters.* Task 8.8 requires that phase 9 not calibrate `k_chg` against the
demand this module already assumed, and phase 8 implements that by grepping its
own source for the number. That keeps phase 8's file clean, which is worth
something, but phase 9 can still read `RECYCLING_DRAIN_MM_PER_S` straight out
of `test/nucleotide_test_models.jl`, which no guard covers. **The obligation
belongs in phase 9's own suite**: assert that phase 9's calibration target is
derived there from 3,484,518 residues over 6,300 s, not imported from a phase 8
constant. Recorded against task 9.6.

**Amendment B — removing the pyrophosphatase strands the phosphate moiety
rather than letting pyrophosphate diverge.**

*Trigger:* phase 8's done-when and task 8.6 both asked for pyrophosphate
"rising without bound" with the enzyme removed. In a model whose phosphate is a
closed moiety — which check 4b asserts Core A′'s is — that cannot happen.
Measured: pyrophosphate rises 128-fold, from 0.1 to 12.7902 mM, and ends holding
72.9% of the 35.0920 mM phosphate budget; free phosphate is then gone, the
substrate-level phosphorylation has nothing to work with, ATP falls below 1% of
its initial value, charging stops for want of ATP, and the rise flattens. The
scoping note's 173 mM is one pyrophosphate per charging event sustained for a
whole cycle — open-pool arithmetic, assuming a phosphate supply the closed
model does not have. An early version of phase 8's glycolytic double did
reproduce 172.7 mM, and only by leaking 345 mM of phosphate into the
composition, which also drove pyrophosphate to 51 mM in the configuration that
*keeps* the enzyme, against the 0.371 mM it settles at once the double takes
its phosphate from the free pool as glycolysis does.

*Change:* the phase 8 done-when, §3 check 4 and task 8.6 now say the enzyme's
removal strands the moiety and stalls the pathway. Both numbers are recorded
with the model each belongs to. This is the stronger demonstration that the
reaction is required, and the one that survives assembly.

**Sections touched:** §3 (the tolerance principle's exception, check 4 and
check 4b), §6 F3 (the slope-zero case now covers the second gate too), §10 R3
(the adenylate tripwire is restated as the per-evaluation gate, since a flat
residual is no longer by itself a warning), §11 phase 8's done-when and tasks
8.6 and 8.7, §12. Tasks 6.7, 7.7 and 14.4 are deliberately **not** touched.
Approved at implementation time, on the ladder and the per-evaluation
measurement above, before the checks were rewritten. Implemented in PR #52.

### 2026-09-10 — what phase 10 needed that the spec did not name: a decay double, two more upstream files, and a predicted band computed from a constant the model never uses

**Trigger:** implementing phase 10. Three of its clauses could not be executed
as written, and one of design D9's numbers turned out to rest on dead code.

**Change:**

**A — task 10.8 gains a transcript-decay double, and it is an acceptance
criterion.** Phase 10 owns transcription; phase 12 owns decay. With no consumer
the transcript counts only accumulate, so each gene's time-averaged count over a
6,300 s cycle lands about an order of magnitude above its measured mean and the
calibration check is vacuous rather than weak. The double runs the published
per-gene law, `(18/452)·88 / n_g` (`MinCell_CMEODE.py:444-460`), which is
derivable from the extract's own length column and introduces no upstream data
the phase did not already vendor. This is the shape phases 8 and 9 already use —
8.6's charging drain double and 9.5's translation-demand double are both
acceptance criteria and say so. Phase 12 supersedes it for composed runs and
keeps it for standalone ones. It also exercises the 2026-09-04 amendment from
the consumer's side: a transcript is not a registry species, so the write is
gated by `written_states` alone.

**B — the Done-when's "each gene" becomes §3's "at least fifteen of
seventeen", because D9's predicted band came from a constant that never
runs.** The archived `add-corea-transcription` design's D9 — the archive's, not §4's — states
predicted steady states of 0.44 to 2.87 copies.
Those are `k_g / rnaDegRate` at `rnaDegRate = 0.00578/2`, which
`MinCell_CMEODE.py` defines at line 295 **and never uses** — it appears on its
own definition line and nowhere else, in that file or in `MinCell_restart.py`.
The reaction the model actually adds is `DegradationRate(rnaMetID,
rnasequence)` (line 444, added at line 729), which is per-gene: one catalytic
constant over transcript length. At the law that runs, predicted counts span
**0.332 to 2.406** against measured **0.292 to 2.178**, and **sixteen of
seventeen** genes agree within a factor of two.

The outlier is **FBA, `JCVISYN3A_0131`** — 3.96× on that deterministic band,
and 3.67× on the simulated eight-replicate time averages the suite reports.
The two are different estimators of the same disagreement, not one number. It is not an extract error:
its 775 copies match the scoping note's own table and its measured 0.3469
matches `mRNA_counts.csv` line 110. FBA is abundant protein with a rare
transcript, which is a property of the two measurements and exactly the kind of
disagreement an out-of-sample check exists to expose. §3's external-validation
row already states the right bound — Spearman ≥ 0.7 and a factor of two for at
least fifteen of seventeen — so the phase Done-when is brought into line with it
rather than the other way round. Measured: Spearman **0.8701** over eight
replicates, **16/17** within a factor of two.

This is the same class of trap as D2's cross-file ambiguity and D3's
first-minute setup constants: a number that looks corroborated because it is
written down, and is not what executes.

**C — the extract takes five upstream files and carries a column D8 omits.**
The archived `add-corea-transcription` design's D8 names three sources. Transcript length and the four base counts do
come from `syn3A.gb` alone, but the protein copy number does not: that record
carries **no AOE protein ids at all** (`grep -c AOE` returns 0, against
`syn2.gb`'s 454). The copy number runs
`JCVISYN3A_xxxx` → `MMSYN1_xxxx` → `FBA/Syn3A_annotation_compilation.xlsx` →
`JCVSYN2_xxxxx` → `syn2.gb` `/protein_id` → `proteomics.xlsx`, which is the
chain `MinCell_CMEODE.py:57-110` and `:146-170` use. The extract also gains a
`!First2` column: the rate law reads the NTP concentrations of the transcript's
first two bases as `C₁` and `C₂` (`:370-372`), so those bases are per-gene data.

**D — a `dev/scripts/Project.toml`.** Reading two `.xlsx` files needs
`XLSX.jl`. It lives in a dev-scripts environment rather than the package's own
`Project.toml`, so Aqua's stale-dependency check stays green and the generator
stays Julia. Build tooling, not framework code, so R15 does not bind.

**E — the deferred counters on CTP and UTP debit pools the registry forbids
anything from owning, so they can never be satisfied. Phase 13 owns it.**

*What happens:* `build_problem([CoreATranscription(), NucleotideRecycling()])`
throws `ArgumentError: Module CoreATranscription declares a
DeferredCounterEdge debiting :M_ctp_c, but no ODE module in this composition
integrates that pool`. The hook's debit path requires an ODE module that
integrates the pool it debits.

*Why no composition fixes it:* the obvious repair — compose the module that
owns CTP — cannot exist. `registry.jl` records `M_ctp_c` and `M_utp_c` as
`:chemostat`, and the resolver refuses any module that integrates a
chemostatted species: *"Module CtpOwner declares dynamics for :M_ctp_c, which
the Core A′ registry holds at a fixed concentration"* (measured, with a stub
owner). So the pool is un-ownable by construction, and a debit that demands an
owner is unsatisfiable by construction. **The two clamps are right** — a clamp
is the correct declaration for a chemostatted pool, and the module is careful
not to clamp ATP or GTP, which recycling really does own. It is the two
*counters* on CTP and UTP that ask for something the registry rules out.

*Why it went unseen:* task 10.6 verifies through `resolve_coupling([m])`, which
is the standalone case and passes, and task 10.8 composes two *jump* modules,
so no hybrid problem is built for this module anywhere in the suite. The phase
is green on its own acceptance criteria and still cannot be assembled.

*Why it is not fixed here:* the repair belongs in the hook's debit check, which
must let a chemostatted pool absorb a debit with no owner — a chemostat is
exactly the thing that can — and §10 R15 freezes framework files against a
module branch. The module-local alternative is to declare no counter on CTP and
UTP, but task 10.6 asks for five counters by name and the cost accounting wants
all four monomers, so that is a change to what the phase delivers rather than a
fix. **Phase 13 owns it**: it composes all seven modules and asserts every
declared edge is executed, so it is where this has to be resolved either way.
Recorded in the shape §12's 2026-09-10 and 2026-09-11 entries use for a
discovered framework gap — a freed initial condition that is sampled and
ignored, a `fixed` parameter with no channel into the composed vector — both
left loud rather than fixed on a module branch.

*Pinned, not merely described:* `test/test_corea_transcription.jl` asserts both
throws and their messages, so the day the framework learns about chemostatted
debits, the test fails and points here.

**A second, smaller thing this surfaced.** The five counters name their products
(`M_adp_c`, `M_pi_c`, `M_ppi_c`) but every edge is declared `direction = :in`.
The driver credits a pool only through an outbound `DeferredCounterEdge`, so a
driven run would debit ATP and return no ADP and no phosphate. `counter_drains`
reports the products, but reporting is not wiring. Same owner, same reason.

**Also recorded, not an amendment:** the archived `add-corea-transcription` design's D4 — again the archive's, since
§4's D4 is the base mapping — `recompute_rate_constants!(m;
atp, ctp, gtp, utp)` with mutable constants on the struct is superseded by
phase 4's task 4.1 — the constants live in the composed parameter vector and
the protocol is `rebuilt_params(m)` plus `rate_constants(p, t, m, pools)`. §11
phase 10's preamble already records the consequence for task 10.5.

**Sections touched:** §3 (external comparison 1, whose predicted band amendment B
falsifies), §5 (the per-gene row's source count, which amendment C raises from
three to five), §11 (the document Status line, phase 10's Done-when, and tasks
10.1 to 10.8 — 10.1 and 10.8 in substance, the rest ticked with their evidence),
§12.

### 2026-09-10 — an exactly conserved moiety cannot satisfy the tolerance principle, and the exception is fenced by a bitwise criterion

*Trigger:* phase 6's redox check was written as the tolerance principle
prescribes — solve at `(1e-10, 1e-8)`, solve again at `(1e-11, 1e-9)`, require
the residual to fall at least fivefold — and it failed on the first Slurm run
(job 16363584): 1.42e-14 mM against 7.99e-15 mM, a ratio of 1.78. A ladder over
eight decades of tolerance, from `(1e-4, 1e-2)` to `(1e-12, 1e-10)`, then showed
the residual **bounded at the floating-point floor and not falling**, wandering
between 8.9e-16 and 2.7e-14 mM — a 30× spread with no downward trend, and if
anything a mild *rise* over the last four decades (3 → 32 → 60 ulps) as the
solver takes more steps. `eps` at NAD⁺ + NADH ≈ 2.2097 mM is 4.4e-16, so every
point on that ladder is 2 to 60 ulps. The residual is linear-algebra roundoff in
the implicit stage solves, accumulated over the step count; it contains no local
truncation error, which is the quantity a tolerance is a knob on.

*Why, and why it was foreseeable:* GAPD and LDH_L are the only reactions that
touch the pair, both are inside this one module, and their stoichiometry
mirrors, so the composed right-hand side returns `du[NAD⁺]` and `du[NADH]` as
bit-for-bit negatives. Their sum is then conserved whatever step the solver
takes. §3 check 3 already said the invariant is "exactly invariant there"; the
done-when clause and task 6.7 asked for a residual that shrinks. The two could
not both hold, and the check that was written to distinguish a structural leak
from a numerical residual has nothing to measure when the residual is neither.

*Change:* the tolerance principle gains a stated exception, fenced three ways so
it cannot become a loophole. **(i) The criterion is mechanical and must be
asserted, not argued:** a check may use the exception only where a test
evaluates the composed right-hand side at sampled states and asserts
`sum(nᵢ · duᵢ) === 0.0` bitwise. Crossing a module boundary is *not* the
criterion — it is the usual reason the criterion fails, and it is neither
necessary nor sufficient, so a phase whose moiety happens to sit inside one
module standalone gains nothing by saying so. **(ii) The assertion carries
numbers:** a ladder over at least six decades, every rung within 100 ulps of the
conserved sum and the largest within 100× of the smallest — as checkable as the
fivefold fall it replaces. **(iii) It is scoped to one continuous solve**,
because a handshake is not the right-hand side: at composition scope growth
dilution rewrites every concentration ~6,300 times, so an assembled restatement
is in particles and carries the `N_restarts` factor. Task 14.4 inherits that
sentence, not a bare module-local bound.

Phase 6's done-when and task 6.7 are annotated in place to assert the ladder
rather than the fall; check 3's row records the measurement and the
composition-scope limit; F3's "slopes required positive" gains the
slope-zero-or-slightly-negative case, labelled rather than counted as a failure.

**The mutation test is unchanged and is what keeps this from being a weakened
check:** a mutated GAPD stoichiometry moves the residual to 2.21 mM — the whole
pool — which is **14.2 orders of magnitude**, and does not shrink with tolerance
either. The test enforces a margin of six orders, well inside what was measured.

Nothing here relaxes a check on an invariant that fails the bitwise criterion,
which is every moiety with a real flux across it: adenylate, guanylate,
phosphate and carbon all keep the fivefold fall. Tasks 7.7, 8.7 and 9.6 are
**not** pre-granted the exception by this entry; each must assert the criterion
in its own composition or keep the fall.

**Sections touched:** §3 (the tolerance principle, check 3), §6 (F3), §11 (phase
6 done-when, task 6.7). Approved at implementation time, on the measurement
above, before the check was rewritten. Landed in PR #50.

### 2026-09-10 — a freed initial condition is sampled and then ignored, so phase 6 refuses to free one

*Trigger:* phase 6 is the first module to expose initial concentrations as
inferable parameters — 13 of its 65 — and its own docstring named `:M_g6p_c0` as
a valid `free` argument. Pre-merge review measured what that does: the freed
parameter appears in the sampled vector, and **nothing reads it.**
`_collect_ic_values` (`src/orchestrator.jl`) builds `u0` from each parameter's
stored `value`, never from the sampled vector, and no rate law indexes a `conc_`
slot. Changing a freed initial condition from 3.7076 to 99.0 left the
right-hand side and a 100 s trajectory bit-identical.

*Why it matters more than it looks:* the failure is silent and it looks like a
result. The likelihood is exactly flat in the parameter, so the posterior
marginal equals the prior, and that reads as "the data does not constrain the
initial pool" — a finding — rather than as a channel that was never wired.
§4 D11's target set is six parameters and the spec does not forbid an initial
concentration among them.

*Change:* **phase 6 refuses `free` on a parameter whose role is
`:initial_condition`**, naming the reason. That converts a wrong posterior into
a build-time error. Wiring `u0` to the sampled vector is framework and is not
done on a module branch (§10 R15), so this is recorded as a known gap rather
than fixed here: **the inference phases 15 to 17 own it**, alongside the
`rebuilt_params` and driver-written-slot problem §12's 2026-09-09 entry assigns
to the same phases — it is the same defect in a different channel, a parameter
that is sampled and then overwritten or ignored. Phases 7, 8, 9 and 10 vendor
initial conditions the same way and inherit the refusal.

**Sections touched:** none — the spec said nothing that is now false, so this
entry records a discovered framework gap and the module-level refusal that keeps
it loud. Landed in PR #50.

### 2026-09-10 — what phase 6's vendored extract falsified: the prior-width column, D11's two controls, and a rounded conversion factor

*Trigger:* vendoring the central balanced table put 65 real numbers in the repo
for the first time, and four statements the spec made about them turned out not
to survive contact. None was caught by a failing test; all four were caught by
checking the prose against the file.

*Change, one per statement.*

**D1's "with its geometric standard deviation" names a column that does not
exist.** The balanced `Parameter` table carries exactly one width,
`UnconstrainedGeometricStd`, so every prior in the project pairs a *balanced*
median with an *unconstrained* width. Measured across phase 6's extract, that is
free for all 13 concentrations and all 32 Michaelis constants — their `Mode`
equals `UnconstrainedGeometricMean` exactly — and not free for **ten of the
twenty catalytic constants**, worst `R_TPI` reverse at mode 4 against an
unconstrained geometric mean of 65,341.7 with a 1.05 width. D1 now says so, and
labels it. **This is the fourth instance of the column trap D2 and D3 record**,
and the first that is not about *which* column but about a column that was
assumed to be there. Phases 7, 8 and 9 vendor from the same tables and inherit
it.

**D11's positive control is not the tightest prior in the core.** `kcat_ENO`'s
gstd is 1.170; eight forward constants are tighter, `kcatF_R_LDH_L` at 1.0512
tightest. And `kcat_FBA` at 1.747 is the loosest *forward* constant only —
`kcatR_R_PFK` (33.03), `kcatR_R_PYK` (11.34) and `kcatR_R_FBA` (8.36) are
looser, and the first two exceed the prior-default width and so enumerate as
uninformed. Both rows are annotated in place. **The target set does not change**:
what the controls need is a tight prior on the capacity-tightest enzyme against a
loose one, and 1.170 against 1.747 delivers that.

**Task 6.4's "copies over 20,180 to six decimals" cannot be satisfied.** The
conversion factor is derived rather than transcribed —
`corea_particles_per_mM()` is 20,180.39 — and the rounded figure disagrees in
the sixth decimal for PFK, GAPD, PGK, ENO and LDH_L. A transcribed constant
would also put this module and the handshake's own conversion into disagreement
the moment translation supersedes these concentrations, so the derived factor
governs and the clause is annotated.

**§9's question about full-cycle checks in the default suite is half answered.**
Phase 6's six full-cycle solves cost about 1m40 to 2m10 of a 4m56 suite, which
is affordable for a 13-state ODE block, so they stay. The assembled hybrid is a
different question and task 13.7 decides it. The question is struck through with
that resolution rather than deleted.

*One thing deliberately not changed.* The archived reference for this phase names
`conc_M_atp_c` among four spot values while fixing the concentration count at
thirteen; ATP is not an owned state, so the two cannot both hold. The count is
what task 6.1 carries, `conc_M_g6p_c` is asserted instead, and
`src/organisms/coreA/data/README.md` records it. That is a correction of record
against an archived document, not an amendment to this spec.

**Sections touched:** §4 (D1, D3, D11), §9, §11 (task 6.4). Landed in PR #50.

### 2026-09-09 — a known framework gap is pulled ahead of the fan-out as phase 5b

*Evidence:* phases 6, 7, 8 and 10 are mutually independent and were about to be
built as four parallel pull requests. R15 makes framework code off-limits on a
module branch, and phase 7 breaks that rule by construction: task 7.2 needs the
radius to reach the lactate exporter's rate law, and phase 5 refuses an inbound
`VolumeEdge` by name in `src/handshake.jl` because the type carries no parameter
slot. The phase 7 agent would meet a live `ArgumentError` in a file it may not
edit. R15's escalation path does not fit, because it is written for a need that
*surfaces* during a phase; this one was already diagnosed and already in the
spec. What was missing was a place in the ordering.

*Change:* a new phase 5b between phases 5 and 6, delivering the inbound half of
the volume channel. Lettered rather than numbered because phases 6 to 17 are
cross-referenced throughout §4, §6, §8, §9, §10, the archived OpenSpec designs
and both progress figures, and renumbering them buys nothing. Three decisions
inside it are worth recording here because they are not derivable from the goal:
the edge carries a **unit-bearing** `quantity` from a fixed vocabulary, since
nothing else in this codebase records a unit anywhere and a `μm`-for-`nm` slip
would be absorbed into the permeability it multiplies rather than caught; the
geometry a rate law receives is derived from the **capped** volume, so the three
quantities describe one sphere, while the reported radius stays the published
uncapped value; the write happens at step 0 of every handshake but *outside*
`_update_volume!`, which returns early for a fixed cell and would therefore skip
exactly the composition phase 7 needs standalone; and an inbound edge's
`species` names **the declaring module's own state whose rate law reads the
geometry**, not the membrane protein. That last one governs how task 7.2 must be
written — the edge goes on `M_lac__L_e`, not on ptsG — and is enforced by a
refusal in `src/resolver.jl`. The geometry is a sum over every flagged state and
so belongs to no module, but every edge must name a registry species, and naming
the protein would be a claim that silently changes meaning once task 7.6 flags
both ptsG phospho-forms.

*What this deliberately does not fix:* a driver-written slot must be a free
parameter and is therefore also sampled. Pre-existing, documented for
`rebuilt_params` and undocumented for the catalytic channel. 5b makes the slots
enumerable through `driver_written_params` and documents the trap on
`VolumeEdge`; the fix stays where `rebuilt_params` already put it, in the
inference phases 15 to 17. It is **not** task 13.3, which is hybrid dispatch:
13.3's test cannot catch this, because the composition that silently samples a
geometry slot is the homogeneous ODE block, not a mixed one.

*Also recorded:* the volume channel now runs in both directions, so §2's and
§3's coupling enumerations say so and §6 F1's volume row carries both gains; and
`driver_declarations` gains a `:capped_rate_law_geometry` label, so §6 T2 and
task 13.6 list it — 13.6 verifies "the count matching the phases that registered
them" and would otherwise fail at phase 13.

*Sections:* §0 (the phase count), §2's model table and G2, §3's coupling bullet,
§4 D6, §6 F1 and T2, §11's opening sentence and parallelism paragraph, §11
phase 5b (new), §11 tasks 5.1, 7.2, 8.4, 8.5, 10.7 and 13.6.
Written with the phase rather than after it, since the phase exists because of
this amendment; landed in PR #48. **Nine of the sections above were added by the
pre-merge review of that PR**, which found the amendment's first draft
assigning the sampled-slot problem to the wrong task, its *Sections:* list
incomplete, and §6 T2 and task 13.6 not carrying the label this phase
registers — the last of which would have failed task 13.6's own count check at
phase 13.


### 2026-09-05 — three things phase 4 contradicted: the continuous cadence, the granularity comparison, and what the clipping census counts

Logged together because one phase found all three, and because the pre-merge
review of PR #45 found them rather than the phase itself — which is the more
useful fact about them.

**Amendment A — a continuous rebuild cadence is refused rather than measured.**
*Evidence:* §3's assumptions table and §9's bullet both promised that the
continuous variant is "declarable, defaults to published, and the deviation is
**measured** rather than argued". Task 4.4 offered either implementing the
cadence or rejecting it by name, and the phase rejected it: the mechanism of
task 3.1 — an outer loop over two stepped integrators — holds a rate constant
between refreshes *by construction*, which is precisely what makes task 4.3's
piecewise-constancy assertion possible, and propensities tracking a pool
continuously would need the `JumpProblem`-over-`ODEProblem` shape task 3.1
rejected. Running it at 1 s and labelling that "continuous" would be a label
that overstates what runs. So the choice is right and the promise is now false:
the deviation can be declared and labelled, but not measured. *Change:* the
deviation that §3, §9 and R11 measure is a shorter **piecewise-constant**
interval; a continuous cadence is labelled and unmeasured, and is refused at
build time with an error naming the module, the species and the mechanism.
*Sections:* §3 assumptions table, §9, annotated in place.

**Amendment B — D10's granularity comparison is a paired ensemble at
drain-aligned instants, not a nominal trajectory.** *Evidence:* two findings
from phase 4's miniature run, and the second was got wrong first. (i)
Coarsening the drain changes how often the hook clears the counters, and the
driver rebuilds the SSA's propensity aggregation whenever it does — which
consumes randomness — so three configurations at one seed give three different
stochastic paths and a single-path difference is Monte Carlo noise wearing the
label of a granularity effect. The same applies to R11's rebuild cadence, since
a rebuild triggers the same aggregation rebuild. (ii) Between drains a coarse
configuration holds up to `drain − interval` seconds of accrued cost in its
counters, so its pools sit high by exactly that unpaid balance. The first
version of phase 4's measurement sampled at every handshake, and its headline
"0.65% granularity cost" was that sawtooth: recomputed in closed form,
59 s × 40 particles/s ÷ 20,180 ÷ 18.13 mM = 0.00645 against the 0.0065
reported — three significant figures, a deterministic bookkeeping offset that
needs no simulation at all. *Change:* D10's blockquote and R11's early warning
now read "nominal trajectory" as a paired ensemble sampled at instants that are
multiples of every drain compared, with the maximum reported alongside the
multiplicity it was selected from and the sawtooth reported separately.
Phase 13 inherits the method with the number. *Sections:* §4 D10, §10 R11,
annotated in place.

**Amendment C — check 7 and K5 are scored at the published 1 s drain.**
*Evidence:* phase 4's `drain_interval` did not exist when check 7 and K5 were
written, and it moves both of their inputs in opposite directions. Only a drain
can clip, so a per-handshake denominator dilutes the clipping fraction by
exactly `steps_per_drain` — at 60 s, a run in which *every* debit clipped would
report 1.7% and pass K5's 5% gate — while each debit is `steps_per_drain` times
larger against the same pool and so clips more readily. Applied to the counter
§3 names as the one to watch, ~553 residues/s against a pool of order 10³
becomes ~33,000 residues in one aggregate debit. *Change:* check 7 and K5 are
scored at the published 1 s drain; `clipping_census`'s fraction is per drain
rather than per handshake, and it reports the accrual still pending in the
counters, because between drains that debt sits in no `deficit` and the ledger
would otherwise read closed. *Sections:* §3 check 7, §8 K5, annotated in place.

**What is *not* amended, recorded because the line matters.** Phase 4 also
added two protocol functions (`rebuilt_params`, `rate_constants`), a
`drain_interval` keyword, `driver_declarations` and two `REDUCTION_CATEGORIES`
members. None of those contradicts anything the spec said: task 4.1 delegated
the mechanism, task 4.6 mandated the comparison the knob serves, and §6 T2
already listed both new label rows. They extend the spec, and extensions are
recorded in the pull request, as tasks 2.2 and 3.1 recorded theirs. Approved
before the fixes landed; landed in PR #45.

### 2026-09-05 — check 7's census on a toy does not decide K5

*Evidence:* phase 3 ran check 7's census on its two-module toy and measured zero
carried deficits at published parameters over 600 handshakes, and zero of 200
prior draws clipping. Read literally against the spec as written — "K5 fires
here or nowhere" (task 3.7), "it fires in phase 3, which is the good news"
(§8 K5) — that would retire the criterion. It does not, because the toy is not
in the regime K5 is about: it drains ~40 particles per second from a pool of
~74,000, a buffer of about half an hour, while the charged-tRNA counter D13 and
D14 put at the centre of K5 debits ~553 residues per second against a pool of
order 10³, a buffer of seconds. The toy's margin is roughly three orders of
magnitude wider, so its census bounds the *mechanism* — that the debit clamps,
carries and repays exactly — and says nothing about the operating regime.
*Change:* the census that scores K5 is the assembled model's, in phase 14, not
phase 3's. Task 3.7's "since K5 fires here or nowhere" is struck, §8 K5's "it
fires in phase 3" is struck, and §10 R2's early warning now names phase 14.
Phase 3's census stands as what it is: evidence the mechanism is correct, and a
lower bound on nothing else. *Consistency note:* the spec was already divided
against itself here — §3's check 7 row scopes the check "Assembled" and task
14.8 already carries "(K5)" — so this records a resolution in favour of §3 and
§14 rather than inventing one. *Sections:* §8 K5, §10 R2, §11 task 3.7.
Approved before the fix landed; landed in PR #44.

### 2026-09-04 — a jump module's peer writes are gated on `written_states`, not on an edge

*Evidence:* phase 2 task 2.5 (and task 12.2, its consumer) said a jump
module's write to a peer's state is "gated on an outbound edge" and must throw
"once the edge is removed". But `_resolve_edges` rejects any edge naming a
species outside the Core A′ registry, and task 10.2 keeps the seventeen
transcripts outside it. Decay decrementing a transcript — the case the task
exists for — cannot carry an edge at all, so the gate as written was
unimplementable. Found at planning time, before the phase's code was written.
*Change:* a new protocol function `written_states(m)`, the jump-side twin of
phase 1's `contributed_states`, names the peer states a jump module's affects
modify. It must be a subset of `inputs(m)`; the write goes through the same
`u_inputs` view as the read, and an undeclared write throws at its first
firing, naming the module and the state. Where the written state *is* a
registry species the resolver additionally requires a mass or currency edge on
it, in either direction, since that write is a boundary crossing: `:in` where
the module consumes the pool, `:out` where it produces into it, and an outbound
edge alone satisfies the inputs contract for that species, because the pool is
an input only as the write channel. The converse is not required, because an
inbound edge may be a pure read. A non-registry state is gated by the
declaration alone. The alternatives — promoting the transcripts to registry
species, which contradicts task 10.2 and reopens the frozen registry; or gating
on `inputs()` membership alone, which makes a read and a write
indistinguishable so no write could ever be refused — were rejected. In the
same phase, `reactions` changed shape: it returns `Reaction(rate, affect!)`
values in local coordinates and the orchestrator wraps them over views of the
composed vectors, mirroring `dynamics`; that is the mechanism task 2.2 asked
for, not an amendment, and the alternative of passing an explicit index
context is recorded in the pull request. *Sections:* §11 tasks 2.2 (annotated),
2.5 and 12.2, annotated in place. Approved at planning time, before
implementation; landed in PR #43.

### 2026-09-03 — the contribution channel executes mass and currency edges in both directions

*Evidence:* phase 1 was specified around *outbound* edges executing. But a
module that draws on a pool it does not own — glycolysis consuming ATP that
recycling owns at PFK, the phosphotransferase cascade consuming the PEP
glycolysis owns — declares an *inbound* mass or currency edge, and that
consumption is a derivative term the owner never sees unless it is executed
too. Executing only the outbound side would leave every cross-module
consumption silently absent from the pools, and phase 6 would have needed a
second mechanism. *Change:* the phase-1 mechanism (`contributed_states` and
`contributions`) carries one signed term per non-owned dynamic registry species
on which an ODE module declares a mass or currency edge, in either direction;
the drift check of task 1.5 holds the two lists to each other on the same set.
A `:jump` module declares its mass and currency edges as before but may not
list contributions: its writes to a peer's state go through its reactions,
which phase 2 builds, so its edges are resolved and not held to the channel.
The phase's goal, done-when and task 1.5 are annotated in place. Registry,
chemostat and ownership rejections, and the both-ways drift check, are what
phases 6 to 9's edge declarations are now written against. *Sections:* §11
phase 1 (goal, done-when, task 1.5). Approved at planning time, before
implementation; landed in PR #42.

### 2026-09-03 — the coupling figure follows amendment 1

*Evidence:* the figure §2 cited, fig 1r, drew the lumped charging step as a
stochastic-block node, contradicting amendment 1 below, the scoping note it
illustrates, and this spec's reaction and module counts. *Change:* a new figure,
fig 1c in `dev/notes/figures/corea-coupling/`, redraws it with charging on the
ODE side, its three boundary crossings greyed and replaced by an intra-ODE
currency flow, a chemostat edge from the medium, and a deferred-counter edge from
the hook for translation's charged-tRNA debit. §2 now cites fig 1c; fig 1r is
kept unaltered as the record of the wave plan's placement. No decision changed.

### 2026-09-03 — four amendments, from the unshared-parameter finding

Logged together because one finding forced all four. Sections touched are named
so a reader can see what moved.

**Amendment 1 — the lumped tRNA charging step moves to the ODE block.**
*Evidence:* the scoping note already placed it there and always did — its
stochastic-side count is 52 reactions with no charging reaction, both tRNA
species appear in its list of dynamic ODE states, and the charged pool sits in
its energy-interface table beside the adenylate and guanylate species. The
frozen registry independently marks both species `:metabolite` regime.
`dev/plans/reduced-syn3a-wave-plan.md` assigned it to the stochastic block and
this spec inherited the error; the precedence rule at the top of this document
says the note wins on model content. Independently, the step fires 553 times a
second, which is 3.49 million events per cycle and 353 times every other
stochastic event combined, against K1's budget. *Sections:* §1, §2, §3, §4 D0,
D6, D7, D9, D13, D14, §5, §6 F1, F5b, T2, §8 K1, K5, §9, §10 R1, R8, R12, and
phases 2, 8, 9, 11, 13, 14. *Costs, not benefits:* the dominant adenylate forward
channel becomes two-stage and buffered by a pool size we assert (D13), and the
non-smooth boundary relocates onto that pool rather than disappearing (check 7).

**Amendment 2 — no parameters are shared across the boundary, accepted as a
first-stage assumption.** *Evidence:* the two namespaces are disjoint; the
derived parameter-coupling figure draws the blocks as separate panels with no
edge between them; metabolism went through parameter balancing and gene
expression did not; and 22 of 23 upstream ODE synthetase reactions are commented
out to avoid duplication. *Sections:* new §4 D13, plus corrections where the
document assumed otherwise — §4 D12 item 3 called the decay constant shared when
it is stochastic-block-local, and phase 17's done-when said the same. §1 now
distinguishes a parameter *shared across blocks* from one whose posterior is
moved by the other block's data. *Extension path:* recorded in
`dev/notes/modular-bayesian-inference-heterogeneous-modules.md`, not here.

**Amendment 3 — the inference architecture becomes a three-block conditional
scheme, with the path update deferred.** *Evidence:* with no shared parameter the
coupling is a fixed point in state space, so conditioning on the path is the
natural factorisation, and it restores gradients to the ODE block that the
survey's pseudo-marginal route forfeits. *Sections:* §4 D10 rewritten, D11 and
§8 K6 gain the observation that the conjugate direction *is* the ridge, §5 gains
the installed-package constraints, phase 16 gains the decision task. *One claim
retracted:* an earlier draft of D10 said conditioning on the path makes the cost
drains deterministic forcing and the ODE block smooth. It does not — the clip
depends on the ODE pool and hence on the ODE parameters — so K5 governs inside
the conditional step.

**Amendment 4 — the joint-versus-cut comparison arm is dropped.** *Evidence:*
`src/boundary.jl` decides what to hand over by matching parameter names
(`boundary_condition`, `:111-135`), and for Core A′ nothing matches, so the cut
degenerates into inferring each block alone. *Sections:* §4 D10's third rung and
D12's fourth check deleted, §6 F11 rebuilt around the composition executing and
auditing itself and F12 deleted, §6's abstract scalars lose the coverage gap, §7
gains the cut beside the iterative protocol, §10 R5 retired, phase 16 task 7 and
phase 17's formulation restated. *Cost:* C5 lost its primary evidence and now
rests on the composition running and enumerating its own departures, which is
weaker as a comparison and is stated as such.

**Nine defects fixed in the same pass**, each wrong independently of the
amendments: translation credited uncharged tRNA without debiting the charged pool
so the composed model had no steady state; the adenylate exhaustion time cannot
survive a mass-action rate law; check 8's summation theorem omitted `k_chg` from
the perturbation set and its frozen-expression variant saturates the tRNA pool;
the charging module's acceptance test was circular against its own calibration;
the reaction count in §2 was wrong twice; phase 8's sole-producer edge rule broke;
the charging module had no consumer standalone; R1 and R12's early warnings slip
later than promised; and three pre-existing contradictions on the reverse-channel
count, the channel total, and 141 against 144 seconds.

**Phases 9 to 12 renumbered** so charging precedes translation: charging 12 to 9,
transcription 9 to 10, translation 10 to 11, decay 11 to 12. Translation's rate
constant reads an ODE state that charging owns, so the old order had the
dependency backwards.

---

Reality will contradict this spec. When it does, the change is recorded here with
its date and its reasoning, never applied silently. An amendment names what
changed, what evidence forced it, and which phase found it. R15 is the most
likely source: a framework need surfacing during a module phase is an interface
bug, and the wave plan's rule stands — that is a conversation on the main branch,
not an edit on a module branch.
