# InferCell: Vision and v0.0.1 Scope

**Status:** Definitive scoping document (brainstorm-validated, 2026-05-29).
**Supersedes:** the framing and module-set discussion in `2026-05-12-v0.0.1-scoping.md`
(which it refines, not replaces — the biology and engineering roadmap there still stand).
**Extends:** `2026-05-29-boundary-representation.md` (the trajectory boundary object) and
`overview.md` (the inference-graph architecture).
**Resolves:** the six strategic questions raised in the 2026-05-29 brainstorm —
framing, module chaining, the role of ML, the ultimate deliverable, the data timeline,
and the place of neural operators.

---

## TL;DR

1. **Framing is sequenced.** The *vision* is a whole-cell model that does inference (a
   scientific instrument). Every *shipped* version through v0.0.x stays a methods
   contribution proven on a deliberately minimal cell. It becomes a biological instrument
   at v0.1+ once real data arrives and the machinery is proven. This reconciles the
   "WCM that does inference" ambition with the "faithful testbed" creed.
2. **North star: put calibrated error bars on a well-mixed reduction of the
   Luthey-Schulten / Thornburg syn3A minimal cell.** Their published simulator becomes our
   synthetic-data generator with known ground truth (their literature parameters *are* the
   injection); their parameter set is our comparison baseline. Complementary to MC4D, not
   competing.
3. **ML's role is settled: a neural *likelihood* (NLE) for the stochastic block, not a
   neural module.** Likelihoods compose; posteriors do not. A learned, differentiable
   `p(x | θ)` for the SSA block replaces ABC-SMC, gives a differentiable boundary message,
   and *is* the learned trajectory-boundary representation §8 of the boundary note calls for.
   It is one rung of a **likelihood ladder**, used only when cheaper exact/approximate
   options are exhausted.
4. **Joint vs. modular is a false dichotomy.** Joint inference = the *oracle* (gold standard
   on small problems, often impossible to assemble at scale). Modular cut→EP over
   trajectories = the *product* (what scales and composes). The gap between them is the
   scientific result. **Central claim: modularity's value flips from cost to benefit as the
   model moves from well-specified to misspecified** — the cut *quarantines* module error
   that the joint lets leak everywhere.
5. **v0.0.1 biology is unchanged: ODE metabolism + ODE bulk gene-expression + SSA bursty
   regulator**, with the SSA block now NLE-backed. "ODE-SSA-ML" read correctly = this, where
   ML is the SSA *backend*. The autocatalytic loop (translation → enzyme → metabolism) is
   retained because it is the testbed for the EP feedback demonstration.
6. **Particle MCMC is killed** (not deferred): never useful at WCM scale. Neural ODEs as
   black-box biology are rejected; they survive only as a v0.1+ *discrepancy* term. Fourier
   neural operators are parked as a v0.2+ spatial-emulator bet.

---

## 1. Framing — what InferCell is, sequenced

Two earlier framings were in tension: the docs commit to "*not* a simulation platform… a
methods contribution" (the testbed creed), while the ambition is "a WCM that does
inference" (a scientific instrument). The resolution is to **sequence them on the same
axis**, not choose:

- **Vision (the pull):** a whole-cell model whose native operation is Bayesian inference
  with honest UQ — the instrument the field lacks.
- **Every shipped version (the discipline):** a methods contribution. The cell is a
  faithful testbed; it need not be biologically *right*, only a faithful *probe of the
  inference machinery*. This holds through all of v0.0.x.
- **The transition:** v0.1+ earns the "instrument" label by conditioning on real data and
  comparing against an external mechanistic ground truth (the Luthey-Schulten cell).

The practical effect: we never over-claim biology before the machinery is proven, and we
never lose the north star that makes the methods worth building.

---

## 2. North star — the well-mixed minimal cell with UQ

The concrete v0.1+ target the methods work walks toward:

> **A well-mixed reduction of the Luthey-Schulten / Thornburg syn3A minimal cell — their
> metabolism + gene-expression + tRNA-charging biology minus the RDME spatial lattice and
> chromosome Brownian dynamics — with full posterior UQ over its kinetic parameters.**

Why this target and not a model of our own:

- **It hands us a ground-truth generator.** Their simulator runs forward from ~1000
  literature parameters. Those parameters are a known injection; their outputs are
  synthetic data we can condition on and recover.
- **It hands us a comparison baseline.** "We put error bars on the Luthey-Schulten minimal
  cell" is a one-line contribution, complementary to MC4D (which structurally cannot do UQ),
  not a competitor to it.
- **The well-mixed limit is honest.** The CME/ODE pieces of their model *are* the
  non-spatial limit; we are not inventing biology, we are removing the spatial machinery
  whose engineering cost dwarfs the inference contribution at this stage.

Spatial dynamics return later (§9), at which point the FNO bet becomes relevant.

---

## 3. The boundary object (point 1, largely pre-decided)

Settled by `2026-05-29-boundary-representation.md` and unchanged here:

- **What crosses a module boundary is a trajectory** `p(a(·))` — a distribution over the
  shared state path — not a parameter posterior. Modules in a dynamical cell couple through
  shared state evolving in time; the trajectory is the Markov blanket between them.
- This is **Axis 1** (the interface object). It is orthogonal to **Axis 2** (the information
  flow: plug-in → cut → EP → Gibbs → joint).
- **Staged scope (option C):** deterministic-ODE modules use the trajectory as a
  representation of the parameter push-forward (zero new degrees of freedom). The object is
  designed so that `process noise = 0` recovers this case, making full state-space inference
  the natural sequel rather than a re-architecture.

What this brainstorm adds: with the SSA block now carrying an **NLE likelihood** (§4), the
ODE↔SSA boundary is differentiable in both directions. The SSA module consumes the upstream
trajectory `a(·)` as input and contributes a learned likelihood — so the trajectory is the
interface object on *both* sides, and the "single differentiable computational graph" claim
stops being aspirational.

---

## 4. ML — a neural likelihood, governed by a ladder (point 2)

### The decision rule

To place a module in one graph — joint inference *or* a coherent cut/EP message — you need
an object you can multiply against the other modules' contributions. **A likelihood composes;
a posterior does not** (a posterior has a prior baked in; combining posteriors means tracking
prior division by hand — that *is* the EP/cut bookkeeping). Therefore:

- **NPE (amortized posterior)** is rejected as the spine: it yields the one object that
  doesn't compose, and it's a posterior over *parameters* — the boundary object the note
  argued against. It remains a fine *optional* amortized readout to bolt on later.
- **Neural ODE (mechanism-free dynamics)** is rejected as a module: its weights are not
  biology, a posterior over a weight matrix answers no scientific question, and it is
  directly off-target for a north star that is a *fully mechanistic* cell. It survives only
  as a discrepancy/closure term (§8).
- **Differentiable surrogate likelihood (NLE)** is the spine. A learned, differentiable
  `p(x | θ)` for the SSA block: replaces ABC-SMC, supplies the differentiable boundary
  message, keeps every parameter mechanistic (the real SSA generates the training data — ML
  only makes its intractable likelihood tractable), and *is* the learned trajectory-boundary
  object. It collapses points 1 and 2 into one mechanism.

### The likelihood ladder (the real design discipline)

ML is not applied everywhere. For each module, in order of preference:

1. **Exact differentiable likelihood** — use it (linear birth-death, Gaussian observation).
2. **Cheap differentiable approximation** — chemical Langevin SDE, moment closure,
   finite-state projection. Mechanistic and differentiable; right for "mildly stochastic"
   modules.
3. **Learned surrogate likelihood (NLE)** — only when 1 and 2 genuinely fail (intractable
   *and* no cheap approximation).
4. **Exact-but-slow (ABC-SMC / particle filter)** — retained as a **validation oracle**,
   never deleted.

The selling point is *"InferCell picks the right inference object per boundary,"* not
*"InferCell is neural."*

### Implementation guidance for v0.0.1

- **Start with a synthetic likelihood, not a normalizing flow.** A conditional Gaussian or
  mixture-density network over summary statistics (Wood 2010) is the minimal NLE: a few
  hundred lines of differentiable Julia, no exotic dependency. This guts the Julia-SBI
  tooling risk (Open Question #7). Upgrade to a flow only if the summary likelihood is too
  crude.
- **The SSA module is the two-state bursty regulator** (scoping-doc Step 5). Its intrinsic,
  discrete, multimodal noise is exactly the load-bearing case with no ODE equivalent.

### Edge cases we commit to handling

| Risk | Why it bites InferCell specifically | Mitigation |
|---|---|---|
| **Overconfidence** | A learned likelihood can be *sharper than the truth*; error hidden in weights. A UQ-first project shipping overconfident posteriors is the worst failure. | **Simulation-based calibration (SBC)** on every surrogate (`2026-05-22-calibration-curve` note); keep the exact oracle. |
| **Information loss in the tails** | Surrogates are worst where bursting/low-copy/extinction biology lives; can't represent atoms (exact-zero counts). | Validate against the SSA in the tail regime; prefer trajectory/flow likelihood over crude summaries where bursting matters. |
| **Off-distribution composability** | NLE trained on a module *in isolation* gets queried in-loop with upstream trajectories from a different distribution. | Condition the NLE on the input trajectory, or refine in-loop; flag where amortization breaks. |
| **Dimensionality** | sim-from-prior coverage is exponential in `dim(θ)`. | Keep each module's free parameters small (~5–8); forbids a monolithic stochastic module. |
| **Surrogate gradients ≠ true gradients** | NUTS exploits network smoothness, mixes confidently to the wrong place. | SBC + exact oracle. |

The throughline: **the simulator and one exact backend stay alive forever — as oracle and as
data source. ML never replaces the mechanism, only its intractable likelihood.**

---

## 5. Joint vs. modular — two objects, two jobs (point 2 cont.)

- **Joint inference = the oracle.** Once NLE makes every block differentiable, one sampler
  over the whole graph is feasible *on a small problem* and is the gold-standard posterior.
  It does **not** scale, and for an assembled WCM it is often *impossible to form* (modules
  from different labs/formalisms; you may hold only a simulator or a learned likelihood). So
  joint is the **validator**, not the product.
- **Modular (cut → EP over trajectories) = the product.** What scales, composes, and
  *is* the contribution. Approximate, so it needs the oracle to measure its cost.
- **The deliverable is the gap between them**, exactly the §9 de-risk experiment of the
  boundary note: joint NUTS [truth] vs. trajectory-boundary graph vs. independent modules
  [naive baseline].

### The central scientific claim

Modularity is not only a tractability concession. **Cut posteriors exist to stop a
misspecified module from corrupting the rest** (Plummer 2015; Jacob et al. 2017). The joint
lets one wrong module's error leak into every parameter; the cut quarantines it. Therefore:

- **Tier A (well-specified):** the joint is exact; modularity *costs* calibration.
- **Tier B (misspecified — richer generator, ultimately the Luthey-Schulten simulator):**
  modularity *quarantines* the error; the cut/EP graph can be **more honest than the joint.**

> **Claim: modular inference's value flips from cost to benefit as the cell moves from a twin
> to a misspecified model.** This directly answers the "severing feedback" worry — severing
> is harmful when the model is right, protective when it is wrong.

### The coherence ladder for v0.0.1

- **L1 — cut (one-way):** the implemented product and the cheap baseline. [DONE — Step 3]
- **L2 — EP / iterative message-passing over trajectories:** demonstrated on the v0.0.1
  autocatalytic loop, where feedback (translation → enzyme → metabolism → translation) makes
  a one-way cut structurally wrong. Restores the coherence the cut discards. EP on a
  linear-Gaussian state-space model *is* the Kalman/RTS smoother, so this is existing
  machinery on the WCM graph. [NEXT — Step 4]
- **Joint NUTS:** the oracle, run only at toy scale. Feasible *because* NLE made the SSA
  likelihood differentiable.
- **L3 — Particle MCMC: KILLED.** Never useful at WCM scale.

---

## 6. The v0.0.1 deliverable

**Biology (unchanged from `2026-05-12-v0.0.1-scoping.md`):** ODE metabolism + ODE bulk
gene-expression (3–5 named syn3A genes, ≥1 metabolic enzyme closing the loop, ≥1
ribosome-component gating translation) + SSA two-state bursty regulator. ~10–15 free
parameters, ~6–8 dynamic states, one closed autocatalytic loop.

**Inference stack:**
- Joint NUTS within the differentiable (ODE) block. [DONE — Step 2]
- **NLE** (synthetic-likelihood / MDN to start) for the SSA block, replacing ABC-SMC as the
  production backend; ABC-SMC retained as oracle.
- L1 cut boundary protocol. [DONE — Step 3]
- L2 EP/message-passing across the ODE↔SSA boundary, closing the loop. [Step 4]
- Joint NUTS over the *whole* graph as the small-scale oracle (enabled by the NLE).
- SBC on the NLE surrogate, mandatory.

**The headline figure:** posterior recovery under three inference schemes — **joint
(oracle) vs. cut vs. EP** — across **Tier A (twin)** and **Tier B (misspecified)**, showing
the cost→benefit flip. Plus architecture schematic, joint ODE posteriors with ground truth,
and posterior-predictive overlays.

**Success criteria:**
- Tier A: 90% CIs cover ground truth for all parameters; chains mix; under ~1 hour on one
  HPC node.
- NLE: passes SBC (calibrated, not overconfident); NLE-based posterior matches the
  ABC-SMC/joint oracle within tolerance on the twin.
- Tier B: the modular (cut/EP) posteriors remain calibrated where the joint degrades —
  demonstrating the central claim.

---

## 7. Data timeline (point 5)

| Stage | Data source | Purpose |
|---|---|---|
| **v0.0.1** | Synthetic from our own toy: Tier A (twin), Tier B (refit under a richer generator) | Correctness, then honesty-under-misspecification. Known injection throughout. |
| **v0.1 bridge** | Synthetic from the *actual Luthey-Schulten well-mixed simulator* | Tier B against a real published cell, with their literature parameters as ground-truth comparison — without wet-lab data yet. |
| **v0.0.2+** | Real experimental: Breuer et al. 2019 (metabolomics → block 1); Hutchison et al. 2016 (bulk RNA-seq → block 2) | The transition from methods contribution to biological instrument. |

Synthetic-first is not timidity: it is the only regime where the injection is known, so it is
the only regime where we can *measure* whether the boundary protocol's UQ is honest.

---

## 8. Parked and rejected

- **Particle MCMC (L3):** killed. Not deferred.
- **Neural ODE as a biological module:** rejected (black-box biology undercuts the
  mechanistic-UQ pitch). **Survives as a discrepancy term:** a small learned correction
  `du = f_mech(u, θ) + g_NN(u)` (a universal differential equation) absorbing
  misspecification — the neural form of the "Kalman process-noise" idea. **v0.1+ only**, tied
  to Tier-B-against-Luthey-Schulten.
- **NPE (amortized posterior):** not the spine; optional amortized readout, later.
- **Fourier neural operators (point 6):** parked as a **v0.2+** bet. When spatial dynamics
  return, FNOs are a candidate differentiable *emulator* of the RDME reaction-diffusion field
  — an instance of likelihood-ladder rung 2/3 for a spatial module — potentially replacing
  Thornburg's explicit lattice discretisation with a faster learned operator. Speculative;
  revisit only after the well-mixed instrument exists.

---

## 9. Relationship to existing roadmap

- The biology and engineering steps in `2026-05-12-v0.0.1-scoping.md` stand, with **one
  change: the SSA production backend becomes NLE, not ABC-SMC** (ABC-SMC demoted to oracle),
  and **one addition: joint NUTS over the whole graph as the oracle**, now feasible.
- Step 4 (Level 2 boundary protocol) is reframed as **EP/message-passing over trajectories**,
  the L2 rung above, demonstrated on the autocatalytic loop.
- The boundary-representation note's §9 de-risk experiment is promoted to **the v0.0.1
  headline result.**
- Open Question #1 (boundary representation) and #9 (misspecification) are answered; #7 (Julia
  SBI tooling) is de-risked by starting NLE at the synthetic-likelihood end of the ladder.
