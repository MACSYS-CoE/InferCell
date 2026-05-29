# Boundary Representation: Trajectories, Not Parameter Summaries

**Status:** Brainstorm / proposed design — extends `2026-03-31-v1-architecture.md`.
**Date:** 2026-05-29
**Scope:** What object crosses a module boundary in the inference graph. Resolves Open
Question #1 of the v1 architecture ("Boundary representation: what exactly crosses the
boundary?") and frames the scientific-identity fork it exposes.

---

## TL;DR

1. Modules in a dynamical whole-cell model are coupled through **shared state evolving over
   time** — a metabolism module hands TX/TL a time course of ATP, not a value of `k_cat`.
   So the natural boundary object is a **distribution over the shared state path**,
   `p(a(·) | D)`, not a distribution over the upstream module's parameters.
2. This is a choice on a different axis from the cut/Gibbs/PMCMC ladder. There are **two
   independent axes**: *what* the interface object is, and *how* information flows across it.
   "Trajectory vs. parameter" is the first; "cut vs. EP vs. Gibbs vs. PMCMC" is the second.
   They compose.
3. Passing trajectories does **not**, by itself, change InferCell from an inverse-problem
   sandbox into a forward whole-cell simulator. That depends on whether the path carries
   degrees of freedom beyond the parameters — and on the *deliverable*, which stays
   inferential.
4. **Proposed v1 scope: option (C), staged.** Deterministic-ODE modules use the trajectory
   boundary as a representation of the parameter push-forward (pure parameter inference,
   zero new degrees of freedom). Design the boundary object so that turning on process
   noise is a strict generalization — `process noise = 0` recovers the v1 case — making
   full state-space inference the natural sequel rather than a re-architecture.

---

## 1. The physical fact

In a dynamical system, modules are not coupled through each other's parameters. They are
coupled through **shared state that evolves over time**. Metabolism does not hand TX/TL its
rate constants; it hands it a trajectory of ATP, NTPs, and amino acids, `a(t)`. Translation's
rates at time `t` depend on `a(t)`. The upstream parameters `θ_M` matter to TX/TL *only
through* `a(·)`.

So `a(·)` — the path of the coupling species over the observation window — is the Markov
blanket between the modules in time. Everything upstream of it (`θ_M`, metabolism's internal
states) is conditionally irrelevant to TX/TL given `a(·)`. That is the precise sense in which
the **trajectory, not the parameter vector, is the natural interface object.**

---

## 2. Three candidate boundary objects

Take metabolism **M** (ODE / NUTS) feeding `ATP(t)` into TX/TL **T**:

| Object | What crosses | How T uses it | Uncertainty carried |
|--------|--------------|---------------|---------------------|
| **Point plug-in** | `θ̂_M` → one path | deterministic `a(t)` | none |
| **Parameter posterior** (current implicit "cut") | `p(θ_M \| D_M)` | **re-simulate M** inside T's loop per draw | upstream parameter uncertainty only |
| **Trajectory** (proposed) | `p(a(·) \| D_M)` directly | condition on the path; **no re-simulation, never sees `θ_M`** | full path uncertainty, incl. intrinsic noise |

The trajectory object is represented as, e.g., a posterior mean function `μ_a(t)` plus
covariance `K_a(t, t′)` (a Gaussian process over the path), or a set of `N` sample paths.

---

## 3. Distinction from the current cut posterior

Honest version, because a *correctly implemented* parameter-cut already induces a trajectory
distribution via re-simulation:

- **Under a one-way cut, a deterministic upstream module, and no compression, the parameter
  posterior and the trajectory encode the *same* posterior.** The trajectory is just the
  push-forward of `p(θ_M | D_M)` through the solver. In that case the change is
  **representational/architectural, not statistical.**

But that representational change buys three things that matter specifically for InferCell:

1. **Cross-formalism composability (the big one).** Modules live in different backends —
   NUTS (Turing) for ODEs, ABC-SMC for SSA. "Pass `θ_M` and re-simulate inside the consumer"
   needs a shared computational graph the two backends do not share. A **trajectory is
   backend-agnostic**: any module consumes `a(t)` regardless of how it was produced. This is
   what actually realizes the project's "interfaces over implementations" principle;
   parameter-passing quietly violates it.
2. **It carries noise parameters structurally cannot.** Fix `θ_G` exactly and the SSA
   gene-expression module's `m(t)` is still random — intrinsic noise is a property of the
   path, not the parameters. A parameter boundary cannot represent it; `p(m(·) | D_G)` can.
   The moment one coupled module is stochastic, the trajectory object is strictly richer.
3. **It makes Open Question #1 well-posed.** "What crosses the boundary?" was vague.
   "What is the right representation of `p(a(·))`?" is a concrete, measurable design axis:
   point path → per-time marginals → full GP with covariance → particle paths.

It becomes a genuinely *different* posterior once (a) a coupled module is stochastic, (b) the
path is **compressed** (GP moments, low rank), or (c) back-flow is allowed (next section).

---

## 4. Distinction from expectation propagation

EP and the trajectory boundary live on **different axes**; conflating them is the usual error.

| | **Axis 1 — interface *object*** | **Axis 2 — information *flow*** |
|---|---|---|
| Options | point → parameter posterior → static summary `z` → **full state path `p(a(·))`** | plug-in → **cut (one-way)** → **EP (iterative, moment-matched)** → Gibbs → full joint PMCMC |
| v1 today | parameter-ish | cut / sequential conditioning (L1); Gibbs (L2), PMCMC (L3) on the roadmap |

- "Trajectory vs. parameter" is **Axis 1** (the object).
- "Cut vs. EP vs. Gibbs vs. PMCMC" is **Axis 2** (the flow). EP is defined by being
  *iterative and bidirectional* (downstream data `D_T` flows back and updates the interface
  message, partially restoring the coherence a cut throws away) and by *moment-matched factor*
  messages (classically low-dim Gaussians).

They compose. The natural synthesis here is **EP whose messages are Gaussian *processes* over
`a(·)`** — moment-match a mean function + covariance kernel and iterate. Key fact: **EP on a
linear-Gaussian state-space model is the Kalman/RTS smoother**, and EP generalizes the smoother
to non-Gaussian module likelihoods. "Trajectory boundary + EP flow" is therefore the existing
Kalman/state-space machinery applied to the WCM graph, not an exotic new method.

The trajectory representation is **orthogonal to** the L1→L2→L3 flow ladder and makes each rung
feasible across formalisms. EP slots in as an Axis-2 option between the cut (L1) and Gibbs (L2).

---

## 5. The identity fork: inverse problem vs. forward WCM

Adopting trajectory boundaries raises a real concern: InferCell was pitched as solving the
**inverse** problem (focus on parameters), not as a forward simulator. Does this drift toward
"a WCM that does forward sim"? The distinction needs to be precise.

**The dividing line is not forward vs. inverse.** Inference always requires a forward map —
NUTS pushes gradients through a forward ODE solve; ABC-SMC forward-simulates the SSA hundreds
of times. InferCell has run forward simulations from day one; that is the engine, not the goal.

The line that actually moved is **what counts as an unknown**:

- **Narrow inverse (original pitch):** unknowns are parameters `θ`; the path `x(·) = f(θ)` is a
  deterministic by-product, not an object of inference.
- **State-space inverse (where trajectory boundaries point):** unknowns are `θ` *and* `x(·)`;
  the target is the joint `p(θ, x(·) | D)`, with the path a first-class latent.

Both are inverse problems. The second is a richer one — and it is exactly what a PTA/Kalman
background is built for (jointly inferring timing parameters, a latent stochastic process, and
a signal, where the hard part is separating them).

**Crucially, trajectory-passing alone does not force the harder regime.** For a deterministic
ODE module, `p(a(·) | D_M)` is just the parameter push-forward — zero new degrees of freedom,
still pure parameter inference. You cross into genuine state-space inference only when the path
gains freedom the parameters lack:

1. a coupled module is **intrinsically stochastic** (the SSA module),
2. you add **process noise** to the path, or
3. you add **model discrepancy** (letting `x(·)` deviate because the mechanism is wrong).

The gravitational pull toward a forward-WCM *product* comes not from the boundary object but
from the **deliverable**. The guardrail is a scoping creed:

> **InferCell infers; it simulates only because inference demands it. State-paths are an
> inference interface and a latent target — not a product. The cell need not be realistic,
> only a faithful testbed.**

"Synthetic data first" already enforces this: the cell does not have to be *right*, it has to
be a faithful *testbed for inference*.

---

## 6. Decision: staged scope for v1

| Option | Path status | Character |
|--------|-------------|-----------|
| (A) Conservative | deterministic ODE; trajectory = representation of the push-forward | pure parameter inference, zero new DOF; fully on the original pitch |
| (B) Full state-space | `x(·)` first-class everywhere, process noise / discrepancy | richer, on-brand for the comparative advantage, harder to identify/validate |
| **(C) Staged (chosen)** | ship (A); design the object so `process noise = 0` recovers (A) | protects the focused first paper; banks state-space as the obvious sequel; no re-architecture |

**Chosen: (C).** Rationale: it keeps the v1 deliverable a focused inference sandbox while making
the unifying abstraction state-space inference — the thing that genuinely turns the two stapled
backends (NUTS, ABC-SMC) into *one* framework (they are the `process noise = 0` and
`intrinsic-noise` special cases of the same object).

---

## 7. Costs and risks (not to be romanticized)

- **Identifiability.** Joint `θ` + `x(·)` + process noise + observation noise is a classic
  trade-off swamp — process noise can masquerade as parameter uncertainty (cf. red-noise /
  white-noise / signal degeneracy in PTA). Tractable, but real work. Largely deferred under
  scope (C), but the representation must not foreclose handling it.
- **Path priors / process-noise models** become modelling choices that must be justified.
- **Validation doubles** under (B): claims about latent *state* recovery, not just parameters.
- **Dimensionality.** A path over `T` timepoints × several species is high-dimensional; a full
  GP is `O(T³)` and naive per-timepoint summaries destroy temporal correlation. Control it with
  a **Markovian / state-space GP** representation — full temporal covariance at `O(T)` via
  filtering/smoothing. Choosing the representation (SSM order, kernel, particle count, number of
  moments) *is* the empirical boundary-format question, now concrete.

---

## 8. Implications for the existing roadmap

- The L1/L2/L3 flow ladder (`2026-03-31-v1-architecture.md`) is unchanged and orthogonal; the
  trajectory object makes each rung work across formalisms.
- Open Question #1 (boundary representation) is answered: **a (Markovian) GP / particle-path
  representation of `p(a(·))`**, with representation fidelity as the dial to sweep.
- Open Question #3 (misspecification) gets a natural home: model discrepancy lives in the
  process-noise term of the state-space form — nowhere to put it under `x(·) = f(θ)`.

---

## 9. Next step: the de-risk experiment

The riskiest assumption is that a low-dimensional boundary is a near-sufficient statistic — that
modular inference over `p(a(·))` recovers a posterior close to the true joint. Test it cheaply
and first, on the two ODE/NUTS modules where the joint is tractable:

- **Modules:** TX/TL ⊕ light metabolism, where the full joint posterior is HMC-computable.
- **Compare:** (a) full joint NUTS [ground truth], (b) trajectory-boundary inference graph,
  (c) independent modules [naive baseline].
- **Measure:** calibration / bias cost of the cut as a function of boundary representation
  (point path → GP moments → full covariance → particle paths) and interface dimensionality.

This both validates the architecture and stands alone as a methods result: *when can you
modularize state-space inference in a coupled cell model, and what boundary representation does
it take?*
