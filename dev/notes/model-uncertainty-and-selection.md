# Model Uncertainty in InferCell: What to Build Now, What to Defer

*Companion to [the-case-for-inference.md](the-case-for-inference.md) and
[modular-bayesian-inference-heterogeneous-modules.md](modular-bayesian-inference-heterogeneous-modules.md).
Both of those documents are about uncertainty over parameters. This one is about
uncertainty over the model itself, and about whether it has to be designed in
early or can safely be left to the modular architecture to absorb later.*

---

## 1. The question

Every parameter in a whole-cell model is uncertain, and that is the argument
`the-case-for-inference.md` makes. But so is every *structural* choice: which
reactions exist, whether transcription is constitutive or bursty, whether a
regulator acts on one gene or several, whether a metabolic step needs an
explicit enzyme pool. Fixing those is a stronger claim than fixing a rate
constant, and it is the claim no whole-cell model in the field currently
examines.

So: bake model selection and model uncertainty into InferCell from the start, or
build parameter inference first and rely on modularity to make the structural
question tractable later?

**Recommendation: defer the machinery, fix the framing now, reserve three cheap
architectural hooks.**

---

## 2. Why deferring the machinery is safe

Sequencing is forced, not chosen. Every model-comparison quantity — a Bayes
factor, a stacking weight, an ELPD difference — is a functional of the
per-model posterior or likelihood. None can be computed before parameter
inference across heterogeneous blocks works and has been validated. Parameter
estimation is not being prioritised *over* structural uncertainty; it is its
prerequisite.

Building a comparison layer now would also mean building it against an inference
layer that is still moving. The ABC-SMC backend is scheduled for replacement by
a learned likelihood; the boundary protocol is scheduled to move from sequential
conditioning to EP. A comparison layer written against today's backends would be
rewritten twice.

---

## 3. Where the bet breaks: normalisation

The one place "modularity will absorb it later" is *not* safe is the
normalising constant.

Parameter inference needs the likelihood only up to a constant. Model comparison
needs the constant. And the modular architecture recommended in §5 of the survey
systematically discards it:

| Component | What happens to the evidence |
|---|---|
| EP sites | Unnormalised moment-matched factors. The EP energy gives an evidence approximation, but a poor one, and it degrades before the posterior moments do. |
| ABC-SMC | The implied "evidence" depends on the tolerance $\varepsilon$ and on the choice of summary statistics. Comparable across structures only under identical summaries *and* identical $\varepsilon$ — and then it compares ABC pseudo-evidences, not evidences. |
| Posterior surrogate (NPE) | Embeds the training prior. Already flagged for messages (§2.3c); the problem is worse for evidence, where there is no cavity to divide it out. |
| Cut posterior / SMI | **No marginal likelihood exists.** §3.6 of the survey says it directly: a cut posterior is not a posterior, and there is no joint distribution of which it is a conditional. Bayes factors under a cut are undefined. |

That last row is a genuine collision between the two existing documents, and
neither notices it. `the-case-for-inference.md` promises Bayes factors; the
survey's recommended architecture puts a per-module $\eta_k$ robustness layer on
top of an EP backbone, and the moment $\eta_k < 1$ anywhere, the object the
Bayes factor is built from has ceased to exist.

The corresponding claim in `docs/positioning.md` — "each candidate is an
`infer()` call; the comparison is a downstream calculation on the resulting
posteriors" — is roughly true for two NUTS runs on a differentiable block and
false for anything touching the ABC path. The two numbers being differenced are
not on the same scale.

---

## 4. The reframe that makes the bet safe

Drop marginal likelihood as the target.

A whole-cell model is irreducibly **M-open**. Every candidate structure is
wrong; the question is which is less wrong and for which predictions. Priors on
kinetic rates are diffuse and semi-arbitrary — scraped from the literature,
widened by judgement. Marginal likelihood is at its most fragile in exactly that
regime: a defensible-but-arbitrary two-fold widening applied across twenty
parameters moves a log-evidence far enough to reverse a conclusion. A Bayes
factor computed under literature-scraped priors is not a number worth defending
to a referee.

Predictive comparison — ELPD, cross-validated, or leave-one-stream-out — is
robust to prior width, and it is *already in the survey*, twice: as the
$\eta$-selection criterion in Stage 3 of §5.3 and as a validation check in
Stage 4. Adopting it as the model-comparison target does not add machinery; it
reuses machinery already scheduled.

### The unification

**$\eta_k$, a structural model index for module $k$, and the copies-and-tie
coupling strength $\lambda_k$ are the same kind of object.**

- $\eta_k \in [0,1]$: a continuous relaxation of "is module $k$'s likelihood
  right?"
- A model index $m_k \in \{1, \dots, M_k\}$: the discrete version of the same
  question.
- $\lambda_k$: the same question posed as "how hard should this module's copy of
  a shared quantity be tied to everyone else's?"

All three are per-module structural uncertainty. All three are selected by
held-out predictive performance. All three need somewhere to live that is not
the kinetic parameter vector. Build that slot once and you have built model
uncertainty, semi-modular robustness, and the surrogate-trust knob together.

### Why modularity genuinely does help

Structural uncertainty in a cell model is naturally *localised*. The question is
almost never "which of $2^K$ global models", it is "which transcription
mechanism, holding the rest of the graph fixed" — a per-slot question the module
graph already isolates. The combinatorial blow-up that makes global model
selection hopeless does not arise, provided the architecture records structural
identity per module rather than per system.

This is the sense in which the original instinct is right: modularity *will*
make structural comparison tractable later. It just will not make it *correct*
later, unless the normalisation semantics are tracked from the start.

---

## 5. The three hooks worth reserving now

Each is cheap now and painful to retrofit, because each changes a type signature
that everything downstream consumes.

**1. A `:structural` parameter role.** `src/parameters.jl` currently declares
`role` as one of `:rate`, `:initial_condition`, `:observation`. Add
`:structural`, plus a `structural_params` filter alongside the existing ones.
This is where $\eta_k$, $\lambda_k$, and mixture weights go. Without it they get
smuggled in as `:rate`, and `model_free_params` — which feeds the ODE/SSA
parameter vector — quietly does the wrong thing with them.

**2. A normalisation tag on the inference interface.** In `src/interface.jl`,
alongside `inference_mode`, something of the shape

```julia
evidence_semantics(m::AbstractSubModel) -> Symbol  # :normalised | :up_to_constant | :incomparable
```

The ABC and NPE routes return `:incomparable`. This costs one method with a
default now, and it is the single thing that prevents a future session
differencing two ABC pseudo-evidences and reporting the result as a Bayes
factor.

**3. Model identity as a slot, not a struct choice.** The comparison pair
already exists in the repository: `src/models/bursty_gene_expression.jl` versus
`src/models/stochastic_gene_expression.jl`. Today, "which one" is a fact about
which struct the caller constructed — invisible to the orchestrator and absent
from anything `infer()` returns. Make it a fact the assembled system *records*:
the orchestrator's identity includes which candidate filled each slot. That is a
field, not a feature, and it is what lets a comparison layer exist later without
rewriting `infer()`'s return type.

None of the three commits the project to computing a model-comparison quantity.

---

## 6. The design discipline, extended

The project's organising question is: *for every module boundary, how does
inference cross this interface?*

The structural version is one clause longer:

> **How does inference cross this interface when the module on the other side
> might be the wrong model?**

It is free to ask and it catches design errors at the point they are cheap. A
boundary that passes a point estimate fails the parameter version. A boundary
that passes a posterior conditioned on one structure, with no record of which
structure, fails the structural version — and that is the failure mode the
current `boundary_condition` KDE path would walk into unremarked.

---

## 7. What this changes in the argument for inference

`the-case-for-inference.md` under-sells itself. Its three-facts argument is
*stronger* for structure than for parameters:

- MC4D does not merely fix ~1000 parameters to point estimates; it fixes which
  reactions exist. That is the larger implicit claim of perfect knowledge.
- Structural error is not a perturbation. Constitutive versus bursty
  transcription changes the *shape* of the mRNA distribution — a qualitative
  change no amount of rate uncertainty reproduces. Parameter uncertainty moves
  a prediction; structural uncertainty can move which predictions are even
  well-posed.
- The asymmetry favours the argument. Parameter uncertainty is at least
  acknowledged-but-unrepresented in the field. Structural uncertainty is not
  acknowledged at all. It is the bigger gap, and it costs nothing to name it
  while declining to solve it yet.

---

## 8. Open problem this exposes

Worth adding to the survey's §6 list, next to SBC for message-passing systems —
same shape of gap, arguably more consequential:

**Evidence and predictive comparison under modular approximation.** Given a
system inferred by EP with heterogeneous site solvers and a per-module
robustness layer, what model-comparison quantity is actually available, and what
does it converge to? Specifically: (a) how badly does the EP energy approximate
the log evidence when sites are ABC- or surrogate-backed rather than analytic;
(b) what is the right leave-one-stream-out predictive when the module supplying
the stream is itself the object under comparison; (c) is there a coherent
comparison criterion under SMI, where the marginal likelihood does not exist,
beyond stacking on held-out predictive density.

---

## 9. Decision summary

| | Now | Later |
|---|---|---|
| Parameter inference across heterogeneous blocks | Build | — |
| Bayes-factor engine | Do not build | Probably never; prefer predictive comparison |
| Predictive comparison (ELPD, LOO, leave-one-stream-out) | — | Build with the Stage 3/4 robustness layer, which needs it anyway |
| `:structural` parameter role | Add | — |
| `evidence_semantics` interface tag | Add | — |
| Model identity recorded per slot | Add | — |
| Structural framing in the positioning and case-for-inference documents | Fix | — |
