# Reframing "Uncertainty Propagation" as "Information Propagation"

Brainstorm notes for the MACSYS Theme 1 talk (May 2026). Captured for later
review — decisions not yet folded into `talk.tex`.

## The question

The talk repeatedly says InferCell *propagates uncertainty* through the
model. Could it equivalently — or additionally — talk about propagating
**information**? Is uncertainty a kind of information?

Short answer: **yes, and the talk already mixes both languages without
unifying them.** That unannounced mixing is the opportunity.

## Where each frame already appears in the talk

- **"Information"** appears structurally:
  - Slide 5: *"No information flows between modules"*
  - Design principles: *"Information flow as design test"*
  - Gibbs panel (How Much Coupling): *"Information flows both ways"*
- **"Uncertainty"** appears representationally:
  - *"propagate uncertainty"*, the boundary-protocol arrow label
  - epistemic vs aleatoric (MC4D backup slide)
  - credible intervals on predictions

So in the current deck: **modules carry information; parameters carry
uncertainty.** Related but never unified.

## Why uncertainty *is* information (the technical case)

1. **Shannon.** Information is *defined* as reduction in uncertainty.
   Entropy $H(\theta)$ quantifies uncertainty; mutual information
   $I(\theta; \text{data}) = H(\theta) - H(\theta \mid \text{data})$ is the
   information the data carries about $\theta$. A posterior's narrowness
   relative to the prior **is** the information extracted.
2. **The posterior is an information object.** Not "uncertainty + a mean"
   — its whole shape (correlations, tails, multimodality) is information.
   A point estimate is the operation that *destroys* almost all of it.
   Standard WCMs aren't merely "no UQ" — they are
   **information-destroying at every module boundary**.
3. **The talk already uses information geometry implicitly.** The
   identifiability slide (tight posterior vs posterior ≈ prior) is a
   Fisher-information / sloppy-models statement. The "What your model
   actually depends on" backup slide is an information-flow diagnostic.
4. **BOED.** Listed as a payoff. Expected Information Gain (EIG) *is*
   the objective. The downstream story is already information-theoretic,
   just unannounced.

## The unifying move: braid forward and backward

Don't swap "uncertainty" for "information." **Braid them** as forward vs
backward operations on the same object — the probability distribution
living on each module boundary:

- **Forward** ($\theta \to$ predictions): a distribution propagates →
  **uncertainty propagation**.
- **Backward** (data $\to \theta$): a distribution sharpens →
  **information flow**.
- Standard WCMs collapse distributions to points at every boundary →
  *simultaneously* uncertainty-discarding and information-destroying.
- InferCell carries full distributions → both directions become possible.

This converts the current ambiguity into a feature: *the same
architectural choice* — distributions at boundaries — gives you both
forward UQ and inverse inference.

**Candidate one-liner for the inference-graph slide:**

> "Every boundary carries a distribution. Run it forward: uncertainty
> propagates. Run it backward against data: information flows."

## Audience-specific hooks (MACSYS Theme 1)

Information framing opens doors uncertainty framing doesn't:

- **Cells as information processors** (Bialek, Tkaçik). Biologists
  already accept that signalling pathways compute and transmit
  information across boundaries. *"Our model boundaries should carry
  information the way the cell's do."* Strong rhetorical bridge.
- **Data processing inequality** (see deep-dive below). A polyglot
  pipeline with point estimates at every boundary is a maximally lossy
  channel. InferCell is the lossless version.
- **Epistemic = missing information; aleatoric = channel noise /
  intrinsic entropy.** The MC4D backup slide already separates these —
  relabelling as "information-recoverable" vs "information-floor" makes
  the distinction land harder.

## Deep dive: the Data Processing Inequality (DPI)

### The result

If $X \to Y \to Z$ is a Markov chain (Z depends on X only through Y):
$$I(X; Z) \leq I(X; Y).$$

Once you've processed $Y$ from $X$, nothing downstream can recover
information about $X$ that $Y$ didn't already carry. Processing can
preserve information; it can never create it.

### A WCM is exactly such a chain

The polyglot pipeline on the "Why Pure Julia?" slide, read as a Markov
chain in $\theta$:
$$\text{data} \to p(\theta_A \mid \text{data}) \to \hat\theta_A \to
  \text{module B output} \to \hat\theta_B \to \ldots$$

Each arrow is a channel. DPI applies at every arrow. Information about
the biology that survives to the final prediction is bounded by the
**worst** channel in the chain.

### Point estimates are (near-)maximally lossy channels

The posterior $p(\theta \mid \text{data})$ is a **sufficient statistic**
— by construction it contains all the information the data carries
about $\theta$. The map
$$p(\theta \mid \text{data}) \longmapsto \hat\theta = \mathbb{E}[\theta \mid \text{data}]$$
is deterministic. By DPI it can only destroy information — and it
destroys nearly everything: variance, shape, tails, and every
off-diagonal entry of the covariance.

The k_tx · γ example (slide 4) is a clean illustration. Data constrains
$k_{tx}/\gamma$ tightly even when it constrains neither alone — an
anti-correlation in the joint posterior. Collapse to $(\hat k_{tx}, \hat\gamma)$
and the correlation is gone. The information the data spent its effort
encoding is the first thing thrown away.

### Why "polyglot" causes this

Module boundaries crossing languages can't easily carry distributions.
Serialising a posterior through a JSON dict or subprocess pipe is awkward;
passing a float is trivial. So the architectural choice (polyglot, IO at
boundaries) **forces** the information-theoretic choice (point
estimates). The DPI loss isn't an oversight — it's baked into the
substrate.

### "Lossless" — with caveats

InferCell passes distributions (or samples) across boundaries, not
points. The channel at each boundary is then the identity on the
relevant sufficient statistic, so DPI saturates:
$$I(\text{data}; \theta_{\text{downstream}}) = I(\text{data}; \theta_{\text{upstream}}).$$

Two caveats:

- **Lossless w.r.t. $\theta$, given the model.** If the model is
  misspecified, you're losslessly carrying the wrong information.
  Architecture doesn't fix that — posterior-predictive checks do.
- **Cut posteriors leak.** Level 1 is one-directional; it preserves
  forward information but discards the feedback information the SSA
  block would give the ODE block. Level 2 (Gibbs) closes the leak.
  DPI actually motivates climbing the three-level boundary ladder —
  each level recovers more of the joint information.

### Possible slide

> **Data Processing Inequality:** no processing step recovers
> information thrown away earlier.
>
> **Standard WCMs:** point estimate at every module boundary ⇒
> information destroyed at every boundary.
>
> **InferCell:** distributions at every boundary ⇒ information
> preserved by construction.

Sells the architecture in one breath: *the polyglot stack isn't just
inconvenient, it's mathematically lossy.*

## Slide 9 walkthrough: how inference delivers each payoff

Slide 9 ("What Becomes Possible") lists six items. Each is a different
operation on the posterior — useful to know which mechanism delivers
which payoff.

### 1. Probabilistic prediction

Posterior predictive distribution:
$$p(y_{\text{pred}} \mid \text{data}) = \int p(y_{\text{pred}} \mid \theta)\, p(\theta \mid \text{data})\, d\theta$$

Recipe: draw $\theta_i \sim p(\theta \mid \text{data})$, push each
through the simulator, collect outputs. Empirical spread = credible
interval. Combines epistemic (from posterior) and aleatoric (intrinsic
stochasticity) uncertainty — the MC4D backup slide already makes this
distinction visible.

### 2. Parameter calibration

Bayes' rule:
$$p(k_{tx} \mid \text{data}) \propto p(\text{data} \mid k_{tx})\, p(k_{tx})$$

The posterior **is** the calibration. Joint inference pays a dividend
when parameters are shared across modules ($k_{tx}$ in both TX/TL and
Metabolism) — both datasets constrain one posterior. Invisible to
isolated module fits.

### 3. Identifiability

Compare posterior to prior:
- Marginal narrower than prior → identifiable
- Marginal ≈ prior → data silent
- Diagonal ridge in joint → only a combination is identifiable
  (the $k_{tx}/\gamma$ case)

See also the FIM-based view below — identifiability is doubly enabled.

### 4. Model criticism

Three operations, all needing a likelihood:
- **Posterior predictive checks.** Simulate $y_{\text{rep},i}$ for each
  $\theta_i$ from posterior; compare distributional summaries (variance,
  bimodality, tails) to real data.
- **Bayes factors.**
  $p(\text{data} \mid M_1) / p(\text{data} \mid M_2)$ via marginal
  likelihood (nested sampling, bridge sampling, VI lower bound).
- **Cross-formalism agreement.** Same data, ODE TX/TL vs SSA TX/TL
  posteriors on shared params. Disagreement ⇒ at least one is
  misspecified. Genuinely novel — needs the multi-formalism architecture.

### 5. Experimental design (BOED)

Expected Information Gain:
$$\text{EIG}(d) = \mathbb{E}_{y \sim p(y \mid d)}\!\left[ \mathrm{KL}\!\left(p(\theta \mid y, d) \,\|\, p(\theta \mid \text{data})\right) \right] = I(\theta; y \mid d)$$

For marginal EIG (e.g. on $\gamma_{prot}$), substitute the marginal
parameter into the formula. Gradient-based optimisation over designs
needs AD through the inner posterior — only tractable on a unified
differentiable substrate.

### 6. Fundamental biology (mechanism discovery)

Two routes:
- **Nested-model trick.** Make the mechanism a parameter $\alpha$ with
  prior centred at 0. If posterior on $\alpha$ pushes far from 0, data
  demands the mechanism. If posterior ≈ prior, data is silent.
- **Bayes factor** between $M_{\text{with}}$ and $M_{\text{without}}$.

Mechanism is **inferred**, not assumed.

### Unifying pattern

| Bullet | Posterior operation |
|---|---|
| Probabilistic prediction | Push forward through model |
| Parameter calibration | Compute from data |
| Identifiability | Compare to prior |
| Model criticism | Compare replicates to data |
| Experimental design | Optimise EIG over hypotheticals |
| Mechanism discovery | Read marginal of structural parameter |

Six operations on one object (the posterior). Candidate slide:
*"Six operations on one object."*

## Identifiability is doubly enabled — gradient + posterior

Identifiability has two complementary views, and InferCell's substrate
gives you both.

### The Fisher information view (gradient-based)

$$\mathcal{I}_{ij}(\theta) = \mathbb{E}\!\left[\frac{\partial \log p(y \mid \theta)}{\partial \theta_i} \cdot \frac{\partial \log p(y \mid \theta)}{\partial \theta_j}\right]$$

The FIM is gradients of the log-likelihood, squared and averaged. AD
through SciMLSensitivity.jl computes it **exactly** at any $\theta$ —
no finite differences, no noise from adaptive ODE solvers, no failure
at language boundaries.

What this buys you:

- **Local identifiability without running NUTS.** Eigendecompose the FIM
  at the MAP. Tiny eigenvalues = sloppy directions; eigenvectors = the
  unidentifiable combinations. The $k_{tx}/\gamma$ ridge falls out
  analytically, cheaply.
- **Sensitivities** $\partial y / \partial \theta$ as a first-class
  object — Sethna-style sloppy-model analysis as a gradient call.
- **NUTS itself only works because of AD.** No gradients ⇒ gradient-free
  samplers (Metropolis, ABC) ⇒ catastrophic scaling in dimension. The
  posterior route to identifiability also rests on the same substrate.

### Complementarity with the posterior view

| | Gives you | Cost | Limitation |
|---|---|---|---|
| FIM via AD | Local identifiability + sloppy directions | One gradient evaluation | Linearisation; misses multimodality, global degeneracy |
| Full posterior | Global identifiability + correlations + tails | Full NUTS run | Convergence-dependent |

**Workflow:** FIM analysis as a free first pass; posterior for what FIM
can't see. Both rest on the same differentiable substrate.

### Connection to information framing

The score $\nabla_\theta \log p(y \mid \theta)$ is the local information
$y$ carries about $\theta$. The FIM is its variance. Cramér–Rao:
$\text{Var}(\hat\theta) \geq \mathcal{I}^{-1}$. So computing gradients
is *literally* computing information. AD is not just a sampler
accelerant — it's the substrate for information measurement, which
dovetails with the DPI / information-flow framing above: the same
architecture that *preserves* information across boundaries lets you
*measure* information about parameters at low cost.

### How this could change the talk

- **"Why Pure Julia?" slide.** AD bullet currently says "HMC/NUTS, VI,
  gradient-based BOED." Add: *"Fisher information & sensitivities —
  sloppy-model analysis as a gradient computation."* Ties substrate to
  a scientific payoff already prominent on slide 9.
- **Identifiability as a both-roots capability.** Currently presented
  as Root-1 only (probabilistic semantics). It's actually the cleanest
  example of a payoff needing probabilistic semantics **and**
  differentiable substrate — gives two views of the same answer for
  free. Strongest candidate to carry the dual-root message on one slide.

## Risks of over-pivoting

- "Information" is colloquially overloaded. Without anchoring to a
  distribution, some listeners hear "stuff we know."
- "Uncertainty propagation" is the standard UQ term — abandoning it
  cuts the talk off from a literature reviewers care about.
- **Recommended discipline:** keep "uncertainty" as the dominant noun
  for what's *carried*; use "information" for what's *gained or lost*
  across boundaries. Two complementary verbs on one noun.

## Open questions to revisit

1. **How quantitative to go?** Pure framing (one liner on the
   inference-graph slide) vs a small info-theory mini-slide (mutual
   information / Fisher / EIG). First is safer for a 20-min talk;
   second is bolder.
2. **Is BOED on the future-work slide?** If yes, "information flow"
   earns its keep by setting up EIG without a vocabulary change later.
3. **Audience appetite for Bialek-style information-in-biology?** If
   yes, cells-as-information-processors is rocket fuel. If no, it's a
   tangent.
4. **Build a TikZ diagram for the DPI slide?** Two parallel chains —
   one losing information at each boundary, one preserving it.
5. ~~**Add Fisher information to the "Why Pure Julia?" AD bullet?**~~
   **Applied** on slide 11 (engineering payoffs): "AD through solvers,
   models, and inference --- HMC/NUTS, gradient-based BOED, Fisher
   information for identifiability". Slide 10's Julia-stack list still
   doesn't mention identifiability — open whether to thread it in there
   too.
6. **Move identifiability to a dual-root showcase?** It's the cleanest
   example of a capability needing both probabilistic semantics and
   differentiable substrate. Slide 11 now advertises the AD/Fisher side;
   slide 9 advertises the posterior side. Verbal "two views, same
   answer" remark is the lightweight option; a dedicated slide is the
   bold option.

## Changes applied to the talk

- **Slide 10:** wording updated to clarify "serialization, not reference"
  for mixed audiences.
- **Slide 11:** End-to-end gradients bullet now lists Fisher information
  for identifiability (VI dropped — wasn't featured elsewhere).
- **Appendix:** two new slides — *What "Crossed by Serialization" Means*
  (reference-vs-serialization deep-dive with the four flavours of
  serialization in a polyglot WCM) and *What Dies at a Module Boundary*
  (the three-objects-die unification table connecting distributions,
  derivatives, and types to the information framing).

## Pending slide edits from the slide-by-slide review

Not yet applied; flag for next pass:

- **Slide 9 #2 (parameter calibration)**: rewrite to lock in the
  joint-inference angle — e.g. *"What's $k_{tx}$, when constrained by
  transcription and metabolism data jointly?"*
- **Slide 9 #4 (model criticism)**: reframe as cross-formalism
  comparison — e.g. *"Does ODE or SSA better explain the data?"*. PPC
  variance question moves to backup.
- **Slide 9 #5 (experimental design)**: comparative-experiments framing
  — e.g. *"Should we measure protein at 30 min, or perturb a gene
  knockdown?"*
- **Slide 11 #3 (hybrid inference)**: option B rewrite — *"Mix
  differentiable and stochastic dynamics in one inference pipeline ---
  the boundary between NUTS and SBI is a function call, not a file"*.
- **Slide 11 #5/#6 (shared infrastructure, cheap prototyping)**: minor
  tightening, optional.
