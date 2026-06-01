# Speaker Script — Brisbane MACSYS, May 2026

**Talk:** InferCell: Inference-First Whole-Cell Modelling
**Slot:** ~20 min. Target ≈ 19 min, leaving ~1 min slack for the click-through.

How to read this: each slide entry has a one-line **take-home** (the message
the audience should leave with), a short list of **talking points** (the
sentences you'd actually say), an **off-slide** section for things that need
to be voiced because they are *not* printed on the slide, and a **transition**
sentence to bridge into the next slide.

Bold = key term or beat worth pausing on.

---

## Section 1: The Problem

### Slide 1 — Title  (≈ 30 s)

**Take-home:** who you are and what the talk is about, in two sentences.

**Talking points:**
- "Thanks — I'm Tom, from Melbourne."

- "Put out a slack message last month..."

- "Present it to the group, get some feedback and see what people think"

---

### Slide 2 — Whole-Cell Models: The State of the Art  (≈ 60 s)

- My representative example of the current SOTA is...

**Talking points:**
- "Reference point: the Luthey-Schulten lab's **4D model**, published in Cell this year."
- "493 genes, four coupled solvers — **RDME, CME, ODE metabolism, and Brownian dynamics** — running on GPU, simulating the full cell cycle."


"This is **genuinely impressive**. Remarkable feats of science & engineering

**Transition →** "But there's a particular *shape* to how this model is built, a particular gap that I want to draw your attention to."

---

### Slide 3 — WCMs are Forward Simulators  (≈ 45 s)

- WCMs are forward simulators
- Params - model - predictions
- Plot
- Uncertainty is a little deceptive. ee-lee-atoric.



---

### Slide 4 — But the parameters themselves are uncertain  (≈ 60 s)

- Uncertainty on the parameters themselves...Pretend (top panel), actaully (bottom panel). Ultimately taken from someone somewhere doing an experiment on somthing, + then doing a statistical analysis to give you a number. Params come from from the literature, often from different organisms, often from experiments measuring something adjacent. 
- uUncertainty propagates **nonlinearly**. e.g. uncertaity on rate parameter into mRNA, into protein, into flux, into growth — and each step compounds."
- Model has no place to store or propagate that uncertainty.** Parameters are floats. There's no slot for a distribution."

*when MC4D predicts a 110-minute doubling time, we cannot distinguish whether the prediction is robust across plausible parameters, or whether it depends fragilely on specific choices.* Say this aloud.



---

### Slide 5 — What Does That Uncertainty Look Like?  (≈ 60 s)

1. Just to give you an example of what this looks like...
2. Smallest possible model..
3. Talk through diagram..

"yes, this assumes independent uniform priors over the literature ranges and reads off the extremes — a real Bayesian propagation would give a tighter interval, but the order of magnitude survives."

**Off-slide — punchline:**
- "MC4D has roughly **two thousand parameters**, feeding coupled, nonlinear dynamics." (let this land)
- "I don't know what the prediction uncertainty actually is for that model. **And neither does anyone else** — because the architecture can't tell you."


---

### Slide 6 — No information flows between modules  (≈ 45 s)

1. That was just one equation one module. A related pioint is how information, hwo uncertainty propagates between modules
2. Talk through diagram...
3. Two things have been thrown away. **Uncertainty** — already gone, because we kept only point estimates. And **coupling** — the parameters that are jointly plausible in the coupled model may never even be sampled during fitting."

**Off-slide:**
- "This is structural — it's not a workflow problem you can fix by being more careful. The architecture doesn't have a place for joint information."



---

## Section 2: InferCell

### Slide 7 — Inverting the Question  (≈ 45 s)

1. "So here's the move I want to propose."
2. "Forward simulation asks: **given θ, what does the cell do?**"
- "Inference asks the inverse: **given data, what θ are consistent with it?**"
- "Same model. Same equations. Different direction through the graph."
- "And critically — the output of inference is a **distribution**, not a point. p of θ given data."



---

## Section 3: Architecture

### Slide 8 — InferCell: Design Principles  (≈ 60 s)

1. What we have been building to answer this question is InferCell (wroking title)
2. - "InferCell is a composable, **inference-first** whole-cell modelling framework in Julia. Four design principles."
- "**Inference-first.** We build in to the architechuter methods for UQ/ probabilistic prediction / inference.
- "**Modular.** Swap a module — say, an ODE TX/TL for an SSA TX/TL — and the inference layer still works without changes."
- "**Information flow as a design test.** At every module boundary, the question is: *how does inference cross this boundary?* If the answer is 'it doesn't,' that's a design bug." The third principle is the one that does most of the work — most of our design arguments reduce to "what does information flow look like here?"
- "**Pure Julia.** One language, one computational graph, one place where autodiff works." Not swapping between different languages, booting up different subprocesses

---

### Slide 9 — What Becomes Possible  (≈ 45 s)

**Take-home:** the scientific payoff is everything Bayesian inference gives you, now available at whole-cell scale.

**Talking points:**
- "I'm not going to read all of these." (gesture across slide)
- "Two clusters. **On the left**: things you do *with* a posterior — credible intervals on predictions, calibration to data, identifiability."
- "**On the right**: things you do *to* the model — posterior predictive checks, experimental design, hypothesis tests on mechanism."
- "**None of these are accessible from forward simulation.** They all need a likelihood and a posterior."


A. Probabilistic prediction
A1. Credible intervals on observables
A2. Robust vs fragile predictions
A3. Prediction risk for downstream use

B. Inference
B1. Formal parameter calibration to data
B2. Joint posteriors & parameter covariance
B3 (related to B2). Parameter bleed through between modules via shared parameters
B4. Principled propagatin of experimental uncertainty into model
B5. Bayesian updatign with new data


**One concrete question per bullet** (pick one when you point at each):

1. **Probabilistic predictions** — *"Is the predicted doubling time 110 ± 5 minutes, or could it plausibly be anywhere from 60 to 200?"*
2. **Inference from data** — *"What is the transcription rate of this gene in syn3A, given our measurements — and how does that estimate sharpen when the next dataset arrives?"*
3. **Identifiability** — *"Does my data actually constrain k_tx on its own, or only the ratio k_tx / γ_mRNA?"*
4. **Model criticism** — *"Does the ODE reproduce the variability in single-cell mRNA counts, or only the mean? Do the ODE and SSA posteriors agree on shared parameters — or does one disagree with the data?"*
5. **Experimental design** — *"Which experiment should I run next to most reduce my uncertainty about γ_protein?"* (your example)
6. **Biology discovery** — *"Does the data require enzyme feedback in metabolism, or can it be explained without that link?"*

**Off-slide:**
- If the room is more biological than statistical, lean on "what mechanism is consistent with the data?" (bullet 6). If it's more computational, lean on identifiability and EIG (bullets 3 and 5).

**Transition →** "Before I show you the architecture, one slide on why this has to be one language."

---

### Slide 10 — Why Pure Julia?  (≈ 60 s)

**Take-home:** standard WCMs are polyglot stacks; gradients and uncertainties don't compose across language boundaries.

**Talking points:**
- "Take MC4D again: **Python plus C++ plus CUDA plus LAMMPS**, glued together by IO."
- "Modules don't talk to each other — they **serialize and send files**. That's the boundary."
- "And critically: **gradients don't compose** across that boundary. You cannot autodiff through a subprocess call."
- "InferCell sits on a single Julia stack: **DifferentialEquations.jl** for solvers, **Turing.jl** for inference, **SciMLSensitivity.jl** for AD through the solvers, **Catalyst.jl** for reaction networks."
- "One language, one computational graph, one autodiff path top to bottom."

**Off-slide:**
- Acknowledge: "Julia is a choice with tradeoffs — I'll come back to that in the backup slides if anyone wants."
- Point at the right diagram briefly: "this is the kind of plumbing you don't want — and right now everyone in the field has it."

**Transition →** "Same idea — what becomes possible — but for the engineering."

---

### Slide 11 — What Becomes Possible (Engineering)  (≈ 45 s)

**Take-home:** a single compiled runtime unlocks performance, end-to-end gradients, and module substitution.

**Talking points:**
- "Symmetric to the science slide — but for the engineering side."
- "**End-to-end gradients**: HMC and NUTS through the whole model — variational inference, gradient-based experimental design."
- "**Hybrid inference**: NUTS and SBI blocks compose, with a boundary protocol that lives in the type system."
- "**Composability**: swap an ODE for an SSA, drop in an ML surrogate, no cross-language wrapping."

**Off-slide:**
- The point isn't speed for its own sake — it's that **inference at whole-cell scale only works if you can autodiff through the simulator.** Without that, you're stuck with gradient-free methods, which are orders of magnitude slower.

**Transition →** "But cells aren't all one kind of dynamics — and that complicates inference."

---

### Slide 12 — Hybrid Inference  (≈ 60 s)

**Take-home:** cells mix differentiable and stochastic dynamics; no single inference method handles both.

**Talking points:**
- "ODE metabolism is **smooth** — gradients exist — so you can use **NUTS or HMC**. Orders of magnitude more efficient than the alternatives."
- "Stochastic gene expression — Gillespie SSA, telegraph promoters — is **discrete and noisy**. There's no gradient path. You need **simulation-based inference** — ABC, neural posterior estimation — methods that treat the simulator as a black box."
- "**No single inference method handles both.** Real cells contain both kinds of dynamics, in the same model."
- "The hybrid graph uses each method **where it's strongest**."

**Off-slide:**
- This is the heart of the architectural argument — most WCMs duck this by picking one formalism. InferCell faces it head-on.

**Transition →** "Which leads to the central architectural picture."

---

### Slide 13 — The Inference Graph  (≈ 75 s)

**Take-home:** modules auto-partition into a DifferentiableBlock and a SimulationBlock; a Boundary Protocol stitches them.

**Talking points:**
- "Here's the architectural payoff."
- "**DifferentiableBlock** on the left — anything smooth lives here. TX/TL, metabolism. Inference is **NUTS or HMC**."
- "**SimulationBlock** on the right — anything stochastic. Gene expression by SSA, for example. Inference is **SBI**."
- "Critically: **modules are auto-partitioned into blocks by their declared inference mode.** A module says 'I'm differentiable' or 'I'm a simulator', and the framework routes it."
- "The two blocks are joined by a **Boundary Protocol** — a contract for how uncertainty propagates across them. That's the red arrow in the middle."
- "Data flows in from the top, posterior flows out the bottom — a joint posterior over the whole graph."

**Off-slide:**
- Pause on "Boundary Protocol" — name it explicitly. It's the most important phrase in the talk.
- "Different choices of boundary protocol give different cost/accuracy tradeoffs — I'll show you that in two slides."

**Transition →** "And what this lets you do — at the model level — is condition the entire thing on data."

---

### Slide 14 — Conditioning the Whole Model  (≈ 60 s)

**Take-home:** the standard workflow throws away uncertainty and coupling at every step; InferCell keeps both.

**Talking points:**
- "Left side — greyed out — the **standard approach**. Each subsystem fit in isolation, point estimates extracted, dropped into the whole-cell model. Prediction comes out with no uncertainty."
- "Two things lost. **Uncertainty discarded at every step.** **Coupling invisible during fitting.**"
- "Right side — the **InferCell approach**. *All* the observed data flows into a **coupled model**, with all modules tied together by their shared parameters. The output is a **joint posterior**, and predictions come with **credible intervals**."

**Off-slide:**
- Point at the red curved arrows between modules: "these are the couplings — they're now part of the inference problem, not assumed away."

**Transition →** "Now — joint inference over the whole model is expensive. So how much coupling do you actually need?"

---

### Slide 15 — How Much Coupling Do You Need?  (≈ 45 s)

**Take-home:** there's a spectrum from cut (cheap, approximate) to full PMCMC (exact, expensive). Pick the level the question deserves.

**Talking points:**
- "Two extremes are both bad. Full joint inference is intractable; isolated fits throw all the coupling away."
- "**Level 1: cut.** Downstream sees upstream's posterior, but there's no feedback. Fast. Approximate."
- "**Level 2: Gibbs-like.** Alternate between blocks. Information flows both ways. Still approximate, but tighter."
- "**Level 3: PMCMC.** Asymptotically correct, much more expensive."
- "**You choose the level that matches the question.** And — usefully — comparing levels gives you a diagnostic: if the posterior doesn't change when you allow feedback, the cut was good enough."

**Off-slide:**
- "This is the kind of decision you can only make *in this framework* — in a forward simulator, the choice doesn't even exist."

**Transition →** "OK. Let me show you what this looks like running on a small example."

---

## Section 4: In Practice

### Slide 16 — Toy example: a 3-module WCM  (≈ 75 s)

**Take-home:** three modules, two formalisms, nine unique parameters with shared linkages — small enough to validate, real enough to show the architecture working.

**Talking points:**
- "Three modules. **TX/TL as an ODE**, **Metabolism as an ODE**, and **Gene Expression as a stochastic SSA** — same biology as TX/TL, just a different formalism."
- "The colour-coded parameters show the **sharing structure**: red is shared across all three modules; green is shared between TX/TL and the SSA; blue is shared within the differentiable block."
- "After deduplication: **about nine unique parameters**, all with priors."
- "The first two modules sit in the **Differentiable Block**, sampled by NUTS. The SSA sits in the **Simulation Block**, sampled by ABC-SMC. The red arrow is the boundary protocol."
- "Crucially: **same biology, two formalisms** — that lets us **compare posteriors across formalisms** on the same data."
- "Validated on synthetic data where ground truth is known."

**Off-slide:**
- "This isn't a whole-cell model — it's a toy. The point is the *architecture*, not the biology. If we can't make the inference graph work on 9 parameters, we have no business scaling to 2000."

**Transition →** "Here's what the inference output actually looks like for the differentiable block."

---

### Slide 17 — Joint parameter estimation: differentiable block  (≈ 60 s)

**Take-home:** NUTS on the differentiable block gives a joint posterior; truth recovered, parameter correlations visible.

**Talking points:**
- "Corner plot. **Diagonals are the marginal posteriors** for each parameter — that's the density we have on that one parameter after seeing the data."
- "**Off-diagonals show pairwise joint posteriors** — they tell you which parameters trade off against each other."
- "The truth is recovered — the posteriors sit on the true values."
- "Where you see **tilted ellipses**, parameters are correlated — individually less identifiable, but jointly constrained."

**Off-slide:**
- "Point estimates would just be dots in the middle of these. Everything else — the width, the shape, the correlations — is gone in a forward simulator."
- If you have time: pick one specific off-diagonal and name what correlation it shows (e.g. transcription rate vs mRNA degradation often trade off because steady-state mRNA only constrains their ratio).

**Transition →** "Same exercise, but for the *stochastic* block — where we don't have a likelihood at all."

---

### Slide 18 — Joint parameter estimation: simulation block  (≈ 60 s)

**Take-home:** ABC-SMC recovers posteriors on the SSA module — without ever writing down a likelihood.

**Talking points:**
- "Same biology — transcription, translation, degradation — but modelled as a **stochastic jump process** instead of an ODE."
- "There's **no closed-form likelihood** for this system at the trajectory level. We can't write down p of data given theta."
- "But we *can* simulate from it. **ABC-SMC** uses simulation + summary statistics to construct a posterior anyway."

**Off-slide — important:**
- "Same biology, swapped formalism — and we get a posterior either way. **That's the architectural promise**: the inference layer doesn't care which formalism a module uses."
- Acknowledge that the SSA posterior is **broader than the ODE one** — that's expected: stochasticity costs you information, and ABC is approximate. The interesting question is *how much* broader, and that's only askable because we have both.

**Transition →** "Final result — what happens when we let the two blocks talk to each other."

---

### Slide 19 — Joint parameter estimation: Gibbs-like  (≈ 75 s)

**Take-home:** Level-2 boundary protocol — the blocks exchange information and the shared-parameter posteriors tighten.

**Talking points:**
- "This is **Level 2** — Gibbs-like alternation. We hold the SSA fixed, sample the differentiable block; then hold the differentiable block fixed, sample the SSA; iterate until the posteriors stop moving."
- "The result: the **shared parameters** — those red and green ones from the toy example — get **jointly informed** by both formalisms' data."
- "Compare with the single-pass cut posterior: the credible intervals on shared parameters **narrow substantially**."
- "This is **information flowing between modules** — the thing that was structurally impossible in the forward-simulator picture."

**Off-slide:**
- This is the slide that demonstrates the central architectural claim — *coupling between modules is now part of the inference problem.* Linger on it.
- The integration-test number from the appendix: Level-2 narrows shared-parameter 90% CIs by more than 50% versus single-pass. Say this if you want a quantitative anchor.

**Transition →** "Let me wrap up."

---

## Section 5: Summary & Outlook

### Slide 20 — Summary  (≈ 30 s)

**Take-home:** three beats — the problem, the inversion, the call to engage.

**Talking points:**
- "**The problem.** Whole-cell models take a thousand parameters from the literature and simulate forward. No UQ, no learning from data, no way to know which predictions to trust."
- "**InferCell inverts the question.** Given observations, infer parameters and their uncertainties. The model becomes something you **condition on data**, not just run forward."
- "**If you're interested — get involved.** Repo is on GitHub at MACSYS-CoE/InferCell. Come talk to me, or open an issue."

**Off-slide:**
- "Thank you. Happy to take questions."
- Have the URL ready to point at; if the room is small, give the verbal handle ("macsys hyphen coe slash infercell on GitHub").

**Transition →** Q&A.

---

## Speaker notes — global reminders

- **Slow down on the architectural words**: *inference-first*, *Boundary Protocol*, *differentiable block*, *simulation block*. These are the load-bearing terms.
- **The 24× punchline (slide 5)** is your "buy-in" moment for the architecture story. Don't rush it.
- **Slides 17–19 are the proof**: corner plots earn you the right to claim the architecture works. Point at the plots, don't read past them.
- If you go long, **drop slide 11** (engineering "what becomes possible") — it mirrors slide 9 and the audience already has the picture.
- If a question hits on **scaling / Julia-vs-JAX / model misspecification**, the backup slides cover all three.
