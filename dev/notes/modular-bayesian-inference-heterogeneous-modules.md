# Modular Bayesian Inference for Heterogeneous Coupled Systems

*A technical survey and recommended architecture for systems combining deterministic ODE modules, discrete stochastic (CTMC/Gillespie) modules, and ML surrogates.*

---

## 1. Problem statement

### 1.1 Setup

A system of $K$ modules $\mathcal{M}_1, \dots, \mathcal{M}_K$. Parameters partition as

$$\theta = (\theta_1, \dots, \theta_K, \varphi)$$

where $\theta_k$ is local to module $k$ and $\varphi$ collects quantities shared across two or more modules. Modules are coupled in three structurally distinct ways:

- **Feed-forward**: output of $\mathcal{M}_A$ is an input to $\mathcal{M}_B$. The coupling variable is a deterministic or stochastic function of $\mathcal{M}_A$'s state.
- **Feedback**: $\mathcal{M}_B \rightleftarrows \mathcal{M}_C$. The coupling is a fixed point, not a composition.
- **Parameter sharing**: $\varphi$ appears in the generative mechanism of both modules with no directional relationship at all.

Data arrives as streams $y_1, \dots, y_J$. The incidence between streams and modules is many-to-many: some $y_j$ informs exactly one module directly; some informs several; and every stream informs *every* module it is connected to through the coupling graph, indirectly.

The target is $p(\theta \mid y_{1:J})$, or some deliberate modification thereof.

### 1.2 Three problems that get conflated

Before choosing a method, decide which of these you are solving. They have different solutions and mixing them up is the main source of confusion in this literature.

| Problem | Statement | Method family |
|---|---|---|
| **Computational** | The joint posterior is the correct target; it is merely intractable at scale. | Message passing, compositional scores, divide-and-conquer SMC, blocked MCMC |
| **Epistemic** | Some modules are misspecified; you do *not* want full information flow. | Cut posteriors, semi-modular inference, generalized Bayes |
| **Constructive** | Submodels were built separately, with their own priors; there is no coherent joint yet. | Markov melding, chained melding, pooling |

A system with ODE, CTMC, and surrogate modules almost always has all three problems at once. The recommended architecture in §6 addresses them in layers.

### 1.3 The second axis: inference interfaces

The literature is written as if modules differ only in their *coupling*. In practice they also differ in what they can computationally provide, and this second axis dominates the engineering.

| Rung | Capability | Typical module |
|---|---|---|
| **L3** | Evaluate $\log p(y_k \mid \cdot)$ and $\nabla_\theta \log p(y_k \mid \cdot)$ | ODE + observation noise, with adjoint sensitivities |
| **L2** | Unbiased *estimate* of $p(y_k \mid \cdot)$; no gradient | Partially observed CTMC via particle filter |
| **L1** | Simulate, and differentiate the simulator | Gillespie via StochasticAD; CLE; tau-leaping |
| **L0** | Simulate only | Exact SSA, black-box legacy code |
| **L?** | Evaluate a learned object with ambiguous probabilistic semantics | ML surrogate |

**Design consequence.** Module boundaries should be redrawn, where possible, so that the tractability partition and the module partition coincide. If a tractability boundary runs through the interior of a module you will be unable to apply a single sampler to that module, and you will pay for it in mixing for the entire run.

---

## 2. Module type profiles

### 2.1 Deterministic ODE modules

**Structure.** $\dot{x} = f(x, \theta_k, u)$ with $u$ the interface input from upstream modules, and an observation model $y_k \sim g(\cdot \mid x(t_i), \sigma)$. The likelihood is available in closed form conditional on solving the IVP.

**Emits.** L3. Full log-density, gradients via forward sensitivities (cheap for few parameters) or the adjoint method (cheap for many), and Hessians if you want a Laplace approximation.

**Best local solver.** NUTS/HMC. For use as an EP site, a Laplace approximation on the tilted distribution is often adequate and dramatically cheaper.

**Pitfalls specific to the modular setting:**

- *Solver tolerance and gradient noise.* Adaptive solvers make $\log p$ a slightly discontinuous function of $\theta$. HMC is unusually sensitive to this — it manifests as divergences or as a step size that adaptation drives to zero. Tighten tolerances well beyond what you would use for forward simulation, and pin the solver's error control so gradients are deterministic given $\theta$.
- *Structural non-identifiability.* ODE modules are frequently non-identifiable from their own stream alone, and identified only through coupling. This is an argument *against* naive cutting: cutting feedback into an ODE module can render it unidentifiable, and the cut posterior will silently collapse onto the prior in some directions. Check identifiability of each module both in isolation and in situ.
- *Stiffness heterogeneity.* If module A is stiff and module B is not, a single global integrator wastes enormous effort. This is an independent argument for modular computation quite apart from statistics.

### 2.2 Discrete stochastic (CTMC / Gillespie) modules

**Structure.** A continuous-time Markov jump process on a discrete state space, with reaction propensities $a_r(x, \theta_k)$. If the jump chain is fully and exactly observed, the likelihood is available. It essentially never is. Under partial or noisy observation, the likelihood requires marginalising over all unobserved event sequences — intractable.

**Emits.** Depends entirely on which approximation you adopt. This is a ladder, and each rung trades bias for tractability:

| Approach | Emits | Bias | Notes |
|---|---|---|---|
| Matrix exponential / finite state projection (FSP) | L3 | Controllable, provable error bound | Only for small, bounded state spaces |
| Particle filter / pseudo-marginal | L2 | **None** — exact target | The gold standard; no gradients |
| Linear noise approximation (LNA) | L3 | Uncontrolled, worst at low copy number | Tractable Gaussian likelihood; differentiable |
| Moment closure | L3 | Uncontrolled; closure can be non-physical | Cheap; verify against SSA |
| Chemical Langevin equation (CLE) | L1/L3 | Fails when species counts are small | Reparameterisable; diffusion-bridge inference |
| Tau-leaping | L1 | Controlled by $\tau$; negativity issues | Between CLE and exact |
| Exact SSA + StochasticAD | L1 | **None** in the derivative estimator | Unbiased gradients through discrete randomness |
| Exact SSA + smoothed AD (DGA) | L1 | Biased by the smoothing | Cheap; bias is uncharacterised for complex nets |
| Exact SSA only | L0 | None | ABC / NPE / synthetic likelihood territory |

**The distinction that matters most.** StochasticAD (Arya, Schauer, Schäfer & Rackauckas, 2022) provides *unbiased* derivative estimates for programs with discrete randomness — the derivative estimator is unbiased for the derivative of the true expectation. The differentiable Gillespie algorithm (Mehta group, eLife 2025) instead replaces the discontinuous operations in the SSA — reaction index selection and species updates — with smooth surrogates (sigmoids, Gaussians), enabling ordinary backpropagation at the cost of a bias that is not characterised for complex reaction networks. If you need gradients for a *variational objective*, unbiasedness matters, because a biased gradient converges to the wrong stationary point rather than merely converging slowly.

**Bias propagation is a modularisation problem.** Every rung you descend introduces a misspecification in module $k$ that will propagate through the coupling into your well-specified modules. An LNA adopted for computational convenience is not a neutral numerical choice — it is a modelling error with system-wide consequences. This is the cleanest practical motivation for semi-modular inference: the influence parameter $\eta_k$ acquires the interpretation *how much do I trust the approximation I made in module $k$*.

### 2.3 ML surrogates

"ML surrogate" denotes at least four distinct objects, which compose in different ways. Getting this taxonomy wrong is, in my experience, the single most common source of silently incorrect modular inference.

#### (a) Forward-map emulator

Approximates the simulator's input→output map: $\hat{f}(u, \theta_k) \approx x$.

- **Composition**: must be *probabilistic*. Propagate the emulator's predictive variance into the likelihood, alongside a Kennedy–O'Hagan discrepancy term for simulator-vs-reality error.
- **For networks of emulated modules**: linked GP emulation (Kyzyurova, Berger & Wolpert 2018; Ming & Guillas 2021) gives closed-form linked predictive distributions for feed-forward systems and has been demonstrated on feedback-coupled systems. Linked deep GPs (Ming & Williamson 2023; `dgpsi`) relax the stationarity requirement by treating the network as a deep GP with partially exposed hidden layers, and outperform both linked GPs and a monolithic DGP fitted to the whole network.
- **Failure mode**: plugging in the emulator *mean* and treating it as the simulator. Produces posteriors that are both biased and overconfident, and the two errors do not announce themselves.
- **Bonus**: the adaptive-design machinery in the linked-GP literature tells you which module to spend simulation budget on — an underused capability.

#### (b) Likelihood surrogate

Approximates $\log p(y_k \mid \theta_k, u)$ directly: neural likelihood estimation, synthetic likelihood (Wood 2010), Bayesian synthetic likelihood.

- **Composition**: trivial. It *promotes an L0 module to L3*, giving you a differentiable log-density that can join a global HMC block.
- **This is usually the most useful surrogate type in a modular system**, precisely because it restores homogeneity of the interface.
- **Caveat**: synthetic likelihood assumes summary statistics are Gaussian; check. Neural likelihoods inherit whatever misspecification the simulator has.

#### (c) Posterior surrogate

An NPE flow or diffusion model returning $\hat{p}(\theta_k \mid y_k)$.

- **Composition trap**: this object *cannot be used directly as a module message*, because it embeds the prior it was trained under. Composing $Q$ such surrogates multiply-counts the prior $Q$ times. You must divide it out:

$$\ell_k(\theta_k) \;\propto\; \frac{\hat{p}(\theta_k \mid y_k)}{p_{\text{train}}(\theta_k)}$$

  which recovers a likelihood factor that composes correctly. This is exactly the prior-double-counting problem that the EP cavity solves structurally, and it recurs throughout the divide-and-conquer NPE literature (where the correction takes the form of dividing by $p(\theta)^{1-Q}$).
- The division is numerically hazardous in the tails, where $p_{\text{train}}$ is small and $\hat{p}$ is least reliable. Truncate, or train with a deliberately diffuse $p_{\text{train}}$, or use (b) instead.
- **NPLE** (jointly learning posterior and likelihood surrogates under a combined loss; Radev et al. 2023 and related) is attractive here because it hands you both objects and you can use whichever composition demands.

#### (d) Black-box point predictor with no UQ

- **Composition**: treat as a deterministic map, fit a discrepancy model on held-out data, and propagate that discrepancy's uncertainty. Be explicit that you now have an error source you have only empirically bounded.
- Strong candidate for $\eta_k < 1$ or an outright cut.

#### Surrogate risk, generally

Surrogates are the dominant misspecification source in a hybrid system, and neural posteriors in particular fail *silently*: NPE can be overconfident or misleading when the observed data falls outside the simulator's support, and robust variants (RNPE, explicit error models) are partial remedies at best. Per-module simulation-based calibration is mandatory, not optional.

### 2.4 Summary

| Module type | Emits | Local solver | Message via | Main risk |
|---|---|---|---|---|
| ODE | L3 | NUTS / Laplace | Analytic moments or short HMC on tilted | Non-identifiability under cutting; solver-induced gradient noise |
| CTMC, small state space | L3 (FSP) | NUTS | Analytic | State-space truncation error |
| CTMC, partially observed | L2 | PMMH | ABC/NPE on tilted (EP-ABC) | No gradients → blocking constraint |
| CTMC via LNA/CLE | L3 | NUTS | Analytic | Uncontrolled approximation bias propagating system-wide |
| CTMC via StochasticAD | L1 | Structured VI | Score/moment estimates | Gradient variance |
| Emulator (a) | L3 | NUTS | Analytic (GP) | Mean-plugging; unpropagated emulator variance |
| Likelihood surrogate (b) | L3 | NUTS | Analytic | Inherited simulator misspecification |
| Posterior surrogate (c) | — | — | Divide out $p_{\text{train}}$ first | Prior double-counting; tail instability |
| Point predictor (d) | L0-ish | — | Discrepancy model | Unquantified error |

---

## 3. Method families

### 3.1 Monolithic full-joint inference

Write the whole thing in one PPL and run NUTS. Always worth attempting on a reduced version of the system, because it is the only unambiguous ground truth you will get. It fails at scale for the obvious reasons, and fails immediately if any module is L0–L2, since HMC requires L3 throughout.

### 3.2 Message passing / Expectation Propagation

**The mechanics.** Factorise the posterior into sites, one per module:

$$p(\theta \mid y) \propto p(\theta) \prod_k f_k(\theta), \qquad f_k(\theta) = p(y_k \mid \theta)$$

Approximate each site by a tractable exponential-family factor, giving $q(\theta) \propto p(\theta)\prod_k \tilde{f}_k(\theta)$. Cycle over sites:

1. **Cavity**: $q_{-k} \propto q / \tilde{f}_k$ — everything the rest of the system currently believes.
2. **Tilted**: $q_{-k}(\theta) f_k(\theta)$ — the only step touching the true module.
3. **Project**: moment-match to the exponential family.
4. **Update**: $\tilde{f}_k \leftarrow \text{proj}(\cdot)/q_{-k}$; refresh $q$.

**Why this is the right shape for heterogeneous systems.** Step 2 requires only *moments of cavity × local likelihood*. That is a low bar, satisfiable by radically different technologies per module:

- ODE site → short NUTS or Laplace on the tilted density.
- Partially observed CTMC site → EP-ABC (Barthelmé & Chopin 2014, JASA). The critical practical point: the site update only needs the mean and covariance of a pseudo-posterior formed from the Gaussian cavity and the local likelihood, and *any* ABC or SBI method can supply it. The local problem is far easier than global ABC because the cavity is much more informative than the true prior, so you rarely propose low-likelihood parameters. **In your setting, the tractable ODE modules and their data act as an informative pseudo-prior for the intractable module — the very information flow you were worried about becomes a computational asset.**
- Surrogate site → if probabilistic, often analytic moments.

**Prior double-counting is solved structurally.** Because you divide out the current site before recomputing it, the prior enters exactly once regardless of $K$. This is the difference between EP and naive "fit A, feed A's posterior into B as a prior" pipelines, which are wrong by a factor of $p(\theta)^{K-1}$.

**EP vs VI.** Mean-field VI minimises $\mathrm{KL}(q\|p)$ globally: mode-seeking, systematically under-disperse. EP minimises $\mathrm{KL}(p\|q)$ *locally at each site*: mass-covering, generally better-calibrated variances. If you care about resolving power rather than point estimates, this asymmetry matters.

**Costs and failure modes.**
- No convergence guarantee. Oscillation is common; damping the updates (partial steps in natural-parameter space) usually fixes it.
- No monotone objective to monitor. You are flying without an ELBO.
- Multimodality in $\varphi$ is fatal to Gaussian sites.
- Requires shared quantities to live somewhere an exponential-family message makes sense. Discrete or hard-constrained $\varphi$ is awkward; transform where possible.
- Power EP / $\alpha$-divergence variants interpolate between EP ($\alpha=1$) and VI ($\alpha\to 0$) and are worth trying when EP is unstable.

**Learned message operators.** The frontier version: amortise the site update itself, learning an operator that maps incoming messages to the outgoing message (Heess, Tarlow & Winn 2013; Jitkrittum et al. 2015). For an expensive CTMC module this converts a costly inner ABC problem into a network evaluation, and the operator is reusable across datasets. See §6.2.

**Implementation.** `RxInfer.jl` / `ReactiveMP.jl` (Bagaev & de Vries) executes hybrid message passing — belief propagation, variational message passing, EP, EM — on Forney-style factor graphs with user-specified local form and factorisation constraints, and scales to state-space models with $10^5$+ variables. Vehtari et al., *Expectation propagation as a way of life* (JMLR 2020) is the definitive treatment of the distributed/modular framing.

### 3.3 Compositional score / diffusion methods

The posterior score is a **sum** of contributions, which means you can mix analytic and learned scores in one sampler:

$$\nabla_\theta \log p(\theta \mid y_{1:J}) = \nabla_\theta \log p(\theta) + \sum_j \nabla_\theta \log p(y_j \mid \theta)$$

Use autodiff for L3 modules; learn the score by neural score estimation for L0 modules. You pay the SBI tax only where you must.

- Foundations: Geffner, Papamakarios & Mnih, *Compositional score modeling for SBI* (ICML 2023).
- Diffusion posterior sampling in tall-data settings: Linhart et al. (TMLR), which replaces annealed Langevin with a Gaussian approximation of the backward kernels.
- Stability at scale: Arruda et al. (arXiv 2505.14429) introduce an error-damping estimator for aggregating large numbers of score contributions, validated to $10^5$ data points. Naive aggregation is unstable; use their correction.
- Time series: Gloeckler, Toyota, Fukumizu & Macke, *Compositional SBI for time series* (ICLR 2025).
- Arbitrary conditioning: Simformer / *All-in-one simulation-based inference* (Gloeckler, Deistler, Weilbach, Wood & Macke) learns the joint over $(\theta, y)$ and conditions on arbitrary subsets — useful when streams are intermittently missing.

**Where this beats EP**: no exponential-family constraint, so multimodal and heavy-tailed shared quantities are fine. **Where it loses**: less mature theory for the mixed analytic/learned case, and the aggregation error is a genuine open problem (see arXiv 2510.15817, 2605.21253 for error analyses).

### 3.4 Blocked Metropolis-within-Gibbs

The exact, boring, essential option.

- ODE block → NUTS.
- CTMC block → pseudo-marginal (Andrieu & Roberts 2009) with a particle-filter likelihood estimate; PMMH (Andrieu, Doucet & Holenstein 2010); for stochastic kinetic models specifically, Golightly & Wilkinson.
- Targets the true joint posterior exactly. Use it as ground truth for a reduced system.

**The hard constraint**: you cannot run HMC on a pseudo-marginal block. The tractability partition therefore *becomes* your blocking, which is why §1.3's advice about aligning module and tractability boundaries matters so much.

**The pain point**: shared parameters $\varphi$ straddling a tractable and an intractable module mix badly. Remedies: non-centred parameterisation; the copies-and-tie device (§3.8); delayed-acceptance schemes using a cheap surrogate as a first-stage filter.

### 3.5 Markov melding (construction)

For submodels built separately with their own priors on shared quantities. Combines $p_m(\psi_m, Y_m \mid \varphi)$ with a *pooled* prior $p_{\text{pool}}(\varphi)$ (linear or logarithmic pooling) to form a coherent joint (Goudie, Presanis, Lunn & De Angelis, JASA 2019).

- **Chained melding** (Manderson & Goudie, *Bayesian Analysis*) handles chains where adjacent submodels share *different* quantities — $\mathcal{M}_1 \cap \mathcal{M}_2 = \varphi_{1\cap2}$, $\mathcal{M}_2 \cap \mathcal{M}_3 = \varphi_{2\cap3}$ — which is the realistic topology.
- **Computation**: multi-stage MCMC using submodel posteriors as proposals; more recently a divide-and-conquer SMC sampler exploiting the tree structure (Liu & Goudie, arXiv 2605.22301), which avoids sampling the full model directly.
- **Landmine**: the marginal prior on $\varphi$ implied by each submodel is usually *implicit* and must be estimated; error in that density estimate destabilises the two-stage sampler, badly in the tails. Manderson & Goudie's weighted-sample estimator of the prior marginal self-density ratios is the fix.

### 3.6 Cut posteriors and semi-modular inference (epistemic)

**Cut**: fix $p(\varphi \mid y_{\text{trusted}})$ from the reliable module; condition the suspect module on it. Characterisation (Yu, Nott & Smith 2023): the cut posterior is the KL-closest approximation to the conventional posterior subject to the constraint that the $\varphi$-marginal equals the trusted one.

- Not the stationary distribution of naive Gibbs (Plummer 2015). Requires dedicated algorithms.
- **The key recent reference for many-module systems**: Liu & Goudie, *A general framework for cutting feedback within modularised Bayesian inference*, JRSS-B 87:1171–1199 (2025). Prior work was almost entirely two-module; this gives a formal definition of "module", methods for identifying modules and determining their order, and construction of the cut distribution on an arbitrary DAG.

**Semi-modular inference (SMI)**: soften the cut with $\eta \in [0,1]$ interpolating between cut and full Bayes (Carmona & Nicholls 2020). Framed as a bias–variance trade-off; formalised by Frazier & Nott (2025); connected to generalized Bayes by Nicholls et al. (2022). $\eta$-selection via amortised variational meta-posteriors (Carmona & Nicholls 2022; Battaglia et al. 2025).

**Generalized Bayes version**: Frazier & Nott, *Cutting Feedback and Modularized Analyses in Generalized Bayesian Inference*, Bayesian Analysis 20(4):1647–1675 (2025). Relevant when some streams are better handled by a targeted loss than a likelihood — and loss scaling is the substantive difficulty there.

**Computation for cut/SMI**: stochastic-approximation adaptive MCMC (SACut; Liu & Goudie); posterior bootstrap (Pompe & Jacob 2021); variational (Yu, Nott & Smith 2023); SMC (Mathews et al. 2025); emulation to increase imputations (Hutchings et al. 2025); **likelihood-free cut** via Gaussian mixture approximations to the joint of parameters and summaries (Chakraborty, Nott, Drovandi, Frazier & Sisson, *Stat. Comput.* 2023) — this last is the one that matches simulator-heavy systems.

**Warning worth taking seriously.** A cut posterior is not a posterior. There is no joint distribution of which it is a conditional; its predictive and decision-theoretic semantics are awkward; and if the "trusted" module is also wrong you have locked in the bias while removing the evidence that would have revealed it. Cut is a targeted robustness intervention, never a default.

**A consequence for model comparison.** Because there is no joint distribution, there is no marginal likelihood either: Bayes factors are undefined under a cut, and approximately so for $\eta < 1$. Any architecture that puts an SMI robustness layer on top has thereby given up evidence-based model selection and must compare structures predictively instead. See [model-uncertainty-and-selection.md](model-uncertainty-and-selection.md).

### 3.7 Linked emulator networks

When modules are expensive simulators, emulate each and link the emulators rather than emulating the composite. Covered in §2.3(a). The adaptive design component — allocating simulation runs to the module contributing most to system-level uncertainty — is the most valuable and least used part.

### 3.8 Copies-and-tie: a unifying device

Give each module its own copy $\varphi^{(k)}$ of each shared quantity, and tie the copies with a coupling prior of strength $\lambda$:

$$p(\varphi^{(1)}, \dots, \varphi^{(K)}) \propto p_0(\bar{\varphi}) \prod_k \exp\!\left(-\lambda \, d(\varphi^{(k)}, \bar\varphi)\right)$$

- $\lambda \to \infty$: recovers full Bayes.
- $\lambda \to 0$: independent modules.
- Asymmetric $\lambda_k$: something cut-like, with $\lambda_k$ playing a role analogous to $\eta_k$.

This turns modularisation from a binary decision into a continuous, checkable modelling choice — and, importantly, a **differentiable** one, so $\lambda$ can be optimised against a predictive criterion rather than chosen by argument. It also decouples the modules for Gibbs sampling (conditional on $\bar\varphi$, modules are independent), which connects to split-and-augment Gibbs schemes.

---

## 4. Decision guide

| Situation | Recommended approach |
|---|---|
| All modules L3, moderate dimension | Monolithic NUTS; don't over-engineer |
| Mixed L3/L2, shared params low-dimensional | Blocked Metropolis-within-Gibbs + PMMH |
| Mixed rungs, many modules, need calibrated variances | **EP-style message passing with heterogeneous site solvers** |
| Mixed rungs, multimodal or heavy-tailed $\varphi$ | Compositional score / diffusion |
| Everything differentiable (incl. via StochasticAD) | Structured VI with per-module factors |
| Submodels built separately, incompatible priors | Markov melding → then any of the above |
| Detected conflict between modules | SMI with $\eta$ chosen by ELPD; hard cut only if conflict is severe |
| Modules are expensive simulators | Linked (deep) GP network + adaptive design |
| A module is a black-box predictor with no UQ | Discrepancy model + $\eta < 1$, or cut |

---

## 5. Recommended architecture

**Summary of the recommendation**: an EP-style message-passing backbone in which the module interface is a distribution over shared quantities, with heterogeneous per-module solvers, and with the stochastic and surrogate modules amortised at the interface. A robustness layer ($\eta$ per module) sits on top, and an exact blocked-MCMC reference implementation validates it on a reduced system.

### 5.1 The interface contract

Define one interface and make every module satisfy it. In Julia this is natural via multiple dispatch:

```julia
abstract type Module end

# --- Required of every module ---
interface_in(m::Module)   # which shared/coupling variables it consumes
interface_out(m::Module)  # which it produces
simulate(m::Module, φ, θ_local, rng)

# --- Optional; presence determines the solver route ---
logdensity(m::Module, φ, θ_local, y)            # L3
∇logdensity(m::Module, φ, θ_local, y)           # L3
loglik_estimate(m::Module, φ, θ_local, y, rng)  # L2, unbiased
message(m::Module, cavity::Distribution, y)     # the EP site update
```

Every module implements `message`. *How* it implements it is private:

| Module | `message` implementation |
|---|---|
| ODE | Laplace, or short NUTS, on cavity × likelihood |
| CTMC (partially observed) | EP-ABC / NPE against the tilted target |
| CTMC (small state space) | FSP likelihood, analytic moments |
| Emulator | Analytic GP moments |
| Likelihood surrogate | Analytic or short NUTS |
| Posterior surrogate | Divide out $p_{\text{train}}$, then multiply by cavity |

**One further field, cheap now and expensive later: normalisation semantics.** Each module should declare whether the density it reports is normalised, correct up to a model-independent constant, or incomparable across model structures. ABC and posterior-surrogate routes are `incomparable` — the ABC pseudo-evidence depends on $\varepsilon$ and on the summaries, and an NPE object carries its training prior. Parameter inference never needs this field; model comparison needs nothing else, and by then the backends have been swapped several times and nobody remembers which returned what.

This is the single most important structural decision. Once the interface is a message rather than a likelihood call, you can swap any module's backend without touching the rest of the system — including swapping a slow exact solver for a fast amortised one, or an LNA for exact SSA, and directly measuring what that substitution costs you.

### 5.2 Amortise at the interface, not at the global posterior

The economics are asymmetric:

- ODE modules are cheap per evaluation; amortisation buys almost nothing and costs training.
- CTMC and surrogate modules are expensive; amortisation is transformative.
- The global posterior is the *wrong* amortisation target — it must be retrained whenever any module, prior, or stream changes.

So: **train, once and offline, a local message/likelihood approximator for each expensive module, conditioned on its data and on its incoming interface variables.** Concretely, for CTMC module $k$, learn either

- $q_\eta(y_k \mid \varphi_{\text{in}}, \theta_k)$ — a neural likelihood, promoting the module to L3; or
- $\mathcal{O}_k: \text{cavity} \mapsto \text{outgoing message}$ — a learned EP message operator.

The first is simpler and gets you into a standard gradient-based sampler. The second is the more ambitious and, I think, more interesting option: it makes each intractable module a fast, reusable inference component, and it is where a genuine methodological contribution is available (§7).

The interface variables are typically low-dimensional even when the module's internal state is not, which is what makes the amortisation tractable. This is the same insight that makes linked emulation work — the network structure exposes low-dimensional bottlenecks that a monolithic surrogate would have to learn around.

### 5.3 Build order

**Stage 0 — Instrument.**
- Draw the module DAG. Classify every shared quantity as feed-forward / feedback / shared-parameter, and every module by tractability rung. Redraw boundaries to align the two partitions.
- Conflict diagnostics *before* modularising: node-splitting statistics, prior-data conflict measures (Presanis et al.; Gåsemyr & Natvig), module-wise posterior predictive checks. Modularise because you detected conflict, not because it's convenient.
- Influence audit: for each stream $y_j$ and each module, measure how much $y_j$ moves that module's parameters (KL with vs without, or sensitivity of posterior moments). This tells you which couplings are actually load-bearing.

**Stage 1 — Ground truth on a reduced system.**
- Two modules, one ODE and one CTMC, sharing one parameter. Small enough for exact blocked MCMC (NUTS + PMMH).
- This is your reference. Every approximation in Stage 2 is validated against it. Without this you will not know whether an EP fixed point is a good approximation or a converged wrong answer.

**Stage 2 — Production message passing.**
- EP backbone over the module graph, heterogeneous site solvers per §5.1.
- Damping from the start. Monitor site-parameter trajectories for oscillation.
- Cross-check against Stage 1 on the reduced system; check against a structured-VI run on the full system (agreement between two different approximations is weak evidence, but disagreement is strong evidence of a problem).

**Stage 3 — Robustness layer.**
- Per-module $\eta_k$ (or coupling strength $\lambda_k$ in the copies-and-tie parameterisation), selected by ELPD on held-out data from streams *not* used to fit that module.
- Expect $\eta_k < 1$ for surrogate modules and for any module using an uncontrolled approximation (LNA, moment closure).

**Stage 4 — Validation.**
- Simulation-based calibration *per module*, then for the full system. SBC on a message-passing system is underexplored and you will likely have to build it: the natural version checks the calibration of each site's outgoing message against a simulated ground truth.
- Coverage checks on the shared quantities specifically — these are where modular approximations degrade first.
- Leave-one-stream-out: refit without $y_j$ and check the predictive. Streams whose removal barely changes anything were never informing the system, and streams whose removal changes everything are load-bearing single points of failure.
- Surrogate error budget: is the emulator/NPE uncertainty a negligible or dominant fraction of the posterior width? If dominant, the whole exercise is measuring your surrogate, not your system.

### 5.4 Failure modes to watch

| Symptom | Likely cause |
|---|---|
| EP oscillates, doesn't converge | Insufficient damping; multimodal $\varphi$; site family too restrictive |
| Posterior variance implausibly small | Prior double-counting (a posterior surrogate used without dividing out $p_{\text{train}}$); emulator mean plugged in |
| A module's posterior collapses to its prior | Cutting removed the feedback it depended on for identifiability |
| Results change a lot with solver tolerance | Gradient noise from adaptive ODE solvers |
| Good fit to each stream, bad joint predictive | Modules absorbing each other's misspecification via shared parameters — the exact case cut/SMI exists for |
| Full Bayes worse than a cut on held-out data | Genuine module misspecification; take the SMI route seriously |

### 5.5 Julia stack

- `DifferentialEquations.jl` + `SciMLSensitivity.jl` — ODE modules, adjoints
- `JumpProcesses.jl` / `Catalyst.jl` — CTMC modules
- `StochasticAD.jl` — unbiased gradients through discrete randomness
- `RxInfer.jl` / `ReactiveMP.jl` — message-passing scaffolding
- `Turing.jl` — site solvers, and the PMMH reference implementation
- `Flux.jl` / `Lux.jl` + a normalising-flow package — amortised message operators / neural likelihoods

---

## 6. Where the open problems are

Flagging these because several sit at the intersection of this survey and existing work, and look like publishable contributions rather than engineering:

1. **Amortised message operators for CTMC modules.** Learned EP message operators exist in the ML literature (Heess et al.; Jitkrittum et al.) but have not, as far as I know, been applied to stochastic biochemical modules inside a scientific multi-module system. The combination of an informative cavity (which makes the local SBI problem easy) with amortisation (which makes it fast) is not obviously explored.

2. **Structured VI over heterogeneous modules with unbiased gradients.** StochasticAD makes the discrete modules differentiable without bias. A per-module factorised variational family, optimised end-to-end with unbiased gradients through the SSA modules, is a natural object and connects directly to existing StochasticAD/VI work.

3. **$\eta$-selection when the untrusted module is a learned surrogate.** SMI's influence parameter has been studied for classically misspecified modules. Its behaviour when the misspecification is *surrogate approximation error* — which is estimable, unlike generic misspecification — should admit sharper theory. There is a natural link to the surrogate's own uncertainty quantification.

4. **Aggregation error for mixed analytic/learned scores.** The error analyses for compositional score methods assume all scores are learned. The hybrid case, where some contributions are exact, should have better rates, and the optimal allocation of simulation budget across the learned components is an open design question.

5. **SBC for message-passing systems.** No established protocol for calibrating a modular approximation module-by-module. This is a methodological gap with immediate practical value.

6. **Evidence and predictive comparison under modular approximation.** Every method family in §3 is designed for the parameter problem and is careless with the normalising constant, which is the one quantity model comparison needs. Concretely: (a) how badly does the EP energy approximate the log evidence when sites are ABC- or surrogate-backed rather than analytic; (b) what is the correct leave-one-stream-out predictive when the module supplying the stream is itself the object under comparison; (c) is there a coherent comparison criterion under SMI, where the marginal likelihood does not exist — stacking on held-out predictive density is the obvious candidate but its behaviour with $\eta$ as a nuisance is unstudied. Same shape of gap as (5), and arguably more consequential.

7. **Resolving power of the module graph.** Given the coupling structure and the streams, what is actually identified about each module? This is the partial-identification question applied to modular systems, and it would tell you *a priori* which modularisations are safe — rather than discovering it after a cut collapses a posterior onto its prior.

---

## 7. Annotated references

**Modular inference, cutting feedback, SMI**
- Liu, Bayarri & Berger (2009), *Modularization in Bayesian analysis*. The origin of the framing.
- Plummer (2015), *Cuts in Bayesian graphical models*, Stat. Comput. Why naive Gibbs fails.
- Jacob, Murray, Holmes & Robert (2017), *Better together? Statistical learning in models made of modules*.
- Carmona & Nicholls (2020, 2022). Semi-modular inference; variational meta-posteriors for $\eta$.
- Yu, Nott & Smith (2023), *Variational inference for cutting feedback in misspecified models*, Statist. Sci.
- Chakraborty, Nott, Drovandi, Frazier & Sisson (2023), *Modularized Bayesian analyses and cutting feedback in likelihood-free inference*, Stat. Comput. **The likelihood-free entry point.**
- Liu & Goudie (2025), JRSS-B 87:1171–1199. **The general multi-module DAG theory. Start here.**
- Frazier & Nott (2025), Bayesian Analysis 20(4):1647–1675. Generalized-Bayes version.
- Pompe & Jacob (2021), asymptotics and posterior bootstrap.

**Melding**
- Poole & Raftery (2000), Bayesian melding for deterministic simulation models.
- Goudie, Presanis, Lunn & De Angelis (2019), JASA. Markov melding.
- Manderson & Goudie, *Combining chains of Bayesian models with Markov melding*, Bayesian Analysis. **The realistic chain topology.**
- Manderson & Goudie (2022), numerically stable two-stage algorithm. Read before implementing.
- Liu & Goudie (arXiv 2605.22301), divide-and-conquer SMC for chained melding.

**Message passing**
- Minka (2001), EP thesis. The mechanics.
- Barthelmé & Chopin (2014), *Expectation propagation for likelihood-free inference*, JASA. **EP-ABC.**
- Vehtari, Gelman, Sivula et al. (2020), *Expectation propagation as a way of life*, JMLR. **The modular framing.**
- Heess, Tarlow & Winn (2013); Jitkrittum et al. (2015). Learned message operators.
- Bagaev & de Vries, *Reactive message passing for scalable Bayesian inference* / RxInfer.jl.

**Compositional / amortised SBI**
- Geffner, Papamakarios & Mnih (2023), ICML. Compositional score modeling.
- Linhart et al., TMLR. Diffusion posterior sampling, tall data.
- Arruda et al. (arXiv 2505.14429). Error-damped compositional aggregation at scale.
- Gloeckler, Toyota, Fukumizu & Macke (ICLR 2025). Compositional SBI for time series.
- Gloeckler, Deistler, Weilbach, Wood & Macke, *All-in-one simulation-based inference* (Simformer).
- Radev et al. (2023). Joint posterior + likelihood surrogates.
- Ward et al. (2022). Robust NPE under misspecification.

**Stochastic kinetics**
- Arya, Schauer, Schäfer & Rackauckas (2022), NeurIPS. StochasticAD — unbiased derivatives through discrete randomness.
- *A differentiable Gillespie algorithm*, eLife (2025) / arXiv 2407.04865. Smoothed, biased, cheap.
- Golightly & Wilkinson. Particle MCMC for stochastic kinetic models.
- Andrieu & Roberts (2009); Andrieu, Doucet & Holenstein (2010). Pseudo-marginal, PMCMC.
- Schnoerr, Sanguinetti & Grima (2017). Review of approximation methods for stochastic chemical kinetics — the LNA/moment-closure landscape.

**Emulation**
- Kennedy & O'Hagan (2001). Calibration with discrepancy.
- Kyzyurova, Berger & Wolpert (2018), SIAM/ASA JUQ. Linked GPs.
- Ming & Guillas (2021), SIAM/ASA JUQ. Matérn linked GPs, adaptive design.
- Ming & Williamson (arXiv 2306.01212). Linked deep GPs for model networks; `dgpsi`.

**Foundations**
- Bissiri, Holmes & Walker (2016), general Bayesian updating. The justification for loss-based module updates.
- Miller & Dunson (2019), coarsened posteriors.
- Talts, Betancourt, Simpson, Vehtari & Gelman (2018), simulation-based calibration.
