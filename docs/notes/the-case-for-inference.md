# The Case for Inference-First Whole-Cell Modelling

## The problem

Whole-cell models simulate a living cell from its molecular parts. One example is the Luthey-Schulten Lab's 4D Whole-Cell Model (MC4D), which wires together 493 genes across 4 coupled solvers with spatial resolution at 10nm voxels.

But MC4D, like every existing whole-cell model, is a forward simulator. It takes ~1000 kinetic parameters from the literature, fixes them, and propagates dynamics forward in time. The output is a single trajectory: one prediction, no error bars, no UQ.

This is the norm across the field. Whole-cell modelling has focused almost exclusively on the forward problem: *"Given these parameters, what does the cell do?"*. An important corollary question, and the one we focus on here, is *"How much should we trust what the model tells us?"*. We argue this is effectively a gap in scientific reasoning.

## The core argument: intellectual honesty in modelling

The case for inference rests on three facts:

**1. Every parameter has uncertainty. This is a fact about the state of biological knowledge, not a modelling choice.**

When a whole-cell model sets `k_transcription = 0.05 s⁻¹`, that value came from somewhere: e.g. a single-molecule experiment in *E. coli* at 37°C in rich media, or a bulk measurement in a different organism under different conditions. The true rate in the organism being modelled, under the specific conditions being simulated, is not known exactly. The uncertainty is real and irreducible at the current state of knowledge.

This applies across the board. Rate constants, binding affinities, Michaelis-Menten constants, degradation rates: every kinetic parameter in a whole-cell model carries measurement uncertainty, cross-organism extrapolation uncertainty, and condition-dependent variation. Fixing 1000 such parameters to point estimates is an implicit claim of perfect knowledge about all of them simultaneously.

**2. That uncertainty propagates through the model into predictions. This is also a fact.**

When uncertain parameters feed into coupled nonlinear dynamics, the uncertainty in inputs becomes uncertainty in outputs. A 2-fold uncertainty in a transcription rate becomes uncertainty in mRNA counts, which becomes uncertainty in protein levels, which becomes uncertainty in metabolic fluxes, which becomes uncertainty in growth rate. The propagation is nonlinear and non-obvious. Small parameter uncertainties can amplify or cancel depending on the network structure.

A forward simulation with fixed parameters produces a single trajectory. That trajectory is one sample from the (unknown) distribution of trajectories consistent with what we actually know about the parameters. Without quantifying the input uncertainty, we have no way of knowing whether the prediction is robust (insensitive to parameter choices) or fragile (a different plausible parameter set gives qualitatively different behaviour).

**3. A framework that cannot represent uncertainty is structurally incapable of answering the most important question: "How confident should we be in this prediction?"**

This is the architectural point. MC4D stores parameters in a `sim_properties` dictionary, a mapping from names to floats. There is no place to put a distribution or mechanism to propagate uncertainty. There is no way to ask "given that I'm uncertain about these 50 rate constants, how uncertain am I about the predicted doubling time?" The architecture does not support UQ. 

This means that when MC4D predicts a 110-minute doubling time, we cannot distinguish between two very different scientific situations: (a) the prediction is robust and would hold across any plausible parameter set, or (b) the prediction is fragile and depends sensitively on particular parameter choices. Both look the same: a single number with no error bar.

## What inference enables

An inference-first architecture does not require full Bayesian inference over all parameters simultaneously. At 500+ parameters, that is computationally dificult and may even be intractable, especially as we scale up to more complex cells. But the architecture must *support* uncertainty quantification at every level, even if it is applied selectively.

Concretely, this means:

- **Parameters carry distributions, not just values.** A parameter can be fixed (point estimate from literature), constrained (informative prior from experimental data), or free (to be inferred). The framework treats all three cases uniformly.

- **Uncertainty propagates across module boundaries.** When one module's output feeds into another, the uncertainty in the first module's parameters becomes uncertainty in the second module's inputs. The boundary protocol handles this explicitly.

- **The architecture has a place to put uncertainty.** Even if you fix 480 parameters and only infer 20, the framework can represent the fact that those 480 are approximations, and can quantify how the 20 inferred parameters interact with the 480 fixed ones via sensitivity analysis.

Beyond uncertainty quantification, inference enables capabilities that forward simulation cannot provide:

- **Model selection.** Formally compare competing biological hypotheses (e.g., constitutive vs. bursty transcription) using Bayes factors, rather than visual trajectory matching.

- **Identifiability analysis.** Discover which parameters the data can actually constrain and which it cannot. When a posterior looks like the prior, the data is uninformative about that parameter. This is scientifically valuable information that forward simulation never reveals.

- **Principled calibration.** Replace subjective hand-tuning ("we adjusted parameters until the growth rate matched") with reproducible, automated parameter estimation from data.

## What we are not claiming

We are not claiming that MC4D's parameter values are wrong, or that forward simulation is not useful. Forward simulation is essential for understanding emergent dynamics, testing whether a model can reproduce observed behaviour, and making predictions.

We are claiming that the *architecture* of existing whole-cell models has a blind spot: it cannot represent or reason about the uncertainty that necessarily exists in its inputs and therefore its outputs. InferCell is designed to fill that gap. It complements forward simulation with the inferential machinery that rigorous science demands.

We are also not claiming that full Bayesian inference over a complete whole-cell model is currently feasible. The argument is that the framework should *support* inference as a first-class operation, applied where it is most needed and most tractable, rather than treating it as an afterthought that the architecture cannot accommodate.
