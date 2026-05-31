# Spatial Dynamics and the RDME Question

## The question

The Luthey-Schulten Lab's 4D Whole-Cell Model (MC4D) is built on a 3D Reaction-Diffusion Master Equation (RDME) running on a lattice via Lattice Microbes. The cell is discretized into ~10nm voxels. Every molecule has a position. Every reaction happens locally within a voxel. Diffusion moves molecules between voxels. The four solvers (RDME, ODE metabolism, Brownian dynamics for the chromosome, CME for gene expression) all hook into this spatial framework.

Space is the organizing principle. Everything else is subordinate to it.

InferCell has no RDME. The `formalism()` method supports `:ode`, `:sde`, `:jump`, and `:ssa`, but not `:rdme`. Spatial dynamics is explicitly listed as deferred scope for v1. Is this a gap?

## Different organizing principles for different questions

In MC4D, the RDME lattice is not a module. It is the world. Space is the architecture, and every solver hooks into it.

In InferCell, inference is the architecture. Modules are grouped not by where they are in the cell, but by how you do inference on them: differentiable (NUTS) vs simulation (SBI). The coupling mechanism is not diffusion between voxels but state variables and shared parameters flowing through the inference graph.

This is a deliberate design choice. MC4D asks "what does the cell look like in 4D?" InferCell asks "How much should we trust what the model tells us?". You do not necessarily need a spatial lattice to answer that question.

## Can RDME be added later?

Yes, at two levels of ambition.

**As a module (straightforward, and probably sufficient).** An RDME sub-model with `formalism() = :rdme` that declares its states (species x voxels), parameters (diffusion coefficients, reaction rates), and enters the simulation block because RDME is not differentiable. Inference on its parameters uses SBI. The sub-model interface already supports this conceptually. The engineering work is adding `:rdme` as a formalism, handling the state space expansion, and providing a solver backend (wrapping Lattice Microbes or using a Julia RDME solver).

**As a foundational layer (deep redesign, probably unnecessary).** Making space foundational the way MC4D does, where all modules exist on the lattice and diffusion couples everything. This would fight the inference-first design. It is hard to see when this would be the right choice for an inference-focused framework, but we should think carefully about this as the explicit spatial structure may prove important. 

## What foundational space gives you that modular space does not

**Emergent spatial organisation across processes.** In MC4D, ribosomes cluster near active transcription sites. This pattern emerges because mRNA is produced at the gene locus, diffuses slowly, and ribosomes bind it before it spreads far. No single module "contains" this phenomenon. It arises because transcription, diffusion, and translation share the same spatial grid and interact through local concentrations.

In a modular approach, if transcription and translation are separate modules coupled through a shared state variable (mRNA concentration), that coupling is a single number with no position. The ribosome clustering pattern cannot emerge because there is no shared space for it to emerge in.

**Physically realistic coupling.** Molecules do not exchange state variables. They diffuse through cytoplasm and encounter each other stochastically. The encounter rate depends on local concentration, which depends on where molecules are produced, how fast they diffuse, and where they are degraded. Foundational space makes this coupling automatic and physical. Modular coupling through shared variables implicitly assumes well-mixed at the interface, even if individual modules are spatial internally.

**Crowding.** The cytoplasm is ~30% macromolecules by volume. Local crowding affects diffusion and reaction rates. This is a global spatial property that depends on what every module is putting into the same physical space. No single module can capture it.

## The modularity argument

In MC4D, if you want to model transcription, you must put it on the lattice even if spatial effects are negligible for the question you are asking. Everything pays the cost of spatial resolution whether it needs to or not.

In InferCell, you can model transcription as a well-mixed ODE, get useful inference results, and later ask: does adding spatial resolution change the posteriors? If it does, you swap in an RDME version of that module. If it does not, you have saved enormous computational cost. The modularity lets you add spatial complexity where the data justifies it, rather than imposing it everywhere by default.

## What each approach can and cannot answer

**Well-mixed coupling is sufficient for:**

- Inferring kinetic rate constants from bulk measurements (mRNA counts, protein levels, growth rate).
- Quantifying parameter uncertainty in ODE/SDE models.
- Comparing competing model structures via Bayes factors (e.g. constitutive vs bursty transcription).
- Propagating uncertainty across module boundaries.

These are InferCell's v1 target questions. Spatial structure is irrelevant to them.

**Spatial modules (without foundational space) are sufficient for:**

- Inferring diffusion coefficients from single-particle tracking or FRAP data, using an internally spatial module.
- Testing whether spatial resolution improves the fit for a specific process (model selection: well-mixed vs RDME for one module).
- Modelling processes that are inherently spatial within a single module (e.g. Min protein oscillations).

**Foundational space is required for:**

- Emergent cross-module spatial phenomena: ribosome clustering near the nucleoid, co-localisation of coupled enzymes, spatial coordination of replication and division.
- Questions about how chromosome position affects transcription rate across the cell.
- Global effects like macromolecular crowding that depend on the spatial contributions of all processes simultaneously.
- Predicting cell geometry and its effect on multiple processes.

The modular approach cannot answer this third category. These are questions about spatial phenomena that emerge from the interaction of multiple processes in shared physical space.

## The well-mixed assumption as a testable hypothesis

This reframing is worth emphasizing. In most whole-cell models, the well-mixed assumption is treated as either a known limitation (acknowledged and accepted) or ignored entirely (implicit in the ODE formulation). In an inference-first framework, it becomes a testable hypothesis.

If you have spatial data (e.g. fluorescence microscopy with subcellular resolution), you can build two versions of a module: one well-mixed, one spatially resolved. Run inference on both. Compare them via Bayes factors. If the spatial model is better supported by the data, the spatial structure matters for that process. If not, the well-mixed model is the appropriate level of description, and adding spatial resolution would introduce unidentifiable parameters without improving the fit.

This is model selection applied to the question of spatial resolution, which is a capability that a forward-simulation-only framework cannot provide.

## Where spatial resolution matters and where it does not

Spatial resolution is not uniformly important across cellular processes.

**Genuinely spatial biology.** Chromosome organisation and gene expression depend on position: gene proximity to the origin of replication affects copy number, nucleoid structure affects accessibility to transcription machinery. Cell division is inherently spatial: Z-ring formation at midcell, Min protein oscillations pole-to-pole, septum placement. Membrane processes (signal transduction, transport, secretion) happen at specific locations. Macromolecular crowding and phase separation create local concentration effects that cannot be captured by a well-mixed model.

**Effectively well-mixed biology.** Small-molecule metabolism in a bacterial cell is fast relative to diffusion. JCVI-syn3A is ~400nm in diameter. Diffusion time for a small metabolite across the entire cell is on the order of milliseconds. Concentration gradients for ATP, amino acids, and nucleotides are negligible. Bulk protein and mRNA counts (how many copies, not where they are) do not require spatial resolution. Steady-state behaviour that is insensitive to spatial organisation gains nothing from voxel-level simulation.

MC4D's spatial resolution has revealed real biology (ribosome localisation patterns, chromosome dynamics). But the RDME is also the most computationally expensive component of the model, requiring GPU acceleration. For many submodels (metabolism, bulk gene expression), it is an open question whether spatial resolution changes the scientific conclusions relative to a well-mixed treatment. MC4D cannot formally test this because it has no model selection framework.

## Spatial resolution in other whole-cell models

MC4D is the exception, not the rule. Almost every whole-cell model ever built is well-mixed.

**Karr et al. 2012 (Covert Lab, Stanford)** is the landmark whole-cell model, for *Mycoplasma genitalium*. 28 submodels covering the full cell cycle: ODEs, FBA, stochastic simulation, boolean logic. No spatial resolution. The submodels are coordinated by a central time-stepping framework that passes state every ~1 second. The contribution was integration (getting 28 processes to run together and produce a coherent cell cycle), not spatial fidelity.

**Macklin et al. 2020 (Covert Lab)** extended this to *E. coli* with the Vivarium framework. Again, primarily well-mixed. The focus was composability and multi-scale coordination.

**E-Cell (Tomita et al.)** was one of the earliest whole-cell modelling efforts, going back to the late 1990s. ODE-based. No spatial resolution.

**Constraint-based models (FBA and variants)** are by definition steady-state and well-mixed. No dynamics at all, let alone spatial dynamics.

The Luthey-Schulten lab is the only group that has made spatial resolution central to whole-cell modelling, first with the 3D minimal cell (Thornburg et al. 2022) and then MC4D. This is their distinctive contribution to the field.

The landscape is clear: the field has pushed toward more spatial resolution (MC4D) and more biological completeness (Karr 2012, Macklin 2020). Nobody has pushed toward inference. That is the gap InferCell fills.

## Why spatial whole-cell models only exist for minimal cells

MC4D models syn3A not because it has few genes (493), but because it is a small, simple cell. The RDME cost scales with the number of molecular species multiplied by the number of voxels, not with gene count.

At ~10nm resolution, a syn3A cell (~400nm diameter) gives roughly 8,000-10,000 voxels. Each voxel tracks copy numbers for every diffusing species. That state space is manageable on GPUs. A spatial *M. genitalium* (525 genes, similar cell size) would be roughly the same cost and is entirely feasible.

Larger cells are a different story. *E. coli* is ~1-2 μm long with ~4,300 genes. The voxel count increases with volume (roughly 10-50x more than syn3A depending on resolution) and the species count increases ~10x. The state space explodes. Eukaryotic cells (10-100 μm) with internal compartments (nucleus, ER, mitochondria) make full RDME essentially impossible at current computational capacity.

This reinforces the inference-first architecture. Spatial simulation only works for the smallest, simplest cells. Well-mixed models scale much further. As organisms get larger and more complex, inference becomes more important (more parameters, more uncertainty) at exactly the point where spatial resolution becomes less feasible.

## The climate analogy

Weather and climate models face the same tension. Space is foundational (the atmosphere is discretised onto a grid), but inference is also essential (data assimilation updates the model state using observations). Climate science resolved this by building spatial simulation first and layering inference on top (EnKF, 4D-Var). The spatial grid came first; data assimilation was added decades later.

InferCell takes the opposite path: inference first, spatial resolution layered on where needed. This is a legitimate alternative ordering. But it does mean that certain spatial questions (emergent cross-module organisation) are deferred until the architecture supports shared spatial coupling.

## Summary

InferCell does not have RDME and does not need it for v1. The architecture can accommodate it as a future module without redesign. The fact that spatial dynamics is not required as a foundation is a feature of the inference-first design: you add spatial complexity where the data justifies it, and you can formally test whether it is justified at all.

The framing for the paper: InferCell addresses the inference axis of whole-cell modelling, complementary to MC4D's spatial axis. The ultimate goal of a spatially resolved cell model with calibrated, uncertainty-quantified parameters is a convergence point that neither project can reach alone.
