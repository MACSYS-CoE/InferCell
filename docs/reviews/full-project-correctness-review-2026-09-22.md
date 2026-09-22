# Full-project correctness review

**Date:** 2026-09-22
**Scope:** Entire reduced-WCM codebase through Phase 10
**Review style:** Correctness-focused; minor style and defensive-programming issues are intentionally omitted.

## Executive summary

The project has a thoughtful forward-model architecture, and the individual Core A′ biological modules are generally internally consistent and well tested. The important problems are at the composition and inference boundaries rather than in the basic rate laws.

The current Core A′ modules cannot be assembled into the intended composite model, and several inference paths either ignore sampled values or do not represent the hybrid model correctly. These issues should be fixed before relying on Phase 13 assembly or posterior results.

## Findings

### 1. Core A′ cannot currently be assembled

**Severity: critical**

The transcription module declares CTP and UTP as deferred-debit resources and uses them as rate inputs ([`transcription.jl`](../../src/organisms/coreA/transcription.jl#L389)). The registry, however, declares those resources as chemostats ([`registry.jl`](../../src/organisms/coreA/registry.jl#L204)). Debit and rebuild lowering require a live ODE-owned state ([`handshake.jl`](../../src/framework/handshake.jl#L794), [`handshake.jl`](../../src/framework/handshake.jl#L1012)), which a chemostat does not provide.

A focused assembly reproduction failed for this reason. Removing the explicit CTP/UTP counter edges did not make the model assemble: rebuilding the transcription rate still requires a live owner for these species.

This is a direct incompatibility between the registry semantics, transcription model, and handshake lowering. One of those contracts must change before full Core A′ assembly is possible.

### 2. Transcription products are declared but never wired into the model

**Severity: critical**

`TRANSCRIPTION_COUNTERS` describes ADP, phosphate, and pyrophosphate products, but the corresponding transcription edges are inbound-only. The driver credits only outbound resource edges ([`handshake.jl`](../../src/framework/handshake.jl#L1582)), so those products are never credited.

The handshake layer also rejects more than one credit pool for an event ([`handshake.jl`](../../src/framework/handshake.jl#L852)). Consequently, a reaction such as ATP consumption producing both ADP and phosphate cannot be represented by the current event-edge contract.

The metadata therefore reports chemistry that the executable composite model does not implement.

### 3. The inference architecture does not support the hybrid model it is meant to infer

**Severity: critical**

Inference dispatches according to `models[1]` ([`inference.jl`](../../src/inference/inference.jl#L78), [`inference.jl`](../../src/inference/inference.jl#L84)). The downstream likelihood paths assume a SciML problem that can be passed through `remake` and `solve`, while a mixed ODE/jump composition returns a `HandshakeDriver`.

The iterative boundary algorithm identifies shared quantities by matching parameter names ([`boundary.jl`](../../src/inference/boundary.jl#L66)). The current Core A′ modules have no such shared parameter names. `iterative_infer` therefore receives an empty shared set, computes a zero KL divergence, and reports convergence after the second iteration ([`boundary.jl`](../../src/inference/boundary.jl#L263)) without exchanging biologically relevant boundary information.

In its current form, the inference layer cannot perform meaningful inference over the intended assembled hybrid model.

### 4. Several sampled parameters do not affect the likelihood

**Severity: high**

Free initial conditions are sampled, but `_collect_ic_values` reads the stored nominal values rather than the sampled values ([`orchestrator.jl`](../../src/inference/orchestrator.jl#L312)). Those samples therefore do not change the simulated trajectory.

Driver-written, rebuilt, catalytic, and geometry parameters can also enter the sampled parameter vector and then be overwritten during model construction or execution. This creates posterior dimensions that appear inferable but have no effect on the likelihood.

One concrete example is PTS transport: `r_cell` is forcibly freed ([`pts_transport.jl`](../../src/organisms/coreA/pts_transport.jl#L169)), so callers cannot pin it even though the public contract describes pinning model parameters.

Inference parameters need a single source of truth, and every free parameter should be verified to perturb the generated likelihood.

### 5. Shared-parameter validation accepts incompatible priors

**Severity: high**

Composition validates shared priors by comparing only `typeof(prior)` ([`orchestrator.jl`](../../src/inference/orchestrator.jl#L81)). Priors with the same distribution family but different parameters are treated as identical.

A focused reproduction composed two modules using `LogNormal(0, 0.1)` and `LogNormal(5, 2)` under the same parameter name. Composition succeeded, and the retained prior depended on module order. This makes the posterior specification order-dependent.

Validation should compare the full prior specification, not just its Julia type.

### 6. The observation layer ignores species identities

**Severity: high**

The likelihood compares observation rows directly with solution rows ([`inference.jl`](../../src/inference/inference.jl#L20)); the observation `species` field is not used to map data onto model states. Reordering observation species silently compares data against the wrong state, while observing a subset of states produces a dimension mismatch.

The same layer uses one additive scalar noise scale for states with potentially very different magnitudes. For stochastic models, ensemble trajectories are averaged before evaluating the likelihood ([`inference.jl`](../../src/inference/inference.jl#L172)), and the summary-statistics path retains means only ([`summary_statistics.jl`](../../src/inference/summary_statistics.jl#L25)). This discards variance and distributional information that is central to jump-process inference.

At minimum, observations need explicit species-to-state indexing. Noise and stochastic summaries should then be defined per observable or on an intentional transformed scale.

### 7. ABC-SMC reports a final tolerance that was not used for acceptance

**Severity: high**

Each population is accepted against the previous tolerance; only afterward is the next tolerance computed from the accepted distances ([`abc_smc.jl`](../../src/inference/abc_smc.jl#L51), [`abc_smc.jl`](../../src/inference/abc_smc.jl#L73), [`abc_smc.jl`](../../src/inference/abc_smc.jl#L94)). The routine can then return without generating a population under the tolerance it advertises as final.

In a deterministic reproduction, the reported final tolerance was approximately `0.349`, while the maximum retained distance was approximately `0.666`; 50 of 100 returned particles exceeded the advertised tolerance.

The returned tolerance should describe the population actually returned, or another population must be sampled after updating the threshold.

### 8. Resolver validation ignores the named peer

**Severity: medium**

Resolver validation checks that a named peer exists ([`resolver.jl`](../../src/framework/resolver.jl#L397)), but it does not verify that the peer owns or otherwise supplies the referenced quantity. Execution later resolves by state ownership instead.

A focused reproduction with an owner, a consumer, and an unrelated bystander successfully resolved an edge that named the bystander as its peer. The edge retained the incorrect peer metadata while execution used the real owner.

This permits invalid wiring declarations to pass validation and makes the declared topology disagree with runtime behavior.

## Smaller Phase 10 issues

These are less urgent but are clearly incorrect:

- The transcription constructor uses its `counters` argument for parameters and edges, but hardcodes the global counter list when creating states ([`transcription.jl`](../../src/organisms/coreA/transcription.jl#L462)). A custom three-counter configuration therefore creates all five counter states.
- `reduction_notes` says the corrected mapping is used even when `base_mapping = :published` ([`transcription.jl`](../../src/organisms/coreA/transcription.jl#L541)). The provenance output is false in that mode.
- The transcription test suite still skips the cross-module protein-coupling case ([`test_corea_transcription.jl`](../../test/test_corea_transcription.jl#L384)), so the most important integration behavior is not exercised.

## Test status observed during review

The default test run completed with **2,513 passing tests, 0 failures, and 2 skipped/broken tests**. That is strong evidence for the isolated modules and covered framework behavior, but it does not contradict the findings above: several tests explicitly preserve current failure behavior, and important integration paths are gated or skipped.

The full integration suite in [`test/runtests.jl`](../../test/runtests.jl#L62) did not complete during the review. In particular, the documented headline conditioned-posterior test remains unverified.

## Recommended order of work

1. Reconcile resource ownership, chemostat semantics, and deferred debit/rebuild behavior so Core A′ assembles.
2. Define executable multi-product event stoichiometry and wire transcription products.
3. Give hybrid inference an explicit model interface rather than dispatching through the first component.
4. Audit free parameters so every sampled value reaches the simulation and likelihood.
5. Add explicit observation-to-state mapping and appropriate stochastic summaries.
6. Correct shared-prior validation, ABC-SMC tolerance reporting, and peer validation.
7. Fix the smaller Phase 10 metadata/configuration issues and replace the skipped integration case with an executable test.

Until items 1–6 are addressed, full-model assembly and posterior results should not be treated as scientifically meaningful.
