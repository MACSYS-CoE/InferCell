## Purpose

A typed contract for the seven distinct ways state crosses a module boundary in
Core A′, so that a module declares not only *which* state it couples to but *how*,
and so that a disagreement between two modules is caught when they are composed
rather than when they are run.

## ADDED Requirements

### Requirement: Modules declare typed coupling edges

The sub-model protocol SHALL be extended so that a module can declare its coupling
as a list of typed edges, each naming a registry species, a direction, and one of
the seven edge kinds below. The declaration SHALL default to empty, so that a
module with no cross-block coupling — and every sub-model written before this
change — remains valid without modification.

The existing untyped input declaration SHALL keep its current meaning and its empty
default. A module MAY declare both; where a species appears in both, the typed
declaration governs.

#### Scenario: Existing sub-models remain valid
- **WHEN** a sub-model written before this change is composed
- **THEN** it composes and runs exactly as before
- **AND** its coupling declaration reads as empty

#### Scenario: A module declares typed coupling
- **WHEN** a module declares a coupling edge naming a species, a direction and a
  kind
- **THEN** the declaration is retrievable from the module through the protocol
- **AND** it survives composition into a resolved form that names the counterpart
  module

### Requirement: The seven edge kinds are distinct and exhaustive

The system SHALL define exactly seven edge kinds, matching the legend of
`fig1r_state_graph_reduced.pdf`, each with the fields its own semantics require:

| Kind | Semantics | Fields beyond species and direction |
|---|---|---|
| mass | shared state, continuous; gradients cross inside the ODE block | none |
| currency | routed via a shared pool node rather than module to module | the pool species, defaulting to the species itself; the debited pool is the edge's own species |
| deferred counter | the stochastic block accrues a cost, the hook debits it against the edge's species one step later | the counter name (required); a clipping policy defaulting to the published clamped drain; the smoothing width, required exactly when the policy is smoothed |
| catalytic | counts enter a rate law as parameters; no mass flows | the rate law's parameter slot (required) |
| rate constant | pools re-enter the stochastic block as recomputed rate constants | a refresh cadence defaulting to the published 60 s piecewise-constant rebuild |
| volume | counts set surface area, which sets volume, which rescales every concentration | none |
| clamped | a real dependence replaced by a constant | whether the clamp is the published model's or ours (required); optionally the held value, which must match the registry where the registry records one |

An edge whose kind is not one of these seven SHALL be rejected. Fields marked
required have no default; the defaulted fields default to the published model's
behaviour, so departing from the published model is something an author has to
write down.

#### Scenario: Each kind carries its own fields
- **WHEN** an edge of a given kind is constructed without a field that kind
  requires
- **THEN** construction fails naming the missing field and the kind
- **AND** the failure is at declaration, not at composition or at run time

#### Scenario: An unknown kind is rejected
- **WHEN** a module declares an edge whose kind is outside the seven
- **THEN** resolving it fails naming the offending kind and listing the seven
  valid ones — kinds are types, so the foreign declaration is caught at the
  first point the vocabulary is consulted

#### Scenario: Catalytic edges move no mass
- **WHEN** the resolver examines a catalytic edge
- **THEN** the edge contributes to no mass or moiety balance
- **AND** an attempt to include a catalytic edge in a conservation check is
  rejected rather than silently counted

### Requirement: Deferred counters declare a clipping policy

A deferred-counter edge SHALL carry a clipping policy stating how a debit larger
than the available pool is handled. The permitted policies are: clamped at zero
with the deficit carried forward to the next step, as the published model does;
unclamped, allowing the pool to go negative; and smoothed, replacing the
discontinuity with a differentiable approximation. The policy defaults to the
published clamped drain, so selecting anything else is an explicit act.

The declaration SHALL record which policy is in force, so that a trajectory
sampled with a gradient-based method can be checked against a policy that admits
gradients.

#### Scenario: The published clipped drain is declarable
- **WHEN** a module declares the expression-cost drain against ATP with the
  clamped-at-zero policy
- **THEN** the edge resolves, and is marked as non-differentiable at the boundary

#### Scenario: A non-differentiable edge under a gradient-based mode is reported
- **WHEN** a composition containing a clamped-at-zero deferred counter is prepared
  for a sub-model whose inference mode is differentiable
- **THEN** the system reports the edge and the species it clips, naming it as an
  obstruction to gradients
- **AND** the report is a diagnosable warning rather than a silent success, so the
  choice to proceed is deliberate

#### Scenario: Smoothing is an explicit choice
- **WHEN** the smoothed policy is selected
- **THEN** the edge records that its behaviour differs from the published model
- **AND** the parameter controlling the smoothing is exposed rather than hidden

### Requirement: Rate-constant edges declare a refresh cadence

A rate-constant edge SHALL declare whether the downstream rate constants are
recomputed at a fixed interval, holding piecewise-constant between refreshes as the
published model's 60 s rebuild does, or track the upstream pools continuously.
A fixed-interval declaration SHALL carry its interval. The cadence defaults to
the published 60 s rebuild, so a continuous edge is an explicit act.

#### Scenario: The published 60 s rebuild is declarable
- **WHEN** a module declares transcription rate constants refreshed from live NTP
  pools at a 60 s interval
- **THEN** the edge resolves with cadence recorded as piecewise-constant at 60 s

#### Scenario: A continuous variant is declarable and marked as a deviation
- **WHEN** a module declares the same edge as continuous
- **THEN** the edge resolves
- **AND** it is marked as deviating from the published model, so any result that
  depends on it can be labelled

### Requirement: Clamped edges record whether the clamp is ours

A clamped edge SHALL record whether replacing the dependence with a constant is
what the published model does, or a simplification introduced by this reduction.
This label SHALL be retrievable from the composed model, so that it can reach any
report of a result that depends on it.

#### Scenario: A reduction-introduced clamp is labelled
- **WHEN** CTP, UTP or the amino-acid pool is declared chemostatted
- **THEN** the edge is labelled as this reduction's simplification rather than the
  published model's
- **AND** enumerating a composed model's reduction-introduced clamps returns them

#### Scenario: Enumerating what is ours
- **WHEN** a composed Core A′ model is queried for every declaration marked as this
  reduction's rather than the published model's
- **THEN** the result includes the chemostatted CTP, UTP and amino-acid pool, the
  lumped tRNA charging step, and any smoothed deferred counter or continuous
  rate-constant edge in force

### Requirement: The resolver validates a composed set of declarations

Given a set of modules, the system SHALL resolve every declared edge and SHALL fail
with a named, actionable error rather than composing a model whose boundary is
inconsistent.

The resolver SHALL detect at least:
- an edge naming a species that is not in the state registry,
- an edge whose counterpart module does not exist in the composition,
- the same species-and-direction pair declared by two modules with different kinds,
- a directional edge with no producer, or none with a consumer.

#### Scenario: Unknown species in an edge
- **WHEN** a module declares an edge naming a species absent from the registry
- **THEN** resolution fails naming the species and the declaring module

#### Scenario: Two modules disagree on an edge kind
- **WHEN** one module declares a species crossing as mass and another declares the
  same species and direction as a deferred counter
- **THEN** resolution fails naming both modules, the species, and the two kinds
- **AND** the error text distinguishes the two semantics rather than reporting a
  generic conflict

#### Scenario: An edge with a missing counterpart
- **WHEN** a module declares an edge to a module that is not in the composition
- **THEN** resolution fails naming the missing counterpart
- **AND** the failure is distinguishable from a deliberate single-module
  composition under test

#### Scenario: A consistent composition resolves
- **WHEN** every declared edge names a registered species and agrees with its
  counterpart
- **THEN** resolution succeeds and returns, for each edge, its kind, its endpoints
  and the registry positions it touches

### Requirement: Every declared cost has a state that pays it and a reaction that returns it

The resolver SHALL report each deferred-counter and mass edge that consumes a
species the composition gives no producer for, and SHALL report the species no
module integrates. Both are reports rather than failures, because a module
validated alone legitimately imports species whose owners and producers are
absent — the standalone-validation constraint the contract exists to serve.
This is the mechanical form of the rule the scoping note derived from two
errors of record: the adenylate dead end that would have exhausted the pool
after 2.3% of the cell cycle, and the stranded GMP.

Because these are reports, **a resolution that succeeds is not evidence that the
boundary is closed.** Asserting completeness — every registry dynamic state
owned, every dead end empty — is out of scope for this capability. A composition
intended to be the whole model therefore needs that assertion made separately,
and MUST NOT infer it from a successful resolution.

#### Scenario: A cost with no paying state
- **WHEN** a module declares a cost debited against a species that no module
  integrates and the registry does not chemostat
- **THEN** resolution succeeds, listing the species among the graph's unowned
  states and — absent a declared producer — as a dead end naming the consumer

#### Scenario: A cost with no return path
- **WHEN** a module declares a consumer of a dynamic species and no module in the
  composition declares a producer of it
- **THEN** resolution reports the species as a dead end, naming the consumer
- **AND** the report identifies which conserved moiety the species belongs to, so
  that AMP with no ADK1 and GMP with no GK1 are recognisable as the same class of
  error

#### Scenario: A chemostatted species is exempt
- **WHEN** a declared cost is debited against a species the registry marks as
  chemostatted
- **THEN** no dead-end error is raised, because the chemostat is the return path
- **AND** the exemption is recorded, so that a later change making the pool live
  reinstates the check
