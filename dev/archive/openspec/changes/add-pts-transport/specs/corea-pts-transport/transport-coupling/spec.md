## Purpose

What this module declares where it meets the rest of Core A′ and the medium —
the clamped external glucose it draws on, the four metabolites it shares with
glycolysis, and the inputs those require — and the requirement that all of it
resolve and integrate with no sibling module present.

## ADDED Requirements

### Requirement: External glucose is a clamped edge, and the clamp is the published model's

`M_glc__D_e` SHALL be declared as a clamped edge holding the medium
concentration the registry records, with its origin recorded as the published
model's rather than this reduction's.

This is a finding rather than a convention: the published model adds every
external species to its equations as a constant parameter, so external glucose is
clamped there too. The clamp SHALL therefore NOT appear in the enumeration of
declarations that are this reduction's — unlike the chemostatted CTP, UTP and
amino-acid pool, which do.

The held value SHALL agree with the registry's. Where the source model's
reconstruction and its initialisation script give different medium
concentrations, the value the model actually runs SHALL govern, and the
discrepancy SHALL be recorded.

#### Scenario: The clamp resolves against the registry
- **WHEN** the module is composed
- **THEN** the clamped edge on `M_glc__D_e` resolves
- **AND** its held value matches the registry's chemostat concentration

#### Scenario: The published clamp is not listed as ours
- **WHEN** the composed model is queried for the declarations that are this
  reduction's rather than the published model's
- **THEN** the external glucose clamp does not appear

#### Scenario: A chemostat is not requested as an input
- **WHEN** the module's untyped inputs are enumerated
- **THEN** `M_glc__D_e` is absent, because nothing integrates a chemostat and an
  input could never deliver it

### Requirement: Shared metabolites are declared as mass edges

The four registry metabolites this module's reactions touch but does not
integrate SHALL each be declared as a mass edge, with the direction mass flows
relative to this module:

| Species | Direction | Reaction |
|---|---|---|
| `M_pep_c` | inbound | GLCpts0 draws phosphoenolpyruvate |
| `M_lac__L_c` | inbound | L_LACt2r draws cytosolic lactate |
| `M_pyr_c` | outbound | GLCpts0 supplies pyruvate |
| `M_g6p_c` | outbound | GLCpts4 supplies glucose-6-phosphate |

Each inbound mass edge names a dynamic state this module does not integrate, so
that species SHALL also appear in the module's untyped inputs — the only channel
by which it reaches the dynamics.

No currency, deferred-counter, rate-constant or volume edge SHALL be declared:
this module moves no energy currency, accrues no stochastic cost, rebuilds no
rate constants, and — although it owns the states a later change reads for
volume — does not itself drive volume.

#### Scenario: The four mass edges and one clamped edge are declared
- **WHEN** the module's coupling declarations are enumerated
- **THEN** exactly the four mass edges above and the clamped glucose edge are
  present

#### Scenario: Inbound metabolites are reachable by the dynamics
- **WHEN** the module's untyped inputs are enumerated
- **THEN** `M_pep_c` and `M_lac__L_c` are both present
- **AND** composing the module raises no disagreement between its typed edges and
  its untyped inputs

#### Scenario: No other edge kind is used
- **WHEN** the module's coupling declarations are enumerated
- **THEN** no currency, deferred-counter, rate-constant or volume edge appears

### Requirement: The module resolves, integrates and is tested alone

The module SHALL be composable, resolvable and integrable by itself, with no
sibling module present. No coupling edge SHALL name a peer module, so that each
edge resolves against whichever module owns its species and a single-module
composition does not fail for a counterpart deliberately absent.

Resolving the module alone SHALL report — not fail on — the states it consumes
with no declared producer and the registry states no module integrates. A
successful resolution SHALL NOT be read as evidence that the Core A′ boundary is
closed.

#### Scenario: The module resolves alone
- **WHEN** the module is composed by itself and its coupling is resolved
- **THEN** resolution succeeds, and no edge fails for a missing counterpart

#### Scenario: Standalone resolution reports rather than fails
- **WHEN** the module is resolved alone
- **THEN** `M_pep_c` and `M_lac__L_c` are reported as unowned states
- **AND** the report distinguishes this deliberate partial composition from an
  incomplete one

#### Scenario: The module integrates alone
- **WHEN** the module is integrated by itself from its initial conditions at
  published parameters, with its inbound metabolites held
- **THEN** a trajectory is produced over the requested interval
- **AND** all nine owned states remain non-negative

### Requirement: Declared production of a state this module does not own is not executed

This module declares outbound mass edges for the glucose-6-phosphate and pyruvate
its reactions produce, but integrates neither, and the composition machinery gives
a module no way to contribute to a derivative it does not own. Those outbound
edges are a declaration of intent that a later change executes.

A trajectory produced by this module alone, or by any composition that does not
itself execute the coupling, SHALL be labelled as one in which the shared
metabolites were held rather than integrated.

#### Scenario: Held metabolites are labelled as held
- **WHEN** a trajectory is produced from a composition that does not execute the
  metabolite coupling
- **THEN** the result records that `M_pep_c` and `M_lac__L_c` were held at fixed
  values

#### Scenario: The declaration survives for the later change to consume
- **WHEN** the module's coupling declarations are enumerated
- **THEN** the outbound mass edges on `M_g6p_c` and `M_pyr_c` are present,
  whether or not any composition executes them
