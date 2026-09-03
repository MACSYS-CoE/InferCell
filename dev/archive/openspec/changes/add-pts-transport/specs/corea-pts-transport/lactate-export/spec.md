## Purpose

The passive-diffusion reaction that carries lactate out of Core A′, the external
lactate state it fills, and the volume ratio that keeps efflux unsaturated —
together with the record of why a one-line rate law is load-bearing, so that a
later reader finds the argument rather than an apparent redundancy.

## ADDED Requirements

### Requirement: Lactate export exists and is required

The sub-model SHALL carry `L_LACt2r`, transporting cytosolic lactate to the
external compartment.

Its necessity SHALL be recorded with the magnitude that establishes it, in the
code rather than only in a design document: glycolysis produces two lactate per
glucose, and at Core A′'s glucose demand that is roughly 13.9 million lactate
over one cell cycle — about 691 mM at initial volume and 345 mM at doubled
volume. Without export, cytosolic lactate exceeds the model's entire measured
phosphate pool by more than an order of magnitude within the first minutes, and
carbon balance cannot be posed at all.

An earlier specification of this model removed the reaction as redundant. That is
the failure mode the record exists to prevent.

#### Scenario: The export reaction is present
- **WHEN** the sub-model's reactions are enumerated
- **THEN** `L_LACt2r` is among them

#### Scenario: The reason it is required is retrievable
- **WHEN** the composed model is queried for this module's recorded rationale
- **THEN** the lactate accumulation figures and their consequence are returned,
  naming both the initial-volume and doubled-volume magnitudes

### Requirement: Export follows the published passive-diffusion law

Export SHALL be driven by the difference between cytosolic and external lactate,
scaled by the membrane permeability and by the surface-to-volume ratio of a
sphere — three over the cell radius — as the published transport model writes it.

No enzyme concentration SHALL scale this rate: the reaction is passive diffusion,
not catalysis, and has no gene-protein-reaction rule.

The cell radius SHALL be a named quantity rather than folded into a lumped
constant, because a later change makes the radius grow, at which point the
export rate constant changes with it.

#### Scenario: Export is proportional to the concentration difference
- **WHEN** cytosolic and external lactate are equal
- **THEN** the export rate is zero

#### Scenario: Export reverses when the gradient reverses
- **WHEN** external lactate exceeds cytosolic lactate
- **THEN** the export rate is negative

#### Scenario: The radius is exposed, not folded in
- **WHEN** the cell radius is changed
- **THEN** the export rate constant changes in inverse proportion
- **AND** the permeability and the radius are separately retrievable

### Requirement: External lactate is dynamic, and its volume ratio is ours

`M_lac__L_e` SHALL be integrated by this sub-model, matching the registry, which
lists it among Core A′'s dynamic states.

The published model instead adds every external species to its equations as a
constant, so published lactate export drains into a pool pinned at zero and never
saturates. Integrating external lactate in a single shared volume would not
reproduce that: export would stall as the two pools equilibrated, and the lactate
would move outside the cell rather than leave the system.

The sub-model SHALL therefore scale external lactate's accumulation by a
medium-to-cell volume ratio, chosen so that external lactate remains negligible
against cytosolic lactate over a full cell cycle. That ratio SHALL be exposed as
a parameter, SHALL be recorded as a prior this project asserts rather than
inherits, and SHALL be registered as a departure from the published model,
retrievable alongside Core A′'s other reduction declarations.

#### Scenario: External lactate is integrated, not held
- **WHEN** the sub-model's states are enumerated
- **THEN** `M_lac__L_e` is among them
- **AND** no clamped edge names it

#### Scenario: Export stays unsaturated over a cycle
- **WHEN** the sub-model is integrated over a full cell cycle at a cytosolic
  lactate concentration representative of Core A′'s production
- **THEN** external lactate remains below one percent of cytosolic lactate
- **AND** the export rate at the end of the trajectory is within one percent of
  its rate at the start for the same cytosolic concentration

#### Scenario: The volume ratio is labelled as ours
- **WHEN** the composed model is queried for the declarations that are this
  reduction's rather than the published model's
- **THEN** the medium-to-cell volume ratio appears, with the behaviour it
  reproduces and the behaviour it replaces

#### Scenario: A shared volume is representable and visibly different
- **WHEN** the volume ratio is set to one, making the two compartments share a
  volume
- **THEN** the total of cytosolic and external lactate is conserved
- **AND** the export rate decays toward zero, demonstrating the saturation the
  default ratio exists to avoid
