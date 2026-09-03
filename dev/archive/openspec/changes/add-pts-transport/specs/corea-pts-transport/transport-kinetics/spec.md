## Purpose

The eleven rate constants this module runs on, every one of them a point value
with no uncertainty anywhere in the source model, and the initial conditions of
the nine states it owns — imported so that the difference between a prior this
project asserted and a prior it inherited is visible in the parameter itself.

## ADDED Requirements

### Requirement: The cascade uses mass-action kinetics

Each of the five cascade steps SHALL be evaluated as a reversible second-order
mass-action rate: the forward constant times the product of its two substrate
concentrations, less the reverse constant times the product of its two product
concentrations.

This differs from the modular saturating rate law Core A′'s intracellular
reactions take, and matches the published transport model, which writes these
five reactions as mass-action and the intracellular reactions as convenience
kinetics.

#### Scenario: A cascade step at equilibrium carries no net flux
- **WHEN** substrate and product concentrations are set so the forward and
  reverse terms are equal
- **THEN** the step's net rate is zero

#### Scenario: No saturation term appears
- **WHEN** every substrate concentration of a cascade step is doubled
- **THEN** its forward term increases fourfold

### Requirement: Every rate constant carries a prior this project asserts

The ten cascade constants and the membrane permeability SHALL be imported with
their informedness recorded as asserted rather than balanced. The source model
provides no balancing distribution for any of them — the transport file carries
no parameter table at all — so any prior placed on them is this project's.

An asserted value SHALL be distinguishable, through the parameter alone, from a
value carrying a wide but genuine inherited prior.

#### Scenario: All eleven constants enumerate as asserted
- **WHEN** the composed model is queried for parameters whose prior this project
  asserted rather than inherited
- **THEN** all ten cascade constants appear, and the membrane permeability appears

#### Scenario: Asserted is distinct from wide-but-inherited
- **WHEN** an asserted constant and a balanced constant with a wide geometric
  standard deviation are compared
- **THEN** the two are distinguishable by informedness, not only by prior width

#### Scenario: A source table with no uncertainty column yields asserted values
- **WHEN** the vendored rate-constant table is read
- **THEN** every row loads as asserted, because the table carries no uncertainty
  column to derive anything else from

### Requirement: Initial conditions carry their derivation

The nine owned states' initial conditions SHALL be imported with provenance, and
the eight phospho-states' values SHALL be derived from two published quantities:
the per-cell copy number of each carrier protein, and the phosphorylated
fraction the source model's own proteomics data records.

Both inputs SHALL be recorded, so that the derived concentration can be traced to
the copy number and the fraction it came from rather than appearing as a bare
number. The resulting per-carrier totals SHALL reproduce the published copy
numbers at the model's initial volume.

The split SHALL be recorded as the published model's, not as this project's: it
is read from live model data, not asserted here. The uniform placeholder the
transport reconstruction gives these species SHALL NOT be used, and the reason
SHALL be recorded — it contradicts the proteomics counts for three of the four
carriers.

#### Scenario: Totals reproduce the published copy numbers
- **WHEN** each carrier's two initial concentrations are summed and converted to
  copies at the model's initial volume
- **THEN** the four totals equal the published per-cell copy numbers

#### Scenario: The derivation is traceable
- **WHEN** a phospho-state's initial condition is read from the composed model
- **THEN** it reports both the copy number and the phosphorylated fraction it was
  derived from, and the files each came from

#### Scenario: The split is published, not asserted
- **WHEN** the composed model is queried for the declarations that are this
  reduction's rather than the published model's
- **THEN** the phospho-state split does not appear

#### Scenario: External lactate starts empty
- **WHEN** `M_lac__L_e`'s initial condition is read
- **THEN** it is zero, matching the value the published model initialises, not
  the transport reconstruction's differing value
- **AND** the discrepancy between the two sources is recorded

### Requirement: The imported values are available without a source-model checkout

Everything this sub-model imports SHALL be available from within the repository,
so that it can be constructed, integrated and tested with no network access and
no checkout of the source model.

The in-repository copies SHALL record which upstream files and which upstream
commit they were derived from, and the transformation SHALL be re-runnable rather
than hand-made. Reshaping SHALL change no value, and a renamed identifier SHALL
remain traceable to its upstream row.

Values coming from different upstream files SHALL be imported as distinct
sources, so that provenance names the file each value actually came from rather
than one file standing in for several.

#### Scenario: The sub-model constructs with no external data
- **WHEN** the sub-model is constructed on a machine with no source-model
  checkout
- **THEN** construction succeeds with every value imported and its provenance
  recorded

#### Scenario: Provenance distinguishes the upstream files
- **WHEN** a rate constant and an initial condition are each queried for their
  source
- **THEN** they report different files, matching where each value actually came
  from

#### Scenario: The in-repository copies are reproducible
- **WHEN** the transformation is re-run against the recorded upstream files at
  the recorded commit
- **THEN** it reproduces the in-repository copies exactly
