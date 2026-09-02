# Core A′ state registry Specification

## Purpose

The single canonical naming and ordering of every species in the Core A′ reduced
syn3A model, so that seven independently developed modules refer to the same state
by the same name in the same position, and none of them invents its own.

## Requirements

### Requirement: Canonical species registry

The system SHALL provide one registry that names and orders every Core A′ species:
32 dynamic states and 5 chemostats, matching the state list in
`dev/notes/reduced-syn3a-scoping.md`. The registry SHALL be the only place these
names are defined.

Each entry SHALL carry:
- a canonical name (a symbol),
- a group — one of glycolytic intermediate, adenylate, guanylate, redox, other,
  tRNA, PTS phospho-state, chemostat,
- a treatment — dynamic or chemostatted,
- a copy-number regime — one of mRNA-scale, protein-scale, metabolite-scale —
  since this is what determines which formalism is defensible for the species.

The registry SHALL be ordered, and its order SHALL be stable across calls and
across processes: a module that indexes position *i* today indexes the same species
tomorrow.

#### Scenario: Registry contents match the scoping note
- **WHEN** the registry is enumerated
- **THEN** it contains exactly 32 dynamic states and 5 chemostats
- **AND** the 11 glycolytic intermediates, 4 adenylate species, 3 guanylate
  species, 2 redox species, 2 other species, 2 tRNA species and 8 PTS phospho-states
  named in the scoping note's state list are all present, each exactly once
- **AND** the 5 chemostats are external glucose, CTP, UTP, the amino-acid pool
  and O2

#### Scenario: Ordering is stable
- **WHEN** the registry order is read twice in the same session, or in two separate
  sessions
- **THEN** both readings give the same species in the same positions

#### Scenario: H2O and H+ are absent
- **WHEN** the registry is queried for water or the hydrogen ion
- **THEN** neither is present, because Core A′ is derived from the `NoH2O` model
  with hydrogen-ion accounting removed

### Requirement: Lookup by name and by position

The system SHALL let a module resolve a species name to its registry position and a
position to its entry, in both directions, without the module holding its own copy
of the ordering.

#### Scenario: Round-trip lookup
- **WHEN** a module resolves a registered species name to a position and then
  resolves that position back to an entry
- **THEN** the entry's name equals the name it started from

#### Scenario: Unknown species is rejected
- **WHEN** a module resolves a name that is not in the registry
- **THEN** the system raises an error naming the unrecognised symbol
- **AND** the error does not silently return a default position

### Requirement: Chemostatted species are marked and cannot be integrated

A chemostatted species SHALL be distinguishable from a dynamic state through the
registry alone, so that composition can reject a module that tries to give a
chemostat its own dynamics.

#### Scenario: Module declares dynamics for a chemostat
- **WHEN** a module declares that it integrates a species the registry marks as
  chemostatted
- **THEN** composition fails with an error naming that species and the declaring
  module

#### Scenario: Chemostats carry their held value
- **WHEN** a chemostatted species is queried
- **THEN** its held concentration is available from the registry, including
  external glucose at 40 mM and the amino-acid pool at 0.1 mM

### Requirement: Each dynamic state is owned by exactly one module

Every dynamic state in the registry SHALL be integrated by exactly one module in a
complete Core A′ composition. A state integrated by two modules is an error in
the composition rather than a runtime surprise; a state integrated by none is
reported rather than failed, so a partial composition — a single module under
test — still resolves.

#### Scenario: Two modules claim the same state
- **WHEN** two modules in one composition both declare they integrate the same
  dynamic state
- **THEN** composition fails with an error naming the state and both modules

#### Scenario: A dynamic state has no owner
- **WHEN** a composition is assembled that leaves a registry dynamic state
  unintegrated
- **THEN** composition reports the unowned state
- **AND** the report distinguishes a genuinely incomplete composition from a
  deliberate partial one, so a single module can still be built and tested alone

### Requirement: Registry provenance for measured values

Where a registry entry carries a nominal initial value, it SHALL also carry the
source file that value was read from, and the entries with no informed value SHALL
be marked as such rather than presented as measurements.

#### Scenario: Nucleotide species resolve to the nucleotide file
- **WHEN** the initial value of GTP, GDP, AMP or GMP is read from the registry
- **THEN** the value is the nucleotide balanced file's value, not the central
  file's 0.1 mM prior default
- **AND** the entry names the nucleotide file as its source

#### Scenario: Uninformed initial conditions are flagged
- **WHEN** the initial value of PPi, cytosolic lactate or G3P is read
- **THEN** the entry is marked as sitting at prior median and prior width rather
  than at a balanced value
