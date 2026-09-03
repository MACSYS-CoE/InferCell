## Purpose

What a module that owns no registry state declares at its boundary: the costs it
debits, the pools its rate constants are rebuilt from, and the two nucleotides it
holds fixed — plus the requirement that it simulate alone, before any module it
couples to exists.

## ADDED Requirements

### Requirement: The module owns no registry state

The sub-model SHALL integrate no Core A′ registry species. Its states are the
seventeen transcripts and the five cost counters, none of which the registry names.

Gene copy number SHALL be carried as a fixed per-gene quantity rather than as a
state, because Core A′ has no replication and nothing in it can change a gene's
copy number. The sub-model SHALL record that the published model carries genes as
species, so a later change adding replication knows what to promote.

#### Scenario: No registry species is integrated
- **WHEN** the sub-model's states are enumerated
- **THEN** none of them is a Core A′ registry species

#### Scenario: Gene copy number is fixed, and the choice is recorded
- **WHEN** the sub-model is simulated
- **THEN** every gene's copy number is unchanged throughout
- **AND** the composed model records that the published model carries genes as
  species

### Requirement: Rate-constant edges declare the rebuild this module does not perform

The sub-model SHALL declare one inbound rate-constant edge per nucleotide
triphosphate its rate constants read — four in total — at the published refresh
cadence, held constant between refreshes.

Direction follows the information: this module's constants are rebuilt from those
pools, so it declares inbound. No mass moves on these edges.

#### Scenario: Four rate-constant edges at the published cadence
- **WHEN** the module's coupling declarations are enumerated
- **THEN** four inbound rate-constant edges are present, one per nucleotide
  triphosphate
- **AND** each records a fixed-interval cadence at the published interval

#### Scenario: Rate-constant edges move no mass
- **WHEN** the resolver examines these edges
- **THEN** none contributes to a mass or moiety balance

### Requirement: Cytidine and uridine triphosphate are clamped, and the clamp is ours

The sub-model SHALL declare a clamped edge for each of the two chemostatted
nucleotides, recording the clamp as this reduction's simplification rather than the
published model's, and holding each at the balanced value the nucleotide source
file records.

This is the opposite finding from external glucose, which the published model also
clamps. Here the sources are in a module Core A′ cut, so the chemostat is ours,
and the two clamps SHALL appear in the enumeration of what is ours.

#### Scenario: Both clamps are labelled as ours
- **WHEN** the composed model is queried for the declarations that are this
  reduction's rather than the published model's
- **THEN** the cytidine and uridine clamps both appear

#### Scenario: The held values carry their provenance
- **WHEN** a clamp's held value is read
- **THEN** it reports the source file it came from and its uncertainty

#### Scenario: A chemostat is not requested as an input
- **WHEN** the module's untyped inputs are enumerated
- **THEN** they are empty, and in particular contain no chemostatted species

### Requirement: The module resolves, simulates and is tested alone

The module SHALL be composable, resolvable and simulable by itself. No coupling
edge SHALL name a peer module.

Resolving alone SHALL report — not fail on — the pools it debits that no module in
the composition produces. Simulating alone SHALL produce a trajectory in which the
cost counters accumulate and are never debited, since no hook is present, and the
result SHALL be labelled as such.

#### Scenario: The module resolves alone
- **WHEN** the module is composed by itself and its coupling is resolved
- **THEN** resolution succeeds, and no edge fails for a missing counterpart

#### Scenario: The module simulates alone
- **WHEN** the module is simulated by itself over a full cell cycle
- **THEN** a trajectory is produced
- **AND** each gene's transcript count stays within the low-copy regime the
  measured data describes

#### Scenario: Undebited counters are labelled
- **WHEN** a trajectory is produced from a composition with no hook to debit the
  costs
- **THEN** the result records that the cost counters accrued without being debited
- **AND** the record names the five counters
