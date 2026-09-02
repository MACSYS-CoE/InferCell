## Purpose

Seventeen genes, one constitutive transcription reaction each, as a discrete
jump process — and the recorded fact that the motif has no promoter switching,
so that the project's existing two-state machinery is not reached for by
resemblance.

## ADDED Requirements

### Requirement: Seventeen genes, one transcription reaction each

The system SHALL provide a Core A′ sub-model carrying exactly seventeen
transcription reactions, one per gene, covering the ten glycolytic loci, the four
phosphotransferase loci and the three nucleotide-recycling loci that Core A′'s
metabolic reactions require.

Each reaction SHALL be catalytic in its gene: the gene appears as both reactant
and product, so transcription does not consume it. Each SHALL produce exactly one
transcript of its own gene.

The seventeen loci SHALL be distinct, and SHALL be exactly the union of the loci
the three ODE modules name. A locus present in the metabolic modules but absent
here, or vice versa, is an inconsistency in the composition rather than a matter
of taste.

#### Scenario: The reaction set is exactly seventeen
- **WHEN** the sub-model's reactions are enumerated
- **THEN** exactly seventeen transcription reactions are present, one per locus
- **AND** the seventeen loci are distinct

#### Scenario: Transcription does not consume the gene
- **WHEN** a transcription reaction fires
- **THEN** the gene's copy number is unchanged
- **AND** the transcript count increases by exactly one

#### Scenario: The gene set matches the metabolic modules
- **WHEN** this module's loci are compared with those the Core A′ metabolic
  modules declare
- **THEN** the two sets are equal
- **AND** a mismatch is reported naming the loci on each side

### Requirement: The formalism is discrete, and inference on it is simulation-based

The sub-model SHALL declare a jump-process formalism and a simulation-based
inference mode.

Transcripts in this organism exist at zero to two copies per cell, so a continuous
state variable would not approximate anything. The measured mean count across
these seventeen genes is below one, and the reactions SHALL be represented as
discrete firings rather than as a rate of change of a concentration.

#### Scenario: The sub-model declares a jump formalism
- **WHEN** the sub-model's formalism is queried
- **THEN** it is a jump process
- **AND** its inference mode is simulation-based rather than gradient-based

#### Scenario: Transcript counts are integers
- **WHEN** the sub-model is simulated
- **THEN** every transcript count is a non-negative integer at every reported time

### Requirement: The promoter is constitutive, and the record says so

The transcription motif SHALL have exactly one rate constant per gene. There SHALL
be no promoter on-state and off-state, and no switching rates between them.

The sub-model SHALL record that this is the published model's motif, and SHALL
record the consequence: identifiability results derived for a two-state promoter
do not transfer to it, because the degeneracy they describe requires two states to
be degenerate between. This record SHALL be retrievable from the composed model.

This module resembles a bursty-transcription model closely enough that reaching
for two-state machinery is the natural mistake, which is why the record exists.

#### Scenario: One rate constant per gene
- **WHEN** a gene's transcription parameters are enumerated
- **THEN** exactly one rate constant governs its transcription
- **AND** no promoter switching rate is present

#### Scenario: The non-transfer of the two-state result is recorded
- **WHEN** the composed model is queried for this module's recorded caveats
- **THEN** the constitutive-motif entry is returned, naming the result that does
  not transfer and why

### Requirement: Initial transcript counts come from the measured distribution

Each gene's initial transcript count SHALL be drawn from the measured mean count
the source model records for that locus, using the same distribution the published
model uses to initialise it.

Initialisation SHALL be reproducible: a run SHALL be repeatable from a recorded
seed, so that a trajectory can be regenerated exactly.

#### Scenario: Initial counts are drawn, not fixed
- **WHEN** the sub-model is initialised twice with different seeds
- **THEN** the initial transcript counts may differ
- **AND** each gene's counts across many seeds have a mean close to its measured
  value

#### Scenario: A seeded run is reproducible
- **WHEN** the sub-model is initialised twice with the same seed
- **THEN** the initial transcript counts are identical
