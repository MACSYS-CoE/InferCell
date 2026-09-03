## Purpose

What a transcript costs, accrued as counters the hook debits one step later:
the energy the polymerase hydrolyses and the four nucleotides it incorporates —
five distinct costs, charged against two live pools and two this reduction
chemostats.

## ADDED Requirements

### Requirement: Transcription accrues five distinct costs

Each transcription firing SHALL increment five counters: one energy counter by the
transcript's length, and four monomer counters by that transcript's count of each
base.

The energy cost and the monomer cost are distinct and both are charged, as the
published model charges them. A reader encountering the energy counter alone would
reasonably assume the nucleotides were free, and one encountering the monomer
counters alone would assume polymerisation cost nothing; the sub-model SHALL record
that both are intended.

#### Scenario: One firing increments five counters
- **WHEN** a gene's transcription reaction fires once
- **THEN** the energy counter increases by that gene's transcript length
- **AND** each monomer counter increases by that gene's count of the corresponding
  base

#### Scenario: The double charge is recorded as intended
- **WHEN** the composed model is queried for this module's recorded caveats
- **THEN** an entry states that adenosine triphosphate is debited twice per
  transcript — once as polymerisation energy and once as an incorporated monomer —
  and that this is the published model's accounting

### Requirement: Each counter declares what it drains into

Every counter SHALL be declared as a deferred cost against a named registry
species, so that the module debiting it and the module owning the pool agree
without either reading the other's source.

The energy counter SHALL be declared against adenosine triphosphate, and the four
monomer counters against their respective nucleotide triphosphates. The sub-model
SHALL record what each drain produces — the energy counter yielding diphosphate
and inorganic phosphate, the monomer counters yielding pyrophosphate — because
those products are another module's states and its conservation checks depend on
them.

Pyrophosphate from transcription SHALL be recorded explicitly as a source the
project's earlier accounting attributed only to amino-acid charging.

#### Scenario: Each counter names its pool
- **WHEN** the sub-model's declared costs are enumerated
- **THEN** five are returned, each naming its counter and a registry species
- **AND** every species named is one the registry carries

#### Scenario: The drain products are recorded
- **WHEN** the composed model is queried for what this module's costs produce
- **THEN** the energy counter reports diphosphate and inorganic phosphate, and the
  monomer counters report pyrophosphate

### Requirement: The drains take the published clamped policy

Every one of the five costs SHALL use the published drain policy: the pool floors
at zero and the shortfall is carried forward to the next step. Selecting any other
policy SHALL be an explicit act and SHALL be recorded as a departure from the
published model.

Because that policy is non-differentiable at the boundary, a composition preparing
these costs for a gradient-based sub-model SHALL report them rather than proceed
silently. This module's own inference mode is simulation-based, so the report
belongs to the composition that joins it to a differentiable block, not to this
module alone.

#### Scenario: The published policy is the default
- **WHEN** the sub-model's declared costs are enumerated
- **THEN** each carries the clamped-at-zero policy with the deficit carried forward
- **AND** none is recorded as departing from the published model

#### Scenario: The gradient obstruction is reported at the join
- **WHEN** this module is composed with a sub-model whose inference mode is
  gradient-based
- **THEN** the composition reports the clamped costs and the species they clip
- **AND** the report is a diagnosable warning rather than a silent success

#### Scenario: This module alone raises no obstruction
- **WHEN** this module is composed by itself
- **THEN** no gradient obstruction is reported, because its own inference mode is
  simulation-based

### Requirement: The chemostatted nucleotide pools are exempt from the return-path check, and the exemption is recorded

Cytidine and uridine triphosphate are held at fixed concentrations by this
reduction, so a cost debited against either has a return path by construction and
SHALL NOT be reported as a dead end.

That exemption SHALL be recorded rather than being silent, so that a later change
making either pool live reinstates the check. Adenosine and guanosine triphosphate
SHALL NOT be exempt: both are integrated elsewhere in Core A′, and a composition
lacking their producers should say so.

#### Scenario: The chemostatted pools are exempt and the exemption is recorded
- **WHEN** the composition is resolved
- **THEN** the cytidine and uridine costs are not reported as dead ends
- **AND** the exemption is recorded against each

#### Scenario: The live pools are not exempt
- **WHEN** this module is resolved without the module owning the adenylate and
  guanylate pools
- **THEN** adenosine and guanosine triphosphate are reported as dead ends naming
  this module as the consumer
- **AND** the report identifies the conserved moiety each belongs to
