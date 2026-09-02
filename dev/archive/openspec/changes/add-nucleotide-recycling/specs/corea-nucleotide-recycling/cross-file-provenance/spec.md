## Purpose

Every value this module imports exists in two balanced files, often with
different numbers, and the wrong choice would run a reaction hundreds of times
too fast or understate a nucleotide pool by an order of magnitude. This capability
is the requirement that the choice be declared, checked, and visible in the
assembled model rather than made once by hand.

## ADDED Requirements

### Requirement: Every imported value declares its governing file

Each of this module's imported values — the five reactions' catalytic and
Michaelis constants and the eight owned states' initial concentrations — SHALL be
imported with an explicit declaration of which source file governs. Importing any
of them without that declaration SHALL fail at load time, naming the identifier,
both files and both values.

A declaration SHALL be checked against the tables rather than trusted: a governing
file that does not hold the identifier, or that is not among the files that do,
SHALL fail.

The governing file for all five reactions' kinetics SHALL be the nucleotide
balanced file. For the initial concentrations, the governing file SHALL be the one
the registry records as that species' source.

#### Scenario: An undeclared ambiguity fails
- **WHEN** any of this module's values is imported without a governing declaration
- **THEN** loading fails naming the identifier, both files and both values

#### Scenario: Kinetics resolve to the nucleotide file
- **WHEN** the five reactions' catalytic constants are read from the composed
  model
- **THEN** each reports the nucleotide balanced file as its source
- **AND** each records the central file as the alternative it was chosen over

#### Scenario: Initial conditions resolve to the registry's source
- **WHEN** each owned state's initial concentration is read
- **THEN** its governing file equals the source file the registry records for that
  species
- **AND** a mismatch between the two fails rather than importing from either

#### Scenario: The guanylate pool is not read from the central file
- **WHEN** `M_gtp_c`, `M_gdp_c`, `M_gmp_c` and `M_amp_c` are imported
- **THEN** each takes the nucleotide file's balanced value, not the central file's
  prior default
- **AND** the resulting guanylate pool is of the order of tens of thousands of
  particles rather than a few thousand

### Requirement: The cross-file disagreements are enumerable and are not only PGK3 and PYK3

The system SHALL be able to report, without running a simulation, every identifier
this module imports that appears in more than one source file, with each file's
value, and SHALL mark whether the values agree.

All five reactions SHALL appear in that report, not only the two the project's
notes previously recorded. Agreeing identifiers SHALL be reported as agreeing
rather than omitted, so that the report distinguishes "checked and consistent"
from "not checked".

#### Scenario: All five reactions appear as disagreements
- **WHEN** the ambiguity report is generated over this module's source tables
- **THEN** the forward catalytic constants of all five reactions appear, each
  marked as a disagreement
- **AND** the largest ratio among them is the pyrophosphatase constant, differing
  by more than two orders of magnitude

#### Scenario: Agreeing initial conditions are reported as agreeing
- **WHEN** the ambiguity report is generated
- **THEN** `M_atp_c`, `M_adp_c`, `M_pi_c` and `M_ppi_c` appear, marked as
  agreeing across the two files
- **AND** `M_amp_c`, `M_gtp_c`, `M_gdp_c` and `M_gmp_c` appear, marked as
  disagreeing

#### Scenario: Every governing choice is enumerable from the composed model
- **WHEN** the composed model is queried for values whose source file was chosen
  in the presence of an alternative
- **THEN** every value this module imports appears, with the file chosen and the
  file rejected

### Requirement: Informedness distinguishes a balanced value from a prior default

An imported value SHALL record whether it came from a balancing distribution or
sits at the prior median at prior width. This SHALL be derived from the source
table's own uncertainty rather than asserted per value.

The central file's values for these five reactions are not a competing estimate
but the signature of balancing with no data: their geometric standard deviations
run to many orders of magnitude, and all of their Michaelis constants sit at
exactly the prior median. The system SHALL classify them accordingly, so that the
reason the nucleotide file governs is visible in the data rather than only
asserted in a document.

`M_ppi_c` SHALL be marked as sitting at prior median and prior width, since
nothing informed it in either file, and its width SHALL be retained rather than
narrowed to look measured.

#### Scenario: The rejected values classify as uninformed
- **WHEN** the central file's rows for these five reactions are read
- **THEN** each classifies as a prior default rather than as a balanced value
- **AND** the nucleotide file's corresponding rows classify as balanced

#### Scenario: Pyrophosphate stays uninformed
- **WHEN** `M_ppi_c`'s initial concentration is read
- **THEN** it is marked as a prior default in whichever file governs
- **AND** its prior width is unchanged

#### Scenario: No parameter of this module is an asserted prior
- **WHEN** the composed model is queried for parameters whose prior this project
  asserted rather than inherited
- **THEN** none of this module's parameters appears

### Requirement: Both files' values are available in the repository

The rival values SHALL both be available from within the repository, so that the
ambiguity is real at load time rather than described in a document, and so that
the module can be constructed and tested with no network access and no checkout of
the source model.

Each in-repository copy SHALL record the upstream file and commit it derives from,
SHALL be regenerable by a re-runnable transformation, and SHALL change no value in
reshaping. Each SHALL be loaded as a distinct source, so that provenance names the
file a value actually came from.

#### Scenario: The module constructs with no external data
- **WHEN** the sub-model is constructed on a machine with no source-model checkout
- **THEN** construction succeeds with every value imported, its governing file
  declared and its alternative recorded

#### Scenario: The ambiguity is exercised, not described
- **WHEN** the two in-repository copies are loaded together
- **THEN** every identifier this module imports is held by both
- **AND** removing the second copy makes the governing declarations fail rather
  than pass unchecked
