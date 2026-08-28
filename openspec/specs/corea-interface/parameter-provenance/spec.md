# Core A′ parameter provenance Specification

## Purpose

A parameter loader that records, for every imported value, the file it was read
from — so that the cross-file inconsistencies in the source model are visible in
the code rather than resolved once by hand and forgotten.

## Requirements

### Requirement: Every imported value carries its source

A parameter declared by a module SHALL be able to carry provenance: the source file
the value was read from, the table or column within it where that is meaningful,
and the identifier under which it appeared there. Provenance SHALL travel with the
parameter through composition, so that a value in an assembled model can be traced
back to its file without re-reading the source.

Provenance SHALL be optional, so that parameters constructed directly — as the
existing sub-models and tests do — remain valid.

#### Scenario: A loaded parameter reports its source
- **WHEN** a parameter is imported through the loader and later read from a
  composed model
- **THEN** it reports the file it came from and the identifier it had there

#### Scenario: Existing parameter construction is unaffected
- **WHEN** a parameter is constructed directly without provenance, as existing
  sub-models do
- **THEN** construction succeeds
- **AND** the parameter reports its provenance as absent rather than as an
  incorrect source

#### Scenario: Provenance survives composition
- **WHEN** modules carrying provenance-tagged parameters are composed, including
  where a parameter name is shared across modules and deduplicated
- **THEN** the surviving parameter still reports a source
- **AND** where two modules supplied the same name from different files, the
  disagreement is reported rather than resolved by taking the first

### Requirement: Values that appear in more than one source file are reported

The loader SHALL be able to report every imported value whose identifier appears in
more than one source file, together with each file's value, whether or not the
values agree. This report SHALL be obtainable without running a simulation.

The cross-file trap this exists for has already produced two errors of record:
kinetics, where PGK3 and PYK3 are parameterised in both the central and nucleotide
balanced files with different modes, and initial concentrations, where GTP reads
1.6627 mM in the nucleotide file and the 0.1 mM prior default in the central file.

#### Scenario: A value present in two files is flagged
- **WHEN** the ambiguity report is generated over a set of imported values
- **THEN** every identifier found in more than one source file appears, with each
  file's value alongside it
- **AND** an identifier found in exactly one file does not appear

#### Scenario: Disagreeing values are distinguished from agreeing ones
- **WHEN** the same identifier carries different values in two files
- **THEN** the report marks it as a disagreement, distinctly from an identifier
  that appears twice with the same value

#### Scenario: The report is empty when it should be
- **WHEN** every imported value comes from exactly one file
- **THEN** the report is empty, and its emptiness is a positive result rather than
  a failure to run

### Requirement: A module states which file governs each value it imports

A module importing a value that exists in more than one source file SHALL state
which file governs. An unresolved ambiguity SHALL be an error at load time, not a
silent choice.

#### Scenario: An ambiguous import without a declared governing file
- **WHEN** a module imports an identifier present in two source files without
  stating which governs
- **THEN** loading fails naming the identifier, both files, and both values

#### Scenario: A declared governing file resolves the ambiguity
- **WHEN** a module imports PGK3's forward kinetic constant and declares the
  nucleotide file as governing
- **THEN** loading succeeds with the nucleotide file's value
- **AND** the parameter's provenance names that file, so the choice is visible in
  the assembled model rather than only in the module's source

#### Scenario: Governing-file declarations are enumerable
- **WHEN** a composed model is queried for every value whose source file was chosen
  in the presence of an alternative
- **THEN** each such value is returned with the file chosen and the file rejected

### Requirement: Values without an informed source are marked, not defaulted silently

Where an imported value is a prior median at prior width rather than a balanced or
measured value, the loader SHALL mark it as uninformed. An uninformed value SHALL
be distinguishable from an informed one through the parameter alone.

#### Scenario: An uninformed value is marked
- **WHEN** a value that sits at the prior default is imported — such as PPi or
  lactate, at 0.1 mM with geometric standard deviation 10
- **THEN** the parameter is marked uninformed
- **AND** its prior width is retained rather than being narrowed to look
  measured

#### Scenario: Values with no prior at all are marked separately
- **WHEN** a value with no quantified uncertainty in the source is imported — such
  as a PTS mass-action constant
- **THEN** the parameter records that any prior on it is asserted by this project
  rather than inherited from the source model
- **AND** this is distinguishable from a value that carries a wide but genuine
  prior

#### Scenario: Enumerating the asserted priors
- **WHEN** a composed model is queried for every parameter whose prior this project
  asserted rather than inherited
- **THEN** the ten PTS mass-action constants appear, and the balanced glycolytic
  parameters do not
