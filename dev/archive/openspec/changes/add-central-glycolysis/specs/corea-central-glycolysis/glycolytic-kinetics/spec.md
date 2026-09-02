## Purpose

The rate law each glycolytic reaction takes, the enzyme concentration that scales
it, and the import of every kinetic constant and initial condition with its source
file, its uncertainty and its informedness intact — so that the best-measured
region of Core A′ carries real inherited priors rather than asserted ones, and so
that the one place this port departs from what the published simulator runs is
visible rather than buried in a column choice.

## ADDED Requirements

### Requirement: Reactions take the published modular rate law

Each of the ten reactions SHALL be evaluated with the reversible modular rate law
the published model builds for central metabolism: a numerator differencing the
forward and reverse catalytic terms over the substrate and product saturation
ratios, divided by the summed saturation denominators less one, and scaled by the
enzyme concentration.

For a reaction with substrates $S_i$ and products $P_j$:

$$v = E \cdot \frac{k_{\text{cat}}^{f} \prod_i (S_i/K_{M,i}) - k_{\text{cat}}^{r} \prod_j (P_j/K_{M,j})}{\prod_i (1 + S_i/K_{M,i}) + \prod_j (1 + P_j/K_{M,j}) - 1}$$

The rate law SHALL be the same form for all ten reactions, differing only in how
many substrate and product terms it carries. The equilibrium constant and Hill
coefficient that appear in the source file's own kinetic-law strings SHALL NOT be
used, because the published model does not use them for these reactions either.

#### Scenario: The rate law form is uniform across the ten reactions
- **WHEN** each reaction's rate is evaluated
- **THEN** every one uses the modular form above, with term counts set by its own
  stoichiometry

#### Scenario: A reaction at equilibrium carries no net flux
- **WHEN** substrate and product concentrations are set so that the forward and
  reverse catalytic terms are equal
- **THEN** the reaction's net rate is zero

#### Scenario: Reverse flux is representable
- **WHEN** a reaction's products are raised far above its substrates
- **THEN** its net rate is negative

### Requirement: Enzyme concentration scales each rate and is replaceable

Each reaction's rate SHALL be scaled by the concentration of the single enzyme its
gene-protein-reaction rule names. The sub-model SHALL carry a nominal enzyme
concentration for each of the ten, derived from the published per-cell copy number
at the model's initial cell volume, and SHALL expose that concentration so that a
module supplying live protein counts can replace it without the rate law changing.

The nominal value SHALL NOT be the published model's default protein
concentration for reactions lacking a gene-protein-reaction rule: every reaction
here has one, so that default does not apply and using it would understate the
fluxes by more than an order of magnitude.

#### Scenario: Nominal enzyme concentrations follow the published copy numbers
- **WHEN** the nominal enzyme concentration of a reaction is read
- **THEN** it equals that enzyme's published copy number converted at the model's
  initial cell volume
- **AND** it is recorded as a nominal standing in for live protein counts, not as
  a measured concentration

#### Scenario: Enzyme concentration is replaceable without touching the rate law
- **WHEN** an enzyme concentration is overridden
- **THEN** every rate that enzyme catalyses scales in proportion
- **AND** no other reaction's rate changes

### Requirement: Every kinetic constant and initial condition is imported with its provenance

Every imported value SHALL be loaded through the Core A′ parameter loader rather
than transcribed into source, and SHALL carry the file it came from, the
identifier it had there, and its informedness.

The imported set SHALL be: a forward and a reverse catalytic constant for each of
the ten reactions, a Michaelis constant for each substrate and product of each
reaction, and an initial concentration for each of the thirteen states the
sub-model owns. Each SHALL carry the geometric standard deviation the source
records, as the width of its prior.

An imported initial condition SHALL agree with the registry's value for the same
species. Where an imported value sits at prior median and prior width rather than
at a balanced value, it SHALL be marked uninformed rather than presented as a
measurement.

#### Scenario: A kinetic constant reports its source
- **WHEN** an imported catalytic or Michaelis constant is read from the composed
  model
- **THEN** it reports the source file it came from and the identifier it had there

#### Scenario: Imported initial conditions agree with the registry
- **WHEN** the thirteen owned states' initial concentrations are imported
- **THEN** each equals the registry's recorded value for that species
- **AND** a disagreement fails the load rather than silently preferring one copy

#### Scenario: Uninformed initial conditions stay uninformed
- **WHEN** the initial concentrations of `M_g3p_c` and `M_lac__L_c` are imported
- **THEN** each is marked as sitting at prior median and prior width
- **AND** its prior width is retained rather than narrowed to look measured

#### Scenario: Priors on this module are inherited, not asserted
- **WHEN** the composed model is queried for parameters whose prior this project
  asserted rather than inherited
- **THEN** none of this sub-model's catalytic or Michaelis constants appears

### Requirement: Michaelis constants come from the balanced column, and the departure is labelled

The thirty-two Michaelis constants SHALL be taken from the source file's balanced
parameter table, which carries a geometric standard deviation for every value, and
NOT from the local-parameter table the published simulator reads, which carries
none.

Twenty-four of the thirty-two agree between the two tables. The eight that do not
SHALL be recorded as a departure from the published model's behaviour — the
largest of them changes a Michaelis constant by more than two orders of magnitude
— and that record SHALL be retrievable from the composed model alongside Core A′'s
other reduction declarations, so that it reaches any result depending on it.

#### Scenario: All thirty-two constants come from one table
- **WHEN** the Michaelis constants are read from the composed model
- **THEN** all thirty-two report the balanced parameter table as their source
- **AND** each carries a geometric standard deviation

#### Scenario: The eight disagreements are enumerable
- **WHEN** the composed model is queried for the declarations that are this
  reduction's rather than the published model's
- **THEN** the Michaelis-constant column choice appears, naming the eight
  constants that differ and both values for each

#### Scenario: The departure is not silent
- **WHEN** the sub-model's reduction notes are enumerated
- **THEN** the column choice is among them, stated as a departure from what the
  published simulator runs

### Requirement: The imported values are available without a source-model checkout

The values this sub-model imports SHALL be available from within the repository,
so that the sub-model can be constructed, integrated and tested on a machine with
no network access and no checkout of the source model.

The in-repository copy SHALL record which upstream file and which upstream commit
it was derived from, and SHALL be reproducible: the transformation from the
upstream file to the in-repository form SHALL be re-runnable rather than
hand-made.

Reshaping for the loader's format SHALL NOT change any value. Where an identifier
is renamed to fit the loader's conventions, the upstream identifier SHALL remain
recoverable.

#### Scenario: The sub-model constructs with no external data
- **WHEN** the sub-model is constructed on a machine with no source-model checkout
- **THEN** construction succeeds with every value imported and its provenance
  recorded

#### Scenario: The in-repository copy is reproducible
- **WHEN** the transformation is re-run against the recorded upstream file at the
  recorded commit
- **THEN** it reproduces the in-repository copy exactly

#### Scenario: Reshaping preserves values
- **WHEN** each in-repository value is compared to the upstream file's value for
  the same quantity
- **THEN** the two agree
