## Purpose

How each of the seventeen rate constants is computed — from transcript length,
base composition, NTP pools and a promoter strength borrowed from proteomics —
and the two things about that computation which must not be silently inherited:
a base-to-NTP mapping error, and the circularity of the promoter proxy.

## ADDED Requirements

### Requirement: Rate constants follow the published formula

Each gene's transcription rate constant SHALL be computed as the published model
computes it: a promoter-scaled polymerase turnover divided by a denominator
combining polymerase binding, per-base nucleotide saturation, and transcript
length.

The nucleotide term SHALL sum, over the four bases, the count of that base in the
transcript times the polymerase dissociation constant divided by that base's NTP
concentration. The turnover SHALL be capped at the published ceiling, and the
sub-model SHALL record that the cap does not bind for any of these seventeen
genes.

#### Scenario: A rate constant depends on transcript length
- **WHEN** two genes with equal promoter strength but different transcript
  lengths are compared
- **THEN** the longer transcript has the smaller rate constant

#### Scenario: A rate constant depends on the NTP pools
- **WHEN** an NTP concentration is lowered
- **THEN** every rate constant whose transcript contains that base decreases

#### Scenario: The turnover cap is recorded as non-binding
- **WHEN** the sub-model reports its rate-constant diagnostics
- **THEN** it records that no gene's promoter-scaled turnover reaches the
  published ceiling

### Requirement: Each base is charged against its own NTP, and the correction is labelled

Each base's count SHALL enter the nucleotide term against the concentration of the
nucleotide that base is polymerised from: adenine against ATP, cytosine against
CTP, guanine against GTP, uracil against UTP.

This differs from the published implementation, which charges cytosine against
UTP, guanine against CTP and uracil against GTP — an ordering error present on
both its code paths, including the periodic rate-constant rebuild. The correction
SHALL be registered as a departure from the published model and SHALL be
retrievable from the composed model, and the published mapping SHALL remain
available so the difference can be measured rather than argued.

The sub-model SHALL record why the correction matters: the effect on a rate
constant's value is about one percent, but the published mapping weights the
guanosine term by the uracil count, making every rate constant roughly twice as
sensitive to the live GTP pool as it should be — and that pool is the model's only
inbound coupling from the metabolic block.

#### Scenario: Each base charges its own nucleotide
- **WHEN** a transcript's cytosine count is increased and all else held
- **THEN** only the CTP term of its rate constant changes

#### Scenario: The correction is labelled as ours
- **WHEN** the composed model is queried for the declarations that are this
  reduction's rather than the published model's
- **THEN** the base-to-nucleotide correction appears, naming both mappings and the
  effect on the coupling

#### Scenario: The published mapping remains measurable
- **WHEN** the published mapping is selected
- **THEN** the rate constants change by a small amount and the sensitivity to the
  guanosine triphosphate pool roughly doubles
- **AND** the selection is recorded, so a trajectory computed under it is
  identifiable as such

### Requirement: Rate constants are recomputable from live pools, not baked in

Every rate constant SHALL be recomputable from a supplied set of NTP
concentrations, without reconstructing the sub-model. The recomputation SHALL be
exposed as an operation a later change can call, and the current values SHALL be
readable.

This module SHALL NOT itself schedule or perform any refresh. It declares the
cadence at which its constants are rebuilt and provides the means; the change that
owns the rebuild performs it.

#### Scenario: Recomputation changes the constants
- **WHEN** recomputation is invoked with NTP concentrations differing from the
  current ones
- **THEN** every rate constant is updated accordingly
- **AND** no other property of the sub-model changes

#### Scenario: The module does not refresh itself
- **WHEN** the sub-model is simulated over an interval longer than its declared
  refresh cadence
- **THEN** its rate constants are unchanged, because nothing in this module
  performs the rebuild

#### Scenario: Constants are readable
- **WHEN** the current rate constants are read
- **THEN** all seventeen are returned, keyed by locus

### Requirement: Promoter strength is a proteomics proxy, and its circularity is recorded

Each gene's promoter strength SHALL be its measured protein copy number scaled by
the published divisor, inherited unchanged as the baseline against which a later
replacement is measured.

The sub-model SHALL record two things about it. First, that it is a proxy: the
seventeen values re-encode proteomics rather than measuring promoter activity, so
recovering them from synthetic data is evidence about the machinery and not about
the biology. Second, that conditioning an inference on measured protein abundance
would be circular, because protein abundance is the input to these constants.

Both records SHALL be retrievable from the composed model, so that the constraint
is visible where an observable is chosen rather than only in a design note.

#### Scenario: Promoter strength follows the protein count
- **WHEN** two genes' promoter strengths are compared
- **THEN** their ratio equals the ratio of their protein copy numbers

#### Scenario: The proxy is flagged as a proxy
- **WHEN** the composed model is queried for this module's recorded caveats
- **THEN** the promoter-strength proxy appears, stating that it re-encodes
  proteomics

#### Scenario: The circularity warning is retrievable
- **WHEN** the same query is made
- **THEN** an entry states that conditioning on measured protein abundance would
  recover this module's own input
- **AND** it names the parameters affected

### Requirement: Per-gene data is available without a source-model checkout

The transcript length, base composition, protein copy number and measured mean
transcript count of each of the seventeen genes SHALL be available from within the
repository, so the sub-model can be constructed and simulated with no network
access and no checkout of the source model.

The in-repository copy SHALL record which upstream files and which upstream commit
it derives from, and SHALL be regenerable by a re-runnable transformation. Base
counts SHALL sum to the transcript length for every gene.

#### Scenario: The sub-model constructs with no external data
- **WHEN** the sub-model is constructed on a machine with no source-model checkout
- **THEN** construction succeeds with all seventeen genes' data present

#### Scenario: Base composition is internally consistent
- **WHEN** each gene's four base counts are summed
- **THEN** the total equals that gene's recorded transcript length

#### Scenario: The copy numbers agree with the metabolic modules
- **WHEN** a gene's protein copy number is compared with the count the metabolic
  module using that locus derives its enzyme concentration from
- **THEN** the two agree
