# Tasks — establish-corea-interface

Groups are internal waves. Tasks within a group are independently implementable in
any order; each group depends on the ones before it. Tests run via
`sbatch test/run_tests.slurm` — Julia is only available on compute nodes.

## 1. Registry

- [x] 1.1 Add `src/corea/registry.jl` defining `SpeciesEntry` (name, group,
  treatment, regime, initial value, source file, informedness) and the ordered
  `const COREA_SPECIES` in scoping-note group order — 11 glycolytic intermediates,
  4 adenylate, 3 guanylate, 2 redox, 2 other, 2 tRNA, 8 PTS phospho-states, then
  5 chemostats. Verify by a test asserting 32 dynamic states, 5 chemostats, and
  each per-group cardinality against the note's table, and asserting that neither
  H2O nor H+ is present.
- [x] 1.2 Build the name→index map once at load and expose lookup in both
  directions. Verify by a round-trip test over every registered name, and a test
  that an unregistered symbol raises an error naming the symbol rather than
  returning a default.
- [x] 1.3 Populate initial values and their source files, taking the nucleotide
  balanced file for GTP (1.6627 mM), GDP (0.2981), AMP (0.0832) and GMP (0.0117),
  and the central file for the glycolytic and adenylate entries. Verify by a test
  asserting GTP resolves to the nucleotide value and names the nucleotide file, not
  the central file's 0.1 mM default.
- [x] 1.4 Mark PPi, cytosolic lactate and G3P as sitting at prior median and prior
  width. Verify by a test that these three report as uninformed and the balanced
  entries do not.
- [x] 1.5 Record chemostat held values — external glucose 40 mM, amino-acid pool
  0.1 mM, and CTP, UTP, O2. Verify by a test reading each held value from the
  registry.
- [x] 1.6 Include `registry.jl` from `src/InferCell.jl` and export the registry
  type, the constant and the lookup functions. Verify the package precompiles and
  the existing suite still passes.

## 2. Edge kinds

- [x] 2.1 Add `src/corea/edges.jl` with `abstract type CouplingEdge` and the seven
  concrete kinds — mass, currency, deferred counter, catalytic, rate constant,
  volume, clamped — each with the fields its own semantics require per the design's
  Decision 2. Verify by a test constructing one of each and a test that omitting a
  required field fails at construction naming the field and the kind.
- [x] 2.2 Give `DeferredCounterEdge` a `clip` field over
  `:clamped_deficit_carried` / `:unclamped` / `:smoothed`, defaulting to
  `:clamped_deficit_carried`. Verify by a test that the default matches the
  published model and that a clamped edge reports itself as non-differentiable at
  the boundary.
- [x] 2.3 Give `RateConstantEdge` a `cadence` field — piecewise-constant with an
  interval, or continuous — defaulting to piecewise-constant at 60 s. Verify by a
  test that the default is 60 s piecewise-constant and that a continuous
  declaration is marked as deviating from the published model.
- [x] 2.4 Give `ClampedEdge` an `origin` field distinguishing the published
  model's clamp from this reduction's. Verify by a test that a chemostatted CTP
  declaration is labelled as this reduction's.
- [x] 2.5 Add `coupling(::AbstractSubModel) = CouplingEdge[]` to
  `src/interface.jl`, leaving `inputs` untouched. Verify by a test that each of the
  five existing sub-models composes unchanged and reports empty coupling.

## 3. Resolver

- [x] 3.1 Add `resolve_coupling(models)` returning a `CouplingGraph` of resolved
  edges with kinds, endpoints and registry positions, callable without building an
  `ODEProblem`. Verify by a test resolving a consistent two-module composition and
  reading back each edge's kind and endpoints.
- [x] 3.2 Reject an edge naming an unregistered species, and an edge whose
  counterpart module is absent from the composition, distinguishing the latter from
  a deliberate single-module composition under test. Verify by one test per case
  asserting the error names the species or the missing counterpart.
- [x] 3.3 Reject the same species-and-direction declared by two modules with
  different kinds, with an error naming both modules, the species and both kinds
  and distinguishing the two semantics. Verify by a test using mass versus deferred
  counter, asserting on the error text.
- [x] 3.4 Reject a module that declares dynamics for a registry chemostat, and a
  dynamic state owned by two modules; report a dynamic state owned by none in a
  form that a single-module composition does not trip over. Verify by one test per
  case.
- [x] 3.5 Implement the dead-end check over declared direction plus registry
  groups: a cost debited against a species no module integrates and the registry
  does not chemostat is an error; a consumed dynamic species with no declared
  producer is reported with its conserved moiety named. Verify by a test that AMP
  consumed with no producer and GMP consumed with no producer both report, naming
  adenylate and guanylate respectively.
- [x] 3.6 Exempt chemostatted species from the dead-end check and record the
  exemption so a later change making the pool live reinstates it. Verify by a test
  that a cost against CTP raises nothing and that the exemption is enumerable.
- [x] 3.7 Exclude catalytic edges from mass and moiety accounting, and reject an
  attempt to include one in a conservation check. Verify by a test asserting the
  rejection rather than a silent zero contribution.
- [x] 3.8 Check that each module's `inputs` is a subset of the species named by its
  inbound coupling edges, so a typed declaration that drifts from `inputs` is
  caught. Verify by a test with a module declaring an inbound edge and an
  inconsistent `inputs`.
- [x] 3.9 Call `resolve_coupling` from the orchestrator's existing
  `_resolve_coupling` path so composition validates the boundary. Verify the full
  existing suite still passes via `sbatch test/run_tests.slurm`.

## 4. Parameter provenance

- [x] 4.1 Add `ParameterSource` (file, table, identifier, informedness,
  alternatives) and add `provenance::Union{ParameterSource,Nothing}` to
  `InferParameter`, retaining the six-positional-argument constructor defaulting it
  to `nothing`. Verify by running the existing `test/test_parameters.jl` unchanged
  and a new test that a directly constructed parameter reports provenance as absent.
- [x] 4.2 Carry provenance through `unique_params` and the orchestrator's parameter
  deduplication, and extend `_validate_shared_params` to report the case where two
  modules give one name equal values from different source files. Verify by a test
  composing two modules sharing a parameter name from two files.
- [x] 4.3 Add `src/corea/loader.jl` reading a balanced table from a caller-supplied
  path and returning provenance-tagged `InferParameter`s. Verify against a fixture
  under `test/fixtures/` by a test asserting each loaded parameter names its file
  and identifier.
- [x] 4.4 Implement the cross-file ambiguity report, listing every identifier found
  in more than one source file with each file's value and marking disagreements
  distinctly from agreements. Verify against a fixture carrying the PGK3 and PYK3
  cases (319.5 vs 140.8; 1874 vs 672.8) and the GTP initial-concentration case,
  asserting all three are reported as disagreements and that the report is empty
  for a single-file load.
- [x] 4.5 Require a governing-file declaration for an identifier present in more
  than one file, failing at load naming the identifier, both files and both values
  when it is absent, and recording the chosen and rejected files when it is
  present. Verify by one test per branch, the second asserting PGK3 loads the
  nucleotide value and its provenance names that file.
- [x] 4.6 Mark loaded values as balanced, prior-default, asserted, or not-imported,
  and make the last two enumerable from a composed model. Verify by a test that a
  fixture PTS mass-action constant enumerates as asserted-by-us and a balanced
  glycolytic parameter does not.
- [x] 4.7 Include `loader.jl` from `src/InferCell.jl` and export
  `ParameterSource`, the loader and the ambiguity report. Verify the package
  precompiles and the suite passes.

## 5. Labelling what is ours

- [x] 5.1 Implement the composed-model query returning every declaration this
  reduction introduced rather than the published model: reduction-origin clamps,
  smoothed deferred counters, continuous rate-constant edges, and parameters whose
  prior is asserted by us. Verify by a test on a composition carrying one of each,
  asserting all four appear.
- [x] 5.2 Add the lumped tRNA charging step to that query's result set as a
  declarable marker, so wave 1's `add-lumped-trna-charging` has somewhere to
  register it. Verify by a test that a module carrying the marker is enumerated.

## 6. Wire-up and hand-off

- [x] 6.1 Add the new test files to `test/runtests.jl` and confirm the whole suite
  passes on a compute node via `sbatch test/run_tests.slurm`, checking the emitted
  log rather than assuming.
- [x] 6.2 Add the two `rules` blocks from `dev/plans/reduced-syn3a-wave-plan.md` —
  the proposal Dependencies-and-wave rule and the tasks numbered-waves-and-balance
  rule — to `openspec/config.yaml`. Verify by running `openspec instructions
  proposal --change establish-corea-interface --json` and confirming the rules
  appear in the response.
- [x] 6.3 Write a short reference page documenting the registry, the seven edge
  kinds and the loader, pointing at `fig1r_state_graph_reduced.pdf` as the source
  of the edge-kind legend. Verify the page builds under `mkdocs` and is linked from
  the docs nav.
- [x] 6.4 Run `openspec validate establish-corea-interface` and confirm it passes.
- [x] 6.5 Update `docs/handoff.md` with what landed and the fact that wave 1 is
  unblocked only after `/opsx:sync establish-corea-interface` and the PR merge to
  `main`.
