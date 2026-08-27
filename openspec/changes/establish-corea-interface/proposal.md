# Establish the Core A′ interface contract

## Wave and dependencies

**Wave 0** of `dev/plans/reduced-syn3a-wave-plan.md`. Serial; blocks everything.

**Dependencies: None.** This is the first change in the Core A′ port. Every wave-1
change (`add-central-glycolysis`, `add-pts-transport`, `add-nucleotide-recycling`,
`add-corea-transcription`, `add-corea-translation`, `add-mrna-decay`,
`add-lumped-trna-charging`) depends on this one, and none of them may be proposed
until this change's delta specs are promoted with `/opsx:sync` and merged to `main`.

## Why

Core A′ is to be built by seven parallel branches that touch each other only
through a coupling boundary. Nothing currently names that boundary. `src/interface.jl`
gives `AbstractSubModel` a protocol — `states` / `parameters` / `dynamics` /
`reactions` / `inputs` / `formalism` / `inference_mode` — but `inputs` is a bare
`Vector{Symbol}`: a module can say *which* state it reads from another module and
nothing about *how*. That is not enough to fan out on. Two branches can both declare
`inputs = [:ATP]` while meaning entirely different things — one meaning shared
continuous mass inside the ODE block, the other meaning a cost counter the hook
debits one step later — and the disagreement will not surface until the join in
wave 2.

`dev/notes/figures/reduced-syn3a-coupling/fig1r_state_graph_reduced.pdf` already
draws the boundary, and its legend already enumerates the seven distinct semantics.
This change turns that legend into a protocol, so that a module's coupling
declaration is checkable before any of it runs.

Two further hazards are properties of the boundary rather than of any module, so
they belong here and nowhere else. The expression-cost drain in the published model
is a clipped subtraction — `max(0, ·)` on the interface — which will obstruct
gradients if the ODE block is sampled with NUTS. And the ODE→CME channel is
piecewise-constant: rate constants recomputed once every 60 s, not propensities
reading pools continuously. Both are choices the contract must let a module state
explicitly rather than assume.

Provenance is the third piece, and it is not bookkeeping neatness. The
central-versus-nucleotide cross-file trap has already produced two errors of record
in the scoping note: PGK3/PYK3 kinetics (319.5 vs 140.8; 1874 vs 672.8) and initial
concentrations (GTP at 1.6627 mM, not the 0.1 mM prior default that the central file
shows). Both were caught by review, not by the code. A loader that carries the source
file for every imported value makes the third instance visible instead of silent.

## What Changes

- **`AbstractSubModel` gains a coupling declaration.** A new `coupling(m)` protocol
  function returning typed edge declarations, alongside the existing `inputs(m)`.
  `inputs` keeps its current meaning and its empty default, so existing sub-models
  (`TranscriptionTranslation`, `StochasticGeneExpression`, `BurstyGeneExpression`,
  `LightMetabolism`, `TierBMetabolism`) continue to work unchanged. Not **BREAKING**.
- **A single Core A′ state registry.** The ~32 dynamic states and 5 chemostats of
  the scoping note's state list, named once and ordered once, with each entry
  carrying its group, its formalism-relevant copy-number regime, and whether it is
  dynamic or chemostatted. Modules index into this registry; no module invents its
  own species names.
- **Seven edge-kind declaration forms**, one per entry in fig1r's legend — `:mass`,
  `:currency`, `:deferred_counter`, `:catalytic`, `:rate_constant`, `:volume`,
  `:clamped` — each carrying the fields its own semantics require. The two
  boundary hazards become fields rather than hard-coded behaviour: a
  `:deferred_counter` declares its clipping policy (clamped-at-zero with deficit
  carried, or unclamped, or smoothed), and a `:rate_constant` declares its refresh
  cadence (piecewise-constant at an interval, or continuous).
- **A validating resolver.** Given a set of sub-models, resolve every declared edge
  against the registry and against the other modules' declarations, and fail with a
  named error when an edge references an unknown species, when two modules declare
  the same edge with different kinds, or when a declared cost has no state that can
  pay it.
- **A provenance-carrying parameter loader.** Every imported value — kinetic
  constant, initial concentration, prior — records the file it came from.
  `InferParameter` gains a provenance field; a query surfaces every value whose
  source file is ambiguous across the central and nucleotide balanced files.

## Capabilities

### New Capabilities

- `corea-interface/state-registry`: the canonical named and ordered Core A′ species
  registry — ~32 dynamic states and 5 chemostats — and the lookup and validation
  behaviour every module uses to index into it.
- `corea-interface/edge-kinds`: the seven coupling-edge declaration forms, the
  `coupling(m)` protocol extension that carries them, and the resolver that
  validates a composed set of modules' declarations against each other and against
  the registry.
- `corea-interface/parameter-provenance`: the parameter loader that records a source
  file for every imported value, the provenance field on `InferParameter`, and the
  cross-file ambiguity report.

### Modified Capabilities

None. `openspec/specs/` is currently empty; this change establishes the first specs
in the project.

## Impact

- **`src/interface.jl`** — extended, not rewritten. `AbstractSubModel` gains
  `coupling(m)` with an empty default; the seven edge-kind types and the registry
  types are added here or in new files included alongside it.
- **`src/parameters.jl`** — `InferParameter` gains a provenance field. Every
  existing construction site must be updated or reach a default; the five sub-models
  in `src/models/` and the tests that build parameters directly are affected.
- **New files** — a state registry, an edge-kind module, and a parameter loader,
  included from `src/InferCell.jl` with their public names exported.
- **`test/`** — new test files added to `test/runtests.jl`, run via
  `sbatch test/run_tests.slurm`.
- **Downstream** — after `/opsx:sync`, seven wave-1 changes consume these specs.
  Nothing else in the repository depends on them yet, so this is the cheapest point
  at which to get the contract wrong and fix it.

## Non-goals

- Importing syn3A parameter values through the loader. This change delivers the
  loader and the provenance mechanism; wave 1 modules import their own values
  through it. (The registry does transcribe the scoping note's initial
  concentrations by hand, each with its source file — the loader checks imports
  against those rows rather than replacing them.)
- Implementing any reaction, rate law, or module. Wave 0 is contract only.
- The orchestrator's use of the resolved edges to actually drive coupling at run
  time. That is wave 2 (`add-hook-1s-coupling`, `add-cme-rebuild-60s`,
  `add-growth-volume-coupling`). Wave 0 declares and validates; wave 2 executes.
