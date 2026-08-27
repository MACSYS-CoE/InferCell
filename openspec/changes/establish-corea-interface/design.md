# Design — Core A′ interface contract

## Context

See `proposal.md` — Why, for motivation. The constraints that shape the approach:

- **`src/interface.jl` is not greenfield.** `AbstractSubModel` already carries
  `states` / `parameters` / `dynamics` / `reactions` / `inputs` / `formalism` /
  `inference_mode`, and five concrete sub-models implement it
  (`TranscriptionTranslation`, `StochasticGeneExpression`, `BurstyGeneExpression`,
  `LightMetabolism`, `TierBMetabolism`), with tests that construct
  `InferParameter` positionally. Wave 0 extends this protocol; it does not replace
  it, and it must not break the ISAB work that runs on it.
- **`_resolve_coupling` in `src/orchestrator.jl` is the existing resolver.** It
  builds a `state_owners` dictionary, errors on a doubly-owned state, and maps each
  module's `inputs` symbols to global state indices. Its shape is the right one;
  it is the semantics that are missing.
- **Seven branches will consume this in parallel.** The contract has to be
  checkable by a module author working alone in a worktree, before any other
  module exists. A validation that only fires when all seven are composed is a
  validation that fires too late to be useful.
- **The figure is the source of truth for the edge kinds**, not this design. The
  seven kinds and their one-line semantics come from the legend built in
  `dev/notes/figures/reduced-syn3a-coupling/make_reduced.py`, rendered in
  `fig1r_state_graph_reduced.pdf`. If the figure and the code disagree, one of
  them is a bug.
- **Tests run through Slurm** (`sbatch test/run_tests.slurm`); Julia is only on
  compute nodes, and compute nodes have no network.

## Goals / Non-Goals

**Goals:**

- A module author can write and validate a coupling declaration with no other
  module present.
- Every failure the specs describe is raised at declaration or at composition,
  never at integration time and never silently.
- The two boundary hazards — the clipped drain, the piecewise-constant rebuild —
  become declared fields with the published model as the default, so a deviation
  from the published model is a thing someone has to write down.
- Additive change only. Every existing sub-model and test compiles and passes
  untouched.

**Non-Goals:**

- Executing coupling. The resolver produces a validated description; wave 2 turns
  it into a callback that actually moves numbers.
- Mixed-formalism composition. `build_problem` still errors on it; this change
  describes the boundary that a later change will teach it to cross.
- Vendoring syn3A data files. See Decision 7.

## Decisions

### 1. A new `coupling(m)` protocol function, not a change to `inputs(m)`

`coupling(::AbstractSubModel) = CouplingEdge[]`, added beside the existing
`inputs(::AbstractSubModel) = Symbol[]`.

*Alternative considered:* widen `inputs` to return typed edges. Rejected — it
changes the return type of a function five sub-models and several test files
already implement, for no benefit that a second function does not give. The two
coexist: `inputs` stays the untyped "which states do I read" used by the current
ODE right-hand-side assembly, and `coupling` carries the semantics. Where a species
appears in both, the typed declaration governs, and the resolver checks that
`inputs` is a subset of the species named by inbound `coupling` edges — so a module
that adds a typed declaration and forgets to keep `inputs` in step is told.

### 2. Seven concrete edge types under one abstract supertype

```
abstract type CouplingEdge end

struct MassEdge            <: CouplingEdge   # species, direction, peer
struct CurrencyEdge        <: CouplingEdge   # + pool
struct DeferredCounterEdge <: CouplingEdge   # + counter, clip, smoothing
struct CatalyticEdge       <: CouplingEdge   # + param_slot
struct RateConstantEdge    <: CouplingEdge   # + cadence, interval
struct VolumeEdge          <: CouplingEdge
struct ClampedEdge         <: CouplingEdge   # + held_value, origin
```

*Alternative considered:* one struct with `kind::Symbol` and a `Dict{Symbol,Any}`
of extras. Rejected on the spec requirement that a missing per-kind field fails at
declaration: with a dictionary that check has to be hand-written and re-run, while
with seven structs Julia's constructor does it and the resolver dispatches on type
rather than branching on a symbol. The cost is seven exported names, which is the
right cost for a contract seven branches read.

The kinds are closed. `CouplingEdge` is not meant to be subtyped outside this file,
and the resolver's dispatch will fail on an unrecognised subtype rather than fall
through to a generic method — which is how an eighth kind announces itself.

### 3. Clipping policy and refresh cadence are fields, not kinds

`clip` is one of `:clamped_deficit_carried` (the published model), `:unclamped`,
`:smoothed`; `cadence` is `:piecewise_constant` with an interval, or `:continuous`.

*Alternative considered:* promote each variant to its own edge kind. Rejected —
it turns seven kinds into eleven and breaks the one-to-one correspondence with
fig1r's legend, which is the thing that makes the contract auditable against a
figure someone drew on purpose. The variants are the same edge with different
numerics, so they belong inside it.

**The default is the published model.** A `DeferredCounterEdge` constructed without
a clip is `:clamped_deficit_carried`; a `RateConstantEdge` without a cadence is
`:piecewise_constant` at 60 s. Deviating is possible, cheap, and requires typing
the deviation — which is exactly the balance the scoping note's open questions ask
for, since neither hazard has an answer yet that is better than "reproduce the
published model and label any departure".

### 4. The registry is a `const` ordered vector plus a name index

A `const COREA_SPECIES::Vector{SpeciesEntry}` in source order, with a
`Dict{Symbol,Int}` built once at load. Source order *is* canonical order, grouped
as the scoping note groups them: glycolytic intermediates, adenylate, guanylate,
redox, other, tRNA, PTS phospho-states, then chemostats last.

*Alternative considered:* sort alphabetically, or derive the order from a data
file. Both rejected. Alphabetical makes groups non-contiguous and makes a rename
silently reorder every module's indices. Deriving from a file adds a load-time
failure mode to something that has to be available for a module author to compile
against.

Chemostats live in the same registry rather than a separate one, marked by a
`treatment` field, because the check that matters — "a module tried to integrate a
chemostat" — needs both in one namespace to be a lookup rather than a
cross-reference.

### 5. `resolve_coupling(models)` is a standalone pass returning a `CouplingGraph`

`_resolve_coupling` inside `_build_ode_problem` keeps doing what it does; the new
pass is callable on its own and returns a description, so a module author can
validate without building an `ODEProblem` and so the specs' "obtainable without
running a simulation" reports are obtainable.

*Alternative considered:* fold everything into `_build_ode_problem`. Rejected —
it couples validation to problem construction, which means a `:jump` module or a
single module under test cannot be validated at all.

### 6. The dead-end check runs on declared direction plus registry groups

Producers and consumers come from edge direction; "which conserved moiety" comes
from the registry's `group` field, which is already adenylate / guanylate / redox
and so on. So the error message for AMP with no ADK1 and the error message for GMP
with no GK1 are the same code path, which is the point — the scoping note found
those two failures a review round apart and they are one class.

This checks *declarations*, not stoichiometry. A module that declares a producer of
AMP and then does not actually produce any will pass. That is a real limit, and it
is why the wave plan puts a conservation check inside each wave-1 change rather
than trusting wave 0 to catch it. Wave 0 catches the structural dead end; the
numerical one is caught by running the module.

### 7. `InferParameter` gains a seventh field with a defaulting outer constructor

```
struct ParameterSource
    file::String
    table::Union{String,Nothing}
    identifier::Union{String,Nothing}
    informedness::Symbol   # :balanced, :prior_default, :asserted, :not_imported
    alternatives::Vector{Pair{String,Float64}}   # same id, other files
end
```

added to `InferParameter` as `provenance::Union{ParameterSource,Nothing}`, with the
existing six-positional-argument constructor retained and defaulting it to
`nothing`. Every current construction site keeps working with no edit.

*Alternatives considered:* a side table keyed by `(name, module_id)` — rejected
because it desynchronises the moment a parameter is copied or deduplicated, and
`unique_params` and `_build_p0` both do exactly that; a wrapper type — rejected
because it forces an unwrap at every read site in `parameters.jl`,
`orchestrator.jl` and `inference.jl`.

No performance concern: `InferParameter` already holds an abstract `Distribution`
field, so it is not `isbits` today, and it never enters the hot path — the
integrator sees a `Vector{Float64}` built by `_build_p0`.

`_validate_shared_params` already errors when two modules give one name
inconsistent values; it gains the case where the values agree but the source files
differ, which is the cross-file trap in its most easily-missed form.

### 8. The loader takes a path; no data files are vendored in wave 0

**Assumption, recorded because it is load-bearing:** the syn3A balanced tables are
not in this repository and compute nodes have no network, so the loader is
format-aware but path-agnostic — it reads a table from a path the caller supplies
and returns provenance-tagged parameters. Wave 0 ships it with a small fixture
under `test/` exercising the cross-file ambiguity path, using the PGK3/PYK3 and GTP
cases as the fixture's content. Vendoring the real tables belongs to the first
wave-1 change that needs them.

If this assumption is wrong — if the intent is for wave 0 to vendor the tables —
the specs do not change, but a task is added and the fixture becomes the real file.

## Risks / Trade-offs

- **The registry hard-codes a count the scoping note owns.** If the note's state
  list changes, the code silently disagrees with it. → A test asserts the totals
  and the per-group cardinalities (11 / 4 / 3 / 2 / 2 / 2 / 8 / 5) against the
  note's table, so drift fails the suite instead of accumulating.
- **Declaration-level validation can be satisfied without being true.** A module
  can declare a producer it does not implement. → Accepted, and mitigated
  elsewhere: each wave-1 change carries its own balance check as a task, per the
  wave plan.
- **Seven parallel branches can still agree on a wrong contract.** Freezing the
  interface before any module is written trades late conflicts for the chance that
  the frozen thing is wrong. → Mitigated by proposing all seven wave-1 changes
  before applying any of them, which is the wave plan's stated reason for doing so;
  and by `/opsx:sync` rather than `/opsx:archive`, so wave 0 stays open for
  amendment while wave 1 reads it.
- **Adding a field to `InferParameter` touches a type used everywhere.** → The
  defaulting constructor keeps every existing call site valid; the test suite is
  the check, and it runs before the PR.
- **The `:smoothed` clip policy is a modelling decision dressed as a field.**
  Offering it makes it easy to select a model that is not the published one. →
  The label is mandatory and enumerable: `Requirement: Clamped edges record
  whether the clamp is ours` and the "enumerating what is ours" scenario exist so
  that a deviation cannot reach a result unlabelled.

## Migration Plan

Additive throughout; there is nothing to migrate. Existing sub-models gain the
empty `coupling` default and are otherwise untouched. The change lands on a feature
branch and merges to `main` by PR, and `/opsx:sync` promotes the delta specs before
any wave-1 change is proposed — the wave plan is explicit that a wave-1 proposal
written before the sync will invent its own interface instead of consuming this
one.

Rollback is deleting the new files and reverting one field on `InferParameter`,
for as long as no wave-1 change has merged. After that, rollback means reworking
seven modules, which is the reason to get the review of this change right.

## Open Questions

Genuinely deferrable — none of these changes the specs, the approach, or the task
breakdown:

- **Which clip policy Core A′ actually runs under for inference.** The contract
  makes all three declarable and defaults to the published one. Whether NUTS on the
  ODE block needs `:smoothed` is answered by trying it, in wave 3.
- **Whether the 60 s rebuild stays piecewise-constant.** Same shape: declarable,
  defaults to published, and the scoping note asks for the deviation to be measured
  rather than assumed. That measurement is a wave-2 task on
  `add-cme-rebuild-60s`.
- **The on-disk format of the vendored balanced tables** (SBtab TSV as shipped, or
  a converted form). The loader's interface is the same either way; only the parser
  behind it differs.
