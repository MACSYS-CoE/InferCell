# Handoff

**Session date:** 2026-08-27
**Branch:** `change/establish-corea-interface` (PR #38, open against `main`)

## What this session did

Implemented OpenSpec change `establish-corea-interface` — wave 0 of the Core A′
port, per `dev/plans/reduced-syn3a-wave-plan.md`. The change was proposed and
applied in the same session, so `openspec/changes/establish-corea-interface/`
carries the proposal, three delta specs, design and tasks alongside the code.

Wave 0 is contract only: it declares and validates a boundary, and executes
nothing. Wave 2 is where the resolved edges start driving actual coupling.

### New source

| File | What it holds |
|---|---|
| `src/corea/registry.jl` | The 32 dynamic states and 5 chemostats, named and ordered once |
| `src/corea/edges.jl` | The seven `CouplingEdge` kinds from fig1r's legend |
| `src/corea/resolver.jl` | `resolve_coupling`, its checks and its reports |
| `src/corea/loader.jl` | Provenance-carrying parameter import, cross-file ambiguity report |
| `src/corea/labels.jl` | `reduction_declarations` — enumerating what is ours, not the model's |

### Changed source

- `src/parameters.jl` — added `ParameterSource` and a seventh `provenance` field
  on `InferParameter`, behind a six-argument outer constructor so all 54
  existing construction sites (47 in `src/models/`, 7 in `test/`) work
  untouched.
- `src/interface.jl` — added `coupling`, `module_id` and `reduction_notes`, all
  with defaults. `inputs` is unchanged.
- `src/orchestrator.jl` — `_resolve_coupling` now calls `resolve_coupling`;
  `_validate_shared_params` warns when one parameter name arrives from two files.
- `src/InferCell.jl` — includes and exports for the above.

### Two decisions worth knowing

**The boundary hazards became fields, not hard-coded behaviour.** The clipped
`max(0, ·)` ATP drain and the 60 s piecewise-constant rebuild are both
declarable, and both default to what the published model does. Departing from
the published model is now something an author has to type, and
`reduction_declarations` enumerates every departure so it can reach a report.
This is how the scoping note's two open questions are handled without answering
them prematurely.

**The loader parses SBtab-shaped TSV with Base only.** `DelimitedFiles` stopped
being a stdlib at Julia 1.9, so using it would mean a new `[deps]` entry, a
`[compat]` bound, and a network fetch. Compute nodes have no network. The format
needs no quoting logic, so a small Base parser was the lower-risk choice.

## Two pre-existing problems fixed in passing

Both were broken on `main` before this branch, and both blocked this change's
own verification steps:

1. **The test suite could not run on this cluster at all.** `Pkg.test()` resolves
   a fresh test environment needing Aqua, Aqua was not in the Julia depot, and
   compute nodes have no network — so every `sbatch test/run_tests.slurm` died in
   ~20 s with `Could not resolve host: pkg.julialang.org`. Fixed by populating
   the depot from the login node, which does have network and does have Julia at
   `/apps/modules/software/Julia/1.10.5/bin/` (the wrapper needs `EBROOTJULIA`
   set and its own `bin` on `PATH`). **If the depot is wiped, this recurs** — see
   `reference_cluster_env` in auto-memory.
2. **`mkdocs build --strict` was failing.** `docs/positioning.md` was deleted in
   `390aa22` but was still referenced from `mkdocs.yml`'s nav *and* from a link
   in `docs/index.md`. Under `strict: true` either one aborts the build, so
   `.github/workflows/docs.yml` was red. Both references removed.

## Status

- `openspec validate establish-corea-interface --strict` — passes.
- `mkdocs build --strict` — passes, with the new `docs/api/corea-interface.md`
  rendering. Verified in a scratch venv, since mkdocs is not installed on the
  cluster.
- Julia test suite — **759 passed, 0 failed, 0 errored** (job 15972827, ~3m
  wall; 701 before the review-driven hardening, 269 before this change). All
  34 tasks in `tasks.md` are ticked.

Two failures were found and fixed along the way, both worth knowing about:
`SpeciesEntry`'s validating constructor originally took eight untyped positional
arguments, which is the signature Julia auto-generates, so it overwrote the
generated one and made precompilation illegal — it is an inner constructor now.
And `test/corea_test_models.jl` extended `states` after a bare `using InferCell`,
which Julia rejects; it needs an explicit `import`. That second one silently took
all five new test files with it, so the first green-looking run had in fact never
executed any of the new tests. Worth remembering as a failure mode: a passing
count that did not go *up* is not a passing run.

## Pre-merge review (2026-08-27, later session)

`/simplify` + `/check-PR` ran against PR #38: nine parallel review agents,
every MAJOR finding adversarially verified. Verdict: merge after fixes, all of
which are applied on this branch. The substance:

- **Hardened the contract** where the review found silent-pass gaps, all
  additive validation: an inbound `MassEdge` must appear in `inputs()` (the
  drift wave-1 authors would hit — only `inputs` wires a state into dynamics);
  producer-with-no-consumer is now a reported dead end (`DeadEnd` gained
  `missing_role`, replacing `consumers` with `modules`); a `ClampedEdge` on a
  chemostat must match the registry's held value; `conc_<species>` imports are
  checked against the registry row's `initial_value` (the GTP 1.6627-vs-0.1
  trap now fails at load); a `governing` declaration binds even with a single
  holder; truncated TSV rows and misspelled `Informedness` cells error instead
  of being skipped; `CurrencyEdge.pool` is registry-checked; positional edge
  construction validates (inner constructors).
- **`DeferredCounterEdge` gained `smoothing`**, required exactly when
  `clip = :smoothed` — the edge-kinds spec's "parameter controlling the
  smoothing is exposed" clause, previously unimplemented.
- **Simplified**: vocab validation consolidated on `_check_vocab` (now in
  `parameters.jl`, shared with the registry's gstd⟺informedness invariant);
  `deviates_from_published` derives from `deviation_reason`; dead code dropped
  (`_module_name`, unused accessors, write-only `SourceTable.path`, unreachable
  `state_index::Nothing` arm); six unused exports removed; corea exports moved
  into their own files so seven wave-1 branches don't all append to one block
  in `src/InferCell.jl`.
- **Spec/docs reconciled with the code**: the edge-kinds table no longer names
  a `debited_pool` field that never existed, marks which fields default to the
  published model's behaviour, and the unknown-kind failure is honestly
  resolution-time; design.md sketches and the informedness vocabulary now match
  the source.

Wave-1 advisories the review surfaced but deliberately did not act on: lower
`ResolvedEdge`s into concrete typed callback state before any hot path iterates
them; resolve chemostat `held_value`s to plain `Float64` at build time if an
RHS ever reads them; the cost-with-no-payer throw means an Expression module
debiting ATP cannot resolve alone (spec-compliant, but the first wave-1 author
will hit it); assert registry-vs-table agreement when wave 1 vendors the real
balanced tables.

## Next steps

1. `/opsx:sync establish-corea-interface` to promote the delta specs into
   `openspec/specs/`, **then** merge the PR to `main`. Sync rather than archive,
   so wave 0 stays open for amendment while wave 1 reads the contract.
2. Only then propose the seven wave-1 changes. A wave-1 proposal written before
   the sync will invent its own interface instead of consuming this one. Propose
   all seven before applying any — that is the last good chance to catch an
   interface assumption two modules disagree about.

## Open questions this change did not settle

Deliberately deferred, and none of them changes the contract:

- Which clip policy Core A′ actually runs under for inference. All three are
  declarable; whether NUTS on the ODE block needs `:smoothed` is answered by
  trying it, in wave 3.
- Whether the 60 s rebuild stays piecewise-constant. Same shape — declarable,
  defaults to published, and the deviation should be measured rather than
  assumed. A wave-2 task on `add-cme-rebuild-60s`.
- The on-disk format and location of the real balanced tables. The loader takes
  a path and a logical file name; only the parser behind it would differ. Wave 1
  vendors them.

## Scope notes for the reviewer

Three things went slightly beyond the literal task list, each flagged rather
than absorbed:

- **`module_id` is a new protocol function** not named in `tasks.md`. The specs
  require coupling errors to name both modules involved; using the type name
  alone would name two instances of one type identically. It defaults to the
  type name, so nothing changes for the five existing sub-models.
- **The two `mkdocs` nav fixes** are unrelated to this change but were blocking
  its docs verification and CI.
- **`test/corea_test_models.jl`** establishes a test-double convention the suite
  did not have — no existing test file subtypes `AbstractSubModel`, because they
  all drive the five real models. The resolver cannot be tested without doubles.
