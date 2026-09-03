# Building Core A′ in waves

> **SUPERSEDED — 2026-09-03.** `spec/spec.md` is now the authoritative plan, and
> it **reorders this file**. Two of the framework changes filed here as "wave 2"
> alter the sub-model protocol itself: letting a module contribute to a state it
> does not own changes what `dynamics` returns, and fixing jump composition
> changes what `reactions` receives. This file freezes those files against module
> branches, so doing them after the seven modules re-opens seven merged PRs. They
> move ahead of the modules, and the 1 s handshake is proven on a two-module toy
> with wall-clock measured before the modules are built — because a mixed
> ODE/jump composition is refused in one line today and no periodic-hook
> machinery exists at all. See `spec/spec.md` §4 D0.
>
> Two further gaps this file does not name: jump-block composition mis-indexes
> state and parameters, and inference dispatch reads only the first module in a
> composition. Its "7-way parallel" wave 1 is at most 5-way and its "3-way
> parallel" wave 2 is serial.
>
> **One content error entered here, not just an ordering one.** The wave-1 table
> below assigns `add-lumped-trna-charging` to the CME block. The scoping note puts
> the charging step and both tRNA pools on the ODE side, and the frozen registry
> marks both species metabolite-scale. `spec/spec.md` inherited the error and
> corrected it in its §12 amendment 1: charging is an ODE module, which also
> removes 3.49 million jump events per trajectory.
>
> **Do not follow this file's ordering.** It is kept for the reasoning behind the
> decomposition — the module boundaries, the file layout, and the rule that each
> module carries its own conservation check — all of which the spec adopts. The
> OpenSpec commands throughout are dead; `openspec/` is retired to
> `dev/archive/openspec/`.

How the Core A′ port is decomposed into parallel work, and which OpenSpec
commands drive each stage.

Read alongside [`../notes/reduced-syn3a-scoping.md`](../notes/reduced-syn3a-scoping.md),
which decides *what* Core A′ is. This file decides *how it gets built*, and
nothing more — it carries no scientific decisions of its own. Where the two
disagree, the scoping note wins and this file is stale.

The decomposition is not invented here either. It is read off
[`fig1r_state_graph_reduced.pdf`](../notes/figures/reduced-syn3a-coupling/fig1r_state_graph_reduced.pdf):
the module boxes are the units of parallel work, the three interface nodes are
the points where parallel work has to rejoin, and the edge-kind legend is
already an interface type system.

---

## The short version

| | |
|---|---|
| **Unit of parallel work** | one module = one OpenSpec change = one branch = one PR |
| **What makes it safe** | the interface contract is a *main spec* on `main` before any module starts |
| **Scope of one `propose`** | one PR-sized change. Never the whole project |
| **Where the DAG lives** | this file. OpenSpec has no project-scope artifact |
| **The gate between waves** | `/opsx:sync` or `/opsx:archive`, which promotes delta specs into `openspec/specs/` |

---

## Do we propose the whole project first?

**No.** `/opsx:propose` generates proposal, delta specs, design and tasks in one
pass, and its `tasks.md` is a flat checklist that `/opsx:apply` walks in order.
A whole-project propose would produce a hundred-item checklist with no branch
boundaries in it — the exact flat-list problem that makes parallel work
impossible. One propose per PR-sized change, every time.

That leaves a real gap: nothing in OpenSpec holds the project-scope dependency
graph. OpenSpec's dependency modelling is artifact-level and *within* a single
change (proposal → specs → design → tasks, the edges you see in
`openspec status --json`). It has no notion of one change depending on another.

This file is that missing layer. It is hand-maintained, it is small enough to
read in one sitting, and it exists because the graph was already drawn before
any tool was asked to infer it.

---

## The dependency graph

```
                    ┌─────────────────────────────┐
  WAVE 0            │  establish-corea-interface  │   serial · blocks everything
                    └──────────────┬──────────────┘
                                   │  /opsx:sync → openspec/specs/
      ┌──────────┬──────────┬──────┴───┬──────────┬──────────┬──────────┐
  W1  │ central  │transport │nucleotide│  transcr │  transl  │ mrna-dec │ trna-chg │
      └──────────┴──────────┴──────┬───┴──────────┴──────────┴──────────┘
                                   │            7-way parallel
                    ┌──────────┬───┴──────┬──────────┐
  WAVE 2            │   hook   │ rebuild  │  growth  │   3-way parallel
                    └──────────┴────┬─────┴──────────┘
                                    │
  WAVE 3                   validation 1→6 → SBC        serial by design
```

Wave 3 is serial on purpose, not for want of trying: the scoping note orders the
balance checks so that each catches errors the next would mask, and the
adenylate dead end only became visible over a full cycle rather than a burst.

---

## Wave 0 — the interface contract

**One change: `establish-corea-interface`. Everything waits on it.**

This is small in code and total in consequence. It is what lets seven branches
run without stepping on each other, because after it lands no module may
unilaterally change how it touches another.

Not greenfield: `src/interface.jl` already defines `AbstractSubModel` with
`states` / `parameters` / `dynamics` / `reactions` / `inputs` / `formalism` /
`inference_mode`. Wave 0 extends that protocol to carry Core A′'s coupling,
rather than inventing a new one.

Three things in scope:

1. **The state registry.** The ~32 dynamic states and 5 chemostats from the
   scoping note's state list, named once, ordered once. Every module indexes
   into this and none of them defines its own names.

2. **The edge-kind contract.** fig1r's legend enumerates seven ways state
   crosses an interface, and they have genuinely different semantics: *mass*,
   *currency*, *deferred counter*, *catalytic*, *rate constant*, *volume*,
   *clamped*. Each needs a declaration form a module can use and the
   orchestrator can resolve.

3. **The parameter loader, carrying provenance.** Every imported value records
   the file it came from. This is not bookkeeping neatness — the
   central-versus-nucleotide cross-file trap has already bitten twice, once in
   kinetics (PGK3/PYK3: 319.5 vs 140.8; 1874 vs 672.8) and once in initial
   concentrations (GTP at 1.6627 mM, not the 0.1 mM prior default). Resolving it
   once by hand is how it bites a third time.

Wave 0 is also the right place to settle the two interface hazards the scoping
note flags, because both are properties of the boundary rather than of any
module: the clipped `max(0, ·)` ATP drain that will obstruct NUTS gradients, and
whether the 60 s rebuild stays piecewise-constant or becomes continuous.

### Commands

```bash
/opsx:propose establish-corea-interface     # proposal + delta specs + design + tasks
# review the artifacts, then:
/opsx:apply establish-corea-interface       # walk tasks.md
openspec validate establish-corea-interface
```

Then the gate that unblocks wave 1:

```bash
/opsx:sync establish-corea-interface        # delta specs → openspec/specs/, change stays active
```

`/opsx:sync` rather than `/opsx:archive` is deliberate. Sync promotes the delta
specs into main specs *without* archiving, so wave 1 can start reading the
contract as an established spec while wave 0's own implementation is still
being polished. Archive it later, once wave 0's PR is merged:

```bash
/opsx:archive establish-corea-interface
```

**Merge wave 0 to `main` before proposing anything in wave 1.** `/opsx:propose`
reads `openspec/specs/` for context, so a wave-1 proposal written before the
sync will invent its own interface instead of consuming the real one.

---

## Wave 1 — the seven modules

Seven independent changes, each one branch and one PR. They touch disjoint files
and read a contract that is already frozen on `main`.

| Change | Block | Contents |
|---|---|---|
| `add-central-glycolysis` | ODE | 10 rxns, PGI→LDH_L. NOX dropped |
| `add-pts-transport` | ODE | PTS cascade ×5 + `L_LACt2r` lactate export |
| `add-nucleotide-recycling` | ODE | PGK3, PYK3, ADK1, PPA, GK1 |
| `add-corea-transcription` | CME | 17 genes, constitutive |
| `add-corea-translation` | CME | 17 genes |
| `add-mrna-decay` | CME | 17 mRNAs, NMP return |
| `add-lumped-trna-charging` | CME | 1 lumped step — ours, not the model's |

Each change carries its own conservation check as part of its own tasks, not
deferred to wave 3: PTS protein conservation belongs to `add-pts-transport`,
redox balance to `add-central-glycolysis`, phosphate closure to
`add-nucleotide-recycling`. A module that cannot pass its own local balance
check is not done, and finding that out on its own branch is much cheaper than
finding it out after the join.

Two of these carry a "ours, not the model's" label that must survive into the
code and any results that depend on them: the lumped tRNA charging step, and the
chemostatting of CTP/UTP/amino acids.

### Where the files go

Seven branches inventing seven locations is the cheapest kind of merge conflict
to avoid, so the layout is fixed here rather than decided per branch.

| | |
|---|---|
| Module source | `src/organisms/coreA/<module>.jl` |
| Module tests | `test/test_corea_<module>.jl`, included from `test/runtests.jl` |
| Wiring | one `include` in `src/InferCell.jl`, after the registry |

`src/organisms/coreA/` is where everything specific to this organism lives — the
registry, the seven modules, and the parameter tables under a `data/` subdirectory
once they are vendored. Anything under a second organism would sit beside it.

The framework files at `src/` root — `edges.jl`, `resolver.jl`, `loader.jl`,
`labels.jl`, plus `interface.jl` and `orchestrator.jl` — are **read-only to a
module branch.** They are the contract wave 0 froze; a module that needs one of
them changed has found an interface bug, and that is a conversation on `main`,
not an edit on a branch that six other branches will merge over.

Keep the `corea_` prefix on the test file. A test is named after the source
file it covers, so a module under `src/organisms/coreA/` keeps the prefix — the
four files this convention was introduced alongside (`test_edges.jl`,
`test_resolver.jl`, `test_loader.jl`, `test_labels.jl`) dropped it only because
the code *they* cover moved to `src/` root.

One caveat on the framework/organism split: `resolver.jl` and `loader.jl` still
reach for the global registry — `resolver.jl` uses nine of its symbols
(`is_registered`, `is_chemostatted`, `is_dynamic`, `dynamic_species`,
`species_group`, `species_index`, `held_value`, `COREA_SPECIES`,
`TRANSCRIPTION_ATOL`) and `loader.jl` four (`is_registered`, `species_entry`,
the `SpeciesEntry` layout, `TRANSCRIPTION_ATOL`). `edges.jl` and `labels.jl`
name none. Every one of those calls sits inside a function body, so the include
order in `src/InferCell.jl` is a reading convention rather than a load-time
constraint. Making resolver and loader take a registry is what a second organism
costs, and it is deliberately not paid yet.

### Commands

Propose all seven up front — it is cheap, and having the seven proposals side by
side is the last good chance to catch an interface assumption that two modules
disagree about:

```bash
/opsx:propose add-central-glycolysis
/opsx:propose add-pts-transport
# … and so on for the remaining five
openspec list                    # confirm seven active changes
openspec view                    # dashboard across specs and changes
```

Then fan out. One worktree per change, one Claude session per worktree:

```bash
git worktree add ../InferCell-central   -b change/add-central-glycolysis
git worktree add ../InferCell-transport -b change/add-pts-transport
# … then in each worktree:  /opsx:apply <change-name>
```

This is conflict-free by construction: each change's artifacts live in their own
`openspec/changes/<name>/` directory, so the `- [ ]` → `- [x]` churn in one
branch's `tasks.md` never meets another's.

**One change per branch, always.** Two branches sharing a change means two
branches editing one `tasks.md`, which conflicts on every task and leaves the
spec-sync state ambiguous. Parallelism inside a single change belongs inside a
single `/opsx:apply` session, where independent tasks can be fanned out to
subagents.

As each PR merges:

```bash
/opsx:archive <change-name>      # on main, after merge
```

Archive per merged change rather than in a batch at the end — the main specs
need to stay current for the changes still in flight.

---

## Wave 2 — the three interfaces

The join. These need wave 1 complete, and each is a node in fig1r rather than a
module: they are where the blocks actually touch.

| Change | Node | Cadence |
|---|---|---|
| `add-hook-1s-coupling` | the hook | 1 s · counts → ODE ICs → integrate → counts back |
| `add-cme-rebuild-60s` | the CME rebuild | 60 s · live pools → k_transcription, k_translation |
| `add-growth-volume-coupling` | growth | 1 s · surface area → V → rescales every count↔mM |

Three-way parallel, same worktree mechanics as wave 1. The rebuild is the only
ODE→CME direction in the model, so it is the change that decides whether
criterion 2 (bidirectional coupling) is actually met; growth is the second
coupling channel, and cheap.

---

## Wave 3 — validation, then inference

Serial, in the scoping note's order, because each check catches errors the next
would mask:

1. non-negativity → 2. carbon balance → 3. redox balance →
4. adenylate/guanylate conservation **over a full cycle** → 4b. phosphate
closure → 5. PTS protein conservation → 6. nominal trajectory

**Wave 3 owns the completeness assertion the interface deliberately does not
make.** Wave 0's `resolve_coupling` reports an unowned state or a dead end
rather than failing on one, because a module author validating a single module
alone imports species whose owners and producers are absent by construction. So
a successful resolve says the declarations are mutually consistent, not that the
boundary is closed. For the assembled model that distinction matters: a cost
with no paying state and a stranded moiety are both errors there, and both would
currently pass. Wave 3's first task is therefore to assert completeness — every
registry dynamic state owned, `dead_ends` empty — and fail on it, whether as a
`complete = true` mode on `resolve_coupling` or a separate assertion over the
graph. Check 4 (adenylate/guanylate conservation) is the numerical form of the
same property; this is the structural form, and it is cheap.

Then, and only then, `add-corea-synthetic-recovery`: simulate at known θ, infer,
check coverage, SBC per module.

Wave 3 is also where the first genuinely unknown number gets measured —
wall-clock per Core A′ trajectory — which is what sizes any inference budget.

---

## Keeping the graph honest

Two `rules` blocks in `openspec/config.yaml` make every future proposal declare
its own position in the graph, so this file does not silently drift out of date:

```yaml
rules:
  proposal:
    - Include a "Dependencies" section naming the other changes this one
      requires, or "None". Name changes, not vague phases.
    - State which wave this change belongs to, per
      dev/plans/reduced-syn3a-wave-plan.md
  tasks:
    - Group tasks into numbered waves. Tasks within a wave must be
      independently implementable, in any order.
    - Every module change carries its own conservation or balance check as a
      task, not deferred to a later validation change.
```

OpenSpec will not enforce these — they are prompt-level, and a "Dependencies:
None" written by an agent that did not check is worth nothing. They exist so
that the information is on the page where a human can see it is wrong.

If you want the contract enforced rather than requested, `openspec schema fork
spec-driven` copies the schema into the project, where a custom artifact with
its own `requires` edges *is* enforced by `openspec status`. Worth doing only if
the convention proves it cannot hold.

---

## Command reference

| Stage | Command |
|---|---|
| Create a change with all artifacts | `/opsx:propose <name>` |
| Revise a change's artifacts | `/opsx:update <name>` |
| Think before proposing | `/opsx:explore` |
| Implement | `/opsx:apply <name>` |
| Promote delta specs, keep change active | `/opsx:sync <name>` |
| Finish after merge | `/opsx:archive <name>` |
| See everything | `openspec view` · `openspec list` |
| Check one change | `openspec status --change <name> --json` |
| Check the whole root | `openspec doctor` · `openspec validate` |
