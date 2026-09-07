# Core A′ progress: what is built, and what each phase still owes

Two figures, added because six pull requests had merged and it had become hard
to see what they delivered. The reason is structural rather than sloppy: **every
merged phase is framework, and the coupling figures draw the model.** A reader
who knows only fig 1c has no way to tell which of its parts exist.

| Figure | What it claims |
|---|---|
| `fig1d_state_graph_progress.{pdf,png}` | **only what the Core A′ reduction keeps**, laid out fresh, with every module and device badged by the spec §11 phase that writes it, and the edge-kind key gaining an **executes since** column |
| `fig2d_phase_roadmap.{pdf,png}` | all eighteen phases in dependency order, including the five that put nothing on a state graph at all |

| File | |
|---|---|
| `make_progress.py` | Builds both; see "Rebuild" |
| `progress_spec.py` | The spec §11 parser, fig 1d's layout and edge table, and the authored phase → figure-element mapping |

## Why fig 1d is not fig 1c's canvas

The chain up to here is cumulative: `../minimal-cell-coupling/` is fig 1, the
published model's state coupling derived from `state_interfaces.csv`;
`../reduced-syn3a-coupling/` is fig 1r, fig 1 with everything Core A′ deletes
greyed and crossed, on fig 1's canvas so the two overlay;
`../corea-coupling/` is fig 1c, fig 1r with tRNA charging moved into the ODE
block.

Every one of those draws the reduction *onto* fig 1, which is the right way to
show **what the reduction did**. It is the wrong way to show **what is built**:
about half the ink is four deleted modules and thirty dead couplings, drawn in
grey, and a reader looking for the live machine has to subtract them by eye. So
fig 1d breaks the chain's one rule and keeps only the survivors —
**seven modules, three devices, two clamped drawing devices and thirty
couplings** — on a layout chosen for that content:

```
   ┌── CME block ──┐                        ┌──── ODE block ────┐
   │  transcription│                        │ Transport  Central│
   │  translation  │ ──── the hook, 1 s ──► │ charging   Nucleo.│
   │  mRNA decay   │                        │   currency pool   │
   └───────────────┘                        └───────────────────┘
            ▲                                         │
            └────── the CME rebuild, 60 s ◄───────────┘
                    growth — on neither clock
```

Two blocks, three coupling devices, one cycle. That is the sentence the figure
is trying to say, and fig 1c cannot say it because fig 1's geometry was not
chosen to.

## Two axes, and why they are drawn differently

Fig 1c spends a red ✗ and a grey wash on *removed by the Core A′ reduction*.
Fig 1d has nothing to mark that way — the removed modules are simply not there —
so its only reduction signal is the amber fill meaning *survives in altered
form*. Progress is a second, independent axis:

| Mark | Axis | Means |
|---|---|---|
| amber fill | the reduction | the module survives the reduction in altered form |
| ✔ green pill | progress | the phase has **merged**; the pill carries its PR |
| ☐ hollow pill | progress | **not written**; the number is the phase in spec §11 |

The two are stated side by side in the "what is built" key, because the one
misreading worth preventing is *amber means unfinished*. Every box in fig 1d is
amber and none of them is written; those are separate facts.

## Nothing here states a status, and nothing here invents a coupling

Two independent derivations, both asserted at build time.

**Status** comes from `progress_spec.phases()`, which parses `spec/spec.md`
§11 — the phase headings, the `**PR:**` line, and the `- [x]` / `- [ ]` ticks.
Every badge, pill, title count and footer count comes from that parse. Only the
mapping from a phase to the thing it puts on the figure is authored:
`ELEMENT_PHASE` (fig 1d's node ids), `EDGE_KIND_PHASE` (which phase turned each
edge kind from a declaration into something that runs) and `DELIVERS` (one line
per phase, for fig 2d's chips).

**The coupling set** comes from fig 1c. `live_couplings()` re-runs fig 1c's own
edge stage, then drops the removed modules and the dead couplings its spec
names; `check_edges()` asserts that fig 1d's authored `EDGES` table is *exactly*
that set, plus only what `ADDED_EDGES` declares with a reason. The table is
authored rather than inherited so the labels can be trimmed and the layout can
breathe — but it cannot lose an arrow by being forgotten, and it cannot gain one
by accident.

Five guards make drift fail the build rather than quietly produce a wrong
figure. All five are exercised in a scratch copy in the PR that added this
directory:

| Guard | Message |
|---|---|
| a coupling dropped from `EDGES` | `fig 1d drops a coupling Core A′ keeps: [('Central', 'Transport')]` |
| a coupling invented and not declared | `fig 1d draws a coupling fig 1c does not, and does not declare it in ADDED_EDGES` |
| a renamed `### Phase N — …` heading | `spec/spec.md §11 has no heading for phase(s) [9]` |
| a merged phase with an unticked task | `phase 4 is merged as #45 but only 6/7 ticked` |
| an `ELEMENT_PHASE` id the figure does not draw | `fig 1d has no node 'HOOKX' to badge` |

## The one coupling fig 1d draws that fig 1c does not

`HOOK → ODE block`, the catalytic channel's second half. Fig 1 draws that
channel once, into Cofactor; Core A′ deletes Cofactor, so fig 1c has to grey the
edge and put the real targets in a footnote (its `CATNOTE`: "the catalytic
channel is not cut … in Core A′ it runs into Transport, Central and
Nucleotide"). With the deleted modules gone there is nothing to grey and no
footnote to hang it on, so fig 1d draws the edge where the channel actually
runs — into the block, once, which is how fig 1 labels it anyway. It is declared
in `ADDED_EDGES` with that reason, and it is the only entry.

## What fig 1d says, in one line

Six of seven edge kinds execute and not one box in the figure exists in code.
The **executes since** column is where a reader sees it: mass and currency since
phase 1, catalytic and deferred counter since phase 3, rate constant since phase
4, volume since phase 5. The clamped edge validates, resolves and appears in
reports, but nothing runs it and no phase schedules it.

The key carries an eighth row below a divider, and it is deliberately **not**
counted among the seven: the two solid intra-CME arrows are not `CouplingEdge`s
at all. Spec §12 (2026-09-04) settled that a jump module's peer writes are gated
on `written_states` rather than on an edge, and phase 2 is what made them work.

## What fig 2d adds

Three things fig 1d structurally cannot show.

- **The dependency structure** §11's parallelism paragraph states: phase 1
  alone, then the ODE modules fanning out beside the framework track, phase 9
  waiting on phase 8, phases 11 and 12 waiting on phase 10, and phase 11 able to
  proceed on a double at a cost R1 records.
- **Phases 13 and 14**, which assert the figure complete and check conservation
  over it — work *about* the figure rather than in it.
- **Phases 15 to 17**, which put nothing on a state graph at all. They are the
  inference layer, and their chips say so in as many words.

## Rebuild

```
python3 make_progress.py     # -> fig1d_state_graph_progress.{pdf,png}
                             #    fig2d_phase_roadmap.{pdf,png}
```

Needs Graphviz (`neato`) and `pdfinfo`. It still builds fig 1c's edge list, so
it also prints fig 1c's note that the rebuilt fig 1 is not byte-identical to the
shipped raster — see `../corea-coupling/README.md` for why that is expected on
the cluster and why it is reported rather than fatal. A clean rebuild reproduces
both PNGs byte for byte.

## Judgement calls

- **Charging sits nearest the hook and Nucleotide furthest.** The other way
  round, the hook's single edge to charging passes about three points above
  whichever box it crosses, and every position for its label lands on that box.
  Putting Central directly over Nucleotide is a free consequence: their two mass
  edges become a clean vertical pair.
- **The ODE block sits further right than the drawing needs.** Five labelled
  couplings leave the hook for it; at closer spacing their labels stack.
- **Every edge label sits on a white table cell with no border.** Graphviz
  centres a label on its own edge, so the line runs through the text. The
  backing is invisible against the page and hides the few millimetres behind it.
- **The two blocks are drawn as ground, not as clusters.** `neato -n` cannot use
  clusters, and a panel node lets an edge terminate on a *block* where the
  coupling is into the block rather than into one module of it — which is what
  the catalytic channel needs. Statement order is paint order: panels, then
  couplings, then module boxes, so each box sits on top of its own arrows.
  Panels are sized from the module positions they contain, so moving a module in
  `POS` resizes its block instead of spilling out of it.
- **Badges straddle the top-right corner, pushed 9 pt out along both axes.**
  Square on the corner they covered the device boxes' own titles.
- **Positions are measured, not guessed.** Badge placement needs the box sizes
  as *drawn*, so the build renders once through `neato -n -Tdot`, reads back the
  true geometry and only then places the badges. `-Tdot` renormalises the origin,
  so the measurement is realigned against the hook, whose position is known.
- **Two phases get a badge that is not on a box.** Phase 1 built the
  contribution channel and phase 2 built jump composition; neither is a box, and
  badging one of the three CME boxes with phase 2 would claim that box exists.
  They sit inside their block panels instead, which is what they are properties
  of.
- **Phase 0 appears on fig 2d and on no element of fig 1d.** It retired
  OpenSpec. Its chip says "— nothing in the figure —", as phases 15 to 17 do.
- **Fig 2d is laid out on a fixed grid, not by `dot`.** The two-row spine —
  framework across the top, assembly and inference across the bottom, the module
  fan-out between — is the shape of §11's argument, and a solver would not find
  it. The cost is that the six edges converging on phase 13 are straight lines
  that clip a chip or two; they are drawn lighter for that reason.
- **Fig 1c is not superseded.** It remains the record of what the reduction cut,
  and it is the only figure that overlays the published model. Fig 1d cannot
  answer that question, and does not try.
