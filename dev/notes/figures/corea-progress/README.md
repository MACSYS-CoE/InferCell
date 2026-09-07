# Core A′ progress: what is built, and what each phase still owes

Two figures, added because six pull requests had merged and it had become hard
to see what they delivered. The reason is structural rather than sloppy: **every
merged phase is framework, and fig 1c draws the model.** A reader who knows only
fig 1c has no way to tell which of its parts exist.

| Figure | What it claims |
|---|---|
| `fig1d_state_graph_progress.{pdf,png}` | fig 1c, with every module and device badged by the spec §11 phase that writes it, and the edge-kind key gaining an **executes since** column |
| `fig2d_phase_roadmap.{pdf,png}` | all eighteen phases in dependency order, including the five that put nothing on a state graph at all |

| File | |
|---|---|
| `make_progress.py` | Builds both; see "Rebuild" |
| `progress_spec.py` | The spec §11 parser, and the authored phase → figure-element mapping |

This is the fourth figure in a chain. `../minimal-cell-coupling/` is fig 1, the
published model's state coupling derived from `state_interfaces.csv`;
`../reduced-syn3a-coupling/` is fig 1r, fig 1 with everything Core A′ deletes
greyed and crossed; `../corea-coupling/` is fig 1c, fig 1r with tRNA charging
moved into the ODE block. **Fig 1d changes no model content.** It re-runs fig
1c's builder and overlays status on the result, so the two cannot disagree about
the model — and it lands on fig 1c's canvas, 3439 × 2845 pt, which the build
prints so a displaced node is caught rather than shipped.

## Two axes, and why they are drawn differently

Fig 1c already spends a red ✗ and a grey wash on *removed by the Core A′
reduction*. Progress is a second, independent axis, so it gets its own
vocabulary and its own key:

| Mark | Axis | Means |
|---|---|---|
| ✗, grey | the reduction | the module is cut from the **model** |
| amber fill | the reduction | the module survives in altered form |
| ✔ green pill | progress | the phase has **merged**; the pill carries its PR |
| ☐ hollow pill | progress | **not written**; the number is the phase in spec §11 |

A crossed-out box says nothing about whether anything is written, and a hollow
badge says nothing about whether the module is in the model. Conflating them is
the one misreading these figures exist to prevent, so the "what is built" key
states both axes side by side.

## Nothing here states a status

A progress figure that drifts from the task list is worse than none. So
`progress_spec.phases()` **parses `spec/spec.md` §11**: the phase headings, the
`**PR:**` line, and the `- [x]` / `- [ ]` ticks. Every badge, every pill, every
count in both titles and both footers comes from that parse.

What is authored is only the mapping from a phase to the thing it puts on the
state graph — `ELEMENT_PHASE` (fig 1c's own node ids), `EDGE_KIND_PHASE` (which
phase turned each edge kind from a declaration into something that runs) and
`DELIVERS` (one line per phase, for fig 2d's chips).

Four guards make drift fail the build rather than quietly produce a wrong
figure:

- a missing or renamed `### Phase N — …` heading;
- a phase recorded as merged with an unticked task;
- a phase recorded as merged with no task list at all;
- a `ELEMENT_PHASE` node id that fig 1c does not draw.

All four are exercised in the PR that added this directory.

## What fig 1d says, in one line

Six of seven edge kinds execute and not one box in the figure exists in code.
That is the whole state of play, and the **executes since** column is where a
reader sees it: mass and currency since phase 1, catalytic and deferred counter
since phase 3, rate constant since phase 4, volume since phase 5; the clamped
edge validates, resolves and reports, but nothing runs it and no phase schedules
it.

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

Needs Graphviz (`neato`) and `pdfinfo`, and inherits fig 1c's build, so it also
prints fig 1c's note that the rebuilt fig 1 is not byte-identical to the shipped
raster — see `../corea-coupling/README.md` for why that is expected on the
cluster and why it is reported rather than fatal. A clean rebuild reproduces
both PNGs byte for byte.

## Judgement calls

- **Badges straddle the top-right corner, pushed 9 pt out along both axes.**
  Square on the corner they covered the device boxes' own titles. Pushed further
  out they stopped reading as belonging to the box.
- **Positions are measured, not guessed.** Three boxes are not the size fig 1
  drew them — charging is re-labelled and re-pinned, and both keys are authored
  here — so the build renders once through `neato -n -Tdot`, reads back the true
  geometry, and only then places the badges. `-Tdot` renormalises the origin, so
  the measurement is realigned against the hook, whose position is known.
- **`MOVENOTE` is dropped from fig 1d.** Fig 1c's note explaining why charging
  moved is settled and recorded in spec §12 amendment 1; on a progress figure it
  is noise, and its slot was the clearest space near the unbuilt ODE boxes.
  `CATNOTE` is kept: the catalytic channel's real targets are still a live fact
  about phases 6 to 8.
- **Two phases get a badge that is not on a box.** Phase 1 built the
  contribution channel and phase 2 built jump composition; neither is a box, and
  badging one of the three CME boxes with phase 2 would claim that box exists.
  They sit under the block labels instead, which is what they are properties of.
- **Phase 0 appears on fig 2d and on no element of fig 1d.** It retired OpenSpec.
  That is honest, and its chip says "— nothing in the figure —" like phases 15
  to 17 do.
- **Fig 2d is laid out on a fixed grid, not by `dot`.** The two-row spine —
  framework across the top, assembly and inference across the bottom, the module
  fan-out between — is the shape of §11's argument, and a solver would not find
  it. The cost is that the six edges converging on phase 13 are straight lines
  that clip a chip or two; they are drawn lighter than the rest for that reason.
