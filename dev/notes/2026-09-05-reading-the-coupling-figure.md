# Reading the coupling figure: what it shows, and the two layers it cannot

**2026-09-05.** Written to reconcile the mental model carried by
[`figures/reduced-syn3a-coupling/fig1r_state_graph_reduced.pdf`](figures/reduced-syn3a-coupling/)
against what phases 0–3 of [`spec/spec.md`](../../spec/spec.md) actually
delivered. It adds no decisions; it is an orientation note. Where it and the
spec disagree, the spec wins.

## First: the figure to navigate by is 1c, not 1r

`fig1r_state_graph_reduced.pdf` is superseded. The live figure is
[`figures/corea-coupling/fig1c_state_graph_corea.pdf`](figures/corea-coupling/),
and the difference is one box: **tRNA charging moved out of the CME block into
the ODE block** (spec §12 amendment 1, §4 D13).

That move greys three couplings fig 1r draws live — `Transport → charging`
(the amino-acid pool is a chemostat), `Central → charging` (charging now draws
ATP inside the ODE block, drawn as an intra-ODE currency flow), and
`charging → translation` (replaced by an orange deferred-counter edge from the
hook, ~553 residues/s against a pool of order 10³, which is the counter check 7
watches). Fig 1r is kept only as the record of what
`dev/plans/reduced-syn3a-wave-plan.md` believed, because moving a node breaks
its one guarantee — that it overlays fig 1 pixel for pixel.

## The figure is one of three layers

**Layer 1 — the state graph.** Boxes are modules, arrows are state couplings.
This is what the figure draws, and for that it is authoritative and complete.

**Layer 2 — the execution layer.** Every arrow in the picture is one of the
**seven `CouplingEdge` kinds on three different clocks** (`src/edges.jl`): mass,
currency, catalytic, deferred counter, rate constant, volume, clamped. The 1 s
arrows, the 60 s arrows and the volume arrow that is on neither clock are drawn
identically and are three different mechanisms. This layer is invisible in the
figure, and it is everything built so far.

**Layer 3 — the inference layer.** Parameters, priors, provenance, observables,
the exact reference posterior, coverage, rank statistics, the shrinkage table.
The figure has no representation of any of it. It is the layer the project
exists for.

A thinner fourth layer is worth holding: **provenance and honesty**
(`src/labels.jl`, outputs T1 and T2). Every amber "changed" box in the figure is
a place where Core A′ departs from the published model, and each departure has
to reach any reported result with its measured cost attached. The amber fill is
the informal version of that; `reduction_declarations` is the checkable one.

## Where we are: arrows are built, boxes are not

The single most useful sentence: **not one box in the figure exists in code.**
`src/organisms/coreA/` holds `registry.jl` and nothing else. Phases 0–3
(PRs #41–#44) are all Layer 2, exercised on two-module toys.

That is D0 on purpose — framework before modules, kill risk on a toy first —
because two framework changes alter the contract the modules are written
against, and the 1 s handshake was the only genuine kill risk in the project.

| In the figure | Phase | Status |
|---|---|---|
| **HOOK** device (the 1 s exchange in the middle) | 3 — *the kill phase* | ✅ merged. Neither kill condition fired |
| Arrows *within* the ODE grid (mass/currency into pools a module does not own) | 1 | ✅ merged |
| The CME block being three coexisting boxes at all | 2 | ✅ merged |
| **REBUILD** device (60 s: the arrows into transcription's and translation's rate constants) | 4 | ⬜ next |
| **GROWTH** device and the volume arrow that dilutes all four ODE modules | 5 | ⬜ |
| Box: **Central** (10 rxns) | 6 | ⬜ |
| Box: **Transport** (6 rxns) | 7 | ⬜ |
| Box: **Nucleotide** (5 rxns) | 8 | ⬜ |
| Box: **tRNA charging** (the box that moved) | 9 | ⬜ |
| Box: **transcription** (17 genes) | 10 | ⬜ |
| Box: **translation** (17, plus ptsG translocation) | 11 | ⬜ |
| Box: **mRNA decay** (17) | 12 | ⬜ |
| *The whole figure asserted complete* — every state owned, every drawn edge executed, no dead ends | 13 | ⬜ |
| Conservation checks over the figure (carbon, redox, adenylate, phosphate, carrier) | 14 | ⬜ |
| — nothing in the figure — | 15–17 | ⬜ |

Four of eighteen phases, all of them plumbing. The greyed-and-crossed part of
the figure is settled, but only as a scoping decision frozen in the registry;
the amber boxes are all still to be written.

## Layer 2, in one line

**Four of seven edge kinds execute.** Mass and currency from phase 1; catalytic
and deferred counter from phase 3. Rate constant, volume and clamped are still
**declare-only** — they validate, resolve and appear in reports, but nothing
runs them. Phase 4 turns on rate constant (REBUILD), phase 5 turns on volume
(GROWTH). That is why two framework phases still precede every module: two of
the figure's devices are, at the moment, drawings.

## Layer 3, and why the figure is load-bearing here

The figure is a forward-model picture, so it cannot show what the claim is
about. Usually that is fine. Here it is not, because of D13: **no parameter is
shared across the ODE/stochastic boundary**, and that is a property of the
published model rather than of our reduction.

So the claim cannot be "one constant appears in both blocks and both data
streams inform it". It has to be: a parameter *local to one block* has its
posterior moved by *the other block's data* — and the only thing that can move
it is the arrows in the figure. The two cells that carry the entire result
(§6 F6) are therefore readable as paths on it:

- the **ptsG promoter strength** shrinking under *metabolite-only* data —
  information travelling CME → HOOK → ODE;
- an **enolase catalytic constant** shrinking under *transcript-only* data —
  information travelling ODE → REBUILD → CME.

Everything in phases 15–17 — observables, the exact reference posterior on a
two-gene toy, coverage curves, rank statistics — exists to make those two cells
mean something, and none of it is drawable on a state graph.

## The short version

The figure is a picture of a *model*. What phases 0–3 built is the *machine that
can run a picture like that*, with the riskiest part of the machine proven on a
two-module toy before a single line of syn3A biology was written. Phases 4–5
finish the machine. Phases 6–12 build the seven boxes, fanning out. Phase 13
asserts the figure complete as an object in code. Then a layer of work begins
that the figure was never able to show.
