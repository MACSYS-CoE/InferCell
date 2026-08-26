# The Core A′ reduction, drawn on fig 1's canvas

One derived-from-a-derived figure for
[`dev/notes/reduced-syn3a-scoping.md`](../../reduced-syn3a-scoping.md).

| File | What it shows |
|---|---|
| `fig1r_state_graph_reduced.{pdf,png}` | Fig 1's state-coupling graph with everything Core A′ deletes greyed and crossed, and everything it keeps in altered form recoloured and annotated `old → new` |

It is the companion to
[`../minimal-cell-coupling/fig1_state_graph.pdf`](../minimal-cell-coupling), and is
drawn on that figure's canvas at that figure's node positions, so **the two overlay**.
Print both and lay one on the other, or flip between them on screen: nothing the
reduction leaves alone moves.

## Rebuild

```
python3 make_reduced.py            # -> fig1r_state_graph_reduced.{pdf,png}
python3 make_reduced.py --verify   # ... and prove it overlays fig 1
```

Needs Graphviz (`neato`) and, for `--verify`, `pdftoppm` and ImageMagick. Nothing is
cloned and nothing is downloaded: the sibling directory's `state_interfaces.csv` and
`draw_state.py` are the only inputs.

## How the alignment is achieved

Fig 1 is not laid out by an algorithm — `draw_state.py` pins every node in `POS` and
renders with `neato -n`. The obvious approach is therefore to reuse fig 1's own DOT
and substitute attributes. That is nearly right, and the trap is worth recording:

**`overlap=false` runs overlap removal, which is content-sensitive.** Add one line of
text to one node label and *every* node in the graph shifts. Measured: lengthening
Cofactor's label by a single line moved Transport by 2.2 pt, TITLE by 3.4 pt and the
canvas by 3.3 pt. Naive substitution produces two figures that look aligned and are
not.

So `make_reduced.py`:

1. rebuilds fig 1's DOT from `state_interfaces.csv`, and **asserts the resulting PNG
   is byte-identical to the shipped `fig1_state_graph.png`** — the layout it is about
   to freeze is provably the shipped figure's;
2. reads back the geometry overlap removal actually produced
   (`neato -n -Tdot`: final `pos`, `width`, `height` per node, and the graph `bb`);
3. re-pins every node at that resolved position and sets `overlap=true`, so no further
   movement is possible;
4. takes fig 1's **edge** statements verbatim and greys or relabels them per
   `reduction_spec.py`; authors the twenty **node** statements afresh;
5. pins two invisible anchors at fig 1's bounding-box corners, so a figure with less
   ink in it cannot shrink the canvas;
6. asserts the rendered page is `3082 × 2534 pt`, the same as fig 1.

`--verify` then rasterises both PDFs at 150 dpi, checks they are the same size
(6421 × 5280 px), checks the hook — a node the reduction does not touch — is identical
*pixel for pixel*, and writes `build/overlay_check.png`, a red/cyan overlay in which
anything black sits in both figures.

**What does differ, and why.** Boxes that changed carry an extra line of text, so they
grow — about their fixed centre — and the edges attached to them shift by a few points.
That is the entire difference. It is confined to the modules the reduction touches,
which is where the reader is meant to be looking.

## Judgement calls

Fig 1 is derived: every edge in it comes from `state_interfaces.csv`, which comes from
the model files. **This figure is not.** It is a hand-written overlay, and the calls
below are ours, not the model's. They live in `reduction_spec.py` so they can be
argued with.

- **Two edges are greyed although both endpoints survive.** `Transport → Nucleotide`
  goes because the salvage network and its transport are cut entirely;
  `Nucleotide → Transport` (pyruvate) goes because no retained transport reaction
  consumes pyruvate — PTS produces it.
- **Both edges into tRNA charging from the CME are greyed.**
  `transcription → tRNA charging` because Core A′'s seventeen genes are all
  protein-coding, so uncharged tRNA becomes a conserved pool rather than a
  transcription product; `translation → tRNA charging` because the lumped charging
  step carries no synthetase species. Both follow from the lumping, which the scoping
  note flags as ours rather than the model's.
- **The catalytic edge is greyed but the channel is not cut.** Fig 1 draws
  `protein counts → enzyme conc.` once, and happens to route it into Cofactor, which
  Core A′ deletes. Rerouting it would have been the only edge in the figure not
  sitting where fig 1 puts it, so the route is kept, greyed, and a note next to it
  says where the channel actually goes. This is a drawing artefact of fig 1's
  "drawn once" device, not a claim about the model.
- **Fifteen surviving edge labels were corrected by hand** from the scoping note —
  `13dpg, 2dr1p +3` becomes `3pg, pyr`, `20 aa pools` becomes `aa pool`, and so on.
  They were not re-derived from a filtered `state_interfaces.csv`. Treat them as
  annotations, not as data.
- **The red cross is a glyph, not two drawn diagonals.** Fig 1 sets
  `outputorder=edgesfirst`, which paints every edge beneath every node, so a
  cross drawn as edges lands underneath the box it is meant to cancel. An image
  overlay was tried instead; this Graphviz build does not honour PNG alpha in
  `image=`. A scaled `✗` in a plaintext node is drawn with the nodes, and works.

## Reaction counts

Fig 1's counts are the model's active reactions per module. Core A′'s are from the
scoping note's Core A′ section, and sum to its twenty-one:

| Module | Fig 1 | Core A′ | Kept |
|---|---|---|---|
| Transport | 87 | 6 | PTS cascade ×5, `L_LACt2r` |
| Central | 26 | 10 | glycolysis through lactate; NOX dropped |
| Nucleotide | 46 | 5 | PGK3, PYK3, ADK1, PPA, GK1 |
| Lipid | 17 | — | cut; membrane lipid supply made exogenous |
| Cofactor | 13 | — | cut |
| Amino acid | 1 | — | cut |
| **ODE total** | **190** | **21** | |

Genes go 455 → 17, and the CME side to 17 × 3 reactions plus one `TranslocRate` for
ptsG. Growth survives but reaches only ~1.07× initial volume, so the 2× clamp is never
approached — report fractional growth, never a doubling time.
