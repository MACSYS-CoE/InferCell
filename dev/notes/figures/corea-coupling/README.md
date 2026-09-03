# Core A′ coupling, with tRNA charging in the ODE block

Fig 1c: the third figure in a chain of three.

| Figure | Where | What it shows |
|---|---|---|
| fig 1 | `../minimal-cell-coupling/` | The published well-stirred model's state coupling, derived from the model files |
| fig 1r | `../reduced-syn3a-coupling/` | Fig 1 with everything Core A′ deletes greyed and crossed, on fig 1's canvas so the two overlay |
| **fig 1c** | here | **Fig 1r with the lumped tRNA charging step moved from the CME to the ODE block**, which is where the scoping note always placed it |

| File | |
|---|---|
| `fig1c_state_graph_corea.{pdf,png}` | The figure |
| `make_corea.py` | Builds it; see "Rebuild" |
| `corea_spec.py` | Only what changes relative to fig 1r's `reduction_spec.py`, which it imports |

## Why a third figure rather than a fix to fig 1r

Fig 1r is a faithful overlay of fig 1, and fig 1 draws charging in the CME because
the published model runs it there (22 of 23 ODE synthetase reactions are commented
out upstream for that reason). Fig 1r inherited the CME placement from
`dev/plans/reduced-syn3a-wave-plan.md`, which mis-assigned the lumped step; the
scoping note, the frozen registry (both tRNA species are `:metabolite` regime) and
`spec/spec.md` §12 amendment 1 all put it in the ODE block. Moving a node breaks
fig 1r's one guarantee, that it overlays fig 1 pixel for pixel, so the corrected
placement is a new figure and fig 1r stands as the record of what the wave plan
believed.

## What moved, and what it did to the edges

The charging box is the only node not at fig 1's position. It sits on the ODE
grid between the two removed modules of the bottom row, with the ODE block's blue
border and the "changed" amber fill. Everything else is fig 1r.

Three couplings fig 1r drew live are greyed, because they were boundary crossings
and there is no boundary there any more:

- `Transport → charging` (amino-acid pools). None of Core A′'s six transport
  reactions moves an amino acid; the pool is a chemostat. Replaced by a red
  clamped edge from the growth medium.
- `Central → charging` (ATP, synced each second). Charging now draws ATP inside
  the ODE block and pays AMP and PP<sub>i</sub> into the pools nucleotide
  recycling owns. Drawn the way fig 1 draws every intra-ODE currency flow: a
  dashed connector to the currency-pool device, which now names four modules.
- `charging → translation` (charged tRNA, intra-CME mass). Replaced by an orange
  deferred-counter edge from the hook: translation debits the charged pool and
  credits the uncharged one, about 553 residues per second against a pool of
  order 10³. This is the counter the spec's check 7 watches.

Two labels change: translation's counter into the hook gains "charged tRNA", and
growth now dilutes four ODE modules. Two greyed edges lose their labels, because
with charging on the ODE grid "20 aa pools" and "atp, synced counts" would land on
top of each other beside Lipid.

The rate-constant channel `charging → rebuild → k_translation` is unchanged. It
was always ODE → CME in substance; it is now drawn from an ODE node.

## Rebuild

```
python3 make_corea.py            # -> fig1c_state_graph_corea.{pdf,png}
```

Needs Graphviz (`neato`) and `pdfinfo`. There is no `--verify` and no ImageMagick
dependency: charging is meant to move, so fig 1c does not overlay fig 1 and no
check claims it does.

**Fig 1r's byte-identity assertion is downgraded here, and the reason is worth
recording.** Fig 1r's build starts by regenerating fig 1 from
`state_interfaces.csv` and asserting the PNG is byte-identical to the shipped one,
so the geometry it freezes is provably the shipped figure's. On the cluster
(Graphviz 2.44, different fonts) the regenerated fig 1 renders at 3253 × 2695 px
against the shipped 2932 × 2384, so that assertion cannot pass and fig 1r itself
does not build there. Fig 1c reports the mismatch and continues: the positions it
freezes are fig 1's hand-set pins after *this* Graphviz's overlap removal. The
layout is the same to the eye; the canvas is 3439 × 2845 pt against fig 1's
3082 × 2534. Rebuilding on the machine that shipped fig 1 would recover fig 1's
canvas exactly.

## Reaction counts

| Module | Fig 1 | Core A′ | Kept |
|---|---|---|---|
| Transport | 87 | 6 | PTS cascade ×5, `L_LACt2r` |
| Central | 26 | 10 | glycolysis through lactate; NOX dropped |
| Nucleotide | 46 | 5 | PGK3, PYK3, ADK1, PPA, GK1 |
| tRNA charging | CME | 1 | the lumped step, moved in |
| Lipid, Cofactor, Amino acid | 31 | — | cut |
| **ODE total** | **190** | **22** | |

The CME side is 17 × 3 reactions plus one `TranslocRate` for ptsG: 52, matching
the scoping note. Fig 1r's README sums the ODE block to 21 because it counted
charging on the other side.

## Judgement calls

Inherited from fig 1r's README, plus:

- **The charging box sits over dead edges.** The bottom-row gap between Cofactor
  and Amino acid carries fig 1's greyed edges between them and the hook's greyed
  fmet-tRNA debit. The box paints over them. They are dead, and the alternative
  was to put the only live ODE module somewhere off the ODE grid.
- **The amino-acid edge comes from the medium, not from Transport.** Fig 1r kept
  fig 1's `Transport → charging` route and relabelled it "aa pool". That was
  already inaccurate for Core A′ and is corrected here rather than carried.
- **Currency is drawn by convention, not by species.** Charging's ATP, AMP and
  PP<sub>i</sub> traffic is the dominant adenylate channel in Core A′ (the spec's
  D13 and F1), and a reader might want it as a labelled blue edge into Nucleotide.
  Fig 1's device routes all intra-ODE currency through the pool node unlabelled,
  and this figure keeps the device so the two read the same way.
