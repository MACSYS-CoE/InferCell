"""Fig 1c spec -- fig 1r with the lumped tRNA charging step in the ODE block.

Everything fig 1r greys, crosses and relabels is inherited from
../reduced-syn3a-coupling/reduction_spec.py. This file lists only what changes
when charging moves, and where the evidence for the move lives:
dev/notes/reduced-syn3a-scoping.md (always placed it there), the frozen registry
(marks both tRNA species metabolite-regime), and spec/spec.md §12 amendment 1.
"""
import os, sys

sys.path.insert(0, os.path.join(os.path.dirname(os.path.abspath(__file__)),
                                "..", "reduced-syn3a-coupling"))
from reduction_spec import (GREY_FILL, GREY_LINE, GREY_TEXT, GREY_EDGE,   # noqa: F401
                            CHANGED_FILL, CROSS, REMOVED_NODES)
import reduction_spec as R

CHG = "CME:tRNAcharging"      # fig 1's node id is kept so its edges still bind

# --- the node that moves ----------------------------------------------------
# Fig 1 and fig 1r draw charging inside the CME. Here it is an ODE module: blue
# border, amber fill (it survives in altered form), placed on the ODE side between
# the two removed modules of the bottom row.
CHANGED_NODES = {k: v for k, v in R.CHANGED_NODES.items() if k != CHG}
CHARGING_LABEL = (
    '<b>tRNA charging</b><br/>'
    '<font point-size="8.5"><b>20 × 5 rxns → 1 lumped step</b></font><br/>'
    '<font point-size="8.5"><b>CME → ODE</b></font><br/>'
    '<font point-size="7.5">owns the charged and uncharged pools<br/>'
    'ATP → AMP + PP<sub>i</sub>, into recycling’s pools<br/>'
    '<i>lumping and formalism ours;<br/>placement the scoping note’s</i></font>')


def charging_pos(g):
    """Midpoint of the removed bottom-row pair, so the fourth ODE module sits on
    the ODE grid without displacing anything fig 1 draws."""
    (xa, ya, *_), (xb, yb, *_) = g["Cofactor"], g["AminoAcid"]
    return (xa + xb) / 2, (ya + yb) / 2


# --- edges ------------------------------------------------------------------
DEAD_ENDPOINT = list(R.DEAD_ENDPOINT)

# Three couplings fig 1r kept live are gone once charging is an ODE module.
DEAD_COUPLING = list(R.DEAD_COUPLING) + [
    (("Transport", CHG, ""),
     "none of Core A′'s six transport reactions moves an amino acid; the aa pool "
     "is a chemostat and is drawn from the medium"),
    (("Central", CHG, ""),
     "ATP is no longer synced across the boundary for charging; the reaction draws "
     "it inside the ODE block, via the currency pool"),
    ((CHG, "CME:Translation", ""),
     "no intra-block transfer remains: translation debits the charged pool and "
     "credits the uncharged one as a deferred counter through the hook"),
]

# Two greyed edges lose their labels: with charging moved onto the ODE grid,
# "20 aa pools" and "atp, synced counts" would land on top of each other at Lipid.
STRIP_LABELS = [("Transport", CHG, ""), ("Central", CHG, "")]

# Greyed edges keep fig 1's labels, as every other dead edge does; the two
# relabels fig 1r applied to them are dropped so each key matches one edge once.
RELABEL = {k: v for k, v in R.RELABEL.items()
           if k not in {("Transport", CHG, ""), (CHG, "CME:Translation", "")}}
RELABEL[("CME:Translation", "HOOK", "fmettrna")] = [
    ("atp, fmettrna, gtp", "atp, gtp · charged tRNA")]
RELABEL[("GROWTH", "HOOK", "")] = [
    ("dilutes all six ODE modules", "dilutes all four ODE modules")]

# Edges fig 1 has no statement for, authored here. Colours are fig 1's kinds.
C_CTR, C_CL = "#b9770e", "#c0392b"
NEW_EDGES = [
    f'  "HOOK" -> "{CHG}" [color="{C_CTR}", style=dashed, penwidth=2.0, '
    f'label=<<font color="{C_CTR}"><b>debit charged tRNA · credit uncharged</b><br/>'
    '~553 residues/s against a pool of order 10³<br/>'
    '<i>a buffer of seconds — check 7’s counter to watch</i></font>>];',
    f'  "Medium" -> "{CHG}" [color="{C_CL}", style=dashed, penwidth=1.7, '
    f'label=<<font color="{C_CL}">aa pool, chemostatted</font>>];',
    f'  "{CHG}" -> "CURRENCY" [color="#b3bfc4", penwidth=0.8, arrowhead=none, style=dashed];',
]
