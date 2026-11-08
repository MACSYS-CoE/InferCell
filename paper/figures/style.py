"""The paper's figure style, shared by make_figures.py and make_schematics.py.

Every figure spans both columns of arxiv.sty, 6.5 in (165.1 mm, measured from \\linewidth), and is
saved at that width so LaTeX draws it at 1:1 and the point sizes here are the printed sizes.
"""
import matplotlib
from matplotlib import font_manager

matplotlib.use("Agg")
import matplotlib.pyplot as plt

WIDTH = 6.5  # inches, the printed width of a figure* in arxiv.sty
SIZE = 7.5  # pt, all text but panel letters: subscripts print at 0.7 of this, above 5 pt
LETTER = 8  # pt, bold panel letters

# Validated together: slot 1 and 2 of the reference categorical palette.
EXACT, MIXED = "#2a78d6", "#eb6834"
INK, MUTED, BAND, FILL = "#0b0b0b", "#52514e", "#e4e3df", "#f0efec"


def _family():
    """Arial where it is installed, else Liberation Sans, which has Arial's metrics."""
    for name in ("Arial", "Liberation Sans", "Helvetica"):
        try:
            font_manager.findfont(name, fallback_to_default=False)
            return name
        except ValueError:
            pass
    return "DejaVu Sans"


FAMILY = _family()

plt.rcParams.update({
    "font.family": "sans-serif", "font.sans-serif": [FAMILY], "font.size": SIZE,
    "axes.titlesize": SIZE, "axes.labelsize": SIZE, "xtick.labelsize": SIZE, "ytick.labelsize": SIZE,
    "legend.fontsize": SIZE, "legend.title_fontsize": SIZE, "figure.titlesize": SIZE,
    # Symbols in the same face as the text, not DejaVu.
    "mathtext.fontset": "custom", "mathtext.rm": FAMILY, "mathtext.it": FAMILY + ":italic",
    "mathtext.bf": FAMILY + ":bold", "mathtext.sf": FAMILY, "mathtext.cal": FAMILY, "mathtext.default": "it",
    "legend.frameon": False, "legend.handlelength": 1.6, "legend.borderaxespad": 0.2,
    "axes.spines.top": False, "axes.spines.right": False, "axes.edgecolor": MUTED,
    "axes.linewidth": 0.6, "axes.titlepad": 3, "axes.labelpad": 2,
    "xtick.color": MUTED, "ytick.color": MUTED, "xtick.labelcolor": INK, "ytick.labelcolor": INK,
    "xtick.major.width": 0.6, "ytick.major.width": 0.6, "xtick.major.size": 2.5, "ytick.major.size": 2.5,
    "xtick.major.pad": 2, "ytick.major.pad": 2,
    "axes.labelcolor": INK, "axes.titlecolor": INK, "lines.linewidth": 1.0, "lines.markersize": 3,
    "pdf.fonttype": 42, "axes.unicode_minus": True,
})


def letter(fig, x, y, text):
    """A panel letter in figure coordinates, so letters align whatever the axes."""
    fig.text(x, y, text, fontsize=LETTER, fontweight="bold", color=INK, ha="left", va="top")
