#!/usr/bin/env python3
"""The paper's three schematic figures, written to this directory as PDFs.

Run from the repository root: python3 paper/figures/make_schematics.py
They read no data. fig_graphs' edge styles follow the paper's boundary section.
"""
import os

import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt
from matplotlib.lines import Line2D
from matplotlib.patches import FancyArrowPatch, FancyBboxPatch
from matplotlib.path import Path

HERE = os.path.dirname(os.path.abspath(__file__))
INK, MUTED, FILL = "#0b0b0b", "#52514e", "#f0efec"
ACCENT, ACCENT2 = "#2a78d6", "#eb6834"
plt.rcParams.update({"font.family": "sans-serif", "font.size": 8, "pdf.fonttype": 42})


def box(ax, x, y, w, h, text, fc=FILL, ec=MUTED, size=8):
    ax.add_patch(FancyBboxPatch((x, y), w, h, boxstyle="round,pad=0.005,rounding_size=0.015",
                                fc=fc, ec=ec, lw=0.8))
    ax.text(x + w / 2, y + h / 2, text, ha="center", va="center", fontsize=size, color=INK)


def arrow(ax, a, b, color=MUTED, ls="-", lw=1.0, rad=0.0):
    ax.add_patch(FancyArrowPatch(a, b, arrowstyle="-|>", mutation_scale=9, color=color, ls=ls, lw=lw,
                                 connectionstyle=f"arc3,rad={rad}", shrinkA=0, shrinkB=0))


def canvas(w, h):
    fig, ax = plt.subplots(figsize=(w, h))
    ax.set_xlim(0, 1)
    ax.set_ylim(0, 1)
    ax.axis("off")
    fig.subplots_adjust(0, 0, 1, 1)
    return fig, ax


def fig_architecture(stem="fig_architecture"):
    # Printed at \linewidth of the two-column template's figure*: 469.75 pt = 165 mm = 6.5 in, 1:1.
    fig, ax = canvas(6.5, 3.0)
    top, low, rh = 0.68, 0.06, 0.24           # the two rows of descriptions and routes
    box(ax, 0.01, 0.62, 0.19, 0.3, "Parameters\nprior, source,\nlabel")
    box(ax, 0.01, 0.10, 0.19, 0.3, "Modules\nstates, formalism,\ntyped edges")
    box(ax, 0.27, 0.40, 0.15, 0.26, "Resolver\nchecks the\ncomposition")
    box(ax, 0.27, 0.06, 0.15, 0.2, "Declared\nreductions", size=7.5)
    box(ax, 0.5, top, 0.2, rh, "Exact hybrid\njump + ODE, stepped")
    box(ax, 0.5, low, 0.2, rh, "Mean field\none ODE problem")
    box(ax, 0.78, top, 0.21, rh, "Exact route\nGibbs, path sampling")
    box(ax, 0.78, low, 0.21, rh, "Mixed route\nNUTS, gradients")
    box(ax, 0.5, 0.425, 0.49, 0.13, "predict: a posterior, a description, an intervention",
        fc="#ffffff", size=7.5)
    arrow(ax, (0.205, 0.77), (0.265, 0.6))
    arrow(ax, (0.205, 0.25), (0.265, 0.46))
    arrow(ax, (0.345, 0.395), (0.345, 0.265))
    arrow(ax, (0.425, 0.6), (0.495, top + rh / 2))
    arrow(ax, (0.425, 0.46), (0.495, low + rh / 2))
    arrow(ax, (0.705, top + rh / 2), (0.775, top + rh / 2))
    arrow(ax, (0.705, low + rh / 2), (0.775, low + rh / 2))
    # Each description and each route feeds predict: from the row above and the row below.
    for x in (0.6, 0.885):
        arrow(ax, (x, top - 0.005), (x, 0.56), lw=1.2)
        arrow(ax, (x, low + rh + 0.005), (x, 0.42), lw=1.2)
    fig.savefig(os.path.join(HERE, stem + ".pdf"), metadata={"CreationDate": None})
    plt.close(fig)


def node(ax, x, y, text, color=INK):
    ax.text(x, y, text, ha="center", va="center", fontsize=8, color=color,
            bbox=dict(boxstyle="round,pad=0.3", fc=FILL, ec=MUTED, lw=0.8))


def fig_graphs():
    fig, axes = plt.subplots(1, 3, figsize=(6.5, 2.1))
    fig.subplots_adjust(0.01, 0.02, 0.99, 0.88, wspace=0.05)
    for ax in axes:
        ax.set_xlim(0, 1)
        ax.set_ylim(0, 1)
        ax.axis("off")
    a, b, c = axes
    a.set_title("A  Coupling: what drives what", fontsize=8, loc="left")
    node(a, 0.22, 0.5, "Gene\nexpression")
    node(a, 0.78, 0.5, "Metabolism")
    arrow(a, (0.38, 0.58), (0.62, 0.58), lw=2.0, rad=-0.2)
    a.text(0.5, 0.78, "gain of order one", ha="center", fontsize=7, color=MUTED)
    arrow(a, (0.62, 0.42), (0.38, 0.42), lw=0.7, rad=-0.2)
    a.text(0.5, 0.2, "gain 0.04–0.09", ha="center", fontsize=7, color=MUTED)

    b.set_title("B  Likelihood, paths integrated out", fontsize=8, loc="left")
    node(b, 0.22, 0.75, r"$\theta_G$")
    node(b, 0.78, 0.75, r"$\theta_M$")
    node(b, 0.22, 0.25, "Transcript\nfactor")
    node(b, 0.78, 0.25, "Bulk\nfactor")
    arrow(b, (0.22, 0.65), (0.22, 0.4), lw=1.5)
    arrow(b, (0.78, 0.65), (0.78, 0.4), lw=1.5)
    arrow(b, (0.3, 0.66), (0.66, 0.38), lw=1.5)
    arrow(b, (0.7, 0.66), (0.36, 0.38), ls=":", lw=0.8)

    c.set_title("C  Information: what moves what", fontsize=8, loc="left")
    node(c, 0.22, 0.75, r"$\theta_G$")
    node(c, 0.78, 0.75, r"$\theta_M$")
    node(c, 0.22, 0.25, "Transcripts")
    node(c, 0.78, 0.25, "Metabolites")
    arrow(c, (0.22, 0.34), (0.22, 0.65), color=ACCENT, lw=1.5)
    arrow(c, (0.78, 0.34), (0.78, 0.65), color=ACCENT, lw=1.5)
    arrow(c, (0.66, 0.34), (0.3, 0.66), color=ACCENT2, ls="--", lw=1.2)
    arrow(c, (0.34, 0.34), (0.7, 0.66), color=MUTED, ls=":", lw=0.8)

    fig.savefig(os.path.join(HERE, "fig_graphs.pdf"), metadata={"CreationDate": None})
    plt.close(fig)


def routed(ax, pts, color, ls="-", lw=1.2):
    """An arrow along a polyline (orthogonal routing through the gap between the blocks)."""
    path = Path(pts, [Path.MOVETO] + [Path.LINETO] * (len(pts) - 1))
    ax.add_patch(FancyArrowPatch(path=path, arrowstyle="-|>", mutation_scale=8, color=color,
                                 ls=ls, lw=lw, shrinkA=0, shrinkB=0, joinstyle="miter"))


def enclosure(ax, x, y, w, h, title):
    ax.add_patch(FancyBboxPatch((x, y), w, h, boxstyle="round,pad=0.004,rounding_size=0.02",
                                fc="none", ec=MUTED, lw=0.8, ls=(0, (4, 2))))
    ax.text(x + w / 2, y + h + 0.015, title, ha="center", va="bottom", fontsize=8,
            fontweight="bold", color=INK)


def fig_model(stem="fig_model"):
    # Printed at \linewidth of figure*: 165 mm = 6.5 in, 1:1. Edges as declared by the modules in
    # src/organisms/coreA/: CatalyticEdge (central_glycolysis, nucleotide_recycling), the
    # DeferredCounterEdge costs of transcription, transcript_decay and translation and translation's
    # credit to the sugar-transport proteins, RateConstantEdge (transcription on the NTPs, translation
    # on charged tRNA). The VolumeEdge lies inside pts_transport and is not drawn (caption). The pools
    # are registry species shared by several modules.
    fig, ax = canvas(6.5, 3.3)
    fig.subplots_adjust(0, 0.11, 1, 1)           # the key sits below the axes, in the figure
    L0, L1, R0, R1 = 0.05, 0.34, 0.66, 0.98      # module columns; the gap holds the lanes
    enclosure(ax, 0.01, 0.15, 0.345, 0.76, "ODE block: metabolism")
    enclosure(ax, 0.645, 0.15, 0.345, 0.76, "Jump block: gene expression")
    # The pools sit at the top, beside transcription, so its two edges run straight.
    box(ax, L0, 0.735, L1 - L0, 0.125, "Shared pools:\nnucleotides, tRNA", fc="#ffffff", size=7.5)
    ode = {"nuc": 0.60, "gly": 0.455, "pts": 0.31, "trna": 0.165}
    for k, t in (("nuc", "Nucleotide recycling"), ("gly", "Glycolysis to lactate"),
                 ("pts", "Sugar transport and lactate export"), ("trna", "Lumped tRNA charging")):
        box(ax, L0, ode[k], L1 - L0, 0.105, t, size=7.5 if k == "pts" else 8)
    box(ax, R0, 0.70, R1 - R0, 0.16, "Transcription")
    box(ax, R0, 0.46, R1 - R0, 0.15, "Transcript decay")
    box(ax, R0, 0.165, R1 - R0, 0.235, "Translation and\nmembrane insertion")
    mid = lambda k: ode[k] + 0.0525
    blue, orange = ACCENT, ACCENT2
    routes = [  # (polyline, colour, line style); arrowheads at the last point
        # Deferred costs: each jump module's drains land on the pools, 1 s
        ([(R0, 0.853), (L1, 0.853)], blue, "--"),
        ([(R0, 0.535), (0.56, 0.535), (0.56, 0.799), (L1, 0.799)], blue, "--"),
        ([(R0, 0.375), (0.62, 0.375), (0.62, 0.826), (L1, 0.826)], blue, "--"),
        # Rate constants: the pools set transcription and translation rates, 60 s
        ([(L1, 0.770), (R0, 0.770)], orange, "-"),
        ([(L1, 0.743), (0.595, 0.743), (0.595, 0.345), (R0, 0.345)], orange, "-"),
        # Catalytic: translated enzyme counts into the rate laws, 1 s
        ([(R0, 0.295), (0.50, 0.295), (0.50, mid("nuc")), (L1, mid("nuc"))], blue, "-"),
        ([(R0, 0.255), (0.47, 0.255), (0.47, mid("gly")), (L1, mid("gly"))], blue, "-"),
        # Deferred credit: translation credits the sugar-transport proteins (PtsG among them), 1 s
        ([(R0, 0.215), (0.44, 0.215), (0.44, mid("pts")), (L1, mid("pts"))], blue, "--"),
    ]
    for pts, col, ls in routes:
        routed(ax, pts, col, ls, lw=1.3 if col == blue else 1.1)
    # Currency: nucleotide recycling, glycolysis and tRNA charging write the pools, inside the ODE solve.
    for k, x in (("nuc", 0.038), ("gly", 0.028), ("trna", 0.018)):
        routed(ax, [(L0, mid(k)), (x, mid(k)), (x, 0.79), (L0, 0.79)], MUTED, "-", lw=0.8)

    # Key, with each edge kind's clock.
    handles = [Line2D([], [], color=blue, lw=1.3, ls="-"), Line2D([], [], color=blue, lw=1.3, ls="--"),
               Line2D([], [], color=orange, lw=1.0), Line2D([], [], color=MUTED, lw=0.8)]
    labels = ["Catalytic: protein counts into rate laws, each 1 s handshake",
              "Deferred costs and credits, each 1 s handshake",
              "Rate constants: pools into jump rates, 60 s",
              "Currency: modules write the pools, in the ODE solve"]
    fig.legend(handles, labels, loc="lower center", bbox_to_anchor=(0.5, 0.0), ncol=2, frameon=False,
              fontsize=7, handlelength=2.6, columnspacing=1.6, borderaxespad=0.2)
    fig.savefig(os.path.join(HERE, stem + ".pdf"), metadata={"CreationDate": None})
    plt.close(fig)


if __name__ == "__main__":
    fig_architecture()
    fig_model()
    fig_graphs()
