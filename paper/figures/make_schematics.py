#!/usr/bin/env python3
"""The paper's two schematic figures, written to this directory as PDFs.

Run from the repository root: python3 paper/figures/make_schematics.py
They read no data. fig_graphs' edge styles follow the paper's boundary section.
"""
import os

import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt
from matplotlib.patches import FancyArrowPatch, FancyBboxPatch

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


def fig_architecture():
    fig, ax = canvas(6.5, 3.0)
    bh = 0.2
    box(ax, 0.01, 0.62, 0.19, 0.3, "Parameters\nprior, source,\nlabel")
    box(ax, 0.01, 0.12, 0.19, 0.3, "Modules\nstates, formalism,\ntyped edges")
    box(ax, 0.27, 0.37, 0.15, 0.3, "Resolver\nchecks the\ncomposition")
    box(ax, 0.27, 0.05, 0.15, 0.2, "Declared\nreductions", size=7.5)
    box(ax, 0.5, 0.62, 0.2, bh + 0.02, "Exact hybrid\njump + ODE, stepped")
    box(ax, 0.5, 0.16, 0.2, bh + 0.02, "Mean field\none ODE problem")
    box(ax, 0.78, 0.62, 0.21, bh + 0.02, "Exact route\nGibbs, path sampling")
    box(ax, 0.78, 0.16, 0.21, bh + 0.02, "Mixed route\nNUTS, gradients")
    box(ax, 0.53, 0.41, 0.46, 0.13, "predict: a posterior, a description, an intervention", fc="#ffffff", size=7.5)
    arrow(ax, (0.205, 0.77), (0.265, 0.56))
    arrow(ax, (0.205, 0.27), (0.265, 0.46))
    arrow(ax, (0.345, 0.365), (0.345, 0.255))
    arrow(ax, (0.425, 0.56), (0.495, 0.72))
    arrow(ax, (0.425, 0.46), (0.495, 0.27))
    arrow(ax, (0.705, 0.73), (0.775, 0.73))
    arrow(ax, (0.705, 0.27), (0.775, 0.27))
    arrow(ax, (0.6, 0.615), (0.6, 0.545))
    arrow(ax, (0.6, 0.38), (0.6, 0.405))
    arrow(ax, (0.885, 0.615), (0.885, 0.545))
    arrow(ax, (0.885, 0.38), (0.885, 0.405))
    fig.savefig(os.path.join(HERE, "fig_architecture.pdf"), metadata={"CreationDate": None})
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


def fig_model():
    fig, ax = canvas(6.5, 3.1)
    ax.text(0.14, 0.95, "ODE block: metabolism", ha="center", fontweight="bold", color=INK)
    ax.text(0.86, 0.95, "Jump block: gene expression", ha="center", fontweight="bold", color=INK)
    for i, t in enumerate(("Glycolysis to lactate", "Sugar transport and\nlactate export",
                           "Nucleotide recycling", "Lumped tRNA charging")):
        box(ax, 0.01, 0.7 - i * 0.215, 0.26, 0.17, t)
    for i, t in enumerate(("Transcription", "Transcript decay", "Translation and\nmembrane insertion")):
        box(ax, 0.73, 0.62 - i * 0.26, 0.26, 0.2, t)
    for y, t, col, lw, a0, a1 in (
            (0.73, "Catalytic edges: protein counts\ninto rate laws, 1 s", ACCENT, 1.5, 0.725, 0.275),
            (0.51, "Deferred costs: drains\non the pools, 1 s", ACCENT, 1.5, 0.725, 0.275),
            (0.29, "Volume edges: membrane proteins\nset the geometry, each step", ACCENT, 1.5, 0.725, 0.275),
            (0.07, "Rate constants: pools\ninto jump rates, 60 s", ACCENT2, 0.8, 0.275, 0.725)):
        arrow(ax, (a0, y), (a1, y), color=col, lw=lw)
        ax.text(0.5, y + 0.02, t, ha="center", va="bottom", fontsize=8, color=INK)
    fig.savefig(os.path.join(HERE, "fig_model.pdf"), metadata={"CreationDate": None})
    plt.close(fig)


if __name__ == "__main__":
    fig_architecture()
    fig_model()
    fig_graphs()
