#!/usr/bin/env python3
"""The paper's three schematic figures, written to this directory as PDFs.

Run from the repository root: python3 paper/figures/make_schematics.py
They read no data. fig_graphs' edge styles follow the paper's boundary section.
"""
import os
import sys

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
from style import EXACT as ACCENT, FILL, INK, LETTER, MIXED as ACCENT2, MUTED, SIZE, WIDTH  # noqa: E402

import matplotlib.pyplot as plt  # noqa: E402
from matplotlib.patches import FancyArrowPatch, FancyBboxPatch  # noqa: E402

HERE = os.path.dirname(os.path.abspath(__file__))


def box(ax, x, y, w, h, text, fc=FILL, ec=MUTED, size=SIZE):
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
    fig, ax = canvas(WIDTH, 2.7)
    bh = 0.2
    box(ax, 0.01, 0.62, 0.19, 0.3, "Parameters\nprior, source,\nlabel")
    box(ax, 0.01, 0.12, 0.19, 0.3, "Modules\nstates, formalism,\ntyped edges")
    box(ax, 0.27, 0.37, 0.15, 0.3, "Resolver\nchecks the\ncomposition")
    box(ax, 0.27, 0.05, 0.15, 0.2, "Declared\nreductions")
    box(ax, 0.5, 0.62, 0.2, bh + 0.02, "Exact hybrid\njump + ODE, stepped")
    box(ax, 0.5, 0.16, 0.2, bh + 0.02, "Mean field\none ODE problem")
    box(ax, 0.78, 0.62, 0.21, bh + 0.02, "Exact route\nGibbs, path sampling")
    box(ax, 0.78, 0.16, 0.21, bh + 0.02, "Mixed route\nNUTS, gradients")
    box(ax, 0.53, 0.435, 0.46, 0.1, "predict: a posterior, a description, an intervention", fc="#ffffff")
    arrow(ax, (0.205, 0.77), (0.265, 0.56))
    arrow(ax, (0.205, 0.27), (0.265, 0.46))
    arrow(ax, (0.345, 0.365), (0.345, 0.255))
    arrow(ax, (0.425, 0.56), (0.495, 0.72))
    arrow(ax, (0.425, 0.46), (0.495, 0.27))
    arrow(ax, (0.705, 0.73), (0.775, 0.73))
    arrow(ax, (0.705, 0.27), (0.775, 0.27))
    arrow(ax, (0.6, 0.615), (0.6, 0.54))
    arrow(ax, (0.6, 0.38), (0.6, 0.43))
    arrow(ax, (0.885, 0.615), (0.885, 0.54))
    arrow(ax, (0.885, 0.38), (0.885, 0.43))
    fig.savefig(os.path.join(HERE, "fig_architecture.pdf"), metadata={"CreationDate": None})
    plt.close(fig)


def node(ax, x, y, text, color=INK):
    ax.text(x, y, text, ha="center", va="center", fontsize=SIZE, color=color,
            bbox=dict(boxstyle="round,pad=0.3", fc=FILL, ec=MUTED, lw=0.8))


def title(ax, letter, text):
    ax.text(0, 1.06, letter, transform=ax.transAxes, fontsize=LETTER, fontweight="bold", va="baseline")
    ax.text(0.07, 1.06, text, transform=ax.transAxes, va="baseline")


def fig_graphs():
    fig, axes = plt.subplots(1, 3, figsize=(WIDTH, 1.9))
    fig.subplots_adjust(0.01, 0.02, 0.99, 0.88, wspace=0.05)
    for ax in axes:
        ax.set_xlim(0, 1)
        ax.set_ylim(0, 1)
        ax.axis("off")
    a, b, c = axes
    title(a, "A", "Coupling: what drives what")
    node(a, 0.22, 0.5, "Gene\nexpression")
    node(a, 0.78, 0.5, "Metabolism")
    arrow(a, (0.38, 0.58), (0.62, 0.58), lw=2.0, rad=-0.2)
    a.text(0.5, 0.71, "gain of order one", ha="center", va="bottom", color=MUTED)
    arrow(a, (0.62, 0.42), (0.38, 0.42), lw=0.7, rad=-0.2)
    a.text(0.5, 0.29, "gain 0.04–0.09", ha="center", va="top", color=MUTED)

    title(b, "B", "Likelihood, paths integrated out")
    node(b, 0.22, 0.75, r"$\theta_G$")
    node(b, 0.78, 0.75, r"$\theta_M$")
    node(b, 0.22, 0.25, "Transcript\nfactor")
    node(b, 0.78, 0.25, "Bulk\nfactor")
    arrow(b, (0.22, 0.65), (0.22, 0.4), lw=1.5)
    arrow(b, (0.78, 0.65), (0.78, 0.4), lw=1.5)
    arrow(b, (0.3, 0.66), (0.66, 0.38), lw=1.5)
    arrow(b, (0.7, 0.66), (0.36, 0.38), ls=":", lw=0.8)

    title(c, "C", "Information: what moves what")
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
    fig, ax = canvas(WIDTH, 2.8)
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
        ax.text(0.5, y + 0.02, t, ha="center", va="bottom", color=INK)
    fig.savefig(os.path.join(HERE, "fig_model.pdf"), metadata={"CreationDate": None})
    plt.close(fig)


if __name__ == "__main__":
    fig_architecture()
    fig_model()
    fig_graphs()
