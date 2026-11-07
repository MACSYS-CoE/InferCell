#!/usr/bin/env python3
"""The paper's data figures, in one style, written to this directory as PDFs.

Run from the repository root:
    python3 paper/figures/make_figures.py [EVIDENCE_ARCHIVE] [B7_RUN_DIR]

Inputs, read and never written:
- the KIM-43 evidence archive: runs/m0-kim40-kim41/score_datasets.tsv and
  runs/budget-kim42/budget_rows.tsv;
- paper/figures/data: av6_replicates.tsv, av6_bands.tsv, av6_v6_draws.tsv, copied from B7's run;
- dev/scripts/f10_16a10_draws.tsv and dev/scripts/v7_16a9b_ranks.tsv.
"""
import csv
import os
import sys

import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np
from scipy.stats import gaussian_kde

HERE = os.path.dirname(os.path.abspath(__file__))
ROOT = os.path.dirname(os.path.dirname(HERE))
ARCHIVE = "/fred/oz022/tkimpson/infercell-predict-archive-2026-11-06"
B7_RUN = os.path.join(HERE, "data")  # copies of three B7 outputs, unchanged

# Validated together: slot 1 and 2 of the reference categorical palette.
EXACT, MIXED = "#2a78d6", "#eb6834"
INK, MUTED, BAND = "#0b0b0b", "#52514e", "#e4e3df"
WIDTH = 6.5  # PLOS text width with \newgeometry, inches

plt.rcParams.update({
    "font.family": "sans-serif", "font.size": 8, "axes.titlesize": 8, "axes.labelsize": 8,
    "xtick.labelsize": 8, "ytick.labelsize": 8, "legend.fontsize": 8, "legend.frameon": False,
    "axes.spines.top": False, "axes.spines.right": False, "axes.edgecolor": MUTED,
    "xtick.color": MUTED, "ytick.color": MUTED, "axes.labelcolor": INK, "lines.linewidth": 1.5,
    "pdf.fonttype": 42,
})

NAMES = ["S_JCVISYN3A_0607", "S_JCVISYN3A_0779", "krnadeg", "kcatF_R_ENO", "kcatF_R_FBA", "sigma_b"]
LABELS = dict(zip(NAMES, ["GAPD promoter", "ptsG promoter", "Transcript decay", "Enolase kcat",
                          "Aldolase kcat", "Assay noise"]))


def rows(path):
    with open(path, newline="") as f:
        return list(csv.DictReader(f, delimiter="\t"))


def save(fig, stem):
    fig.savefig(os.path.join(HERE, stem + ".pdf"), metadata={"CreationDate": None})
    plt.close(fig)


def panel(fig, letter):
    fig.suptitle(letter, x=0.005, ha="left", fontweight="bold", fontsize=9, color=INK)


def legend_row(fig, axes):
    h, l = axes[0].get_legend_handles_labels()
    fig.legend(h, l, loc="outside lower center", ncol=len(h))


def ecdf_diff(ranks, L):
    z = np.arange(1, L + 1) / (L + 1)
    return z, np.array([(ranks < i).sum() / len(ranks) for i in range(1, L + 1)]) - z


def bands(b7):
    b = np.array([[float(r[k]) for k in r] for r in rows(os.path.join(b7, "av6_bands.tsv"))])
    return b[:, 0], b[:, 1], b[:, 2]


def exact_ranks():
    v = rows(os.path.join(ROOT, "dev/scripts/v7_16a9b_ranks.tsv"))
    return np.array([[int(r[n]) for n in NAMES] for r in v])


def mixed_ranks(b7):
    v = rows(os.path.join(b7, "av6_replicates.tsv"))
    key = lambda n: "rank_log_" + n
    return np.array([[int(r[key(n)]) for n in NAMES] for r in v])


def reference_draws():
    d = rows(os.path.join(ROOT, "dev/scripts/f10_16a10_draws.tsv"))
    ref = np.array([[float(r[n]) for n in NAMES] for r in d
                    if r["chain"] in "1234" and int(r["sweep"]) > 100])
    prod = np.array([[float(r[n]) for n in NAMES] for r in d
                     if r["chain"] == "5" and int(r["sweep"]) > 100])
    return ref, prod


def density_row(axes, ref, other, label, colour):
    for p, ax in enumerate(axes):
        x = np.linspace(min(ref[:, p].min(), other[:, p].min()),
                        max(ref[:, p].max(), other[:, p].max()), 300)
        ax.fill_between(x, gaussian_kde(ref[:, p])(x), color=EXACT, alpha=0.25, lw=0,
                        label="Exact reference")
        ax.plot(x, gaussian_kde(other[:, p])(x), color=colour, label=label)
        ax.set_title(LABELS[NAMES[p]])
        ax.set_yticks([])
        ax.spines["left"].set_visible(False)
        ax.set_xlabel("ln(value)")


def ecdf_row(axes, z, lo, hi, curves):
    for p, ax in enumerate(axes):
        R = len(curves[0][0])
        ax.fill_between(z, lo / R - z, hi / R - z, color=BAND, lw=0, step="mid",
                        label="Band, joint 95% over six parameters")
        for ranks, colour, label in curves:
            zz, d = ecdf_diff(ranks[:, p], len(z))
            ax.step(zz, d, where="mid", color=colour, label=label)
        ax.axhline(0, color=MUTED, lw=0.5)
        ax.set_title(LABELS[NAMES[p]])
        ax.set_xlabel("Fractional rank")
    axes[0].set_ylabel("ECDF − fractional rank")


def fig_calibration(b7):
    z, lo, hi = bands(b7)
    ref, prod = reference_draws()
    fig = plt.figure(figsize=(WIDTH, 4.2), layout="constrained")
    top, bottom = fig.subfigures(2, 1)
    a = top.subplots(1, 6, sharey=True)
    ecdf_row(a, z, lo, hi, [(exact_ranks(), EXACT, "Exact sampler")])
    legend_row(top, a)
    panel(top, "A")
    b = bottom.subplots(1, 6)
    density_row(b, ref, prod, "Production-length chain", MIXED)
    legend_row(bottom, b)
    panel(bottom, "B")
    save(fig, "fig_calibration")


def fig_mixed(b7):
    z, lo, hi = bands(b7)
    ref, _ = reference_draws()
    m = rows(os.path.join(b7, "av6_v6_draws.tsv"))
    mixed = np.array([[float(r["log_" + n]) for n in NAMES] for r in m])
    fig = plt.figure(figsize=(WIDTH, 4.2), layout="constrained")
    top, bottom = fig.subfigures(2, 1)
    a = top.subplots(1, 6, sharey=True)
    ecdf_row(a, z, lo, hi, [(exact_ranks(), EXACT, "Exact"), (mixed_ranks(b7), MIXED, "Mixed")])
    legend_row(top, a)
    panel(top, "A")
    b = bottom.subplots(1, 6)
    density_row(b, ref, mixed, "Mixed", MIXED)
    legend_row(bottom, b)
    panel(bottom, "B")
    save(fig, "fig_mixed")


def fig_budget(archive):
    d = rows(os.path.join(archive, "runs/budget-kim42/budget_rows.tsv"))
    get = lambda w, q: next((float(r["estimate"]), float(r["se"])) for r in d
                            if float(r["width"]) == w and r["target"] == "q"
                            and r["law"] == "per_theta_I" and r["quantity"] == q)
    widths = [0.5, 1.0, 2.0]
    tick = ["×½", "×1", "×2"]
    fig, (a, b) = plt.subplots(1, 2, figsize=(WIDTH, 2.4), layout="constrained")
    x = np.arange(3)
    for q, off, marker, colour, label in (("share_i", -0.08, "o", EXACT, r"$\theta_I$ first"),
                                          ("share_h", 0.08, "s", MIXED, r"$\theta_H$ first")):
        v = [get(w, q) for w in widths]
        a.errorbar(x + off, [e[0] for e in v], yerr=[e[1] for e in v], fmt=marker, ms=5,
                   color=colour, capsize=0, label=label)
    a.set_xticks(x, tick)
    a.set_xlim(-0.5, 2.5)
    a.set_ylim(0, 1.05)
    a.set_xlabel(r"Prior width of $\theta_H$")
    a.set_ylabel("Held share of Var(μ)")
    a.legend(loc="lower right")
    a.set_title("A", loc="left", fontweight="bold", fontsize=9)
    a.set_title("UNSCORED", loc="right", color=MUTED, fontweight="bold")
    mean = [get(w, "mean") for w in widths]
    total = [np.sqrt(get(w, "total")[0]) for w in widths]
    b.errorbar(x, [m[0] / 1e6 for m in mean], yerr=[t / 1e6 for t in total], fmt="o", ms=5,
               color=EXACT)
    b.set_xticks(x, tick)
    b.set_xlim(-0.5, 2.5)
    b.set_ylim(bottom=0)
    b.set_xlabel(r"Prior width of $\theta_H$")
    b.set_ylabel("Lactate exported, millions of particles")
    b.set_title("B", loc="left", fontweight="bold", fontsize=9)
    b.set_title("UNSCORED", loc="right", color=MUTED, fontweight="bold")
    save(fig, "fig_budget")


ARMS = [("v7", "hybrid", "Exact, hybrid"), ("v7", "mean_field", "Exact, mean field"),
        ("v7", "mean_field_delta", "Exact, mean field + δ"), ("b7", "hybrid", "Mixed, hybrid"),
        ("b7", "mean_field", "Mixed, mean field"), ("b7", "mean_field_delta", "Mixed, mean field + δ")]


def fig_crossed(archive):
    d = rows(os.path.join(archive, "runs/m0-kim40-kim41/score_datasets.tsv"))
    fig, axes = plt.subplots(2, 3, figsize=(WIDTH, 3.6), layout="constrained", sharey=True)
    yy = np.arange(len(ARMS))[::-1]
    for i, (target, name) in enumerate((("q", "Q, lactate export"), ("dq", "ΔQ, knockdown effect"))):
        for k, (post, arm, _) in enumerate(ARMS):
            sel = [r for r in d if r["target"] == target and r["posterior"] == post and r["arm"] == arm]
            assert len(sel) == 50, (target, post, arm, len(sel))
            colour = EXACT if post == "v7" else MIXED
            marker = "o" if post == "v7" else "s"
            for j, key in enumerate(("cover50", "cover90")):
                p = np.mean([r[key] == "true" for r in sel])
                se = np.sqrt(p * (1 - p) / 50)
                axes[i, j].errorbar(100 * p, yy[k], xerr=100 * se, fmt=marker, ms=4, color=colour)
            w = np.array([float(r["width_ratio"]) for r in sel])
            q = np.percentile(w, [25, 50, 75])
            axes[i, 2].errorbar(q[1], yy[k], xerr=[[q[1] - q[0]], [q[2] - q[1]]], fmt=marker, ms=4,
                                color=colour)
        for j, (ref, title) in enumerate(((49, "Covered at nominal 50%"), (88.2, "Covered at nominal 90%"))):
            axes[i, j].axvline(ref, color=MUTED, lw=0.8, ls="--")
            axes[i, j].set_title(title)
        axes[i, 0].set_xlim(10, 90)
        axes[i, 1].set_xlim(60, 105)
        axes[i, 2].axvline(1, color=MUTED, lw=0.8, ls="--")
        axes[i, 2].set_title("Width of 90% interval")
        axes[i, 0].set_yticks(yy, [a[2] for a in ARMS])
        axes[i, 0].set_ylabel(name)
        axes[i, 0].text(-0.75, 1.04, "AB"[i], transform=axes[i, 0].transAxes, fontweight="bold", fontsize=9)
    for ax in axes[1, :2]:
        ax.set_xlabel("Percent of 50 datasets")
    axes[1, 2].set_xlabel("Median ratio and IQR")
    fig.suptitle("SCORED against the truth", x=0.995, ha="right", color=MUTED, fontweight="bold", fontsize=8)
    save(fig, "fig_crossed")


if __name__ == "__main__":
    archive = sys.argv[1] if len(sys.argv) > 1 else ARCHIVE
    b7 = sys.argv[2] if len(sys.argv) > 2 else B7_RUN
    fig_calibration(b7)
    fig_mixed(b7)
    fig_budget(archive)
    fig_crossed(archive)
