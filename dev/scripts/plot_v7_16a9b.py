"""F8 for 16a.9b's V7: rank-ECDF difference plots with the merge's simultaneous
bands, one per parameter, and rank histograms as secondary evidence.

Reads the TSVs that v7_16a9b_export.jl writes from the merge's ranks.jls, so the
ECDF and bands are the ones the verdict was computed on (job 18689602).

    python3 dev/scripts/plot_v7_16a9b.py   # needs numpy, scipy, matplotlib
"""
import os

import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np
from scipy.stats import binom

HERE = os.path.dirname(os.path.abspath(__file__))
OUT = os.path.join(HERE, "..", "notes", "figures", "v7-16a9b")
LABELS = {
    "S_JCVISYN3A_0607": "GAPD promoter (S_0607)",
    "S_JCVISYN3A_0779": "ptsG promoter (S_0779)",
    "krnadeg": "krnadeg",
    "kcatF_R_ENO": "kcatF ENO",
    "kcatF_R_FBA": "kcatF FBA",
    "sigma_b": "σ_b",
}
NBINS = 10


def load():
    path = os.path.join(HERE, "v7_16a9b_ranks.tsv")
    with open(path) as f:
        names = f.readline().rstrip("\n").split("\t")[2:]
    data = np.loadtxt(path, skiprows=1)
    bpath = os.path.join(HERE, "v7_16a9b_bands.tsv")
    with open(bpath) as f:
        header = f.readline()   # "# R = .., L = .., α = <per parameter> per parameter, .."
    alpha = float(header.split("α = ")[1].split()[0])
    bands = np.loadtxt(bpath, skiprows=2)
    return names, data[:, 2:].astype(int), bands, alpha


def main():
    names, ranks, bands, alpha = load()
    R, P = ranks.shape
    z, lo, hi = bands[:, 0], bands[:, 1], bands[:, 2]
    L = len(z)
    os.makedirs(OUT, exist_ok=True)

    # Rank-ECDF difference: the merge's counts are #(rank < i) at i = 1..L, against z_i = i/(L+1).
    within = 0
    fig, axes = plt.subplots(2, 3, figsize=(10, 5.6), sharex=True, sharey=True)
    for p, ax in enumerate(axes.flat):
        counts = np.array([(ranks[:, p] < i).sum() for i in range(1, L + 1)])
        within += bool(np.all((lo <= counts) & (counts <= hi)))
        ax.fill_between(z, lo / R - z, hi / R - z, step="mid", color="0.85", lw=0,
                        label=f"simultaneous band, 1 − {alpha:.4f}")
        ax.axhline(0, color="0.5", lw=0.6)
        ax.step(z, counts / R - z, where="mid", color="C0", lw=1.4, label="ECDF − uniform")
        dev = np.max(np.abs(counts / R - z))
        ax.set_title(f"{LABELS.get(names[p], names[p])}\nmax |ECDF − z| = {dev:.3f}", fontsize=10)
        ax.set_xlim(0, 1)
    for ax in axes[1]:
        ax.set_xlabel("fractional rank z")
    for ax in axes[:, 0]:
        ax.set_ylabel("ECDF − z")
    axes.flat[0].legend(fontsize=7.5, loc="lower left", frameon=False)
    fig.suptitle(f"V7 on M0: rank-ECDF differences, R = {R} replicates, L = {L} draws each: "
                 f"{within} of {P} within their bands", fontsize=11)
    fig.tight_layout()
    for ext in ("pdf", "png"):
        fig.savefig(os.path.join(OUT, f"f8_rank_ecdf.{ext}"), dpi=150, metadata={"CreationDate": None} if ext == "pdf" else None)
    plt.close(fig)

    # Histograms, secondary: ranks 0..L in NBINS bins of whole ranks. L + 1 need not
    # divide evenly, so each bin's expectation and pointwise 99% interval use its own
    # share of the L + 1 values.
    edges = np.round(np.linspace(0, L + 1, NBINS + 1)) - 0.5
    pbin = np.diff(edges) / (L + 1)
    blo, bhi = binom.ppf(0.005, R, pbin), binom.ppf(0.995, R, pbin)
    fig, axes = plt.subplots(2, 3, figsize=(10, 5.6), sharex=True, sharey=True)
    for p, ax in enumerate(axes.flat):
        ax.stairs(bhi, edges, baseline=blo, fill=True, color="0.88", lw=0,
                  label="pointwise 99% under uniform")
        ax.hist(ranks[:, p], bins=edges, color="C0", edgecolor="white")
        ax.stairs(R * pbin, edges, color="0.5", lw=0.6)
        ax.set_title(LABELS.get(names[p], names[p]), fontsize=10)
    for ax in axes[1]:
        ax.set_xlabel(f"rank of the truth (of {L})")
    for ax in axes[:, 0]:
        ax.set_ylabel("replicates")
    axes.flat[0].legend(fontsize=7.5, loc="upper right", frameon=False)
    fig.suptitle(f"V7 on M0: rank histograms, R = {R}, {NBINS} bins (secondary to the ECDF plot)",
                 fontsize=11)
    fig.tight_layout()
    for ext in ("pdf", "png"):
        fig.savefig(os.path.join(OUT, f"f8_rank_hist.{ext}"), dpi=150, metadata={"CreationDate": None} if ext == "pdf" else None)
    plt.close(fig)


if __name__ == "__main__":
    main()
