"""Figures for #89 B7 (AV6), from the TSVs av6_89_b7.jl's merge writes.

Figure 4: the mixed posterior on m0_v6 against V6's (overlaid marginals, the
ptsG-ENO contour) and a cost bar per posterior. Figure 5, middle panel: the
rank-ECDF of the mixed method over V7's datasets, beside the exact sampler's (V7),
against the same simultaneous bands.

    /fred/oz022/tkimpson/miniconda3/bin/python3 dev/scripts/av6_89_b7.py <out>
"""
import os
import sys

import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np
from scipy.stats import gaussian_kde

HERE = os.path.dirname(os.path.abspath(__file__))
FIG = os.path.join(HERE, "..", "notes", "figures", "av6-89b7")
LABELS = {
    "log_S_JCVISYN3A_0607": "GAPD promoter (S_0607)",
    "log_S_JCVISYN3A_0779": "ptsG promoter (S_0779)",
    "log_krnadeg": "krnadeg",
    "log_kcatF_R_ENO": "kcatF ENO",
    "log_kcatF_R_FBA": "kcatF FBA",
    "log_sigma_b": "σ_b",
}
PAIR = ("log_S_JCVISYN3A_0779", "log_kcatF_R_ENO")
F10_CPUH = 507.5     # F10's production posterior (job 18751635)
V6_CPUH = 6809.7     # V6 in all (v6_16a9b_result.md)


def tsv(path):
    with open(path) as f:
        head = f.readline().rstrip("\n").split("\t")
        rows = [l.rstrip("\n").split("\t") for l in f if l.strip()]
    return head, rows


def v6_reference(names):
    head, rows = tsv(os.path.join(HERE, "f10_16a10_draws.tsv"))
    d = np.array(rows, dtype=float)
    keep = np.isin(d[:, 0], [1, 2, 3, 4]) & (d[:, 1] > 100)
    col = {h: j for j, h in enumerate(head)}
    key = lambda n: "sigma_b" if n == "log_sigma_b" else n[4:]
    return np.column_stack([d[keep, col[key(n)]] for n in names])


def save(fig, stem):
    os.makedirs(FIG, exist_ok=True)
    for ext in ("pdf", "png"):
        fig.savefig(os.path.join(FIG, f"{stem}.{ext}"), dpi=150,
                    metadata={"CreationDate": None} if ext == "pdf" else None)
    plt.close(fig)


def figure4(out):
    names, rows = tsv(os.path.join(out, "av6_v6_draws.tsv"))
    mixed = np.array(rows, dtype=float)
    ref = v6_reference(names)
    rh, rrows = tsv(os.path.join(out, "av6_replicates.tsv"))
    cpuh = np.array([float(r[rh.index("cpuh")]) for r in rrows])

    fig = plt.figure(figsize=(12, 6.6))
    gs = fig.add_gridspec(2, 4)
    for p, n in enumerate(names):
        ax = fig.add_subplot(gs[p // 3, p % 3])
        lo = min(ref[:, p].min(), mixed[:, p].min())
        hi = max(ref[:, p].max(), mixed[:, p].max())
        x = np.linspace(lo, hi, 300)
        ax.fill_between(x, gaussian_kde(ref[:, p])(x), color="C0", alpha=0.3, label="exact (V6)")
        ax.plot(x, gaussian_kde(mixed[:, p])(x), color="C3", lw=1.8, label="mixed")
        ax.set_title(LABELS.get(n, n), fontsize=9)
        ax.set_yticks([])
        ax.set_xlabel("log value", fontsize=8)
    fig.axes[0].legend(fontsize=7, loc="upper left")

    ax = fig.add_subplot(gs[0, 3])
    i, j = names.index(PAIR[0]), names.index(PAIR[1])
    for data, colour, lab in ((ref, "C0", "exact (V6)"), (mixed, "C3", "mixed")):
        k = gaussian_kde(data[:, [i, j]].T)
        xs = np.linspace(min(ref[:, i].min(), mixed[:, i].min()), max(ref[:, i].max(), mixed[:, i].max()), 80)
        ys = np.linspace(min(ref[:, j].min(), mixed[:, j].min()), max(ref[:, j].max(), mixed[:, j].max()), 80)
        X, Y = np.meshgrid(xs, ys)
        Z = k(np.vstack([X.ravel(), Y.ravel()])).reshape(X.shape)
        lv = np.sort(Z.ravel())[::-1]
        c = np.cumsum(lv) / lv.sum()
        levels = sorted(lv[np.searchsorted(c, q)] for q in (0.9, 0.5))
        ax.contour(X, Y, Z, levels=levels, colors=colour)
        ax.plot([], [], color=colour, label=lab)
    ax.set_xlabel(LABELS[PAIR[0]], fontsize=8)
    ax.set_ylabel(LABELS[PAIR[1]], fontsize=8)
    ax.set_title("ptsG–ENO, 50% and 90%", fontsize=9)
    ax.legend(fontsize=7)

    ax = fig.add_subplot(gs[1, 3])
    bars = [("exact reference\n(V6)", V6_CPUH), ("exact production\n(F10)", F10_CPUH),
            ("mixed, median\nover V7's 50", float(np.median(cpuh)))]
    ax.bar(range(3), [b[1] for b in bars], color=["C0", "C0", "C3"])
    ax.set_yscale("log")
    ax.set_xticks(range(3), [b[0] for b in bars], fontsize=7)
    ax.set_ylabel("CPU-h per posterior", fontsize=8)
    ax.axhline(F10_CPUH / 10, color="k", ls=":", lw=0.8)
    ax.text(2.4, F10_CPUH / 10, " KC", fontsize=7, va="center")
    fig.suptitle("Figure 4 (AV6): the mixed posterior against the exact sampler's on m0_v6", fontsize=10)
    fig.tight_layout()
    save(fig, "f4_mixed_vs_exact")


def ecdf_diff(ranks, L):
    z = np.arange(1, L + 1) / (L + 1)
    R = len(ranks)
    return z, np.array([(ranks < i).sum() / R for i in range(1, L + 1)]) - z


def figure5(out):
    rh, rrows = tsv(os.path.join(out, "av6_replicates.tsv"))
    names = [h[5:] for h in rh if h.startswith("rank_")]
    mixed = np.array([[int(r[rh.index("rank_" + n)]) for n in names] for r in rrows])
    bh, brows = tsv(os.path.join(out, "av6_bands.tsv"))
    b = np.array(brows, dtype=float)
    z, lo, hi = b[:, 0], b[:, 1], b[:, 2]
    R, L = mixed.shape[0], len(z)
    vh, vrows = tsv(os.path.join(HERE, "v7_16a9b_ranks.tsv"))
    vkey = lambda n: "sigma_b" if n == "log_sigma_b" else n[4:]
    exact = np.array([[int(r[vh.index(vkey(n))]) for n in names] for r in vrows])

    fig, axes = plt.subplots(2, 3, figsize=(11, 6.2), sharex=True, sharey=True)
    for p, ax in enumerate(axes.flat):
        ax.fill_between(z, lo / R - z, hi / R - z, color="0.85", step="mid",
                        label="1 − 0.05/6 band" if p == 0 else None)
        if exact.shape[0] == R:
            zz, d = ecdf_diff(exact[:, p], L)
            ax.step(zz, d, where="mid", color="C0", lw=1, label="exact (V7)" if p == 0 else None)
        zz, d = ecdf_diff(mixed[:, p], L)
        ax.step(zz, d, where="mid", color="C3", lw=1.5, label="mixed" if p == 0 else None)
        ax.axhline(0, color="k", lw=0.5)
        ax.set_title(LABELS.get(names[p], names[p]), fontsize=9)
    axes[0, 0].legend(fontsize=7)
    for ax in axes[1]:
        ax.set_xlabel("fractional rank z", fontsize=8)
    for ax in axes[:, 0]:
        ax.set_ylabel("ECDF − z", fontsize=8)
    fig.suptitle(f"Figure 5, middle: rank-ECDF over V7's {R} datasets, mixed beside exact", fontsize=10)
    fig.tight_layout()
    save(fig, "f5_rank_ecdf_m0")


if __name__ == "__main__":
    out = sys.argv[1]
    figure4(out)
    figure5(out)
    print(f"wrote figures to {os.path.normpath(FIG)}")
