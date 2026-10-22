"""F10 for 16a.10: the production chain against the M0 reference, under the rule of
spec/phases/16-recovery.md §12 2026-10-22.

Reads the TSVs that f10_16a10_export.jl writes. Per marginal, on the log scale:
the Wasserstein-1 distance between production and reference in reference SDs,
beside V8's split-half floor and four null draws matched to production's length
(each V6 chain's sweeps 101-656 against the other three chains). Prints the
result tables as markdown and draws F10.

    python3 dev/scripts/f10_16a10.py   # needs numpy, scipy, matplotlib
"""
import os

import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np
from scipy.stats import gaussian_kde, wasserstein_distance

HERE = os.path.dirname(os.path.abspath(__file__))
OUT = os.path.join(HERE, "..", "notes", "figures", "f10-16a10")
LABELS = {
    "S_JCVISYN3A_0607": "GAPD promoter (S_0607)",
    "S_JCVISYN3A_0779": "ptsG promoter (S_0779)",
    "krnadeg": "krnadeg",
    "kcatF_R_ENO": "kcatF ENO",
    "kcatF_R_FBA": "kcatF FBA",
    "sigma_b": "σ_b",
}
BURN = 100        # V6's merge and production alike
PROD_SWEEPS = 656
REF_CHAINS = (1, 2, 3, 4)
PROD_CHAIN = 5
PAIR = ("S_JCVISYN3A_0779", "kcatF_R_ENO")


def load():
    path = os.path.join(HERE, "f10_16a10_draws.tsv")
    with open(path) as f:
        names = f.readline().rstrip("\n").split("\t")[2:]
    d = np.loadtxt(path, skiprows=1)
    tpath = os.path.join(HERE, "f10_16a10_truth.tsv")
    truth = np.loadtxt(tpath, skiprows=1)
    return names, d, truth


def draws(d, chains, lo, hi):
    """Sweeps lo..hi (1-based, inclusive) of the given chains, pooled."""
    m = np.isin(d[:, 0], chains) & (d[:, 1] >= lo) & (d[:, 1] <= hi)
    return d[m, 2:]


def ess(x):
    """Single-chain ESS by Geyer's initial monotone sequence (not rank-normalised)."""
    x = x - x.mean()
    n = len(x)
    f = np.fft.rfft(x, 2 * n)
    acf = np.fft.irfft(f * np.conj(f))[:n] / (n * x.var())
    pairs = acf[0:n - 1:2] + acf[1:n:2]
    s, prev = 0.0, np.inf
    for p in pairs:
        if p <= 0:
            break
        p = min(p, prev)
        s += p
        prev = p
    return n / (2 * s - 1)


def w1(a, b, sd):
    return wasserstein_distance(a, b) / sd


def main():
    names, d, truth = load()
    P = len(names)
    last = int(d[d[:, 0] == 1, 1].max())
    ref = draws(d, REF_CHAINS, BURN + 1, last)
    prod = draws(d, [PROD_CHAIN], BURN + 1, PROD_SWEEPS)
    sd = ref.std(axis=0, ddof=1)

    # The export reproduces V6's merge: its means and SDs, to compare by eye
    # with v6_16a9b_result.md.
    print(f"Reference: chains {REF_CHAINS}, sweeps {BURN + 1}-{last}, {len(ref)} draws. "
          f"Production: chain {PROD_CHAIN}, sweeps {BURN + 1}-{PROD_SWEEPS}, {len(prod)} draws.\n")

    prod_w = np.array([w1(prod[:, p], ref[:, p], sd[p]) for p in range(P)])
    split = np.array([w1(draws(d, [1, 2], BURN + 1, last)[:, p],
                         draws(d, [3, 4], BURN + 1, last)[:, p], sd[p]) for p in range(P)])
    null = np.zeros((4, P))
    for i, k in enumerate(REF_CHAINS):
        seg = draws(d, [k], BURN + 1, PROD_SWEEPS)
        rest = draws(d, [c for c in REF_CHAINS if c != k], BURN + 1, last)
        null[i] = [w1(seg[:, p], rest[:, p], sd[p]) for p in range(P)]
    rank = 1 + (null < prod_w).sum(axis=0)
    prod_ess = [ess(prod[:, p]) for p in range(P)]
    null_ess = [[ess(draws(d, [k], BURN + 1, PROD_SWEEPS)[:, p]) for k in REF_CHAINS] for p in range(P)]

    print("| Parameter (log) | Reference mean ± SD | Production mean ± SD | Δmean / SD_ref | SD ratio | "
          "Production ESS | Null chains' ESS |")
    print("|---|---|---|---|---|---|---|")
    for p in range(P):
        print(f"| {names[p]} | {ref[:, p].mean():.3f} ± {sd[p]:.3f} | {prod[:, p].mean():.3f} ± "
              f"{prod[:, p].std(ddof=1):.3f} | {(prod[:, p].mean() - ref[:, p].mean()) / sd[p]:+.2f} | "
              f"{prod[:, p].std(ddof=1) / sd[p]:.2f} | {prod_ess[p]:.0f} | "
              f"{', '.join(f'{e:.0f}' for e in null_ess[p])} |")

    print("\n| Parameter (log) | W1 / SD_ref, production | Matched nulls (V6 chains 1-4, sweeps 101-656) | "
          "V8 split-half floor | Production's rank of 5 | Reading |")
    print("|---|---|---|---|---|---|")
    for p in range(P):
        reading = "not resolved" if rank[p] < 5 else f"above all four nulls, by {prod_w[p] - null[:, p].max():.3f}"
        print(f"| {names[p]} | {prod_w[p]:.3f} | {', '.join(f'{v:.3f}' for v in null[:, p])} | "
              f"{split[p]:.3f} | {rank[p]} | {reading} |")
    print(f"\nParameters ranked 5 of 5: {int((rank == 5).sum())} of {P} "
          f"(about {P / 5:.1f} expected by chance if chain 5 is exchangeable with the nulls).")

    os.makedirs(OUT, exist_ok=True)
    plot(names, d, ref, prod, sd, truth, prod_w, null, split, last)


def plot(names, d, ref, prod, sd, truth, prod_w, null, split, last):
    P = len(names)
    fig, axes = plt.subplots(2, 3, figsize=(11, 6.4))
    for p, ax in enumerate(axes.flat):
        lo = min(ref[:, p].min(), prod[:, p].min())
        hi = max(ref[:, p].max(), prod[:, p].max())
        x = np.linspace(lo, hi, 300)
        for k in REF_CHAINS:
            seg = draws(d, [k], BURN + 1, PROD_SWEEPS)[:, p]
            ax.plot(x, gaussian_kde(seg)(x), color="0.7", lw=0.8,
                    label="V6 chain, sweeps 101–656" if k == 1 else None)
        ax.fill_between(x, gaussian_kde(ref[:, p])(x), color="C0", alpha=0.3, label="reference (V6)")
        ax.plot(x, gaussian_kde(prod[:, p])(x), color="C3", lw=1.8, label="production (chain 5)")
        ax.axvline(truth[p], color="k", ls="--", lw=0.8, label="truth")
        ax.set_title(f"{LABELS.get(names[p], names[p])}\nW1 = {prod_w[p]:.2f} SD "
                     f"(nulls ≤ {null[:, p].max():.2f})", fontsize=9)
        ax.set_yticks([])
        ax.set_xlabel("log value", fontsize=8)
    axes.flat[0].legend(fontsize=7, loc="upper left")
    fig.suptitle("F10: production-length chain against the M0 reference (PG, N = 20 for both; D16.5's PG case)",
                 fontsize=10)
    fig.tight_layout()
    for ext in ("pdf", "png"):
        fig.savefig(os.path.join(OUT, f"f10_marginals.{ext}"), dpi=150)
    plt.close(fig)

    i, j = names.index(PAIR[0]), names.index(PAIR[1])
    fig, ax = plt.subplots(figsize=(5, 4.5))
    for data, color, label in ((ref, "C0", "reference (V6)"), (prod, "C3", "production (chain 5)")):
        kde = gaussian_kde(data[:, [i, j]].T)
        gx = np.linspace(ref[:, i].min(), ref[:, i].max(), 120)
        gy = np.linspace(ref[:, j].min(), ref[:, j].max(), 120)
        X, Y = np.meshgrid(gx, gy)
        Z = kde(np.vstack([X.ravel(), Y.ravel()])).reshape(X.shape)
        # Contours at the density levels enclosing about 50% and 90% of the draws.
        dens = np.sort(kde(data[:, [i, j]].T))
        levels = sorted({dens[int(0.10 * len(dens))], dens[int(0.50 * len(dens))]})
        ax.contour(X, Y, Z, levels=levels, colors=color, linewidths=[1.0, 1.8])
        ax.plot([], [], color=color, label=label)
    ax.plot(truth[i], truth[j], "k+", ms=10, label="truth")
    ax.set_xlabel(f"log {LABELS[PAIR[0]]}")
    ax.set_ylabel(f"log {LABELS[PAIR[1]]}")
    ax.set_title("F10: ptsG × ENO, 50% and 90% contours", fontsize=10)
    ax.legend(fontsize=8)
    fig.tight_layout()
    for ext in ("pdf", "png"):
        fig.savefig(os.path.join(OUT, f"f10_contour.{ext}"), dpi=150)
    plt.close(fig)

    fig, ax = plt.subplots(figsize=(7, 3.6))
    y = np.arange(P)
    for k in range(4):
        ax.plot(null[k], y, "o", color="0.6", ms=5, label="matched nulls" if k == 0 else None)
    ax.plot(split, y, "s", mfc="none", color="C0", ms=6, label="V8 split-half floor")
    ax.plot(prod_w, y, "D", color="C3", ms=6, label="production (chain 5)")
    ax.set_yticks(y)
    ax.set_yticklabels([LABELS.get(n, n) for n in names], fontsize=8)
    ax.set_xlabel("Wasserstein-1 against the reference, in reference SDs")
    ax.legend(fontsize=7, loc="upper center", bbox_to_anchor=(0.5, -0.2), ncol=3)
    ax.set_title("F10: divergence per marginal, beside its floors", fontsize=10)
    fig.tight_layout()
    for ext in ("pdf", "png"):
        fig.savefig(os.path.join(OUT, f"f10_divergence.{ext}"), dpi=150)
    plt.close(fig)


if __name__ == "__main__":
    main()
