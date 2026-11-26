#!/usr/bin/env python3
"""F4 (fig:reducedpost): the joint posterior over the reduced model's 21 parameters, full width.

Run from anywhere:
    python3 paper/figures/fig_reduced_posterior.py

Inputs, read and never written, both in paper/figures/data/ and copied unchanged from the KIM-92
archive /fred/oz022/tkimpson/infercell-kim92-archive-2026-11-17 (copy of the run's output
directory; result note dev/scripts/corea_pilot_kim92_result.md; driver dev/scripts/corea_pilot_kim92.jl
at de24b27):
- kim92_pilot_coordinates.tsv: pilot_coordinates.tsv, written by merge job 21052425 from 1,600 pooled
  draws (4 chains x 400; chains 20660933-36, extensions 20968974-77); SHA-256
  c7288568a73702805e8769528bfcfc05eb2528fb5d66a102f99fe160c826a8ed, as in the archive's SHA256SUMS.
- kim92_pilot_ppc.tsv: pilot_ppc.tsv, the PPC of PPC job 21052424 (200 draws x 10 replicates, central
  99% band), written by merge job 21052425; SHA-256
  8b1e3ff676b719fc93f0af49b2974480c7200d7a529aa2ddfc796c8942b8b140, as in the archive's SHA256SUMS.
"""
import os

import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np

from make_figures import BAND, INK, LABELS, MIXED, MUTED, WIDTH, rows, save  # also applies rcParams

HERE = os.path.dirname(os.path.abspath(__file__))
DATA = os.path.join(HERE, "data")

OTHER = ["log_krnadeg", "log_kcatF_R_ENO", "log_kcatF_R_FBA", "log_sigma_b"]
OTHER_LABELS = [LABELS["krnadeg"], LABELS["kcatF_R_ENO"], LABELS["kcatF_R_FBA"], LABELS["sigma_b"]]


def promoter_label(name):
    locus = name.removeprefix("log_S_JCVISYN3A_")
    named = {"0607": "GAPD", "0779": "ptsG"}
    return f"{locus} ({named[locus]})" if locus in named else locus


def main():
    c = {r["coordinate"]: r for r in rows(os.path.join(DATA, "kim92_pilot_coordinates.tsv"))}
    prom = [k for k in c if k.startswith("log_S_")]
    assert len(c) == 21 and len(prom) == 17 and set(c) == set(prom) | set(OTHER)
    order = prom + OTHER  # file order within the promoters
    f = lambda k, key: np.array([float(c[n][key]) for n in k])
    mean, sd, truth = f(order, "post_mean"), f(order, "post_sd"), f(order, "truth")
    lo, hi, psd, contr = f(order, "lo90"), f(order, "hi90"), f(order, "prior_sd"), f(order, "contraction")
    z, zlo, zhi = (mean - truth) / sd, (lo - truth) / sd, (hi - truth) / sd
    ratio = sd / psd
    assert np.allclose(contr, 1 - ratio ** 2, atol=1e-9)  # contraction = 1 - (post sd / prior sd)^2
    inside = (lo <= truth) & (truth <= hi)
    miss = [n for n, i in zip(order, inside) if not i]
    zp = z[:17]
    assert inside.sum() == 20 and miss == ["log_S_JCVISYN3A_0234"]
    assert round(contr.min(), 3) == 0.919 and order[int(np.argmin(contr))] == "log_kcatF_R_ENO"
    assert (zp < 0).sum() == 14 and round(zp.mean(), 2) == -0.81

    p = rows(os.path.join(DATA, "kim92_pilot_ppc.tsv"))
    row = np.array([int(r["row"]) for r in p])
    logz, blo, bhi = (np.array([float(r[k]) for r in p]) for k in ("logz", "lo", "hi"))
    out = (logz < blo) | (logz > bhi)
    assert (out == np.array([r["outside"] == "1" for r in p])).all()
    n, limit = len(p), int(np.floor(0.05 * len(p)))
    assert n == 884 and out.sum() == 11 and limit == 44
    u = (logz - blo) / (bhi - blo)  # position within the band: 0 = lower edge, 1 = upper edge

    print(f"contraction min {contr.min():.3f} ({order[int(np.argmin(contr))]}), "
          f"sorted lowest: {np.round(np.sort(contr)[:3], 3)}; promoters min {contr[:17].min():.3f}")
    print(f"truth inside 90%: {inside.sum()} of 21, miss {miss}; promoters below truth "
          f"{(zp < 0).sum()} of 17, mean z {zp.mean():.2f} (sd {zp.std(ddof=1):.2f})")
    print(f"PPC: {out.sum()} of {n} outside (limit {limit}); rows {list(row[out])}; "
          f"band position of outside rows {np.round(u[out], 3)}")

    # y positions: 17 promoters on top, a gap, then the four others.
    y = np.r_[np.arange(17)[::-1] + 4.8, np.arange(4)[::-1]].astype(float)
    ylabels = [promoter_label(k) for k in prom] + OTHER_LABELS
    pcol = np.where(np.arange(21) < 17, MIXED, INK)
    top = 23.3

    fig, (a, b, ax) = plt.subplots(1, 3, figsize=(WIDTH, 3.3), layout="constrained",
                                   gridspec_kw={"width_ratios": [0.7, 1.0, 1.25]})
    a.axvline(1, color=MUTED, lw=0.8, ls="--")
    a.scatter(ratio, y, s=10, color=pcol, lw=0, zorder=3)
    a.set_xscale("log")
    a.set_xlim(0.005, 1.6)
    a.set_xticks([0.01, 0.1, 1], ["0.01", "0.1", "1"])
    a.minorticks_off()
    a.set_xlabel("Posterior SD /\nprior SD")
    a.text(1, 21.2, "prior", ha="right", va="bottom", color=MUTED)
    a.set_yticks(y, ylabels)
    a.set_ylim(-0.8, top)
    a.tick_params(axis="y", length=0)
    # Bracket the 17 promoter strengths, left of their tick labels.
    xb = -0.95
    a.plot([xb, xb], [y[16] - 0.3, y[0] + 0.3], color=MUTED, lw=0.8, transform=a.get_yaxis_transform(),
           clip_on=False)
    a.text(xb - 0.04, (y[0] + y[16]) / 2, "17 promoter strengths", rotation=90, ha="right", va="center",
           color=INK, transform=a.get_yaxis_transform())
    a.text(xb - 0.2, 1.0, "A", transform=a.transAxes, ha="left", va="bottom", fontweight="bold", fontsize=9)

    b.axvspan(-1, 1, color=BAND, lw=0)  # within one posterior SD of the truth
    b.axvline(0, color=MUTED, lw=0.8)
    for yi, l, h, zi, col, ok in zip(y, zlo, zhi, z, pcol, inside):
        b.plot([l, h], [yi, yi], color=col, lw=1.2, solid_capstyle="butt")
        b.scatter(zi, yi, s=12, color=col if ok else "white", edgecolor=col, lw=1.0, zorder=3)
    b.plot([zp.mean()] * 2, [4.3, 22.6], color=MIXED, lw=0.8, ls=":")
    b.text(0.15, 21.2, "truth", ha="left", va="bottom", color=MUTED)
    # In the gap between the promoters and the other four, clear of the spine and of panel A.
    b.text(zp.mean() + 0.12, 3.9, f"promoter mean {zp.mean():.2f}".replace("-", "−"),
           ha="left", va="center", color=MIXED, fontsize=7, bbox=dict(fc="white", ec="none", pad=0.4))
    b.set_ylim(-0.8, top)
    b.set_yticks(y, [])
    b.set_xlim(-4.2, 4.2)
    b.set_xticks([-4, -2, 0, 2, 4])
    b.set_xlabel("(Posterior mean − truth) /\nposterior SD")
    b.tick_params(axis="y", length=0)
    b.set_title("B", loc="left", x=-0.04, fontweight="bold", fontsize=9)

    ax.axhspan(0, 1, color=BAND, lw=0)
    ax.scatter(row[~out], u[~out], s=3, color="#8d8c88", lw=0, label="inside")
    ax.scatter(row[out], u[out], s=14, color=INK, lw=0, zorder=3,
               label="outside")
    ax.set_xlim(-10, n + 11)
    ax.set_xticks([1, 200, 400, 600, 884])
    ax.set_ylim(-0.25, 1.3)
    ax.set_yticks([0, 0.5, 1], ["Lower\nedge", "Middle", "Upper\nedge"])
    ax.set_xlabel(f"Bulk row (of {n})")
    ax.set_ylabel("Observed value within\nthe central 99% band")
    ax.text(0.02, 0.985, f"{out.sum()} of {n} outside (limit {limit})", transform=ax.transAxes,
            ha="left", va="top", color=INK)
    ax.set_title("C", loc="left", x=-0.12, y=1.0, fontweight="bold", fontsize=9)
    save(fig, "fig_reduced_posterior")


if __name__ == "__main__":
    main()
