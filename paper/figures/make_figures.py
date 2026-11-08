#!/usr/bin/env python3
"""The paper's data figures, in one style, written to this directory as PDFs.

Run from the repository root:
    python3 paper/figures/make_figures.py [EVIDENCE_ARCHIVE] [B7_RUN_DIR]

Inputs, read and never written:
- the KIM-43 evidence archive: runs/m0-kim40-kim41/score_datasets.tsv and
  runs/budget-kim42/budget_rows.tsv;
- paper/figures/data: av6_replicates.tsv, av6_bands.tsv, av6_v6_draws.tsv, copied from B7's run,
  budget_mu.tsv, from the archive by paper/figures/budget_dump.jl, and traj_*.tsv, copied unchanged
  from the trajectory runs' tsv/ (infercell-trajectories-2026-11-08; dev/scripts/trajectories_predict_result.md on main, #159);
- dev/scripts/f10_16a10_draws.tsv and dev/scripts/v7_16a9b_ranks.tsv.

The style, sizes and colours are in style.py. Layouts are placed in inches, so the columns of
stacked rows line up and the printed sizes are the ones set here.
"""
import csv
import os
import sys

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
from style import BAND, EXACT, INK, LETTER, MIXED, MUTED, WIDTH, letter  # noqa: E402  (sets rcParams)

import matplotlib.pyplot as plt  # noqa: E402
import numpy as np  # noqa: E402
from scipy.stats import gaussian_kde  # noqa: E402

HERE = os.path.dirname(os.path.abspath(__file__))
ROOT = os.path.dirname(os.path.dirname(HERE))
ARCHIVE = "/fred/oz022/tkimpson/infercell-predict-archive-2026-11-06"
B7_RUN = os.path.join(HERE, "data")  # copies of three B7 outputs, unchanged

NAMES = ["S_JCVISYN3A_0607", "S_JCVISYN3A_0779", "krnadeg", "kcatF_R_ENO", "kcatF_R_FBA", "sigma_b"]
LABELS = dict(zip(NAMES, ["GAPD promoter", "ptsG promoter", "Transcript decay", r"Enolase $k_\mathrm{cat}$",
                          r"Aldolase $k_\mathrm{cat}$", "Assay noise"]))
DASH = (0, (3, 1.5))  # the truth and equality lines


def rows(path):
    with open(path, newline="") as f:
        return list(csv.DictReader(f, delimiter="\t"))


def save(fig, stem):
    fig.savefig(os.path.join(HERE, stem + ".pdf"), metadata={"CreationDate": None})
    plt.close(fig)


def grid_row(fig, H, top, height, n, left=0.5, right=0.06, gap=0.12, ratios=None, sharey=False):
    """One row of n axes, placed in inches from the top of a figure H inches high."""
    gs = fig.add_gridspec(1, n, left=left / WIDTH, right=1 - right / WIDTH, top=1 - top / H,
                          bottom=1 - (top + height) / H, width_ratios=ratios,
                          wspace=gap * n / (WIDTH - left - right - gap * (n - 1)))
    axes = [fig.add_subplot(gs[0, 0])]
    for i in range(1, n):
        axes.append(fig.add_subplot(gs[0, i], sharey=axes[0] if sharey else None))
    return axes


def letters_on(fig, pairs, dx=-0.02, dy=0.0):
    """Panel letters level with the top of each axes, after the layout has run."""
    fig.canvas.draw()
    for ax, text in pairs:
        p = ax.get_position()
        letter(fig, max(p.x0 + dx, 0.004), p.y1 + dy, text)


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
        ax.set_yticks([])
        ax.spines["left"].set_visible(False)
        ax.set_ylim(bottom=0)
        ax.locator_params(axis="x", nbins=3)
    return "Natural log of the parameter"


def ecdf_row(axes, z, lo, hi, curves):
    R = len(curves[0][0])
    for p, ax in enumerate(axes):
        ax.fill_between(z, lo / R - z, hi / R - z, color=BAND, lw=0, step="mid",
                        label="Band, joint 95% over six parameters")
        for ranks, colour, label in curves:
            zz, d = ecdf_diff(ranks[:, p], len(z))
            ax.step(zz, d, where="mid", color=colour, lw=0.9, label=label)
        ax.axhline(0, color=MUTED, lw=0.5)
        ax.set_xticks([0, 0.5, 1], ["0", "0.5", "1"])
        if p:
            ax.tick_params(labelleft=False)
    axes[0].set_ylabel("ECDF − fractional rank")
    return "Fractional rank"


def six_by_two(draw_first, draw_second, stem):
    """Two rows of the six parameters, columns aligned, titles once, one x label and key per row."""
    H = 3.75
    fig = plt.figure(figsize=(WIDTH, H))
    top = grid_row(fig, H, 0.32, 0.95, 6, left=0.52, right=0.14, gap=0.18)
    bottom = grid_row(fig, H, 2.08, 1.0, 6, left=0.52, right=0.14, gap=0.18)
    for ax, name in zip(top, NAMES):
        ax.set_title(LABELS[name])
    for axes, draw, y0, text in ((top, draw_first, 0.32 + 0.95, "A"), (bottom, draw_second, 2.08 + 1.0, "B")):
        xlabel = draw(axes)
        fig.text(0.5 * (axes[0].get_position().x0 + axes[-1].get_position().x1), 1 - (y0 + 0.2) / H, xlabel,
                 ha="center", va="top")
        h, l = axes[0].get_legend_handles_labels()
        fig.legend(h, l, loc="center", ncol=len(h), bbox_to_anchor=(0.5, 1 - (y0 + 0.44) / H))
    letter(fig, 0.004, 1 - 0.02 / H, "A")
    letter(fig, 0.004, 1 - 1.93 / H, "B")
    save(fig, stem)


def fig_calibration(b7):
    z, lo, hi = bands(b7)
    ref, prod = reference_draws()
    six_by_two(lambda a: ecdf_row(a, z, lo, hi, [(exact_ranks(), EXACT, "Exact sampler")]),
               lambda a: density_row(a, ref, prod, "Production-length chain", MIXED),
               "fig_calibration")


def fig_mixed(b7):
    z, lo, hi = bands(b7)
    ref, _ = reference_draws()
    m = rows(os.path.join(b7, "av6_v6_draws.tsv"))
    mixed = np.array([[float(r["log_" + n]) for n in NAMES] for r in m])
    six_by_two(lambda a: density_row(a, ref, mixed, "Mixed", MIXED),
               lambda a: ecdf_row(a, z, lo, hi, [(exact_ranks(), EXACT, "Exact"), (mixed_ranks(b7), MIXED, "Mixed")]),
               "fig_mixed")


def unscored(ax):
    ax.set_title("UNSCORED", loc="right", color=MUTED, fontweight="bold")


def fig_budget(archive):
    d = rows(os.path.join(archive, "runs/budget-kim42/budget_rows.tsv"))
    get = lambda w, q: next((float(r["estimate"]), float(r["se"])) for r in d
                            if float(r["width"]) == w and r["target"] == "q"
                            and r["law"] == "per_theta_I" and r["quantity"] == q)
    widths = [0.5, 1.0, 2.0]
    tick = ["×½", "×1", "×2"]
    fig, (a, b) = plt.subplots(1, 2, figsize=(WIDTH, 2.2), layout="constrained", width_ratios=[1, 1.4])
    fig.get_layout_engine().set(w_pad=0.06, wspace=0.08)
    x = np.arange(3)
    for q, off, marker, colour, label in (("share_i", -0.08, "o", EXACT, r"$\theta_I$ first"),
                                          ("share_h", 0.08, "s", MIXED, r"$\theta_H$ first")):
        v = [get(w, q) for w in widths]
        a.errorbar(x + off, [e[0] for e in v], yerr=[e[1] for e in v], fmt=marker, ms=4,
                   color=colour, capsize=0, lw=1, label=label)
    a.set_xticks(x, tick)
    a.set_xlim(-0.5, 2.5)
    a.set_ylim(0, 1.05)
    a.set_xlabel(r"Prior width of $\theta_H$")
    a.set_ylabel(r"Held share of Var($\mu$)")
    a.legend(loc="lower right")
    unscored(a)
    # Each point averages many pairs, so it shows one parameter block's reach, with the other
    # block averaged out. Pairs the viability rule rejected are left out.
    mu = rows(os.path.join(HERE, "data", "budget_mu.tsv"))
    rng = np.random.default_rng(1)
    for k, w in enumerate(widths):
        kept = [r for r in mu if float(r["width"]) == w and r["viable"] == "true"]
        for key, off, marker, colour, label in (("row", -0.16, "o", EXACT, r"One draw of $\theta_I$"),
                                                ("column", 0.16, "s", MIXED, r"One draw of $\theta_H$")):
            groups = {}
            for r in kept:
                groups.setdefault(r[key], []).append(float(r["mu"]))
            y = np.array([np.mean(v) for v in groups.values()]) / 1e6
            b.scatter(k + off + rng.uniform(-0.09, 0.09, len(y)), y, s=6, marker=marker, color=colour,
                      alpha=0.6, lw=0, label=label if k == 0 else None)
    b.set_xticks(x, tick)
    b.set_xlim(-0.5, 2.5)
    b.set_ylim(bottom=0)
    b.set_xlabel(r"Prior width of $\theta_H$")
    b.set_ylabel("Lactate exported (millions)")
    b.legend(loc="upper left", handletextpad=0.2, markerscale=1.5)
    unscored(b)
    letters_on(fig, [(a, "A"), (b, "B")], dx=-0.06, dy=0.055)
    save(fig, "fig_budget")


def fig_predict_m0(archive):
    d = rows(os.path.join(archive, "runs/m0-kim40-kim41/score_datasets.tsv"))
    # Four square panels in two pairs; the second of each pair shares the first's y axis.
    left, right, gap, pair_gap, top, bottom = 0.5, 0.06, 0.24, 0.55, 0.22, 0.36
    s = (WIDTH - left - right - 2 * gap - pair_gap) / 4
    H = s + top + bottom
    x0 = [left, left + s + gap, left + 2 * s + gap + pair_gap, left + 3 * s + 2 * gap + pair_gap]
    fig = plt.figure(figsize=(WIDTH, H))
    axes = [fig.add_axes([x / WIDTH, bottom / H, s / WIDTH, s / H]) for x in x0]
    panels = (("q", "v7", 1e6, "Lactate export, exact", EXACT, "o", r"$Q$", "millions"),
              ("q", "b7", 1e6, "Lactate export, mixed", MIXED, "s", r"$Q$", "millions"),
              ("dq", "v7", 1e3, "Knockdown effect, exact", EXACT, "o", r"$\Delta Q$", "thousands"),
              ("dq", "b7", 1e3, "Knockdown effect, mixed", MIXED, "s", r"$\Delta Q$", "thousands"))
    for i, (ax, (target, post, scale, title, colour, marker, sym, unit)) in enumerate(zip(axes, panels)):
        sel = [r for r in d if r["target"] == target and r["posterior"] == post and r["arm"] == "hybrid"]
        assert len(sel) == 50, (target, post, len(sel))
        t, m, lo, hi = (np.array([float(r[k]) for r in sel]) / scale
                        for k in ("truth", "mean", "lo90", "hi90"))
        # One range per quantity, over both posteriors, so panels of a quantity compare directly.
        both = [r for r in d if r["target"] == target and r["arm"] == "hybrid"]
        low = min(min(float(r["truth"]), float(r["lo90"])) for r in both) / scale
        high = max(max(float(r["truth"]), float(r["hi90"])) for r in both) / scale
        pad = 0.05 * (high - low)
        lim = [low - pad, high + pad]
        ax.plot(lim, lim, color=MUTED, lw=0.7, ls=DASH)
        ax.errorbar(t, m, yerr=[m - lo, hi - m], fmt=marker, ms=2.5, lw=0.7, color=colour)
        ax.set_xlim(lim)
        ax.set_ylim(lim)
        ax.locator_params(nbins=4)
        ax.set_title(title, loc="left")
        ax.set_xlabel(f"True {sym} ({unit})")
        if i % 2 == 0:
            ax.set_ylabel(f"Predicted {sym} ({unit})")
        else:
            ax.tick_params(labelleft=False)
        letter(fig, (x0[i] - (0.42 if i % 2 == 0 else 0.2)) / WIDTH, 1 - 0.01 / H, "ABCD"[i])
    save(fig, "fig_predict_m0")


ARMS = [("v7", "hybrid", "Exact, hybrid"), ("v7", "mean_field", "Exact, mean field"),
        ("v7", "mean_field_delta", "Exact, mean field + δ"), ("b7", "hybrid", "Mixed, hybrid"),
        ("b7", "mean_field", "Mixed, mean field"), ("b7", "mean_field_delta", "Mixed, mean field + δ")]


def fig_crossed(archive):
    d = rows(os.path.join(archive, "runs/m0-kim40-kim41/score_datasets.tsv"))
    fig, axes = plt.subplots(2, 3, figsize=(WIDTH, 3.3), layout="constrained", sharey=True, sharex="col")
    fig.get_layout_engine().set(h_pad=0.08)
    yy = np.arange(len(ARMS))[::-1]
    for i, (target, name) in enumerate((("q", r"$Q$, lactate export"), ("dq", r"$\Delta Q$, knockdown effect"))):
        for k, (post, arm, _) in enumerate(ARMS):
            sel = [r for r in d if r["target"] == target and r["posterior"] == post and r["arm"] == arm]
            assert len(sel) == 50, (target, post, arm, len(sel))
            colour = EXACT if post == "v7" else MIXED
            marker = "o" if post == "v7" else "s"
            for j, key in enumerate(("cover50", "cover90")):
                p = np.mean([r[key] == "true" for r in sel])
                se = np.sqrt(p * (1 - p) / 50)
                axes[i, j].errorbar(100 * p, yy[k], xerr=100 * se, fmt=marker, ms=3.5, lw=1, color=colour)
            w = np.array([float(r["width_ratio"]) for r in sel])
            q = np.percentile(w, [25, 50, 75])
            axes[i, 2].errorbar(q[1], yy[k], xerr=[[q[1] - q[0]], [q[2] - q[1]]], fmt=marker, ms=3.5, lw=1,
                                color=colour)
        for j, ref in enumerate((49, 88.2)):
            axes[i, j].axvline(ref, color=MUTED, lw=0.7, ls=DASH)
        axes[i, 2].axvline(1, color=MUTED, lw=0.7, ls=DASH)
        axes[i, 0].set_yticks(yy, [a[2] for a in ARMS])
        axes[i, 0].set_ylabel(name)
        axes[i, 0].set_ylim(-0.6, len(ARMS) - 0.4)
    for j, title in enumerate(("Covered at nominal 50%", "Covered at nominal 90%", "Width of 90% interval")):
        axes[0, j].set_title(title)
    axes[0, 0].set_xlim(10, 90)
    axes[0, 1].set_xlim(60, 105)
    for ax in axes[1, :2]:
        ax.set_xlabel("Percent of 50 datasets")
    axes[1, 2].set_xlabel("Ratio to exact, hybrid (median, IQR)")
    fig.suptitle("SCORED", x=0.995, ha="right", color=MUTED, fontweight="bold")
    letters_on(fig, [(axes[0, 0], "A"), (axes[1, 0], "B")], dx=-1, dy=0.05)  # at the left edge
    save(fig, "fig_crossed")


# Exported lactate over the change in M_lac__L_e, from the run's replay.md (job 20271501). It puts the
# bulk lactate data on the export's axis; its largest deviation over time and cells was 0.245%.
RATIO = 99920.98568131763
# GTP's K_d in particles, a lower and an upper bound: reduced_budget_kim42_result.md, "K_d diagnostic".
KD = (2018, 2160)
LIGHT = "#a9c8ee"  # lighter steps of EXACT's hue, for outer bands and layers
LAYER2, LAYER3 = "#7aaae6", "#cfe0f6"
PALE = "#f8c9b4"  # a light step of MIXED's hue


def traj(name):
    return rows(os.path.join(HERE, "data", name + ".tsv"))


def band(d, species, layer, weights="uniform", arm=None):
    sel = sorted((r for r in d if r["species"] == species and r["layer"] == layer and r["weights"] == weights
                  and (arm is None or r["arm"] == arm)), key=lambda r: float(r["t"]))
    assert sel, (species, layer, weights, arm)
    return {k: np.array([float(r[k]) for r in sel]) for k in ("t", "mean", "lo50", "hi50", "lo90", "hi90")}


def shade(ax, b, scale, colour, alpha=1.0, inner=None, label=None, t=None):
    t = b["t"] / 60 if t is None else t
    ax.fill_between(t, b["lo90"] / scale, b["hi90"] / scale, color=colour, alpha=alpha, lw=0, label=label)
    if inner:
        ax.fill_between(t, b["lo50"] / scale, b["hi50"] / scale, color=inner, lw=0)


def kd(ax):
    ax.axhspan(KD[0] / 1e3, KD[1] / 1e3, color=MUTED, alpha=0.35, lw=0)
    ax.text(1.01, sum(KD) / 2e3, r"$K_d$", transform=ax.get_yaxis_transform(), ha="left", va="center",
            color=MUTED)


def fig_trajectories():
    a, truth, data, b, c, split = (traj("traj_" + k) for k in ("a", "a_truth", "a_data", "b", "c", "split"))
    tru = lambda s: band_truth(truth, s)
    obs = lambda s: (np.array([float(r["t"]) for r in data if r["species"] == s]),
                     np.array([float(r["value"]) for r in data if r["species"] == s]))
    # The truth's own ratio should sit within the replay's 1% bound of RATIO.
    tt, lac = tru("M_lac__L_e")
    _, meter = tru("lac_export_meter")
    assert np.all(np.abs(meter[1:] / lac[1:] / RATIO - 1) < 0.01)

    # One column geometry for all three rows, so the columns line up; the fourth holds the key.
    H, height, step = 6.0, 1.2, 2.0
    fig = plt.figure(figsize=(WIDTH, H))
    grid = [grid_row(fig, H, 0.38 + k * step, height, 4, left=0.48, right=0.02, gap=0.6,
                     ratios=[1, 1, 1, 0.95]) for k in range(3)]
    ax = np.array([g[:3] for g in grid])
    key = [g[3] for g in grid]
    for k in key:
        k.axis("off")
    ax[2, 2].axis("off")
    data_style = dict(ls="none", marker="o", ms=3, mfc="white", mec=INK, mew=0.6)

    # A: M0 over 600 s, the exact posterior on V7 replicate 1, against its truth.
    row = ax[0]
    nb = band(a, "M_lac__L_e", "new_bulk")
    nb = {k: (v if k == "t" else v * RATIO) for k, v in nb.items()}
    shade(row[0], nb, 1e3, LIGHT, label="New bulk\nmeasurement")
    shade(row[0], band(a, "lac_export_meter", "rows"), 1e3, EXACT, label="Posterior")
    t, v = obs("M_lac__L_e")
    row[0].plot(t / 60, v * RATIO / 1e3, label="Data", **data_style)
    row[0].plot(tt / 60, meter / 1e3, color=INK, lw=0.8, ls=DASH, label="Truth")
    row[0].set_ylabel("Particles (thousands)")
    row[0].set_title("Lactate exported (data: lactate out)")
    shade(row[1], band(a, "M_gtp_c", "new_bulk"), 1e3, LIGHT)
    shade(row[1], band(a, "M_gtp_c", "rows"), 1e3, EXACT)
    t, v = obs("M_gtp_c")
    row[1].plot(t / 60, v / 1e3, **data_style)
    tg, g = tru("M_gtp_c")
    row[1].plot(tg / 60, g / 1e3, color=INK, lw=0.8, ls=DASH)
    row[1].set_ylabel("Particles (thousands)")
    row[1].set_title("GTP (observed)")
    # Pyruvate rises by a third in the first minute and then holds, so its band shows only against the truth.
    p = band(a, "M_pyr_c", "rows")
    tp, ptrue = tru("M_pyr_c")
    assert np.array_equal(tp, p["t"])
    later = tp > 0
    off = {k: (v[later] if k == "t" else v[later] - ptrue[later]) for k, v in p.items()}
    shade(row[2], off, 1, EXACT)
    row[2].axhline(0, color=INK, lw=0.8, ls=DASH)
    row[2].set_ylabel("Deviation from truth\n(particles)")
    row[2].set_title("Pyruvate (not observed)")
    for x in row:
        x.set_xlim(0, 10)
        x.set_xlabel("Time (min)")
    key[0].legend(*row[0].get_legend_handles_labels(), loc="center left", title="90% bands",
                  alignment="left", bbox_to_anchor=(-0.3, 0.5))

    # B: the reduced model over one cycle, 20 parameters drawn, in three nested layers.
    row = ax[1]
    layers = (("cells", LAYER3, "A new cell"), ("pairs", LAYER2, r"$\theta_I$ and $\theta_H$"),
              ("rows", EXACT, r"$\theta_I$ only"))
    for k, (species, scale) in enumerate((("lac_export_meter", 1e6), ("M_gtp_c", 1e3))):
        for layer, colour, label in layers:
            shade(row[k], band(b, species, layer), scale, colour, label=label)
        m = band(b, species, "rows")
        row[k].plot(m["t"] / 60, m["mean"] / scale, color=INK, lw=0.9, label="Predictive mean")
    end = band(b, "lac_export_meter", "rows")["mean"][-1]
    assert round(end, -2) == 9576500, end  # KIM-42's mean Q at width 1
    row[0].set_ylabel("Particles (millions)")
    row[0].set_title("Lactate exported")
    row[1].set_ylabel("Particles (thousands)")
    row[1].set_title("GTP")
    kd(row[1])
    shown = [r for r in split if r["flagged"] == "false" and r["mc_undefined"] == "false"
             and float(r["total"]) > float(r["mc"])]
    ts = np.array([float(r["t"]) for r in shown]) / 60
    for name, colour, label in (("held_share_i_first", EXACT, r"$\theta_I$ first"),
                                ("held_share_h_first", MIXED, r"$\theta_H$ first")):
        row[2].plot(ts, [float(r[name]) for r in shown], color=colour, label=label)
    assert float(split[-1]["t"]) == 6300.0
    final = float(split[-1]["held_share_i_first"])
    assert round(final, 3) == 0.908, final  # KIM-42's headline at width 1
    row[2].set_ylim(0, 1)
    row[2].set_ylabel(r"Share of Var($\mu$)")
    row[2].set_title(r"Export: share from $\theta_H$")
    row[2].legend(loc="lower right")
    for x in row:
        x.set_xlim(0, 105)
        x.set_xlabel("Time (min)")
    h, l = row[0].get_legend_handles_labels()
    key[1].legend(h, l, loc="center left", title="90%, over", alignment="left", bbox_to_anchor=(-0.3, 0.5))

    # C: 62 more constants released to their priors, against B's 20: 90% over parameter pairs.
    row = ax[2]
    for k, (species, scale) in enumerate((("lac_export_meter", 1e6), ("M_gtp_c", 1e3))):
        twenty = band(b, species, "pairs")
        cond = band(c, species, "pairs", "per_theta_I")
        unc = band(c, species, "pairs", "unconditioned")
        row[k].fill_between(cond["t"] / 60, cond["lo90"] / scale, cond["hi90"] / scale, facecolor=PALE,
                            edgecolor=MIXED, lw=0.9, label="20 + 62 drawn, viable pairs")
        shade(row[k], twenty, scale, LAYER2, label="20 drawn, as in B")
        for y in (cond["lo90"], cond["hi90"]):  # the released band's edges, again above B's band
            row[k].plot(cond["t"] / 60, y / scale, color=MIXED, lw=0.9)
        for i, y in enumerate((unc["lo90"], unc["hi90"])):
            row[k].plot(unc["t"] / 60, y / scale, color=MIXED, lw=0.9, ls=(0, (3, 2)),
                        label="20 + 62 drawn, all pairs" if i == 0 else None)
        row[k].set_xlim(0, 105)
        row[k].set_xlabel("Time (min)")
    row[0].set_ylabel("Particles (millions)")
    row[0].set_title("Lactate exported")
    row[1].set_ylabel("Particles (thousands)")
    row[1].set_title("GTP")
    kd(row[1])
    # The key sits where B's third panel is, in the space this row leaves empty.
    ax[2, 2].legend(*row[0].get_legend_handles_labels(), loc="center left", alignment="left",
                    title=r"90% over pairs of $\theta_I$, $\theta_H$", bbox_to_anchor=(-0.12, 0.5))

    for k, title in enumerate(("M0 over 10 min, one dataset (scored)",
                               "The reduced model over one cycle, 20 parameters drawn (unscored)",
                               "62 more constants released to their priors (unscored)")):
        y = 1 - (0.02 + k * step) / H
        letter(fig, 0.004, y, "ABC"[k])
        fig.text(0.004 + 0.2 / WIDTH, y, title, fontsize=LETTER, color=INK, ha="left", va="top")
    save(fig, "fig_trajectories")


def band_truth(d, species):
    sel = sorted((r for r in d if r["species"] == species), key=lambda r: float(r["t"]))
    return np.array([float(r["t"]) for r in sel]), np.array([float(r["mean"]) for r in sel])


if __name__ == "__main__":
    archive = sys.argv[1] if len(sys.argv) > 1 else ARCHIVE
    b7 = sys.argv[2] if len(sys.argv) > 2 else B7_RUN
    fig_calibration(b7)
    fig_mixed(b7)
    fig_budget(archive)
    fig_predict_m0(archive)
    fig_crossed(archive)
    fig_trajectories()
