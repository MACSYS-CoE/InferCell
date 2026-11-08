#!/usr/bin/env python3
"""The paper's data figures, in one style, written to this directory as PDFs.

Run from the repository root:
    python3 paper/figures/make_figures.py [EVIDENCE_ARCHIVE] [B7_RUN_DIR]

Inputs, read and never written:
- the KIM-43 evidence archive: runs/m0-kim40-kim41/score_datasets.tsv and
  runs/budget-kim42/budget_rows.tsv;
- paper/figures/data: av6_replicates.tsv, av6_bands.tsv, av6_v6_draws.tsv, copied from B7's run,
  budget_mu.tsv, from the archive by paper/figures/budget_dump.jl, and traj_*.tsv, copied unchanged
  from the trajectory runs' tsv/ (infercell-trajectories-2026-11-08, dev/scripts/trajectories_predict_result.md);
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
WIDTH = 6.5  # arxiv.sty text width, inches: these figures span both columns

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
    a = top.subplots(1, 6)
    density_row(a, ref, mixed, "Mixed", MIXED)
    legend_row(top, a)
    panel(top, "A")
    b = bottom.subplots(1, 6, sharey=True)
    ecdf_row(b, z, lo, hi, [(exact_ranks(), EXACT, "Exact"), (mixed_ranks(b7), MIXED, "Mixed")])
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
    # Each point averages many pairs, so it shows one parameter block's reach, with the other
    # block averaged out. Pairs the viability rule rejected are left out.
    mu = rows(os.path.join(HERE, "data", "budget_mu.tsv"))
    rng = np.random.default_rng(1)
    for k, w in enumerate(widths):
        kept = [r for r in mu if float(r["width"]) == w and r["viable"] == "true"]
        for key, off, marker, colour, label in (("row", -0.16, "o", EXACT, r"one draw of $\theta_I$"),
                                                ("column", 0.16, "s", MIXED, r"one draw of $\theta_H$")):
            groups = {}
            for r in kept:
                groups.setdefault(r[key], []).append(float(r["mu"]))
            y = np.array([np.mean(v) for v in groups.values()]) / 1e6
            b.scatter(k + off + rng.uniform(-0.09, 0.09, len(y)), y, s=7, marker=marker, color=colour,
                      alpha=0.55, lw=0, label=label if k == 0 else None)
    b.set_xticks(x, tick)
    b.set_xlim(-0.5, 2.5)
    b.set_ylim(bottom=0)
    b.set_xlabel(r"Prior width of $\theta_H$")
    b.set_ylabel("Lactate exported, millions")
    b.legend(loc="upper left", handletextpad=0.2)
    b.set_title("B", loc="left", fontweight="bold", fontsize=9)
    b.set_title("UNSCORED", loc="right", color=MUTED, fontweight="bold")
    save(fig, "fig_budget")


def fig_predict_m0(archive):
    d = rows(os.path.join(archive, "runs/m0-kim40-kim41/score_datasets.tsv"))
    fig, axes = plt.subplots(1, 4, figsize=(WIDTH, 2.05), layout="constrained")
    panels = (("q", "v7", 1e6, "Export, exact", EXACT, "o", "millions"),
              ("q", "b7", 1e6, "Export, mixed", MIXED, "s", "millions"),
              ("dq", "v7", 1e3, "Knockdown, exact", EXACT, "o", "thousands"),
              ("dq", "b7", 1e3, "Knockdown, mixed", MIXED, "s", "thousands"))
    for ax, letter, (target, post, scale, title, colour, marker, unit) in zip(axes, "ABCD", panels):
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
        ax.plot(lim, lim, color=MUTED, lw=0.8, ls="--")
        ax.errorbar(t, m, yerr=[m - lo, hi - m], fmt=marker, ms=2.5, lw=0.8, color=colour)
        ax.set_xlim(lim)
        ax.set_ylim(lim)
        ax.set_aspect("equal")
        ax.set_title(title, loc="left", x=0.02)
        ax.set_xlabel(f"Truth, {unit}")
        ax.set_ylabel(f"Predicted, {unit}")
        ax.text(-0.04, 1.045, letter, transform=ax.transAxes, fontweight="bold", fontsize=9, ha="right")
    save(fig, "fig_predict_m0")


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

    fig = plt.figure(figsize=(WIDTH, 6.4), layout="constrained")
    sub = fig.subfigures(3, 1)
    # The same column ratios in every row, with a fourth column for the row's key, keep the columns aligned.
    grid = [s.subplots(1, 4, width_ratios=[1, 1, 1, 0.7]) for s in sub]
    ax = np.array([g[:3] for g in grid])
    key = [g[3] for g in grid]
    for k in key:
        k.axis("off")
    ax[2, 2].axis("off")

    # A: M0 over 600 s, the exact posterior on V7 replicate 1, against its truth.
    row = ax[0]
    nb = band(a, "M_lac__L_e", "new_bulk")
    nb = {k: (v if k == "t" else v * RATIO) for k, v in nb.items()}
    shade(row[0], nb, 1e3, LIGHT, label="New bulk\nmeasurement")
    shade(row[0], band(a, "lac_export_meter", "rows"), 1e3, EXACT, label="Posterior")
    t, v = obs("M_lac__L_e")
    row[0].plot(t / 60, v * RATIO / 1e3, "o", ms=3, mfc="white", mec=INK, mew=0.7, label="Data")
    row[0].plot(tt / 60, meter / 1e3, color=INK, lw=0.8, ls=(0, (3, 1.5)), label="Truth")
    row[0].set_ylabel("Thousands")
    row[0].set_title("Lactate exported, observed")
    shade(row[1], band(a, "M_gtp_c", "new_bulk"), 1e3, LIGHT)
    shade(row[1], band(a, "M_gtp_c", "rows"), 1e3, EXACT)
    t, v = obs("M_gtp_c")
    row[1].plot(t / 60, v / 1e3, "o", ms=3, mfc="white", mec=INK, mew=0.7)
    tg, g = tru("M_gtp_c")
    row[1].plot(tg / 60, g / 1e3, color=INK, lw=0.8, ls=(0, (3, 1.5)))
    row[1].set_ylabel("Thousands")
    row[1].set_title("GTP, observed")
    # Pyruvate rises by a third in the first minute and then holds, so its band shows only against the truth.
    p = band(a, "M_pyr_c", "rows")
    tp, ptrue = tru("M_pyr_c")
    assert np.array_equal(tp, p["t"])
    later = tp > 0
    off = {k: (v[later] if k == "t" else v[later] - ptrue[later]) for k, v in p.items()}
    shade(row[2], off, 1, EXACT)
    row[2].axhline(0, color=INK, lw=0.8, ls=(0, (3, 1.5)))
    row[2].set_xlim(0, 10)
    row[2].set_ylabel("Deviation from truth, particles")
    row[2].set_title("Pyruvate, not observed")
    for x in row:
        x.set_xlim(0, 10)
        x.set_xlabel("Time, min")
    key[0].legend(*row[0].get_legend_handles_labels(), loc="center left", title="90% bands")

    # B: the reduced model over one cycle, 20 parameters drawn, in three nested layers.
    row = ax[1]
    layers = (("cells", LAYER3, "A new cell"), ("pairs", LAYER2, r"$\theta_I$ and $\theta_H$"),
              ("rows", EXACT, r"$\theta_I$ only"))
    for k, (species, scale) in enumerate((("lac_export_meter", 1e6), ("M_gtp_c", 1e3))):
        for layer, colour, label in layers:
            shade(row[k], band(b, species, layer), scale, colour, label=label)
        m = band(b, species, "rows")
        row[k].plot(m["t"] / 60, m["mean"] / scale, color=INK, lw=1, label="Predictive mean")
    end = band(b, "lac_export_meter", "rows")["mean"][-1]
    assert abs(end / 9.5765e6 - 1) < 1e-4, end  # KIM-42's mean Q at width 1
    row[0].set_ylabel("Millions")
    row[0].set_title("Lactate exported")
    row[1].set_ylabel("Thousands")
    row[1].set_title("GTP")
    kd(row[1])
    shown = [r for r in split if r["flagged"] == "false" and r["mc_undefined"] == "false"
             and float(r["total"]) > float(r["mc"])]
    ts = np.array([float(r["t"]) for r in shown]) / 60
    for name, colour, label in (("held_share_i_first", EXACT, r"$\theta_I$ first"),
                                ("held_share_h_first", MIXED, r"$\theta_H$ first")):
        row[2].plot(ts, [float(r[name]) for r in shown], color=colour, label=label)
    final = float(split[-1]["held_share_i_first"])
    assert round(final, 3) == 0.908, final  # KIM-42's headline at width 1
    row[2].set_ylim(0, 1)
    row[2].set_ylabel(r"Share of Var($\mu$)")
    row[2].set_title(r"Export: share from $\theta_H$")
    row[2].legend(loc="lower right")
    for x in row:
        x.set_xlim(0, 105)
        x.set_xlabel("Time, min")
    h, l = row[0].get_legend_handles_labels()
    key[1].legend(h, l, loc="center left", title="90%, over")

    # C: 62 more constants released to their priors, against B's 20: 90% over parameter pairs.
    row = ax[2]
    for k, (species, scale) in enumerate((("lac_export_meter", 1e6), ("M_gtp_c", 1e3))):
        twenty = band(b, species, "pairs")
        cond = band(c, species, "pairs", "per_theta_I")
        unc = band(c, species, "pairs", "unconditioned")
        row[k].fill_between(cond["t"] / 60, cond["lo90"] / scale, cond["hi90"] / scale, facecolor=PALE,
                            edgecolor=MIXED, lw=1, label="20 + 62 drawn, viable pairs")
        shade(row[k], twenty, scale, LAYER2, label="20 drawn, as in B")
        for y in (cond["lo90"], cond["hi90"]):  # the released band's edges, again above B's band
            row[k].plot(cond["t"] / 60, y / scale, color=MIXED, lw=1)
        for i, y in enumerate((unc["lo90"], unc["hi90"])):
            row[k].plot(unc["t"] / 60, y / scale, color=MIXED, lw=1, ls=(0, (3, 2)),
                        label="20 + 62 drawn, all pairs" if i == 0 else None)
        row[k].set_xlim(0, 105)
        row[k].set_xlabel("Time, min")
    row[0].set_ylabel("Millions")
    row[0].set_title("Lactate exported")
    row[1].set_ylabel("Thousands")
    row[1].set_title("GTP")
    kd(row[1])
    row[2].legend(*row[0].get_legend_handles_labels(), loc="center left",
                  title=r"90% over pairs of $\theta_I$, $\theta_H$")

    for s, title in zip(sub, ("A   M0 over 10 min, one dataset: scored",
                              "B   The reduced model over one cycle, 20 parameters drawn: unscored",
                              "C   62 more constants released to their priors: unscored")):
        s.suptitle(title, x=0.005, ha="left", fontweight="bold", fontsize=9, color=INK)
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
