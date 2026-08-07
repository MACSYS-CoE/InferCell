#!/usr/bin/env python3
"""ISAB results figure: information crossing a module boundary.

Overlays two ABC-SMC posteriors over the bursty regulator's six parameters:
unconditioned (its own priors) and conditioned (KDE priors fitted to the NUTS
posterior of the differentiable block).

The parameter order matters and is not cosmetic. The three parameters the
simulation block does NOT share with the differentiable block (k_on, k_off,
k_tx_burst -- promoter switching, which nothing in the bulk data sees) come
first; the three shared ones (k_tl, gamma_mRNA, gamma_prot -- the same
ribosomes, the same degradation machinery) come last. So the tightening is
confined to the bottom-right block of the corner, and the boxed region is
literally where the coupling is.

--mode alone draws the SAME figure with the conditioned posterior suppressed,
for the setup slide that precedes the result. The two are a build: identical
axes, identical shared-parameter box, identical legend geometry, so advancing
the slide makes the green appear and moves nothing else. That is why the axis
ranges below are always computed from BOTH posteriors even when only one is
drawn -- ranging the alone figure on its own particles would shift every panel
between the two slides and destroy the effect.

Usage:
    python3 dev/talks/ISAB/figures/plot_bursty_boundary.py \
        --results examples/results \
        --out dev/talks/ISAB/figures/isab_ssa_posteriors

    python3 dev/talks/ISAB/figures/plot_bursty_boundary.py \
        --results examples/results --mode alone
"""

import argparse
import csv
from pathlib import Path

import corner
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np
from matplotlib.lines import Line2D
from matplotlib.patches import Rectangle

plt.rcParams.update({
    "font.family": "serif",
    "mathtext.fontset": "cm",
})


def set_fonts(ndim):
    """Size type so it is still legible after the slide downscales the figure.

    corner draws roughly 2 inches per panel, so a 6-parameter grid is ~13 in
    wide and gets reduced ~5.5x to reach the 6.4 cm it occupies on the results
    slide. At PR #36's 26 pt that lands near 4 pt on the projector. Sizes scale
    with ndim because the 4-parameter backup corner is reduced far less.

    Returns the label size, which corner also needs passed via label_kwargs.
    """
    label = 26 + 6.5 * max(0, ndim - 4)
    plt.rcParams.update({
        "axes.labelsize": label,
        "xtick.labelsize": 0.66 * label,
        "ytick.labelsize": 0.66 * label,
    })
    return label

# One colour code across the whole deck, matching the slide-10 schematic:
# blue is the differentiable block, orange is the stochastic one. Before Aug
# 2026 both corners were drawn in the same blue, which was survivable while they
# lived on separate slides and stopped being so once "One Way" put them side by
# side -- two adjacent same-coloured corners read as the same quantity.
#
# Truth moved off orange to make room. Black is the corner-plot convention
# anyway, and it stops the truth lines competing with a posterior for the eye.
COL_DIFF = "#3B7EA1"     # blue   -- the differentiable block's own posterior
COL_ALONE = "#D2762B"    # orange -- the regulator inferred alone
COL_COND = "#3D9970"     # green  -- conditioned across the boundary
COL_TRUTH = "#1A1A1A"    # black  -- ground truth

LABELS = {
    "k_on": r"$k_{\mathrm{on}}$",
    "k_off": r"$k_{\mathrm{off}}$",
    "k_tx_burst": r"$k_{\mathrm{tx}}^{\mathrm{burst}}$",
    "k_tl": r"$k_{tl}$",
    "gamma_mRNA": r"$\gamma_{mRNA}$",
    "gamma_protein": r"$\gamma_{prot}$",
    "k_tx": r"$k_{tx}$",
    "k_atp": r"$k_{ATP}$",
    "k_ntp": r"$k_{NTP}$",
    "k_aa": r"$k_{AA}$",
    "sigma_obs": r"$\sigma_{obs}$",
}

# Which of the differentiable block's parameters the backup corner shows. The
# full 8x8 grid puts the axis labels below reading size at slide height, and the
# rates are the part that carries the identifiability claim -- sigma_obs and the
# initial conditions are nuisance. Set to None to plot all of them.
ODE_SUBSET = ["k_tx", "k_tl", "gamma_mRNA", "gamma_protein"]


def read_particles(path):
    with open(path) as fh:
        reader = csv.reader(fh)
        names = next(reader)
        rows = [[float(x) for x in row] for row in reader if row]
    return names, np.array(rows)


def read_weights(path):
    with open(path) as fh:
        next(fh)
        return np.array([float(line) for line in fh if line.strip()])


def weighted_quantile(values, weights, quantiles):
    """Quantiles of a weighted sample, by interpolating the weighted CDF."""
    order = np.argsort(values)
    v, w = np.asarray(values)[order], np.asarray(weights)[order]
    total = w.sum()
    if total <= 0:
        return np.quantile(values, quantiles)
    # Midpoint convention: cumulative weight at each point less half its own,
    # so a single atom lands at the median rather than at 1.0.
    cdf = (np.cumsum(w) - 0.5 * w) / total
    return np.interp(quantiles, cdf, v)


def read_truth(path):
    truth, shared = {}, {}
    with open(path) as fh:
        for row in csv.DictReader(fh):
            truth[row["param"]] = float(row["truth"])
            shared[row["param"]] = row["conditioned"].strip().lower() == "true"
    return truth, shared


def plot_ode_corner(res, out):
    """Backup slide: the differentiable block's own posterior.

    One posterior, not two -- this is the identifiability check promised on
    "What Becomes Possible (Science)", not a boundary result. Tilted ellipses
    are the claim that the data constrains combinations of rates.
    """
    names, chain = read_particles(res / "isab_ode_chain.csv")
    truth = {}
    with open(res / "isab_ode_truth.csv") as fh:
        for row in csv.DictReader(fh):
            truth[row["param"]] = float(row["truth"])

    keep = [j for j, n in enumerate(names)
            if ODE_SUBSET is None or n in ODE_SUBSET]
    missing = set(ODE_SUBSET or []) - set(names)
    if missing:
        raise SystemExit(f"ODE_SUBSET names not in chain: {sorted(missing)}")

    sub = chain[:, keep]
    labels = [LABELS.get(names[j], names[j]) for j in keep]
    truths = [truth.get(names[j]) for j in keep]
    fs = set_fonts(len(keep))

    fig = corner.corner(
        sub, labels=labels, color=COL_DIFF,
        truths=truths, truth_color=COL_TRUTH,
        plot_datapoints=False, fill_contours=True, smooth=1.0,
        levels=(0.68, 0.95), hist_kwargs={"linewidth": 2.4},
        contour_kwargs={"linewidths": 1.6},
        label_kwargs={"fontsize": fs}, max_n_ticks=3,
    )
    for ext in ("png", "pdf"):
        fig.savefig(f"{out}.{ext}", dpi=200, bbox_inches="tight")
        print(f"wrote {out}.{ext}")


def weighted_kde(samples, weights, grid, bw_scale=1.0):
    """Gaussian KDE of a weighted sample. Hand-rolled to avoid a scipy dep."""
    s, w = np.asarray(samples), np.asarray(weights, dtype=float)
    w = w / w.sum()
    mean = (w * s).sum()
    sd = np.sqrt((w * (s - mean) ** 2).sum())
    neff = 1.0 / (w ** 2).sum()                       # Kish effective sample size
    bw = bw_scale * 1.06 * sd * neff ** (-0.2)        # Silverman, weighted
    if not bw > 0:
        bw = 1e-6
    z = (grid[:, None] - s[None, :]) / bw
    return (w[None, :] * np.exp(-0.5 * z ** 2)).sum(axis=1) / (bw * np.sqrt(2 * np.pi))


def plot_burst_rate(res, out):
    """Backup slide: why the promoter rates miss truth.

    A telegraph promoter enters a replicate-MEAN trajectory only through the
    effective burst rate k_tx_burst * k_on/(k_on + k_off). The individual rates
    are therefore unidentifiable from these summary statistics, but the
    combination is -- and the boundary sharpens it too. This plot is that claim
    shown rather than asserted, so it can answer the question the results slide
    invites: "why do three of your parameters miss the orange line?"
    """
    names, cond = read_particles(res / "isab_cond_particles.csv")
    names_u, uncond = read_particles(res / "isab_uncond_particles.csv")
    if names != names_u:
        raise SystemExit(f"parameter order differs: {names} vs {names_u}")
    w_cond = read_weights(res / "isab_cond_weights.csv")
    w_uncond = read_weights(res / "isab_uncond_weights.csv")
    truth, _ = read_truth(res / "isab_ssa_truth.csv")

    i = {k: names.index(k) for k in ("k_on", "k_off", "k_tx_burst")}

    def burst(P):
        return P[:, i["k_tx_burst"]] * P[:, i["k_on"]] / (P[:, i["k_on"]] + P[:, i["k_off"]])

    b_cond, b_uncond = burst(cond), burst(uncond)
    b_truth = (truth["k_tx_burst"] * truth["k_on"]
               / (truth["k_on"] + truth["k_off"]))

    lo = 0.0
    hi = weighted_quantile(np.concatenate([b_uncond, b_cond]),
                           np.concatenate([w_uncond, w_cond]), [0.99])[0]
    hi = max(hi, b_truth * 1.6)
    grid = np.linspace(lo, hi, 400)

    set_fonts(4)
    fig, ax = plt.subplots(figsize=(7.4, 4.6))
    for vals, wts, col, lab in (
        (b_uncond, w_uncond, COL_ALONE, "Regulator alone"),
        (b_cond, w_cond, COL_COND, "Conditioned"),
    ):
        d = weighted_kde(vals, wts, grid)
        ax.plot(grid, d, color=col, lw=3.0, label=lab)
        ax.fill_between(grid, d, color=col, alpha=0.22)
    ax.axvline(b_truth, color=COL_TRUTH, lw=2.5, label="Truth")

    ax.set_xlabel(r"burst rate  $k_{\mathrm{tx}}^{\mathrm{burst}}\,"
                  r"k_{\mathrm{on}}/(k_{\mathrm{on}}+k_{\mathrm{off}})$")
    ax.set_yticks([])
    ax.set_xlim(lo, hi)
    ax.spines[["top", "right", "left"]].set_visible(False)
    ax.legend(frameon=False, fontsize=20)
    fig.tight_layout()
    for ext in ("png", "pdf"):
        fig.savefig(f"{out}.{ext}", dpi=200, bbox_inches="tight")
        print(f"wrote {out}.{ext}")


def plot_ssa(res, out, overlay=True):
    """The boundary figure. overlay=False draws the unconditioned posterior only.

    Everything except the green contours is identical between the two modes --
    see the module docstring for why that matters.
    """
    names, cond = read_particles(res / "isab_cond_particles.csv")
    names_u, uncond = read_particles(res / "isab_uncond_particles.csv")
    if names != names_u:
        raise SystemExit(f"parameter order differs: {names} vs {names_u}")

    w_cond = read_weights(res / "isab_cond_weights.csv")
    w_uncond = read_weights(res / "isab_uncond_weights.csv")
    truth, shared = read_truth(res / "isab_ssa_truth.csv")

    labels = [LABELS.get(n, n) for n in names]
    truths = [truth[n] for n in names]
    ndim = len(names)
    fs = set_fonts(ndim)

    # Shared ranges across both posteriors, so the two are directly comparable
    # panel by panel. Without this, corner picks per-call limits and the
    # conditioned overlay silently rescales the axes it is drawn onto.
    #
    # WEIGHTED quantiles, not plain ones: these are ABC-SMC particles, and a
    # handful of near-zero-weight survivors sit far out in the tails. Ranging on
    # unweighted quantiles let those set the limits and squashed every posterior
    # into the bottom-left corner of its panel.
    ranges = []
    for j in range(ndim):
        vals = np.concatenate([uncond[:, j], cond[:, j]])
        wts = np.concatenate([w_uncond, w_cond])
        lo, hi = weighted_quantile(vals, wts, [0.005, 0.995])
        if not hi > lo:                      # degenerate: a near-delta posterior
            lo, hi = lo - 0.05 * abs(lo) - 1e-9, hi + 0.05 * abs(hi) + 1e-9
        pad = 0.08 * (hi - lo)
        ranges.append((min(lo - pad, truths[j] - pad), max(hi + pad, truths[j] + pad)))

    fig = corner.corner(
        uncond, weights=w_uncond, labels=labels, range=ranges,
        color=COL_ALONE, truths=truths, truth_color=COL_TRUTH,
        plot_datapoints=False, fill_contours=True, smooth=1.0,
        levels=(0.68, 0.95), hist_kwargs={"linewidth": 2.4},
        contour_kwargs={"linewidths": 1.6},
        label_kwargs={"fontsize": fs}, max_n_ticks=3,
    )
    if overlay:
        corner.corner(
            cond, weights=w_cond, labels=labels, range=ranges, fig=fig,
            color=COL_COND, truths=truths, truth_color=COL_TRUTH,
            plot_datapoints=False, fill_contours=True, smooth=1.0,
            levels=(0.68, 0.95), hist_kwargs={"linewidth": 2.4},
            contour_kwargs={"linewidths": 1.6},
            label_kwargs={"fontsize": fs}, max_n_ticks=3,
        )

    axes = np.array(fig.axes).reshape((ndim, ndim))

    # Box the shared-parameter block. These are contiguous by construction (the
    # Julia export orders unshared first), so one rectangle covers them; bail
    # out rather than draw a misleading box if that ever stops being true.
    #
    # Drawn in BOTH modes. It marks which parameters the two modules share --
    # a fact about the model, not about the inference -- so on the alone figure
    # it says "we already know these three from elsewhere and are refitting
    # them anyway", and on the overlay the green then lands inside it.
    shared_idx = [j for j, n in enumerate(names) if shared[n]]
    if shared_idx and shared_idx == list(range(shared_idx[0], shared_idx[-1] + 1)):
        lo, hi = shared_idx[0], shared_idx[-1]
        p0 = axes[hi, lo].get_position()
        p1 = axes[lo, hi].get_position()
        fig.add_artist(Rectangle(
            (p0.x0 - 0.012, p0.y0 - 0.012),
            p1.x1 - p0.x0 + 0.024, p1.y1 - p0.y0 + 0.024,
            transform=fig.transFigure, fill=False,
            edgecolor=COL_COND, linewidth=2.5, linestyle=(0, (6, 4)), zorder=10,
        ))
        fig.text(
            p1.x1, p1.y1 + 0.030,
            "shared with the\ndifferentiable block",
            ha="right", va="bottom", fontsize=0.62 * fs, color=COL_COND,
            linespacing=1.25,
        )
    else:
        print(f"WARNING: shared parameters are not contiguous ({shared_idx}); "
              "skipping the annotation box")

    # The conditioned row is RESERVED rather than dropped when overlay=False: an
    # invisible handle of the same height holds the slot, so "Truth" sits at the
    # same y on both slides and advancing the build only adds the green line
    # instead of shuffling the legend under the audience.
    cond_handle = (
        Line2D([], [], color=COL_COND, lw=6, label="Conditioned across the boundary")
        if overlay else
        Line2D([], [], color="none", lw=6, label=" ")
    )
    handles = [
        Line2D([], [], color=COL_ALONE, lw=6, label="Regulator inferred alone"),
        cond_handle,
        Line2D([], [], color=COL_TRUTH, lw=3, label="Truth"),
    ]
    fig.legend(handles=handles, loc="upper right",
               bbox_to_anchor=(0.995, 0.995), frameon=False,
               fontsize=0.72 * fs, handlelength=1.4, labelspacing=0.35,
               borderaxespad=0.0)

    out = Path(out)
    out.parent.mkdir(parents=True, exist_ok=True)
    for ext in ("png", "pdf"):
        fig.savefig(f"{out}.{ext}", dpi=200, bbox_inches="tight")
        print(f"wrote {out}.{ext}")


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--results", default="examples/results")
    ap.add_argument("--out", default=None,
                    help="output stem; defaults per --mode")
    ap.add_argument("--mode", choices=("ssa", "alone", "ode", "burst"),
                    default="ssa",
                    help="ssa: the boundary figure. alone: the same figure with "
                         "the conditioned posterior suppressed, for the setup "
                         "slide. ode / burst: backup slides.")
    args = ap.parse_args()

    res = Path(args.results)
    figdir = Path("dev/talks/ISAB/figures")
    default = {
        "ssa": "isab_ssa_posteriors",
        "alone": "isab_ssa_alone",
        "ode": "isab_ode_corner",
        "burst": "isab_burst_rate",
    }[args.mode]
    out = Path(args.out) if args.out else figdir / default
    out.parent.mkdir(parents=True, exist_ok=True)

    if args.mode == "ode":
        plot_ode_corner(res, out)
    elif args.mode == "burst":
        plot_burst_rate(res, out)
    else:
        plot_ssa(res, out, overlay=(args.mode == "ssa"))


if __name__ == "__main__":
    main()
