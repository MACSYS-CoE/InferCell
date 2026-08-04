"""Step 4 (Level-2 boundary protocol) publication-style figures.

Reads CSV artefacts written by `test/run_integ_iterative.jl` (export section)
and produces three figures:

1. Corner plot — single-pass vs iterative ODE posteriors with ground truth.
2. Per-parameter 90% CI width per iteration vs single-pass.
3. KL trace across iterations (log-scale).
"""

import numpy as np
import matplotlib.pyplot as plt
import corner
import os
import glob

plt.rcParams.update({
    "font.family": "serif",
    "font.serif": ["Computer Modern Roman"],
    "text.usetex": True,            # match step3 plots; LaTeX is available in talk env
    "axes.spines.top": False,
    "axes.spines.right": False,
    "axes.grid": False,
    "figure.dpi": 150,
})

RESULTS = "examples/results"
FIGS = "examples/figures"
os.makedirs(FIGS, exist_ok=True)

# Wong-palette colours (colourblind-safe)
SINGLE_C = "#0072B2"   # blue
ITER_C   = "#D55E00"   # vermillion
TRUTH_C  = "#E69F00"   # orange
ITER_FADE = ["#D55E00", "#E69F00", "#F0E442", "#009E73"]


def load_chain(path):
    """Load CSV chain (n_samples x n_params) with header row of param names."""
    with open(path) as f:
        header = f.readline().strip().split(",")
    data = np.genfromtxt(path, delimiter=",", skip_header=1)
    return header, data


# --- Load all artefacts ---
single_names, single_chain = load_chain(f"{RESULTS}/step4_single_chain.csv")
iter_names,   iter_chain   = load_chain(f"{RESULTS}/step4_iter_chain.csv")

# Per-iter histories: step4_iter_history_1.csv, _2.csv, ...
history_paths = sorted(
    glob.glob(f"{RESULTS}/step4_iter_history_*.csv"),
    key=lambda p: int(os.path.basename(p).split("_")[-1].split(".")[0]),
)
history_chains = [load_chain(p)[1] for p in history_paths]

# Ground truth
truth_arr = np.genfromtxt(f"{RESULTS}/step4_truth.csv", delimiter=",",
                           skip_header=1, dtype=None, encoding="utf-8")
truth = {row[0]: float(row[1]) for row in truth_arr}
truth_vec = [truth[n] for n in single_names]

# KL trace
kl_data = np.genfromtxt(f"{RESULTS}/step4_kl_trace.csv", delimiter=",",
                         skip_header=1, dtype=None, encoding="utf-8")
# kl_data rows are (iter_pair, kl); coerce to list of floats
if kl_data.ndim == 0:
    kl_vals = [float(kl_data["f1"])]
else:
    kl_vals = [float(r[1]) for r in kl_data]

# Pretty parameter labels for corner plot
LABEL_MAP = {
    "k_tx": r"$k_{tx}$", "k_tl": r"$k_{tl}$",
    "gamma_mRNA": r"$\gamma_{mRNA}$", "gamma_protein": r"$\gamma_{prot}$",
    "k_atp": r"$k_{atp}$", "k_ntp": r"$k_{ntp}$", "k_aa": r"$k_{aa}$",
    "sigma_obs": r"$\sigma$",
}
labels = [LABEL_MAP.get(n, n) for n in single_names]

# Shared parameters between ODE Block 2 and SSA Block 3 (where the protocol
# does work). These are the rows highlighted in the bar chart.
SHARED = ["k_tl", "gamma_mRNA", "gamma_protein"]

# ============================================================
# Figure 1: corner plot — single-pass vs iterative + ground truth
# ============================================================
print("Plotting corner ...")
levels = (1 - np.exp(-0.5), 1 - np.exp(-2))  # 1,2 sigma

fig1 = corner.corner(
    single_chain,
    labels=labels,
    truths=truth_vec,
    truth_color=TRUTH_C,
    color=SINGLE_C,
    levels=levels,
    fill_contours=True,
    smooth=1.2, smooth1d=1.0,
    quantiles=[0.05, 0.5, 0.95],
    show_titles=False,
    plot_datapoints=False,
    plot_density=False,
)
corner.corner(
    iter_chain,
    labels=labels,
    color=ITER_C,
    levels=levels,
    fill_contours=False,
    smooth=1.2, smooth1d=1.0,
    quantiles=[0.05, 0.5, 0.95],
    show_titles=False,
    plot_datapoints=False,
    plot_density=False,
    fig=fig1,
)

# Manual legend in the upper-right
ax_leg = fig1.axes[len(labels) - 1]  # top-right cell
handles = [
    plt.Line2D([], [], color=SINGLE_C, lw=2, label="single-pass"),
    plt.Line2D([], [], color=ITER_C, lw=2, label="iterative"),
    plt.Line2D([], [], color=TRUTH_C, lw=2, label="truth"),
]
ax_leg.legend(handles=handles, loc="upper right", frameon=False, fontsize=10)

fig1.savefig(f"{FIGS}/step4_corner.png", bbox_inches="tight")
plt.close(fig1)

# ============================================================
# Figure 2: CI width per iteration vs single-pass (shared params only)
# ============================================================
print("Plotting CI widths ...")
def ci_width(arr, q_lo=0.05, q_hi=0.95):
    lo, hi = np.quantile(arr, [q_lo, q_hi])
    return hi - lo

fig2, ax = plt.subplots(figsize=(7, 4))
n_iter = len(history_chains)
x = np.arange(len(SHARED))
bar_w = 0.18

# Bar group: single-pass + each iter
single_widths = [ci_width(single_chain[:, single_names.index(p)]) for p in SHARED]
ax.bar(x - bar_w * (n_iter / 2), single_widths, width=bar_w,
       color=SINGLE_C, label="single-pass")

for k, ch in enumerate(history_chains):
    widths = [ci_width(ch[:, single_names.index(p)]) for p in SHARED]
    ax.bar(x + bar_w * (k - n_iter / 2 + 1), widths, width=bar_w,
           color=ITER_FADE[k % len(ITER_FADE)], label=f"iter {k+1}")

ax.set_xticks(x)
ax.set_xticklabels([LABEL_MAP[p] for p in SHARED])
ax.set_ylabel("90% CI width")
ax.set_title("Posterior tightening across iterations")
ax.legend(frameon=False, fontsize=9, loc="upper right")
fig2.tight_layout()
fig2.savefig(f"{FIGS}/step4_ci_widths.png", bbox_inches="tight")
plt.close(fig2)

# ============================================================
# Figure 3: KL trace
# ============================================================
print("Plotting KL trace ...")
fig3, ax = plt.subplots(figsize=(6, 4))
iters = np.arange(2, 2 + len(kl_vals))
ax.semilogy(iters, kl_vals, "o-", color=ITER_C, lw=2, ms=8,
            label="KL(iter_k || iter_{k-1})")
ax.axhline(0.05, color="grey", ls="--", lw=1, label="kl_tol = 0.05")
ax.set_xlabel("iteration k")
ax.set_ylabel("KL divergence (summed over shared params)")
ax.set_title("Level-2 boundary protocol convergence")
ax.set_xticks(iters)
ax.legend(frameon=False)
fig3.tight_layout()
fig3.savefig(f"{FIGS}/step4_kl_trace.png", bbox_inches="tight")
plt.close(fig3)

# ============================================================
# Talk-style corners (Figures 4 + 5)
# ------------------------------------------------------------
# Aesthetic conventions for both:
#   - Translucent filled 1D marginals with crisp outlines
#   - Line-only 2D contours at 1, 2 sigma with subtle drop shadow
#   - Bold black truth crosshairs drawn after corner (not via corner's
#     thin truths= path) so we can control linewidth/alpha
#   - Large fonts, rotated tick labels, fewer ticks per axis
# ============================================================
from matplotlib.patheffects import Stroke, Normal
from matplotlib.patches import Patch
from scipy.stats import gaussian_kde

TRUTH_TALK = "#1a1a1a"
LEVELS_TALK = (1 - np.exp(-0.5), 1 - np.exp(-2))  # 1, 2 sigma
HIST_ALPHA = 0.30


def _draw_kde_marginal(ax, samples, color, alpha=HIST_ALPHA, lw=3.0):
    """Replace the diagonal histogram with a smooth filled KDE."""
    kde = gaussian_kde(samples, bw_method="scott")
    lo, hi = np.quantile(samples, [0.001, 0.999])
    pad = 0.08 * (hi - lo)
    grid = np.linspace(lo - pad, hi + pad, 400)
    density = kde(grid)
    ax.fill_between(grid, density, color=color, alpha=alpha,
                    linewidth=0, zorder=2)
    ax.plot(grid, density, color=color, lw=lw, alpha=0.95, zorder=3)


def _talk_corner(chain_single, chain_iter, names, label_list, truth_list,
                 figsize, label_size, tick_size, legend_size,
                 line_w=3.0, savepath=None):
    """Render a talk-grade overlaid corner plot.

    Both posteriors share style: filled translucent KDE marginals, crisp
    contour lines at 1 + 2 sigma with a white halo for readability, and
    bold dark truth crosshairs on top.
    """
    fig = plt.figure(figsize=figsize)

    # Corner needs hist_kwargs that fit a Line2D when smooth1d is set;
    # we overwrite the diagonals manually below, so just keep these neutral.
    line_only_kw_s = {"color": SINGLE_C, "lw": 0.0, "alpha": 0.0}
    line_only_kw_i = {"color": ITER_C, "lw": 0.0, "alpha": 0.0}

    contour_kw = {"linewidths": line_w, "alpha": 0.95}
    halo = [Stroke(linewidth=line_w + 2.4, foreground="white",
                   alpha=0.85), Normal()]

    corner.corner(
        chain_single,
        labels=label_list,
        color=SINGLE_C,
        levels=LEVELS_TALK,
        fill_contours=False,
        no_fill_contours=True,
        plot_contours=True,
        smooth=1.4, smooth1d=1.2,
        quantiles=None,
        show_titles=False,
        plot_datapoints=False,
        plot_density=False,
        label_kwargs={"fontsize": label_size},
        hist_kwargs=line_only_kw_s,
        contour_kwargs=dict(contour_kw, path_effects=halo),
        fig=fig,
    )
    corner.corner(
        chain_iter,
        labels=label_list,
        color=ITER_C,
        levels=LEVELS_TALK,
        fill_contours=False,
        no_fill_contours=True,
        plot_contours=True,
        smooth=1.4, smooth1d=1.2,
        quantiles=None,
        show_titles=False,
        plot_datapoints=False,
        plot_density=False,
        hist_kwargs=line_only_kw_i,
        contour_kwargs=dict(contour_kw, path_effects=halo),
        fig=fig,
    )

    ndim = len(label_list)
    axes = np.array(fig.axes).reshape(ndim, ndim)

    # Manual KDE marginals on the diagonal (overwrites corner's invisible
    # placeholder histograms).
    for k in range(ndim):
        ax = axes[k, k]
        xlim = ax.get_xlim()
        ax.clear()
        _draw_kde_marginal(ax, chain_single[:, k], SINGLE_C, lw=line_w)
        _draw_kde_marginal(ax, chain_iter[:, k],   ITER_C,   lw=line_w)
        ax.axvline(truth_list[k], color=TRUTH_TALK, lw=2.4,
                   alpha=0.85, zorder=5)
        ax.set_xlim(xlim)
        ax.set_yticks([])
        # Restore corner's label conventions for the diagonal: x-label only
        # on bottom row, no y-label.
        if k == ndim - 1:
            ax.set_xlabel(label_list[k], fontsize=label_size)
        ax.spines["left"].set_visible(False)

    # Truth crosshairs on the off-diagonal 2D panels.
    for i in range(ndim):
        for j in range(i):
            ax = axes[i, j]
            ax.axvline(truth_list[j], color=TRUTH_TALK, lw=2.4,
                       alpha=0.85, zorder=5)
            ax.axhline(truth_list[i], color=TRUTH_TALK, lw=2.4,
                       alpha=0.85, zorder=5)

    # Tick + label styling: large rotated x-ticks, at most ~3 ticks
    # per axis so labels never collide.
    for ax in fig.axes:
        ax.tick_params(axis="both", which="major",
                       labelsize=tick_size, length=6, width=1.4, pad=4)
        ax.locator_params(nbins=3)
        for lbl in ax.get_xticklabels():
            lbl.set_rotation(30)
            lbl.set_horizontalalignment("right")

    # Legend: filled patches matching marginal fills + truth line.
    leg_ax = axes[0, -1] if ndim > 1 else axes[0, 0]
    handles = [
        Patch(facecolor=SINGLE_C, edgecolor=SINGLE_C,
              alpha=HIST_ALPHA + 0.2, lw=line_w, label="single-pass"),
        Patch(facecolor=ITER_C, edgecolor=ITER_C,
              alpha=HIST_ALPHA + 0.2, lw=line_w, label="iterative"),
        plt.Line2D([], [], color=TRUTH_TALK, lw=2.4,
                   alpha=0.85, label="truth"),
    ]
    leg_ax.legend(handles=handles, loc="upper right", frameon=False,
                  fontsize=legend_size, handlelength=2.0,
                  handleheight=1.4, borderaxespad=0.6,
                  labelspacing=0.6)

    if savepath:
        fig.savefig(savepath, bbox_inches="tight", dpi=200)
    plt.close(fig)
    return fig


# ----- Figure 4: full 8-param corner, talk style ------------------
print("Plotting corner (talk style, 8 params) ...")
_talk_corner(
    single_chain, iter_chain, single_names, labels, truth_vec,
    figsize=(14, 14),
    label_size=22, tick_size=15, legend_size=24,
    line_w=2.6,
    savepath=f"{FIGS}/step4_corner_talk.png",
)

# ----- Figure 5: focused 3-param corner (shared params) -----------
print("Plotting corner (talk style, 3 shared params) ...")
shared_idx = [single_names.index(p) for p in SHARED]
shared_labels = [LABEL_MAP[p] for p in SHARED]
shared_truth = [truth[p] for p in SHARED]
_talk_corner(
    single_chain[:, shared_idx], iter_chain[:, shared_idx],
    SHARED, shared_labels, shared_truth,
    figsize=(9, 9),
    label_size=28, tick_size=18, legend_size=24,
    line_w=3.2,
    savepath=f"{FIGS}/step4_corner_talk_shared.png",
)

# ============================================================
# Figure 5b: shared 3-param corner — talk style matched to step3
# ------------------------------------------------------------
# Overwrites step4_corner_talk_shared.png produced above so slide 19 uses
# the same visual language as slides 17/18: filled corner.corner contours
# with default histograms, blue base posterior, reddish-purple overlay,
# orange truth crosshairs.
# ============================================================
print("Plotting shared 3-param corner (step3-matched style) ...")

ITER_TALK_C = "#CC79A7"  # Wong reddish purple — distinct from blue base + orange truth

shared_levels = (1 - np.exp(-0.5), 1 - np.exp(-2))
fig5b = corner.corner(
    single_chain[:, shared_idx],
    labels=shared_labels,
    truths=shared_truth,
    truth_color=TRUTH_C,
    color=SINGLE_C,
    levels=shared_levels,
    fill_contours=True,
    plot_datapoints=False,
    plot_density=False,
    smooth=1.5, smooth1d=1.3,
    quantiles=None,
    show_titles=False,
    label_kwargs={"fontsize": 15},
    hist_kwargs={"linewidth": 1.5},
)
corner.corner(
    iter_chain[:, shared_idx],
    labels=shared_labels,
    truths=shared_truth,
    truth_color=TRUTH_C,
    color=ITER_TALK_C,
    levels=shared_levels,
    fill_contours=True,
    plot_datapoints=False,
    plot_density=False,
    smooth=1.5, smooth1d=1.3,
    quantiles=None,
    show_titles=False,
    label_kwargs={"fontsize": 15},
    hist_kwargs={"linewidth": 1.5},
    fig=fig5b,
)

for ax in fig5b.axes:
    ax.tick_params(axis="both", labelsize=11)

legend_elements = [
    Patch(facecolor=SINGLE_C, alpha=0.4, label="single-pass"),
    Patch(facecolor=ITER_TALK_C, alpha=0.4, label="iterative"),
    plt.Line2D([0], [0], color=TRUTH_C, linewidth=1.5, label="truth"),
]
fig5b.legend(handles=legend_elements, loc="upper right", fontsize=16,
             frameon=False, bbox_to_anchor=(0.95, 0.95),
             handlelength=2.0, handleheight=1.4, labelspacing=0.5)

fig5b.savefig(f"{FIGS}/step4_corner_talk_shared.png",
              dpi=300, bbox_inches="tight")
plt.close(fig5b)

print(f"\nWrote figures to {FIGS}/:")
print(f"  step4_corner.png")
print(f"  step4_ci_widths.png")
print(f"  step4_kl_trace.png")
print(f"  step4_corner_talk.png")
print(f"  step4_corner_talk_shared.png")
