"""Step 3 publication-quality figures — Python/matplotlib version."""

import numpy as np
import matplotlib.pyplot as plt
import corner

plt.rcParams.update({
    "font.family": "serif",
    "font.serif": ["Computer Modern Roman"],
    "text.usetex": True,
    "axes.spines.top": False,
    "axes.spines.right": False,
    "axes.grid": False,
    "figure.dpi": 150,
})

# Font sizes are set for the on-slide size, not the on-screen size. The 4-panel
# corner exports at ~9.4in wide and is shown 6.2cm tall in the talk, i.e. shrunk
# ~4x, so anything below ~25pt here is unreadable on a projector.
FS_LABEL = 38   # axis labels (parameter names)
FS_TICK = 28    # tick numbers
FS_TITLE = 26   # per-panel quantile summaries
FS_LEGEND = 30   # kept below FS_LABEL so it clears the diagonal panel titles
MAX_TICKS = 3   # thin the ticks; 5 (corner's default) crowds at this scale
LABELPAD = 0.28  # push axis labels clear of the enlarged tick numbers

# corner's default titles stack the asymmetric errors as sub/superscripts, which
# render at ~0.7x the title size (~5pt on a 6.2cm slide) and make the string ~2x
# wider than its panel. Median-only titles fit the panel and stay legible.
# Set False to restore corner's stacked-quantile titles.
TITLE_MEDIAN_ONLY = True


def set_median_titles(fig, samples, labels, fontsize):
    """Retitle the diagonal panels with 'label = median' only."""
    k = samples.shape[1]
    for i in range(k):
        med = np.median(samples[:, i])
        base = labels[i].strip("$")
        fig.axes[i * k + i].set_title(
            rf"${base} = {med:.2f}$", fontsize=fontsize, pad=12
        )

# --- Load data ---
ode_chain = np.genfromtxt("examples/results/ode_chain.csv", delimiter=",", skip_header=1)
cond_particles = np.genfromtxt("examples/results/cond_particles.csv", delimiter=",", skip_header=1)
cond_weights = np.genfromtxt("examples/results/cond_weights.csv", delimiter=",")
uncond_particles = np.genfromtxt("examples/results/uncond_particles.csv", delimiter=",", skip_header=1)
uncond_weights = np.genfromtxt("examples/results/uncond_weights.csv", delimiter=",")

# ============================================================
# Figure 1: Corner plot — joint ODE posteriors (8 parameters)
# ============================================================

param_labels = [
    r"$k_{tx}$", r"$k_{tl}$", r"$\gamma_{mRNA}$", r"$\gamma_{prot}$",
    r"$k_{atp}$", r"$k_{ntp}$", r"$k_{aa}$", r"$\sigma$",
]
truths = [1.0, 2.0, 0.5, 0.1, 1.0, 0.2, 0.4, 0.3]

# The full 8x8 is shown at only 3.35cm tall in the talk, where no per-panel
# label can be legible — it reads as texture (the shape of the degeneracies).
# Fonts are bumped so it survives being shown larger, but the legible cut is
# the reduced 4x4 exported below.
fig1 = corner.corner(
    ode_chain,
    labels=param_labels,
    truths=truths,
    truth_color="#E69F00",       # wong orange
    quantiles=[0.16, 0.5, 0.84],
    show_titles=True,
    title_kwargs={"fontsize": 22},
    label_kwargs={"fontsize": 30},
    max_n_ticks=MAX_TICKS,
    labelpad=LABELPAD,
    smooth=1.5,
    smooth1d=1.3,
    levels=(1 - np.exp(-0.5), 1 - np.exp(-2), 1 - np.exp(-4.5)),  # 1,2,3 sigma
    fill_contours=True,
    plot_datapoints=True,
    plot_density=False,
    data_kwargs={"alpha": 0.05, "ms": 1},
    color="#0072B2",             # wong blue
)
for ax in fig1.axes:
    ax.tick_params(axis="both", labelsize=20)
fig1.savefig("examples/figures/step3_corner.png", dpi=300, bbox_inches="tight")
fig1.savefig("examples/figures/step3_corner.pdf", bbox_inches="tight")
print("Saved examples/figures/step3_corner.{png,pdf}")
plt.close(fig1)

# ============================================================
# Figure 1b: Reduced corner — the degenerate block only
#   k_tx, k_tl, gamma_mRNA, k_ntp are mutually correlated at
#   |r| = 0.68-0.88; sigma_obs is uncorrelated with everything.
#   A 4x4 cut carries the degeneracy story at a legible size.
# ============================================================

reduced_idx = [0, 1, 2, 5]   # k_tx, k_tl, gamma_mRNA, k_ntp
fig1b = corner.corner(
    ode_chain[:, reduced_idx],
    labels=[param_labels[i] for i in reduced_idx],
    truths=[truths[i] for i in reduced_idx],
    truth_color="#E69F00",
    quantiles=[0.16, 0.5, 0.84],
    show_titles=not TITLE_MEDIAN_ONLY,
    title_kwargs={"fontsize": FS_TITLE},
    label_kwargs={"fontsize": FS_LABEL},
    max_n_ticks=MAX_TICKS,
    labelpad=LABELPAD,
    smooth=1.5,
    smooth1d=1.3,
    levels=(1 - np.exp(-0.5), 1 - np.exp(-2), 1 - np.exp(-4.5)),
    fill_contours=True,
    plot_datapoints=True,
    plot_density=False,
    data_kwargs={"alpha": 0.05, "ms": 1},
    color="#0072B2",
)
if TITLE_MEDIAN_ONLY:
    set_median_titles(fig1b, ode_chain[:, reduced_idx],
                      [param_labels[i] for i in reduced_idx], FS_TITLE)
for ax in fig1b.axes:
    ax.tick_params(axis="both", labelsize=FS_TICK)
fig1b.savefig("examples/figures/step3_corner_reduced.png", dpi=300, bbox_inches="tight")
fig1b.savefig("examples/figures/step3_corner_reduced.pdf", bbox_inches="tight")
print("Saved examples/figures/step3_corner_reduced.{png,pdf}")
plt.close(fig1b)

# ============================================================
# Figure 2: SSA posteriors — conditioned vs unconditioned
#   Corner plot with two overlaid series
# ============================================================

ssa_labels = [r"$k_{tx}$", r"$k_{tl}$", r"$\gamma_{mRNA}$", r"$\gamma_{prot}$"]
ssa_truths = [1.0, 2.0, 0.5, 0.1]

# Weighted resampling
rng = np.random.default_rng(42)
n_resample = 5000

def resample(particles, weights, n, rng):
    w = weights / weights.sum()
    idx = rng.choice(len(w), size=n, p=w)
    return particles[idx]

cond_samples = resample(cond_particles, cond_weights, n_resample, rng)
uncond_samples = resample(uncond_particles, uncond_weights, n_resample, rng)

# Axis ranges centered on the posteriors (not dominated by unconditioned tails)
# Each range: ~0 to ~3x the true value
ssa_range = [(0, 4.0), (0, 8.0), (0, 2.0), (0, 0.5)]

# Unconditioned first (background), then conditioned on top.
# Convention across the talk: blue = base posterior; comparison overlay differs.
color_uncond = "#0072B2"  # wong blue (base)
color_cond = "#009E73"    # wong green (conditioned overlay)

fig2 = corner.corner(
    uncond_samples,
    range=ssa_range,
    labels=ssa_labels,
    truths=ssa_truths,
    truth_color="#E69F00",
    quantiles=[0.16, 0.5, 0.84],
    show_titles=False,
    label_kwargs={"fontsize": FS_LABEL},
    max_n_ticks=MAX_TICKS,
    labelpad=LABELPAD,
    smooth=1.4,
    smooth1d=1.1,
    levels=(1 - np.exp(-0.5), 1 - np.exp(-2)),  # 1,2 sigma
    fill_contours=True,
    plot_datapoints=False,
    color=color_uncond,
    hist_kwargs={"linewidth": 3.0},
)

# Overlay conditioned
corner.corner(
    cond_samples,
    range=ssa_range,
    labels=ssa_labels,
    truths=ssa_truths,
    truth_color="#E69F00",
    quantiles=[0.16, 0.5, 0.84],
    show_titles=not TITLE_MEDIAN_ONLY,
    title_kwargs={"fontsize": FS_TITLE},
    label_kwargs={"fontsize": FS_LABEL},
    max_n_ticks=MAX_TICKS,
    labelpad=LABELPAD,
    smooth=1.4,
    smooth1d=1.1,
    levels=(1 - np.exp(-0.5), 1 - np.exp(-2)),
    fill_contours=True,
    plot_datapoints=False,
    color=color_cond,
    hist_kwargs={"linewidth": 3.0},
    fig=fig2,
)

if TITLE_MEDIAN_ONLY:
    # titles summarise the conditioned posterior (the green series)
    set_median_titles(fig2, cond_samples, ssa_labels, FS_TITLE)
for ax in fig2.axes:
    ax.tick_params(axis="both", labelsize=FS_TICK, width=1.6, length=8)

# Manual legend
from matplotlib.patches import Patch
from matplotlib.lines import Line2D
legend_elements = [
    Patch(facecolor=color_uncond, alpha=0.4, label="Unconditioned (ABC-SMC)"),
    Patch(facecolor=color_cond, alpha=0.4, label="Conditioned (ABC-SMC)"),
    Line2D([0], [0], color="#E69F00", linewidth=3.0, label="True value"),
]
fig2.legend(handles=legend_elements, loc="upper right", fontsize=FS_LEGEND,
            frameon=False, bbox_to_anchor=(0.99, 0.99),
            handlelength=1.6, handleheight=1.2, labelspacing=0.35)

fig2.savefig("examples/figures/step3_ssa_posteriors.png", dpi=300, bbox_inches="tight")
fig2.savefig("examples/figures/step3_ssa_posteriors.pdf", bbox_inches="tight")
print("Saved examples/figures/step3_ssa_posteriors.{png,pdf}")
plt.close(fig2)

print("Done.")
