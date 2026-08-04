"""Calibration curve for the iterative boundary protocol (issue #11).

Reads shard CSVs written by `examples/calibration_sweep.jl`, computes empirical
coverage of ground truth across seeds at each nominal CI level, and plots
empirical vs nominal for the single-pass and iterative protocols side-by-side.
"""

import glob
import os

import numpy as np
import pandas as pd
import matplotlib.pyplot as plt

plt.rcParams.update({
    "font.family": "serif",
    "font.serif": ["Computer Modern Roman"],
    "text.usetex": True,
    "axes.spines.top": False,
    "axes.spines.right": False,
    "axes.grid": False,
    "figure.dpi": 150,
})

RESULTS = "examples/results"
FIGS = "examples/figures"
os.makedirs(FIGS, exist_ok=True)

# Wong-palette colours, one per ODE parameter.
PARAM_COLOURS = {
    "k_tx":          "#0072B2",  # blue
    "k_tl":          "#D55E00",  # vermillion
    "gamma_mRNA":    "#009E73",  # green
    "gamma_protein": "#E69F00",  # orange
}
PARAM_ORDER = ["k_tx", "k_tl", "gamma_mRNA", "gamma_protein"]
PARAM_LABELS = {
    "k_tx":          r"$k_{\mathrm{tx}}$",
    "k_tl":          r"$k_{\mathrm{tl}}$",
    "gamma_mRNA":    r"$\gamma_{\mathrm{mRNA}}$",
    "gamma_protein": r"$\gamma_{\mathrm{protein}}$",
}

# --- Load + concatenate shard CSVs ---
shard_paths = sorted(glob.glob(f"{RESULTS}/step4_calibration_shard_*.csv"))
if not shard_paths:
    raise SystemExit(f"No shard CSVs found in {RESULTS}/step4_calibration_shard_*.csv")
df = pd.concat([pd.read_csv(p) for p in shard_paths], ignore_index=True)

n_seeds = df["seed"].nunique()
print(f"Loaded {len(df)} rows from {len(shard_paths)} shard(s); {n_seeds} unique seeds")

# Empirical coverage per (protocol, param, level): fraction of seeds where covered=1.
coverage = (df.groupby(["protocol", "param", "level"])["covered"]
              .mean()
              .reset_index())

# --- Plot ---
fig, axes = plt.subplots(1, 2, figsize=(10, 4.5), sharey=True)
protocols = [("single", "Single-pass"), ("iter", "Iterative")]

for ax, (proto_key, proto_label) in zip(axes, protocols):
    ax.plot([0, 1], [0, 1], color="k", linestyle="--", alpha=0.3, linewidth=1,
            label="nominal = empirical")
    sub = coverage[coverage["protocol"] == proto_key]
    for param in PARAM_ORDER:
        rows = sub[sub["param"] == param].sort_values("level")
        ax.plot(rows["level"], rows["covered"],
                marker="o", color=PARAM_COLOURS[param],
                label=PARAM_LABELS[param], linewidth=1.5)
    ax.set_xlim(0.75, 1.01)
    ax.set_ylim(0.0, 1.05)
    ax.set_xticks([0.80, 0.90, 0.95, 0.99])
    ax.set_xlabel("Nominal CI level")
    ax.set_title(proto_label)

axes[0].set_ylabel(f"Empirical coverage (n={n_seeds} seeds)")
axes[1].legend(loc="lower right", frameon=False)

out_png = f"{FIGS}/step4_calibration_curve.png"
out_pdf = f"{FIGS}/step4_calibration_curve.pdf"
fig.savefig(out_png, bbox_inches="tight", dpi=150)
fig.savefig(out_pdf, bbox_inches="tight", dpi=300)
plt.close(fig)

print(f"Wrote {out_png}")
print(f"Wrote {out_pdf}")

# Print the underlying table for quick inspection / paste into the writeup.
print("\nCoverage table:")
pivot = coverage.pivot_table(index=["protocol", "param"], columns="level",
                              values="covered")
print(pivot.to_string(float_format=lambda x: f"{x:.2f}"))
