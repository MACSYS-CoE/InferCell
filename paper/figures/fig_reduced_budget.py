#!/usr/bin/env python3
"""F6 (fig:reducedbudget): Q's parameter variance in three states, one column wide.

Run from anywhere:
    python3 paper/figures/fig_reduced_budget.py

Inputs, read and never written, both in paper/figures/data/ and copied unchanged:
- kim96_budget_rows.tsv: budget_rows.tsv from /fred/oz022/tkimpson/infercell-kim96-archive-2026-11-19,
  written by Slurm job 21188187 at c507315 (dev/scripts/corea_budget_kim96.sh); SHA-256
  8328257b79149649be905d6cef2c0c3c3c68b8b8534cd0a5c303b438ea0a04eb, as in that archive's SHA256SUMS.
  Its three states: prior (KIM-95 prior arm, jobs 21154499/21154500 at baf5331), KIM-42 (width 1),
  posterior (KIM-92 joint posterior, merge job 21052425; KIM-95 posterior arm, job 21154498).
- kim42_budget_rows.tsv: budget_rows.tsv from /fred/oz022/tkimpson/infercell-m5-kim42-2026-11-05,
  jobs 20041038-40 at f23fd06; SHA-256 48b37c9e145505863203db01a0b6756dbfff4952697326530fe90583ab8b9319,
  as in that archive's SHA256SUMS and in the KIM-96 archive's inputs.tsv. Used for the width-1/2 and
  width-2 totals only (target q, law per_theta_I, quantity total).
The result note is dev/scripts/corea_budget_kim96_result.md.
"""
import os
import sys

import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np

sys.dont_write_bytecode = True
from make_figures import EXACT, INK, MIXED, MUTED, rows  # importing also applies the house rcParams

HERE = os.path.dirname(os.path.abspath(__file__))
COLUMN = 229.9 / 72.27  # \columnwidth in inches: printed at 1:1


def main():
    d = rows(os.path.join(HERE, "data", "kim96_budget_rows.tsv"))
    get = lambda s, q: next((float(r["estimate"]), float(r["se"])) for r in d
                            if r["state"] == s and r["quantity"] == q)
    k42 = rows(os.path.join(HERE, "data", "kim42_budget_rows.tsv"))
    width_total = {float(r["width"]): float(r["estimate"]) for r in k42
                   if r["target"] == "q" and r["law"] == "per_theta_I" and r["quantity"] == "total"}

    states = ["prior", "KIM-42", "posterior"]
    var = np.array([get(s, "var_theta")[0] for s in states])
    se = np.array([get(s, "var_theta")[1] for s in states])
    cv = np.array([get(s, "sd_theta")[0] / get(s, "mean")[0] for s in states])
    held, _ = get("KIM-42", "held_share_I_first")
    assert np.isclose(width_total[1.0], var[1])  # the two archives agree on the middle state
    ratio = var[1] / var[2]
    # Log-scale interval, r * exp(±2 sqrt(sum of squared relative SEs)), as in the result note.
    half = 2 * np.hypot(se[1] / var[1], se[2] / var[2])
    lo, hi = ratio * np.exp(-half), ratio * np.exp(half)
    cmin, _ = get("contraction", "minimum")
    reached, _ = get("contraction", "reached_count")
    rest = (1 - held) * var[1]  # what M0's five carry in the middle state (the note's i_between)
    print(f"Var: {var}; SE: {se}; SD/mean: {cv}; held share {held:.3f}; rest {rest:.3g}")
    print(f"contraction: {reached:.0f} of 21 reach 0.5, minimum {cmin:.4f}")
    print(f"middle/posterior {ratio:.0f} ({lo:.0f}-{hi:.0f}); prior/middle {var[0] / var[1]:.1f}")
    print("width ½, 1, 2 totals:", [f"{width_total[w]:.3g}" for w in (0.5, 1.0, 2.0)],
          "factors:", [f"{width_total[w] / var[2]:.0f}" for w in (0.5, 1.0, 2.0)])

    fig, ax = plt.subplots(figsize=(COLUMN, 2.45), layout="constrained")
    x = np.array([0.0, 1.0, 2.0])
    ax.set_yscale("log")
    # The middle state split: the 15 held promoters' share, and the five's remainder.
    ax.plot([x[1], x[1]], [rest, var[1]], color=MIXED, lw=3, solid_capstyle="butt", zorder=2)
    ax.plot(x[1], rest, "o", ms=4, mfc="white", mec=EXACT, mew=1, zorder=3)
    ax.text(x[1] - 0.07, np.sqrt(rest * var[1]), f"{held:.0%} from the 15\npromoters",
            ha="right", va="center", color=MIXED, fontsize=7)
    ax.text(x[1] - 0.07, rest, f"{1 - held:.0%} from the five", ha="right", va="center",
            color=EXACT, fontsize=7)
    # Width dependence of the middle state: its total at prior widths ½ and 2.
    xw = x[1] + 0.13
    ws = [width_total[w] for w in (0.5, 1.0, 2.0)]
    # A dashed grey range, unlike the states' solid black marks: the middle state's total at three widths.
    ax.plot([xw, xw], [ws[0], ws[2]], color=MUTED, lw=0.8, ls=(0, (2, 1.5)), zorder=1)
    ax.plot([xw] * 3, ws, "_", ms=5, mew=0.8, color=MUTED)
    ax.text(xw + 0.04, width_total[2.0], "prior width ×2", ha="left", va="center", color=MUTED,
            fontsize=7)
    ax.text(xw, width_total[0.5] / 1.15, "×½", ha="center", va="top", color=MUTED, fontsize=7)
    # The three states, with jackknife ±1 SE.
    ax.errorbar(x, var, yerr=se, fmt="o", ms=3.2, color=INK, capsize=2.5, capthick=0.8, lw=1, zorder=4)
    for xi, v, c in zip(x, var, cv):
        label = f"SD {100 * c:.0f}% of $Q$" if c > 0.05 else f"SD {100 * c:.1f}% of $Q$"
        if xi == x[2]:
            label += f"\nevery contraction ≥ {np.floor(cmin * 1000) / 1000:.3f}"
        ax.text(xi + (0.12 if xi == 0 else -0.12), v, label,
                ha="left" if xi == 0 else "right", va="center", fontsize=7, color=INK)
    # The step the claim is about, rounded to two figures (relative SE about 30%).
    ax.text(1.6, 1.8e11, f"≈{float(f'{ratio:.2g}'):.0f}× lower\n({lo:.0f}–{hi:.0f})",
            ha="left", va="bottom", fontsize=7, color=INK)

    ax.set_xticks(x, ["All at\nprior", "Five inferred,\n15 at prior", "All 21\ninferred"])
    ax.tick_params(axis="x", length=0)
    ax.set_xlim(-0.35, 2.35)
    ax.set_ylim(4e9, 4e13)
    ax.set_ylabel(r"Var$_\theta(Q)$ (particles$^2$)")
    fig.savefig(os.path.join(HERE, "fig_reduced_budget.pdf"), metadata={"CreationDate": None})


if __name__ == "__main__":
    main()
