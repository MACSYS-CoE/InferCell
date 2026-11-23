#!/usr/bin/env python3
"""F5 (fig:reducedpredict): the joint posterior's predictions on the reduced model, scored against the
true simulator, with the prior predictive beside them. Full text width.

Run from anywhere:
    python3 -B paper/figures/fig_reduced_predict.py

Input, read and never written, in paper/figures/data/ and copied unchanged:
- kim95_score_targets.tsv: score_targets.tsv from /fred/oz022/tkimpson/infercell-kim95-archive-2026-11-19,
  SHA-256 2931a9f97e40d8794294840f4cd56817c582c658626dd189353a56983f458047, as in that archive's SHA256SUMS.
  Runs at baf5331: truth job 21154497 (10,810 cells at theta*), posterior job 21154498 (100 draws x 973
  cells, KIM-92's joint posterior), prior job 21154499 (100 draws x 200 cells); collect and score job
  21154500, with score_targets.tsv rewritten on the login node by corea_predict_kim95_score.jl at 2333ae2.
The result note is dev/scripts/corea_predict_kim95_result.md. A species' offset is the median over its 53
later-horizon saves of z = (truth - posterior median) / posterior SD, with posterior SD = se / se_over_sd.
post_lo and post_hi are the 5% and 95% quantiles over draws; pools and transcripts are at baseline.
In B, a pool's whisker is the 5-95% range of z over its 53 saves, which are strongly correlated: spread, not
uncertainty. Filled marks |median z| > 0.9, which picks out the note's ten.
"""
import os
import sys

sys.dont_write_bytecode = True
import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np

from make_figures import EXACT, INK, MUTED, WIDTH, rows, save  # importing also applies the house rcParams

HERE = os.path.dirname(os.path.abspath(__file__))
POST, PRIOR, PRIOR_EDGE = EXACT, "#d3d2cd", "#8f8e89"  # posterior in the paper's blue, the prior in grey

POOLS = [("M_g6p_c", "G6P"), ("M_f6p_c", "F6P"), ("M_fdp_c", "FDP"), ("M_g3p_c", "G3P"), ("M_dhap_c", "DHAP"),
         ("M_3pg_c", "3PG"), ("M_2pg_c", "2PG"), ("M_lac__L_c", "Lactate in"), ("M_lac__L_e", "Lactate out"),
         ("M_atp_c", "ATP"), ("M_adp_c", "ADP"), ("M_amp_c", "AMP"), ("M_gtp_c", "GTP"), ("M_gdp_c", "GDP"),
         ("M_gmp_c", "GMP"), ("M_pi_c", "Pi"), ("M_trna_chg_c", "Charged tRNA")]
GENES = {"0607": "0607 (GAPD)", "0779": "0779 (ptsG)"}  # names the paper already uses; locus numbers otherwise
# The note's offsets for the ten shifted pools (corea_predict_kim95_result.md).
NOTE = {"M_g3p_c": 1.39, "M_dhap_c": 1.23, "M_atp_c": -1.33, "M_adp_c": -1.00, "M_gdp_c": -1.40,
        "M_amp_c": 1.14, "M_gmp_c": 0.97, "M_g6p_c": 0.94, "M_f6p_c": 1.00, "M_fdp_c": 0.96}


def col(sel, k):
    return np.array([float(r[k]) for r in sel])


def offsets(d, species):
    sel = [r for r in d if r["species"] == species]
    assert len(sel) == 53, (species, len(sel))
    z = (col(sel, "y") - col(sel, "post_med")) / (col(sel, "se") / col(sel, "se_over_sd"))
    return np.median(z), np.percentile(z, [5, 95])


def interval_rows(ax, items, scale, prior=True):
    for y, r in items:
        if prior:
            ax.plot([float(r["prior_lo"]) / scale, float(r["prior_hi"]) / scale], [y + 0.17] * 2, color=PRIOR,
                    lw=4.5, solid_capstyle="butt")
            ax.plot(float(r["prior_med"]) / scale, y + 0.17, "|", ms=5, mew=1.2, color=PRIOR_EDGE)
        ax.plot([float(r["post_lo"]) / scale, float(r["post_hi"]) / scale], [y - 0.17] * 2, color=POST, lw=4.5,
                solid_capstyle="butt")
        ax.plot(float(r["post_med"]) / scale, y - 0.17, "|", ms=4.5, mew=0.8, color="white")
        t, se = float(r["y"]) / scale, float(r["se"]) / scale
        ax.fill_betweenx([y - 0.42, y + 0.42], t - se, t + se, color=INK, alpha=0.2, lw=0)
        ax.plot([t, t], [y - 0.42, y + 0.42], color=INK, lw=0.9)


def main():
    d = rows(os.path.join(HERE, "data", "kim95_score_targets.tsv"))
    by = {r["id"]: r for r in d}
    q, qkd, dq = by["Q"], by["Q_kd"], by["dQ"]

    # Checks against the result note.
    for r in (q, qkd):
        half = (float(r["post_hi"]) - float(r["post_lo"])) / 2 / float(r["post_med"])
        low = 1 - float(r["prior_med"]) / float(r["y"])
        span = float(r["prior_hi"]) / float(r["prior_lo"])
        print(f"{r['id']}: posterior 90% ±{100 * half:.2f}%, prior median {100 * low:.1f}% low, "
              f"prior 90% spans ×{span:.2f}, skill {float(r['skill']):.3f}")
        assert round(100 * half, 1) == 1.4 and 0.24 <= low <= 0.265 and round(span, 1) == 2.8
        assert r["post90"] == "true" and 0.98 <= float(r["skill"]) < 0.99  # the note: "0.98 on both"
    assert dq["post90"] == "true" and float(dq["post_lo"]) < 0 < float(dq["post_hi"])
    assert round(float(dq["mc_post"]), 2) == 1.09
    pool = [r for r in d if r["family"] == "pool"]
    tx = [r for r in d if r["family"] == "transcript"]
    assert len(pool) == len(tx) == 901
    assert sum(r["post90"] == "true" for r in tx) == 900
    h90 = np.mean([r["post90"] == "true" for r in pool])
    h50 = np.mean([r["post50"] == "true" for r in pool])
    assert round(100 * h90) == 92 and round(100 * h50) == 34, (h90, h50)
    zp = {s: offsets(d, s) for s, _ in POOLS}
    tx_species = sorted({r["species"] for r in tx})
    zt = {s: offsets(d, s) for s in tx_species}
    assert len(zp) == len(zt) == 17
    shifted = {s for s in zp if abs(zp[s][0]) > 0.9}
    assert shifted == set(NOTE), shifted
    for s, v in NOTE.items():
        assert round(zp[s][0], 2) == v, (s, zp[s][0], v)
    assert max(abs(zp[s][0]) for s in zp if s not in shifted) < 0.8
    assert max(abs(v[0]) for v in zt.values()) < 0.64
    for s, _ in POOLS:
        print(f"  {s:14s} {zp[s][0]:+.3f}")
    print(f"pools: {100 * h90:.1f}% in 90%, {100 * h50:.1f}% in 50%; shifted "
          f"{min(abs(zp[s][0]) for s in shifted):.3f}-{max(abs(zp[s][0]) for s in shifted):.3f}; "
          f"transcripts max |z| {max(abs(v[0]) for v in zt.values()):.3f}")

    fig = plt.figure(figsize=(WIDTH, 3.7), layout="constrained")
    left, right = fig.subfigures(1, 2, width_ratios=[1.55, 1], wspace=0.02)
    top, bottom = left.subfigures(2, 1, height_ratios=[1, 1.15], hspace=0.04)

    # A: Q and Q_kd on the prior's scale, the same on the posterior's scale, and the knockdown effect.
    a_full, a_zoom, a_dq = top.subplots(1, 3, width_ratios=[1.15, 1, 0.95])
    items = [(1, q), (0, qkd)]
    interval_rows(a_full, items, 1e6)
    a_full.set_yticks([1, 0], [r"$Q$", r"$Q_{\mathrm{kd}}$"])
    a_full.set_xlim(5, 18.5)
    a_full.set_xlabel("Millions")
    a_full.set_title("Full range", loc="left", color=MUTED)
    interval_rows(a_zoom, items, 1e6, prior=False)
    a_zoom.set_yticks([1, 0], [])
    a_zoom.set_xlim(14.75, 15.4)
    a_zoom.set_xticks([14.8, 15.0, 15.2])  # fits once B is one axis
    a_zoom.set_xlabel("Millions")
    a_zoom.set_title("Zoomed", loc="left", color=MUTED)
    lo, hi = a_zoom.get_xlim()
    a_full.axvspan(lo, hi, color=POST, alpha=0.12, lw=0)
    interval_rows(a_dq, [(0, dq)], 1e3)
    a_dq.axvline(0, color=MUTED, lw=0.5, ls=(0, (2, 2)))
    a_dq.set_yticks([0], [r"$\Delta Q$"])
    a_dq.set_xlabel("Thousands")
    a_dq.set_title("Knockdown", loc="left", color=MUTED)
    for ax in (a_full, a_zoom):
        ax.set_ylim(-0.6, 1.6)
    a_dq.set_ylim(-0.9, 0.9)
    for ax in (a_full, a_zoom, a_dq):
        ax.tick_params(axis="y", length=0)
        ax.spines["left"].set_visible(False)
    a_full.text(float(q["prior_lo"]) / 1e6, 1.45, "Prior", color=PRIOR_EDGE, ha="left", va="bottom", fontsize=7)
    a_full.text(float(q["post_lo"]) / 1e6 - 0.3, 0.83, "Posterior", color=POST, ha="right", va="center",
                fontsize=7)
    a_zoom.text(float(q["y"]) / 1e6, 1.45, "Truth ± SE", color=INK, ha="center", va="bottom", fontsize=7)

    # C: one shifted pool and the knocked-down gene's transcript over the later horizon.
    c = bottom.subplots(1, 2)
    for ax, species, scale, title, ylab in ((c[0], "M_g3p_c", 1e3, f"G3P, offset {zp['M_g3p_c'][0]:+.2f} SD", "Thousands"),
                                            (c[1], "mRNA_JCVISYN3A_0779", 1, "Transcript 0779 (ptsG), covered",
                                             "Mean per cell")):
        sel = sorted((r for r in d if r["species"] == species), key=lambda r: float(r["t"]))
        t = col(sel, "t") / 60
        ax.fill_between(t, col(sel, "prior_lo") / scale, col(sel, "prior_hi") / scale, color=PRIOR, lw=0,
                        label="Prior, 90%")
        ax.fill_between(t, col(sel, "post_lo") / scale, col(sel, "post_hi") / scale, color=POST, alpha=0.75,
                        lw=0, label="Posterior, 90%")
        ax.plot(t, col(sel, "y") / scale, color=INK, lw=0.9, label="Truth")
        ax.set_xlim(t[0], t[-1])
        ax.set_xticks([60, 80, 100])
        ax.set_xlabel("Time, min")
        ax.set_ylabel(ylab)
        ax.set_title(title, loc="left", color=MUTED)
    c[1].set_yscale("log")
    c[1].set_yticks([0.2, 0.5, 1, 2, 5], ["0.2", "0.5", "1", "2", "5"])
    c[1].minorticks_off()
    c[0].legend(loc="upper left", fontsize=7, handlelength=1.2, borderaxespad=0.2, labelspacing=0.25)

    # B: each pool's offset over the later horizon, and the 17 transcripts' on one row, on one z scale.
    ax = right.subplots()
    y = np.arange(len(POOLS))[::-1] + 1.4
    ax.axvspan(-1, 1, color="#f0efeb", lw=0, zorder=0)
    ax.axvline(0, color=MUTED, lw=0.5)
    for yi, (s, _) in zip(y, POOLS):
        med, (lo, hi) = zp[s]
        ax.plot([lo, hi], [yi, yi], color=POST, lw=0.8, alpha=0.6)
        ax.plot(med, yi, "o", ms=4.2, mew=1, mec=POST, mfc=POST if s in shifted else "white", zorder=3)
    tz = np.array([zt[s][0] for s in tx_species])
    ax.plot(tz, np.zeros_like(tz), "o", ms=4.2, mew=0.8, mec=POST, mfc="white", alpha=0.8, zorder=3)
    ax.axhline(0.7, color=MUTED, lw=0.5)
    ax.text(-0.95, y[0] + 0.9, "±1 SD", ha="left", va="center", fontsize=7, color=MUTED)
    ax.set_yticks(list(y) + [0], [n for _, n in POOLS] + [f"17 transcripts,\nall within {max(abs(tz)):.2f}"])
    ax.set_ylim(-0.9, y[0] + 1.3)
    ax.tick_params(axis="y", length=0)
    ax.spines["left"].set_visible(False)
    ax.set_xlim(-2.1, 2.1)
    ax.set_xticks([-2, -1, 0, 1, 2], ["−2", "−1", "0", "+1", "+2"])
    ax.set_xlabel("Truth − posterior median,\nin posterior SDs")
    lo10, hi10 = min(abs(zp[s][0]) for s in shifted), max(abs(zp[s][0]) for s in shifted)
    ax.set_title(f"Filled: {len(shifted)} pools offset {lo10:.2f}–{hi10:.2f} SD", loc="right", color=POST)

    for sub, letter in ((top, "A"), (right, "B"), (bottom, "C")):
        sub.text(0.0, 1.0, letter, transform=sub.transSubfigure, ha="left", va="top", fontweight="bold",
                 fontsize=9, color=INK)
    fig.suptitle("SCORED against the truth", x=0.995, ha="right", color=MUTED, fontweight="bold", fontsize=8)
    save(fig, "fig_reduced_predict")


if __name__ == "__main__":
    main()
