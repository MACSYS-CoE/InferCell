"""16b.3's interval coverage, read off V7's ranks (16a.11's proposal, spec/phases/
16-recovery.md §12 2026-10-22, "16a.11"). No new compute.

Rank r of L = 100 thinned draws gives the truth's quantile (r + 0.5)/(L + 1). A
central interval at level q covers the truth when that quantile lies strictly
inside ((1 - q)/2, (1 + q)/2). The band is the binomial 2.5% and 97.5% quantiles
at R replicates. With L = 100, the 99% interval sits next to the extreme ranks.

    python3 dev/scripts/v7_coverage_16a11.py   # needs numpy, scipy
"""
import os

import numpy as np
from scipy.stats import binom

HERE = os.path.dirname(os.path.abspath(__file__))
L = 100
LEVELS = (0.5, 0.8, 0.9, 0.95, 0.99)


def main():
    path = os.path.join(HERE, "v7_16a9b_ranks.tsv")
    with open(path) as f:
        names = f.readline().rstrip("\n").split("\t")[2:]
    ranks = np.loadtxt(path, skiprows=1)[:, 2:]
    R = len(ranks)
    u = (ranks + 0.5) / (L + 1)
    print(f"V7: R = {R}, L = {L}. Covered of {R}, with the binomial 95% band.\n")
    print("| Level | " + " | ".join(names) + " | Band |")
    print("|---" * (len(names) + 2) + "|")
    for q in LEVELS:
        cov = ((u > (1 - q) / 2) & (u < (1 + q) / 2)).sum(axis=0)
        lo, hi = int(binom.ppf(0.025, R, q)), int(binom.ppf(0.975, R, q))
        cells = [f"{c}" if lo <= c <= hi else f"**{c}**" for c in cov]
        print(f"| {q:.0%} | " + " | ".join(cells) + f" | {lo}–{hi} |")


if __name__ == "__main__":
    main()
