"""KIM-161: the promoter shift and the ENO/FBA marginals, from the KIM-92 pilot archive.

Read-only and light (two TSVs, no model evaluation):

    python3 dev/scripts/pilot_promoter_shift.py /fred/oz022/tkimpson/infercell-kim92-archive-2026-11-17

It reads pilot_draws.tsv (the 1,600 pooled kept draws, 4 chains x (200 + 200), one column per
coordinate, natural-log scale) and pilot_coordinates.tsv (per coordinate: posterior mean and
SD, and the truth), both written by dev/scripts/corea_pilot_kim92.jl at de24b27, and prints:

  1. the promoters' standardised errors z_i = (posterior mean - truth) / posterior SD, their
     mean, and the standard error of that mean two ways: independent (1/sqrt(17)) and with
     the posterior's correlations among the 17 promoters (sqrt(1' R 1) / 17), with the mean
     over each SE;
  2. the promoters' mean pairwise posterior correlation (off-diagonal of R);
  3. the promoters' mean shift on the natural log, mean(posterior mean - truth), and the
     corresponding factor exp(shift) - 1;
  4. the sample skewness of the ENO and FBA marginals (log kcatF), from the pooled draws.

The z-vector's covariance is taken as R, the posterior correlation of the promoters: under
calibration each (mean - truth)/sd has unit variance and the errors share the posterior's
correlation structure. That is an approximation, stated as such.
"""
import csv
import math
import os
import sys

import numpy as np


def main():
    if len(sys.argv) != 2:
        sys.exit(__doc__)
    out = sys.argv[1]
    with open(os.path.join(out, "pilot_draws.tsv")) as fh:
        rd = csv.reader(fh, delimiter="\t")
        head = next(rd)
        draws = np.array([[float(x) for x in row] for row in rd])
    coords = {}
    with open(os.path.join(out, "pilot_coordinates.tsv")) as fh:
        for row in csv.DictReader(fh, delimiter="\t"):
            coords[row["coordinate"]] = row
    assert draws.shape == (1600, len(head)), draws.shape
    assert set(coords) == set(head), "coordinate names differ between the two files"

    prom = [c for c in head if c.startswith("log_S_")]
    assert len(prom) == 17, len(prom)
    idx = [head.index(c) for c in prom]
    pm = np.array([float(coords[c]["post_mean"]) for c in prom])
    ps = np.array([float(coords[c]["post_sd"]) for c in prom])
    tr = np.array([float(coords[c]["truth"]) for c in prom])

    # the TSV's means and SDs agree with the draws
    d = draws[:, idx]
    assert np.allclose(d.mean(0), pm, atol=1e-6), "post_mean disagrees with the draws"
    assert np.allclose(d.std(0, ddof=1), ps, rtol=1e-3), "post_sd disagrees with the draws"

    z = (pm - tr) / ps
    R = np.corrcoef(d, rowvar=False)
    n = len(prom)
    se_ind = 1 / math.sqrt(n)
    se_cor = math.sqrt(R.sum()) / n
    off = R[~np.eye(n, dtype=bool)]
    shift = (pm - tr).mean()

    print("## 1. Promoters' standardised errors (n = %d)" % n)
    print("mean z %.3f, SD %.3f, below truth %d of %d" % (z.mean(), z.std(ddof=1), (z < 0).sum(), n))
    print("SE independent %.3f -> mean/SE %.2f" % (se_ind, z.mean() / se_ind))
    print("SE with the posterior's correlations %.3f -> mean/SE %.2f" % (se_cor, z.mean() / se_cor))
    print("## 2. Mean pairwise posterior correlation among the promoters: %.3f (min %.3f, max %.3f)"
          % (off.mean(), off.min(), off.max()))
    print("## 3. Mean shift in log rate: %.4f, factor exp(shift) - 1 = %.2f%%" % (shift, 100 * (math.exp(shift) - 1)))
    print("## 4. Skewness of the marginals (pooled draws, log scale)")
    for c in ("log_kcatF_R_ENO", "log_kcatF_R_FBA"):
        x = draws[:, head.index(c)]
        m, s = x.mean(), x.std(ddof=0)
        g1 = ((x - m) ** 3).mean() / s ** 3
        print("  %-18s skewness %.3f" % (c, g1))


if __name__ == "__main__":
    main()
