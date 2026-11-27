"""KIM-180: apply the doubling-time definition of spec/prereg/kim-180.md §1 to the default run.

    python3 -I dev/scripts/fullcell_doubling_kim180/doubling_reads.py <outdir>

Reads two default ensembles and prints, for each, the per-draw doubling time and the
check against Thornburg et al. 2022 (Cell 185:345, Fig. 3F: surface area doubles at
88-112 min, mean 97):

- syn3a-ensemble's 10,000 default cells (results/cycle/cells.csv, 120 recorded minutes,
  syn3a-ensemble 2075809), which this definition did not produce;
- this issue's 100 default cells to minute 180 (<outdir>/ensemble/ensemble_cells.csv,
  KIM-183's driver).

The definition (prereg §1):
- t_SA2 is the first recorded minute with CellSA >= 2 x its minute-0 value, +inf if it is
  not reached by the end of the window (minute 180).
- A cell is unhealthy if it collapses (fewer than 59 hooks in a minute, t_short) or ends in
  an ODE failure or an error at or before its t_SA2, or anywhere in the window if t_SA2 is
  not reached. Unhealthy cells are left out.
- The draw's doubling time is the ceil(n/2)-th smallest t_SA2 over its n healthy cells.
"""
import csv
import math
import sys

LO, HI, PAPER_MEAN = 88, 112, 97
SYN3A_CELLS = "/fred/oz022/tkimpson/syn3a-ensemble/results/cycle/cells.csv"


def num(x, missing=math.inf):
    """A CSV event minute: blank (syn3a-ensemble) or -1 (KIM-183's driver) is not reached."""
    if x in ("", None):
        return missing
    v = float(x)
    return missing if v < 0 else v


def median_upper(v):
    v = sorted(v)
    return v[math.ceil(len(v) / 2) - 1] if v else math.nan


def doubling(cells):
    """cells: (t_SA2, t_end) per cell, t_end the minute it collapsed or failed (inf if
    neither). Returns the draw's doubling time and the counts."""
    healthy = [t for t, end in cells if not (math.isfinite(end) and end <= t)]
    return median_upper(healthy), len(healthy), sum(math.isfinite(t) for t in healthy)


def syn3a_cells():
    out = []
    for r in csv.DictReader(open(SYN3A_CELLS)):
        t = num(r["t_SA2"])
        end = 0.0 if r["status"] != "ok" else num(r["t_short"])
        out.append((t, end))
    return out


def ensemble_cells(path, window):
    out = []
    for r in csv.DictReader(open(path)):
        if r["status"] == "error":
            out.append((math.inf, 0.0))
            continue
        t = num(r["t_SA2"])
        t = t if t <= window else math.inf
        end = num(r["t_short"])
        if r["status"] == "ode_failure":
            end = min(end, float(r["failure"].split(":")[0]))
        out.append((t, end))
    return out


def report(name, cells):
    d, n, reached = doubling(cells)
    ok = LO <= d <= HI
    print(f"{name}: {len(cells)} cells, {n} healthy, {reached} of them reach t_SA2; "
          f"doubling time {d:g} min, against {LO}-{HI} (mean {PAPER_MEAN}): {'holds' if ok else 'FAILS'}")
    return ok


def main(outdir):
    ok = report("syn3a-ensemble default (120 min)", syn3a_cells())
    ens = ensemble_cells(f"{outdir}/ensemble/ensemble_cells.csv", 180)
    ok &= report("KIM-180 default (180 min)", ens)
    cut = [(t if t <= 120 else math.inf, end if end <= 120 else math.inf) for t, end in ens]
    report("KIM-180 default, cut at 120 min", cut)
    print("check:", "holds" if ok else "FAILS")


if __name__ == "__main__":
    main(sys.argv[1])
