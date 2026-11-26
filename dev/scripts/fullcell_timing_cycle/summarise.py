"""KIM-183: the numbers in fullcell_timing_cycle_result.md, from timing_cycle.jl's tables.

    python3 -I summarise.py serial_cells.csv ensemble_cells.csv serial_minutes.csv TRACE_PROCID

Standard library only. Prints the serial set's walls (all cells, and the cells that
reach t_SA2), the ensemble's per-cell walls and thread use, and one cell's per-minute
CellSA and CellV every 15 minutes, as Markdown rows.
"""
import csv
import statistics
import sys

THRESHOLD_S = 120.0


def rows(path):
    with open(path) as f:
        return list(csv.DictReader(f))


def num(r, k):
    return float(r[k]) if r[k] not in ("", None) else None


def stats(xs):
    xs = [x for x in xs if x is not None]
    if not xs:
        return "n = 0"
    return f"n = {len(xs)}: median {statistics.median(xs):.1f} s, range {min(xs):.1f}-{max(xs):.1f} s"


def main(serial, ensemble, minutes, trace_procid):
    s = rows(serial)
    print(f"## serial ({len(s)} cells)")
    print("status:", {k: sum(r["status"] == k for r in s) for k in sorted({r["status"] for r in s})})
    healthy = [r for r in s if int(r["t_SA2"]) >= 0]
    collapsed = [r for r in s if int(r["t_short"]) >= 0]
    print(f"reach t_SA2: {len(healthy)}; collapse (t_short >= 0): {len(collapsed)}")
    for k in ("wall_to_105_s", "wall_to_SA2_s", "wall_to_120_s", "wall_outer_s"):
        print(f"{k}, all: {stats([num(r, k) for r in s])}")
        print(f"{k}, reach t_SA2: {stats([num(r, k) for r in healthy])}")
    over = [r["procid"] for r in s if (num(r, "wall_to_105_s") or 0) > THRESHOLD_S]
    print(f"cells over {THRESHOLD_S:.0f} s to minute 105: {over or 'none'}")
    print("t_SA2 of the cells that reach it:", sorted(int(r["t_SA2"]) for r in healthy))
    print("t_Vcap:", sorted(int(r["t_Vcap"]) for r in s if int(r["t_Vcap"]) >= 0))
    print("t_short:", sorted(int(r["t_short"]) for r in collapsed))
    per_min = [num(r, "wall_to_120_s") / int(r["minutes_run"]) for r in s if r["minutes_run"]]
    print(f"s per simulated minute, all: median {statistics.median(per_min):.3f}")
    comp = ("cme_s", "ode_s", "ode_setup_s", "hook_other_s", "rebuild_s", "io_s")
    tot = sum(num(r, "wall_total_s") for r in s)
    print("component shares:", {c: round(sum(num(r, c) for r in s) / tot, 3) for c in comp})
    print("| procid | status | t_SA2 | t_Vcap | t_short | to 105 min (s) | to t_SA2 (s) | to 120 min (s) |")
    print("|---|---|---|---|---|---|---|---|")
    dash = lambda v: "-" if v in ("-1", "") else v
    for r in s:
        print(f"| {r['procid']} | {r['status']} | {dash(r['t_SA2'])} | {dash(r['t_Vcap'])} | {dash(r['t_short'])} | "
              f"{dash(r['wall_to_105_s'])} | {dash(r['wall_to_SA2_s'])} | {dash(r['wall_to_120_s'])} |")

    e = rows(ensemble)
    print(f"\n## ensemble ({len(e)} cells)")
    print("status:", {k: sum(r["status"] == k for r in e) for k in sorted({r["status"] for r in e})})
    wall = max(num(r, "t_end") for r in e)
    busy = sum(num(r, "t_end") - num(r, "t_start") for r in e)
    nthreads = len({r["thread"] for r in e})
    print(f"wall (first start to last end) {wall:.1f} s; threads used {nthreads}; cells per node-hour "
          f"{len(e) * 3600 / wall:.0f}")
    print(f"thread-seconds busy {busy:.0f} of {64 * wall:.0f} (64 threads): {100 * busy / (64 * wall):.1f}%")
    print(f"per-cell wall_outer, all: {stats([num(r, 'wall_outer_s') for r in e])}")
    print(f"per-cell to 105 min, all: {stats([num(r, 'wall_to_105_s') for r in e])}")
    he = [r for r in e if int(r["t_SA2"]) >= 0]
    print(f"reach t_SA2: {len(he)}; collapse: {sum(int(r['t_short']) >= 0 for r in e)}")
    print(f"per-cell to t_SA2, reach t_SA2: {stats([num(r, 'wall_to_SA2_s') for r in he])}")
    starts = sorted(num(r, "t_start") for r in e)
    print(f"last cell started at {starts[-1]:.1f} s")
    # The same procids, serial against threaded (each cell's output depends only on procid and seed).
    se = {r["procid"]: r for r in s}
    same = [r["procid"] for r in e if r["procid"] in se and r["nevents"] == se[r["procid"]]["nevents"]]
    print(f"procids with the same event count serial and threaded: {len(same)} of {len(se)}")
    ratio = [num(r, "wall_to_120_s") / num(se[r["procid"]], "wall_to_120_s") for r in e if r["procid"] in se]
    print(f"threaded / serial per-cell wall (same procids): median {statistics.median(ratio):.2f}")

    print(f"\n## trace, procid {trace_procid}")
    m = [r for r in rows(minutes) if r["procid"] == trace_procid]
    sa0 = int(m[0]["CellSA"])
    print("| minute | CellSA (nm²) | CellSA / initial | CellV | hooks | wall (s) |")
    print("|---|---|---|---|---|---|")
    for r in m:
        k = int(r["minute"])
        if k % 15 == 0 or k == len(m) - 1 or r is next((x for x in m if int(x["CellSA"]) >= 2 * sa0), None) \
                or r is next((x for x in m if int(x["CellV"]) >= 670), None):
            print(f"| {k} | {r['CellSA']} | {int(r['CellSA']) / sa0:.3f} | {r['CellV']} | {r['hooks']} | {r['wall_s']} |")


if __name__ == "__main__":
    main(*sys.argv[1:5])
