"""KIM-183: the numbers in fullcell_timing_cycle_result.md, from timing_cycle.jl's tables.

    python3 -I summarise.py serial_cells.csv ensemble_cells.csv serial_minutes.csv TRACE_PROCID

Standard library only. Prints the serial set's walls (all cells, and the cells that
reach t_SA2), the ensemble's per-cell walls and thread use (by wave), one cell's per-minute
CellSA and CellV every 15 minutes as Markdown rows, the lines of each run's *_meta.txt
that the note quotes (model build, warm-up, wall, GC), and the note's arithmetic on them.
A cell that errored or whose ODE failed has empty fields; it is counted, and left out of
the statistics it has no value for.
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


def ival(r, k):
    """An integer field, -1 where empty (an error row) or absent."""
    return int(r[k]) if r.get(k) not in ("", None) else -1


def meta_lines(cells_path, keys):
    """The lines of the run's *_meta.txt (next to its *_cells.csv) that start with one of keys."""
    path = cells_path.replace("_cells.csv", "_meta.txt")
    try:
        with open(path) as f:
            return [l.rstrip() for l in f if l.startswith(keys)]
    except OSError:
        return [f"(no {path})"]


def stats(xs):
    xs = [x for x in xs if x is not None]
    if not xs:
        return "n = 0"
    return f"n = {len(xs)}: median {statistics.median(xs):.1f} s, range {min(xs):.1f}-{max(xs):.1f} s"


def main(serial, ensemble, minutes, trace_procid):
    s = rows(serial)
    print(f"## serial ({len(s)} cells)")
    print("status:", {k: sum(r["status"] == k for r in s) for k in sorted({r["status"] for r in s})})
    healthy = [r for r in s if ival(r, "t_SA2") >= 0]
    collapsed = [r for r in s if ival(r, "t_short") >= 0]
    for line in meta_lines(serial, ("CellModel", "serial:")):
        print("meta:", line)
    print(f"reach t_SA2: {len(healthy)}; collapse (t_short >= 0): {len(collapsed)}")
    for k in ("wall_to_105_s", "wall_to_SA2_s", "wall_to_120_s", "wall_outer_s"):
        print(f"{k}, all: {stats([num(r, k) for r in s])}")
        print(f"{k}, reach t_SA2: {stats([num(r, k) for r in healthy])}")
    by105 = [r for r in healthy if ival(r, "t_SA2") <= 105]
    print(f"wall_to_105_s, t_SA2 <= 105: {stats([num(r, 'wall_to_105_s') for r in by105])}")
    print("reach t_SA2 after minute 105:", [(r["procid"], r["t_SA2"]) for r in healthy if ival(r, "t_SA2") > 105])
    over = [r["procid"] for r in s if (num(r, "wall_to_105_s") or 0) > THRESHOLD_S]
    print(f"cells over {THRESHOLD_S:.0f} s to minute 105: {over or 'none'}")
    tsa = sorted(ival(r, "t_SA2") for r in healthy)
    print("t_SA2 of the cells that reach it:", tsa, f"(median {statistics.median(tsa) if tsa else '-'})")
    print("t_Vcap:", sorted(ival(r, "t_Vcap") for r in s if ival(r, "t_Vcap") >= 0))
    print("t_short:", sorted(ival(r, "t_short") for r in collapsed))
    print("never reach t_SA2:", [(r["procid"], r["t_short"]) for r in s if ival(r, "t_SA2") < 0])
    print("collapse after t_SA2:", [(r["procid"], r["t_short"]) for r in healthy if ival(r, "t_short") >= 0])
    per_min = [num(r, "wall_to_120_s") / int(r["minutes_run"]) for r in s
               if r["minutes_run"] and num(r, "wall_to_120_s") is not None]
    print(f"s per simulated minute, all: median {statistics.median(per_min):.3f}")
    comp = ("cme_s", "ode_s", "ode_setup_s", "hook_other_s", "rebuild_s", "io_s")
    ok = [r for r in s if num(r, "wall_total_s") is not None]
    tot = sum(num(r, "wall_total_s") for r in ok)
    print("component shares:", {c: round(sum(num(r, c) for r in ok) / tot, 3) for c in comp})
    print("| procid | status | t_SA2 | t_Vcap | t_short | to 105 min (s) | to t_SA2 (s) | to 120 min (s) |")
    print("|---|---|---|---|---|---|---|---|")
    dash = lambda v: "-" if v in ("-1", "") else v
    for r in s:
        print(f"| {r['procid']} | {r['status']} | {dash(r['t_SA2'])} | {dash(r['t_Vcap'])} | {dash(r['t_short'])} | "
              f"{dash(r['wall_to_105_s'])} | {dash(r['wall_to_SA2_s'])} | {dash(r['wall_to_120_s'])} |")

    e = rows(ensemble)
    print(f"\n## ensemble ({len(e)} cells)")
    print("status:", {k: sum(r["status"] == k for r in e) for k in sorted({r["status"] for r in e})})
    meta = meta_lines(ensemble, ("threads", "CellModel", "ensemble:", "cells per node-hour"))
    for line in meta:
        print("meta:", line)
    nthreads = next((int(l.split()[1]) for l in meta if l.startswith("threads ")), None)
    wall = max(num(r, "t_end") for r in e)
    busy = sum(num(r, "t_end") - num(r, "t_start") for r in e)
    used = len({r["thread"] for r in e})
    nthreads = nthreads or used
    print(f"wall (first start to last end) {wall:.1f} s; threads used {used} of {nthreads}; cells per node-hour "
          f"{len(e) * 3600 / wall:.0f}")
    print(f"thread-seconds busy {busy:.0f} of {nthreads * wall:.0f} ({nthreads} threads): "
          f"{100 * busy / (nthreads * wall):.1f}%")
    # Waves: each thread's cells in start order. Wave 1 is every thread's first cell, wave 2 its second.
    per_thread = {}
    for r in sorted(e, key=lambda r: num(r, "t_start")):
        per_thread.setdefault(r["thread"], []).append(r)
    counts = {}
    for cs in per_thread.values():
        counts[len(cs)] = counts.get(len(cs), 0) + 1
    print("threads by cells run:", dict(sorted(counts.items())))
    wave = {r["procid"]: i + 1 for cs in per_thread.values() for i, r in enumerate(cs)}
    print(f"per-cell wall_outer, all: {stats([num(r, 'wall_outer_s') for r in e])}")
    print(f"per-cell to 105 min, all: {stats([num(r, 'wall_to_105_s') for r in e])}")
    he = [r for r in e if int(r["t_SA2"]) >= 0]
    print(f"reach t_SA2: {len(he)}; collapse: {sum(ival(r, 't_short') >= 0 for r in e)}; "
          f"collapse after t_SA2: {sum(ival(r, 't_short') >= 0 for r in he)}; "
          f"collapse without reaching it: {sum(ival(r, 't_short') >= 0 and ival(r, 't_SA2') < 0 for r in e)}")
    print(f"per-cell to t_SA2, reach t_SA2: {stats([num(r, 'wall_to_SA2_s') for r in he])}")
    starts = sorted(num(r, "t_start") for r in e)
    print(f"last cell started at {starts[-1]:.1f} s")
    # The same procids, serial against threaded (each cell's output depends only on procid and seed).
    se = {r["procid"]: r for r in s}
    same = [r["procid"] for r in e if r["procid"] in se and r["nevents"] == se[r["procid"]]["nevents"]]
    print(f"procids with the same event count serial and threaded: {len(same)} of {len(se)}")
    def ratio(rs):
        xs = [num(r, "wall_to_120_s") / num(se[r["procid"]], "wall_to_120_s") for r in rs
              if r["procid"] in se and num(r, "wall_to_120_s") and num(se[r["procid"]], "wall_to_120_s")]
        return f"median {statistics.median(xs):.2f} (n = {len(xs)})" if xs else "n = 0"
    print(f"threaded / serial per-cell wall (same procids): {ratio(e)}")
    for w in sorted(set(wave.values())):
        rw = [r for r in e if wave[r["procid"]] == w]
        print(f"wave {w}: {len(rw)} cells; per-cell wall_outer {stats([num(r, 'wall_outer_s') for r in rw])}; "
              f"threaded / serial {ratio(rw)}")
    w1 = statistics.median(num(r, "wall_outer_s") for r in e if wave[r["procid"]] == 1)
    print(f"arithmetic: {nthreads} threads at the wave-1 median {w1:.1f} s a cell: "
          f"{nthreads * 3600 / w1:.0f} cells per node-hour")

    # The note's arithmetic on the one-core time (to minute 120, all serial cells).
    one = statistics.median(x for x in (num(r, "wall_to_120_s") for r in ok) if x is not None)
    workers = nthreads * 3600 / one
    print(f"\n## arithmetic (one-core median to minute 120: {one:.1f} s)")
    print(f"{nthreads} one-core workers at that speed: {workers:.0f} cells per node-hour")
    for label, n in (("one width", 1e5), ("three widths", 3e5)):
        print(f"KIM-58 {label} ({n:.0e} cells): {n * one / 3600:.0f} core-h at one-core speed; "
              f"{n / (len(e) * 3600 / wall):.0f} node-h ({nthreads * n / (len(e) * 3600 / wall):.0f} CPU-h) threaded; "
              f"{n / workers:.0f} node-h on workers")

    print(f"\n## trace, procid {trace_procid}")
    m = [r for r in rows(minutes) if r["procid"] == trace_procid]
    sa0 = int(m[0]["CellSA"])
    tr = next(r for r in s if r["procid"] == trace_procid)
    tsa2, tvcap = ival(tr, "t_SA2"), ival(tr, "t_Vcap")
    print("| minute | CellSA (nm²) | CellSA / initial | CellV | hooks | wall (s) |")
    print("|---|---|---|---|---|---|")
    for r in m:
        k = int(r["minute"])
        if k % 15 == 0 or k == len(m) - 1 or k in (tsa2, tvcap):
            print(f"| {k} | {r['CellSA']} | {int(r['CellSA']) / sa0:.3f} | {r['CellV']} | {r['hooks']} | {r['wall_s']} |")
    w = [float(r["wall_s"]) for r in m if int(r["minute"]) >= 1]
    print(f"minute wall: minute 1 {w[0]:.3f} s, mean of minutes 1-10 {statistics.mean(w[:10]):.3f} s, "
          f"minute {len(w)} {w[-1]:.3f} s")


if __name__ == "__main__":
    main(*sys.argv[1:5])
