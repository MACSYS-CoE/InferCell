"""KIM-183: this run's per-minute record against syn3a-ensemble's own (Julia 1.11.2) run, and the
collapse counts of note section 3.

    python3 -I bitexact_check.py SYN3A_RESULTS MINUTES_CSV...

SYN3A_RESULTS is syn3a-ensemble's results/ (port_minutes_cycle.csv is seed 1 and
port_minutes_cycle_seed2.csv is seed 2, procids 1-200, minutes 0-120). Each MINUTES_CSV is a
timing_cycle.jl *_minutes.csv. For every (procid, seed, minute) row, CellSA, CellV, hooks and
events must be equal as strings. A cell counts as collapsed if any minute from 1 on has fewer
than 59 hooks. Standard library only; exits 1 on any mismatch or missing row.
"""
import csv
import sys

COLS = ("CellSA", "CellV", "hooks", "events")


def load(path):
    with open(path) as f:
        return {(int(r["procid"]), int(r["seed"]), int(r["minute"])): tuple(r[c] for c in COLS)
                for r in csv.DictReader(f)}


def collapsed(d, procids):
    return sorted({p for (p, _, k), v in d.items() if p in procids and k >= 1 and int(v[2]) < 59})


def main(results, *minutes):
    ref = {1: load(f"{results}/port_minutes_cycle.csv"), 2: load(f"{results}/port_minutes_cycle_seed2.csv")}
    bad = 0
    for path in minutes:
        m = load(path)
        miss = sum(ref.get(s, {}).get((p, s, k)) is None for (p, s, k) in m)
        mis = sum(ref.get(s, {}).get((p, s, k)) not in (None, v) for (p, s, k), v in m.items())
        print(f"{path}: {len(m)} rows, {mis} mismatches, {miss} missing from the reference")
        bad += mis + miss
    for seed, d in ref.items():
        c20, c200 = collapsed(d, range(1, 21)), collapsed(d, range(1, 201))
        print(f"seed {seed}: collapse by minute 120 in procids 1-20: {len(c20)} {c20}; in 1-200: {len(c200)}")
    sys.exit(1 if bad else 0)


if __name__ == "__main__":
    main(*sys.argv[1:])
