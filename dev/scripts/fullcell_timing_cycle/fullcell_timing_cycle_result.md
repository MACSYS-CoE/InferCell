# KIM-183: one SpeedyCell cell cycle on a Milan core

Linear: [KIM-183](https://linear.app/kimpson/issue/KIM-183/full-cell-timing-1-one-speedycell-cell-cycle-on-a-milan-core).
Threshold, frozen in the Decision log entry of 2026-11-26 that closed KIM-176: the forward
simulator is usable if one full cell cycle takes **at most 120 s on one core**.

## 1. Which copy of the port was timed

There are two copies of the SpeedyCell port, and they are not the same code.

| | standalone `SpeedyCell.jl` | `syn3a-ensemble/code/SpeedyCell` |
|---|---|---|
| Tree | `/fred/oz022/tkimpson/SpeedyCell.jl` | `/fred/oz022/tkimpson/syn3a-ensemble/code/SpeedyCell` |
| Commit | `0823b53144e2` (HEAD, 2026-11-09) | repo HEAD `2075809cbb82` (2026-11-18); `code/SpeedyCell` last changed at `db40e6d` (2026-11-15), `src/` at `72fedc8` |
| SSA | vendored as `src/SSA/` (6 files, 470 lines) | the `LowPrecSSA` package, pinned from `tomkimpson/wcm-rdme-lowprec` at `cb0abb6` |
| Run entry point | `run_cell` in `src/Cell/run.jl` | `run_cell` in `src/Cell/cell.jl` (the same file; the standalone renamed it for macOS's case-insensitive file system, `6c69c31`) |
| Parameterised path | none | `src/{Cell,CME,ODE}/params.jl`: named kinetic constants with a multiplier each (`params=` in `run_cell`) |
| Early stop | none | `stop=` callback and `ode_failure=:record` (`562396e`, 2026-11-15) |
| Manifest Julia | 1.13.3 | 1.11.2 |
| Model tables | `data/model/` in its own tree | `syn3a-ensemble/data/model/` |

`diff -r` of the two `src/` trees is 164 lines (90 changed lines). Seven files differ, and six
entries exist on one side only (`diff -rq`):
- **On one side only:** `Cell/params.jl`, `CME/params.jl`, `ODE/params.jl` and `Cell/cell.jl`
  (syn3a-ensemble); `Cell/run.jl` and the directory `SSA/` (standalone). `cell.jl` is `run.jl`
  renamed and then extended. `diff run.jl cell.jl` is 173 lines, 127 of them changed.
- **`Cell/cell.jl` against `Cell/run.jl`:**
  - `LowPrecSSA.Model` in place of `SSA.Model`.
  - The `params`, `stop` and `ode_failure` keywords, and the `stopped`/`failure` fields of
    `CellRun`, which the `fnv1a` output digest now mixes in.
  - A `cell_rules(::RateRules)` method, so the parameterised path can build its own rules.
  - `MIN_HEALTHY_HOOKS` and `collapse_stop`.
  - Docstrings.
  - **The timed path calls `ode_integrate` and checks the retcode** (`ode_failure=:record`)
    instead of calling `ode_solve`. `ode_solve` is that same call and check plus a throw
    (`hook.jl`), so the computation and its cost are the same.
- **`CME/segment.jl`** (22 lines each way): `LowPrecSSA.` in place of `SSA.` and comment
  wording; `MODEL_DIR` points at the syn3a-ensemble repo's `data/model/`.
- **`ODE/ODE.jl`** (5 out, 7 in): `rhs_body(exprs=FLUX_EXPRS)`, so the parameterised path can
  generate its own right-hand side; `include("params.jl")`; a comment.
- **`Cell/Cell.jl`, `Cell/hook.jl`, `CME/CME.jl`, `CME/rates.jl`, `SpeedyCell.jl`** (2 to 5
  lines each): exports and includes for the above, and the module docstring.

On the published path (`params=nothing`), the code each copy runs computes the same thing; it
differs only in the SSA's module name and the bookkeeping above. The tables it reads are byte-identical: every file under `data/model/cme/`
and `data/model/ode/`, and `data/model/nodes.csv`, has the same sha256 in both trees (the
CME `reactions.csv` is `63fbeca1…`). Only the docs (`README.md`, `inventory.md`,
`export_meta.json`) differ, and `data/model/params/` exists in syn3a-ensemble only.

**The syn3a-ensemble copy was timed.** KIM-55 and KIM-56 build on it, it is the likely
forward engine for the paper, and it already has the override hook (`params=`) and the early
stop that a prior predictive needs. It has a full-cycle entry point (`run_cell` with
`minutes=121`, as its own `scripts/run_cycle.jl` uses), so no switch to the standalone was
needed. Neither tree was written to. This directory's environment develops the
syn3a-ensemble copy by path (section 2).

## 2. How SpeedyCell runs on the cluster

SpeedyCell needs Julia 1.11 (`julia = "1.11"` in its `Project.toml`); InferCell's
julia/1.10.5 module cannot load it. syn3a-ensemble itself runs it with the Julia 1.11.2
tarball at `/fred/oz022/tkimpson/julia-1.11.2` and the depot `/fred/oz022/tkimpson/julia_depot`
(`scripts/reference/env.sh` there). Here it ran under **Julia 1.11.9**
(`/fred/oz022/tkimpson/julia-1.11.9`), as `versioninfo()` records in every Slurm log of the
runs below:

```
Julia Version 1.11.9
Commit 53a02c0720c (2026-02-06 00:27 UTC)
Platform Info:
  OS: Linux (x86_64-linux-gnu)
  CPU: 64 × AMD EPYC 7543 32-Core Processor
  LLVM: libLLVM-16.0.6 (ORCJIT, znver3)
```

The environment is this directory's `Project.toml` and `Manifest.toml`:
- It develops `/fred/oz022/tkimpson/syn3a-ensemble/code/SpeedyCell` by path. A path is needed:
  `MODEL_DIR` is resolved relative to the source tree, so the tables come from
  `syn3a-ensemble/data/model/`.
- Its Manifest started as a copy of `code/SpeedyCell/Manifest.toml`, and `Pkg.develop` and
  `Pkg.add` ran with `PRESERVE_ALL`. Against the package's own Manifest, only `julia_version`
  (1.11.2 to 1.11.9) and the stdlib `OpenLibm_jll` (0.8.1+2 to 0.8.5+0, bundled with Julia)
  changed. Every SciML, Sundials and LowPrecSSA version is the one syn3a-ensemble pins.
- `LowPrecSSA` is unregistered and pinned over SSH (`tomkimpson/wcm-rdme-lowprec` at
  `cb0abb6`). The jobs stack a depot of their own on the shared one that already holds it:
  `JULIA_DEPOT_PATH=/fred/oz022/tkimpson/infercell-kim183-2026-11-26/depot:/fred/oz022/tkimpson/julia_depot`.
  The environment was instantiated on the login node with `JULIA_PKG_OFFLINE=true` and no
  precompilation. Precompilation ran in one milan job (21663322, dave23): 230 packages in
  1,193 s on one core. Its last step, `versioninfo()`, failed because `-e` does not load
  InteractiveUtils (fixed in `2789c21`). The cache was complete, so the job was not rerun.

All jobs ran on milan nodes with an AMD EPYC 7543. `timing_cycle.slurm` refuses any other
CPU (lscpu guard). The node list excludes `gina[1-18,301-305]` and `dave[148-149]`, the
exclusion syn3a-ensemble uses for its EPYC 7543 timings; it contains the turin nodes that
`launcher_guard.sh` excludes (dave148, dave149, gina18).

## 3. A full cycle runs, and what marks its end

Each cell ran with `run_cell(model, procid; seed=1, minutes=121, ode_failure=:record)`:
121 recorded minutes (blocks 0 to 120), so 120 simulated restart minutes, syn3a-ensemble's
T_end. A run of any shorter cycle is a prefix of this one, bit for bit, so one run gives the
wall time to minute 105, to t_SA2 and to minute 120. All 20 serial cells and all 100
ensemble cells ran the full 120 minutes. None threw and none had a failed ODE step.

**What marks the end of a cycle.** SpeedyCell has no doubling-time or division readout. It
records counts each minute, and the cycle is read off them. The definitions are those of
syn3a-ensemble's `code/cell_cycle/events.py`:
- **t_SA2, the cycle end:** the first minute with `CellSA` ≥ 2 × its minute-0 value. This is
  the paper's criterion (Thornburg et al. 2022, as `events.py` quotes it). Division is assumed
  to begin when the volume has doubled, at t_Vcap below. The surface area goes on growing
  until two cells have formed, which is when it has doubled.
- **t_Vcap, the volume cap:** the first minute with `CellV` = 670, i.e. a volume of
  6.70e-17 L, twice the initial 335. `calcCellVolume` stops the volume growing there
  (`V_CAP` in `src/Cell/writeback.jl`), while the surface area goes on growing to t_SA2.
- **t_short, collapse:** the first minute with fewer than 59 hooks. The cell has run out of
  PEP and its CME goes quiet ("unhealthy cells" in the paper). A collapsed cell reaches t_SA2
  late or never.

The issue's "105-minute cycle" is about the paper's mean doubling time (97 min, range 88 to
112). In the serial set, 15 of 20 cells reach t_SA2, at minutes 92 to 114 (median 100). The
other 5 collapse first, at minutes 55 to 105. Two cells collapse after t_SA2 (procids 9 and 11,
at 112 and 116). Every cell reaches the volume cap, at minutes 57 to 81.

**These collapses are the model's behaviour, not this harness's or this Julia's.**
- **Bit-identical to syn3a-ensemble's own run.** syn3a-ensemble ran the same cells under
  Julia 1.11.2 (`results/port_minutes_cycle.csv`: seed 1, procids 1–200, minutes 0–120). This
  run's per-minute CellSA, CellV, hooks and CME events match it in every row. That holds for
  all 20 serial cells (2,420 rows) and all 100 ensemble cells (12,100 rows), with 0 mismatches.
  So the move to Julia 1.11.9 and OpenLibm_jll 0.8.5 changed nothing in these four columns.
  The reference's own provenance records its Julia: `"julia": "1.11.2"` in
  syn3a-ensemble's `results/port_observables_cycle.json` (and in the `_seed2` file).
  `bitexact_check.py` in this directory runs the comparison: a row-wise string equality of
  the four columns, keyed on (procid, seed, minute). Its output: 2,420 and 12,100 rows, 0
  mismatches, 0 missing. It covers seed 1 and those four columns only.
- **The port against the reference** is settled in syn3a-ensemble's README (C4-cycle): 13
  against 19 of 99 cells collapse by minute 120 (Fisher p 0.33).
- **Procids 1–20 at seed 1 are a high-collapse draw.** Counting a cell as collapsed if any
  minute from 1 to 120 has fewer than 59 hooks (`bitexact_check.py` prints these counts from
  syn3a-ensemble's `results/port_minutes_cycle*.csv`):

  | | procids 1–20 | procids 1–200 |
  |---|---|---|
  | seed 1 | 7 | 33 (16.5%) |
  | seed 2 (`port_minutes_cycle_seed2.csv`) | 0 | 30 |

  So the 20-cell picture is not typical of a batch. Collapsed cells are cheap, so the 15-cell
  rows of section 4 are closer to a typical batch: 50.6 s to minute 105, and 46.8 s to t_SA2,
  which is the exact cost of a cycle.

**One full-cycle trace: procid 1, seed 1** (Slurm 21664418, dave23). It is a healthy cell.
`CellV` hits the cap at minute 62, and `CellSA` doubles at minute 102. The full per-minute
record of all 20 cells is `results/serial_minutes.csv`.

| minute | CellSA (nm²) | CellSA / initial | CellV | hooks | wall (s) |
|---|---|---|---|---|---|
| 0 | 502831 | 1.000 | 335 | 0 | 0.0000 |
| 15 | 568491 | 1.131 | 403 | 59 | 0.3940 |
| 30 | 622320 | 1.238 | 462 | 59 | 0.4247 |
| 45 | 700079 | 1.392 | 551 | 59 | 0.5681 |
| 60 | 785446 | 1.562 | 654 | 59 | 0.5507 |
| 62 (t_Vcap) | 797556 | 1.586 | 670 | 59 | 0.5833 |
| 75 | 862587 | 1.715 | 670 | 59 | 0.5859 |
| 90 | 943919 | 1.877 | 670 | 59 | 0.5124 |
| 102 (t_SA2) | 1011936 | 2.012 | 670 | 59 | 0.6272 |
| 105 | 1032576 | 2.054 | 670 | 59 | 0.6220 |
| 120 | 1142235 | 2.272 | 670 | 59 | 0.7476 |

"wall" is that simulated minute's wall time. It grows from 0.34 s in minute 1 (0.37 s on average
over minutes 1–10) to 0.75 s at minute 120, as the cell and its event rate grow. So a cycle's cost is not 105 × the
20-minute rate that the syn3a-ensemble README quotes (0.425 s per simulated minute).

## 4. Twenty full cycles on one core

**Job:** Slurm **21664418**, at commit `917bcb2`, on dave23 (AMD EPYC 7543 32-Core
Processor), with 1 CPU, 1 Julia thread and 1 BLAS thread. One process ran the cells one after
another:
- `CellModel()` took 1.53 s.
- A warm-up cell followed (procid 1000, 3 recorded minutes, with the timed path's keywords),
  so that compilation fell before the clock: 16.58 s.
- Then procids 1 to 20 ran, seed 1, each to minute 120.

The wall time to minute k is the run's setup (about 0.01 s) plus its first k minutes, from
`run_cell`'s own per-minute timers. The 120-minute wall is `run_cell`'s own total
(`timings.total`). The outer wall clock around `run_cell` agrees with it to the printed 0.001 s here,
and to 0.034 s in the ensemble. dave23 was shared with other jobs.

| per cell, one core | cells | median | range |
|---|---|---|---|
| **to simulated minute 105, all cells (the verdict number, named before the run)** | 20 | **49.5 s** | 28.3–52.9 s |
| to t_SA2 (the exact cost of a full cycle), cells that reach it | 15 | 46.8 s | 40.8–51.0 s |
| to minute 105, cells that go on to reach t_SA2 (added after the run was read) | 15 | 50.6 s | 43.1–52.9 s |
| to minute 105, cells with t_SA2 ≤ 105 | 12 | 51.7 s | 46.4–52.9 s |
| to minute 120 (the whole run) | 20 | 58.1 s | 28.4–62.3 s |
| to minute 120, cells that reach t_SA2 | 15 | 59.6 s | 49.4–62.3 s |

The five cells under 40 s to minute 105 are the five that never reach t_SA2 (procids 2, 4, 7,
12 and 13, which collapse at minutes 55 to 105). A collapsed cell fires few CME events, so its
later minutes are cheap, and its run is not a full cycle. Not all 15 cells that reach t_SA2
have done so by minute 105: procids 8, 11 and 16 reach it at minutes 111, 108 and 114, so their
time to minute 105 stops before their cycle ends. The bottom of that row's range, 43.1 s, is
procid 16. The time to t_SA2 is the exact cost of a full cycle.

| procid | t_SA2 | t_Vcap | t_short | to 105 min (s) | to t_SA2 (s) | to 120 min (s) |
|---|---|---|---|---|---|---|
| 1 | 102 | 62 | – | 51.7 | 49.8 | 61.8 |
| 2 | – | 78 | 55 | 28.7 | – | 28.7 |
| 3 | 96 | 62 | – | 50.4 | 45.5 | 59.6 |
| 4 | – | 66 | 63 | 33.3 | – | 34.0 |
| 5 | 92 | 57 | – | 48.3 | 40.8 | 57.0 |
| 6 | 100 | 65 | – | 49.0 | 46.0 | 59.2 |
| 7 | – | 73 | 105 | 39.3 | – | 41.4 |
| 8 | 111 | 70 | – | 44.0 | 46.8 | 51.1 |
| 9 | 104 | 66 | 112 | 46.4 | 46.0 | 50.7 |
| 10 | 103 | 68 | – | 52.2 | 51.0 | 61.3 |
| 11 | 108 | 67 | 116 | 49.9 | 50.6 | 52.5 |
| 12 | – | 67 | 90 | 39.8 | – | 41.3 |
| 13 | – | 81 | 57 | 28.3 | – | 28.4 |
| 14 | 97 | 65 | – | 52.9 | 48.0 | 62.1 |
| 15 | 94 | 62 | – | 52.8 | 46.6 | 61.7 |
| 16 | 114 | 69 | – | 43.1 | 46.8 | 49.4 |
| 17 | 102 | 71 | – | 50.6 | 48.7 | 59.6 |
| 18 | 97 | 63 | – | 51.6 | 46.5 | 62.3 |
| 19 | 96 | 60 | – | 52.3 | 46.9 | 60.0 |
| 20 | 99 | 63 | – | 51.9 | 48.7 | 59.3 |

Where the time goes (all 20 cells): the CME event loop takes 59.8% and the ODE integration
38.9%; the rest is about 1%. The median over cells is 0.484 s per simulated minute over the
whole 120-minute run.

### Verdict

**PASS against the 2026-11-26 threshold** (at most 120 s per cell cycle on one core).

- **The verdict number passes.** It was named on the issue before the run: the wall time to
  simulated minute 105 over procids 1–20, with no cell replaced. Its median is **49.5 s**, and
  the slowest cell takes 52.9 s.
- **Every other definition passes too.** A full cycle, the time to t_SA2 for the 15 cells that
  reach it, takes a median of 46.8 s, at most 51.0 s. The whole 120-minute run takes at most
  62.3 s.
- **The completing-cell view was added after the run was read.** It is the time to minute 105
  for the 15 cells that go on to reach t_SA2: median 50.6 s, 43.1 to 52.9 s. It passes as well.
- **Five of the 20 cells are not full cycles.** They collapse first. The set was fixed before
  the run and was not topped up.

This is the time on one core, with nothing else in the process. Threaded, a cell is slower
(section 5).

## 5. A 100-cell threaded ensemble on one node

**Job:** Slurm **21664422**, at commit `917bcb2`, on dave104 (AMD EPYC 7543 32-Core
Processor, 64 cores). It requested `--cpus-per-task=64`, the whole node; the site refuses
`--exclusive`. The configuration was set before the run:
- one process with 64 Julia threads, `--gcthreads=64,1` (64 mark threads and 1 sweep thread,
  so `Threads.ngcthreads()` = 65), and `OPENBLAS_NUM_THREADS=1` (`BLAS.get_num_threads()` = 1);
- one shared `CellModel`, with threads not pinned;
- `CellModel()` (1.92 s) and a warm-up cell (16.95 s) ran before the clock;
- procids 1 to 100, seed 1, to minute 120 under `Threads.@threads :dynamic`, each cell with
  `ode_failure=:record`.

**Result: 100 cells in 188.3 s of wall time, 1,912 cells per node-hour.** All 100 cells ran
the full 120 minutes, with no errors and no ODE failures. 84 reach t_SA2, 3 of which collapse
later, and 16 collapse without reaching it.
- **Each cell is slower threaded, and more so on a full node.** Per cell, the whole run took a
  median of 90.6 s (range 35.4 to 113.0 s), or 80.9 s to minute 105. For procids 1 to 20 that
  is a median of 1.60× the one-core time. That median mixes two loads (see the chunking bullet below):

  | | cells | median per-cell wall | ratio to one core, procids 1–20 |
  |---|---|---|---|
  | wave 1 (all 64 threads busy) | 64 | 98.7 s | 1.71× (n = 10) |
  | wave 2 (36 threads busy) | 36 | 77.6 s | 1.35× (n = 10) |
  | both | 100 | 90.6 s | 1.60× (n = 20) |
- **The cells' output is the same.** All 20 of those procids fire the same number of CME
  events threaded as on one core.
- **GC takes 50.8 s, 27.0% of the wall time.** A collection stops every thread.
- **100 cells do not fill 64 threads twice, and the loop's split is fixed.** On Julia 1.11,
  `Threads.@threads :dynamic` over 100 items on 64 threads makes 64 contiguous chunks. Procids
  (1,2), (3,4) … (71,72) went to 36 threads, and procids 73–100 one each to the other 28. Each
  second cell started the moment its thread's first cell ended. So the 36/28 split is the
  chunking, not load balancing. The last cell started at 113.1 s, and the threads were busy
  71.8% of the node's thread-seconds.
- **The full-node figure is lower than the mixed one.** At the wave-1 median of 98.7 s a cell,
  64 busy threads would give about 2,330 cells per node-hour. syn3a-ensemble's 2.16× slowdown
  at 64 threads (F-B3) would give about 1,840. Both are arithmetic, not measurements; the
  measured figure is the 1,912 above.
- **Not rerun here.** `:greedy`, or one task per cell, would balance a batch that is not a
  multiple of the thread count. KIM-58 and KIM-187 should use one of them, or one worker per
  core (next paragraph).

**Context, not measured here.** At the one-core time of section 4 (58.1 s to minute 120), 64
independent one-core cells would make about 3,970 cells per node-hour. This needs one worker
process per core, and assumes they keep their one-core speed. syn3a-ensemble measured the gap on
20-minute cells (README, Phase 6 lane 6a, F-B3). There, one worker per core gave 1.79× the cells
per node-hour of one 64-thread process. GC took 44% of the wall time at 64 threads, and workers
scaled to 85% of linear at 64 cores. So KIM-58 should be sized on workers, not on the threaded
figure above.

For scale only: KIM-58 runs 10³ draws × 100 cells at each of three widths (the default, halved
and doubled), so 10⁵ cells a width and 3 × 10⁵ in all.

| KIM-58 | at 58.1 s a cell on one core | threaded, 1,912 cells per node-hour | workers, 3,970 cells per node-hour if they keep one-core speed |
|---|---|---|---|
| one width (10⁵ cells) | 1,610 core-h | 52 node-h (3,350 CPU-h) | 25 node-h (1,610 CPU-h) |
| three widths (3 × 10⁵) | 4,840 core-h | 157 node-h (10,040 CPU-h) | 76 node-h (4,840 CPU-h) |

At one-core speed that is about 4,840 CPU-h for the three widths. The Decision log entry of
2026-11-26 brackets the prior predictive at 5,000 to 6,600 CPU-h for three widths, at 60 to
80 s a cell; this estimate sits just under it. Two caveats:
- **The `params=` path was not timed.** All cells here ran the published path.
- **Draws will change the cost.** Parameter draws change the event rate and the collapse rate,
  and so the cost per cell.

The per-minute record of the 100 cells stays in the log directory, outside the repo (445 KB,
not committed; its `wall_s` column cannot be regenerated):
`/fred/oz022/tkimpson/infercell-kim183-2026-11-26/ensemble/ensemble_minutes.csv`, sha256
`e6c416ce00217a3f4ba6e5de35438999079012be512dcdda71f125fac06ffdb5`.

## Reproduce

Run these commands from the worktree root, at a pushed commit.
- **Pick a new output directory `L`.** The runs of record are in
  `/fred/oz022/tkimpson/infercell-kim183-2026-11-26/`. A rerun there would overwrite
  `ensemble/ensemble_minutes.csv`, whose sha256 is quoted above.
- **Set the depot.** `timing_cycle.slurm` uses this issue's depot,
  `/fred/oz022/tkimpson/infercell-kim183-2026-11-26/depot`, stacked on the shared one. It is
  only read and extended, never overwritten. To use another depot, edit the script's
  `JULIA_DEPOT_PATH`.
- **Set the environment up once, on the login node.** The committed Project and Manifest
  give it as it ran.

```bash
L=/fred/oz022/tkimpson/infercell-kim183-rerun-$(date +%F)    # new; never the directory of record
mkdir -p $L/logs
export JULIA_DEPOT_PATH=/fred/oz022/tkimpson/infercell-kim183-2026-11-26/depot:/fred/oz022/tkimpson/julia_depot
export JULIA_PKG_OFFLINE=true JULIA_PKG_PRECOMPILE_AUTO=0
/fred/oz022/tkimpson/julia-1.11.9/bin/julia --project=dev/scripts/fullcell_timing_cycle -e 'using Pkg; Pkg.instantiate()'
X="gina[1-18,301-305],dave148,dave149"
S=dev/scripts/fullcell_timing_cycle/timing_cycle.slurm
sbatch --partition=milan --exclude="$X" --cpus-per-task=1 --mem=8G --time=01:00:00 \
  --output=$L/logs/%x_%j.log --export=ALL,MODE=precompile $S
sbatch --partition=milan --exclude="$X" --cpus-per-task=1 --mem=8G --time=01:30:00 --job-name=kim183-serial \
  --output=$L/logs/%x_%j.log --export=ALL,MODE=serial,PROCIDS=1-20,MINUTES=121,OUT=$L/serial/serial $S
sbatch --partition=milan --exclude="$X" --nodes=1 --cpus-per-task=64 --mem=200G --time=01:00:00 --job-name=kim183-ensemble \
  --output=$L/logs/%x_%j.log --export=ALL,MODE=ensemble,PROCIDS=1-100,MINUTES=121,OUT=$L/ensemble/ensemble $S
/fred/oz022/tkimpson/miniconda3/bin/python3 -I dev/scripts/fullcell_timing_cycle/summarise.py \
  $L/serial/serial_cells.csv $L/ensemble/ensemble_cells.csv $L/serial/serial_minutes.csv 1
/fred/oz022/tkimpson/miniconda3/bin/python3 -I dev/scripts/fullcell_timing_cycle/bitexact_check.py \
  /fred/oz022/tkimpson/syn3a-ensemble/results $L/serial/serial_minutes.csv $L/ensemble/ensemble_minutes.csv
```

Each cell's output depends only on (procid, seed, minutes), and is bit-exact only on an
EPYC 7543, as syn3a-ensemble's README notes. The wall times are not reproducible to the digit.
They vary with the node and with what else runs on it.

`results/` holds the run tables: `serial_cells.csv`, `serial_minutes.csv`,
`ensemble_cells.csv`, and each run's `*_meta.txt`.
- **The meta files** record the host, CPU, thread counts, both git commits, `versioninfo()`
  and the sha256 of every file under `syn3a-ensemble/data/model/`.
- **What `summarise.py` prints.** Run on these files, it prints:
  - the trace table (section 3) and the tables of sections 4 and 5;
  - the meta-file lines quoted here (model build, warm-up, wall, GC);
  - the waves and the thread split;
  - the arithmetic it lists: about 2,330 cells per node-hour from the wave-1 median, 3,970 for
    workers, and the KIM-58 core-h, node-h and threaded CPU-h. The workers' CPU-h column equals
    the core-h line.
- **What comes from elsewhere:**
  - The syn3a-ensemble figures (1.79×, 2.16×, 44%, 85%, C4-cycle) come from its README.
  - The bit-exact match and the seed/procid collapse counts (7, 33, 0, 30) come from
    `bitexact_check.py`, run on syn3a-ensemble's `results/port_minutes_cycle*.csv`
    (section 3).
  - The 1,840 cells per node-hour is arithmetic on F-B3's 2.16×: 64 × 3600 / (2.16 × 58.1 s).
  - The section 1 comparison comes from git and `diff` of the two trees, and the sha256
    table from the meta files.
- **Failed cells.** `summarise.py` takes rows with an error or a failed ODE step: it counts
  them and leaves them out of the statistics they have no value for.

**The committed `timing_cycle.jl` is the one that ran.** It is unchanged since `917bcb2`; the
`.slurm` changed only in comments. Two latent limits are left as they are, since neither fired
here:
- **A failed ODE step skews that cell's walls.** `run_cell` stops before it records that
  minute's wall time, so the partial minute lands in `setup_s` and every `wall_to` value of
  that cell. There were 0 ODE failures.
- **The warm-up's result is not checked.** If it threw, compilation would fall in the first
  timed cell without warning. Here the warm-up took 16.6 s, and procid 1's setup was 0.009 s.

The tables the published path reads are:

| table | sha256 |
|---|---|
| `nodes.csv` | `1bd4b35e1c33f3547fd5683213b5b4c32d7211dd9046d8a4b7ee628ee1f39b89` |
| `cme/initial_counts.csv` | `1149da1f0ac53d53d1eb8d63c302e4d638e5e15666bddbb66cfd1f08276c67b5` |
| `cme/rate_rules.csv` | `a95619a05e47b32e7e97e8f65128e31a25403c386e80b5ccfb298d905aa3af0e` |
| `cme/reactions.csv` | `63fbeca180fe8bdebec646df0ddbffb0149958fc2ff9fc3e6ec8803cd24c6aaf` |
| `cme/species.csv` | `c1acb80878e8a500a86e0e1d9e195a491b46a01c337788c0a6554df19cbe6ec9` |
| `ode/parameters.csv` | `60e2ac5ae374daa4496268a0c01ff000043cea02da131355dcd1144be6e9a58c` |
| `ode/reactions.csv` | `4da4941402a0e176ad6adad2c35b6ad0f67f992c023a306ea43f5edb3b8b9a27` |
| `ode/species.csv` | `77d003243ec52770287ba76658f4deeb00cb3e8ad2e51172d0d3a7517b02bd55` |
| `ode/stoichiometry.csv` | `be4293c536cac2ef1a73cabf0144e5a64b696a7bff5f07ab15c871c8a710bf2f` |

These were read at syn3a-ensemble commit `2075809cbb825140cb31a7bec8bc3529770c8e13`, with
`code/SpeedyCell` and `data/model` clean.

**Compute used:** 5.15 CPU-h (sacct). That is 21664418 (0.30), 21664422 (64 × 4 min, 4.27)
and the precompile and smoke jobs (0.59).
