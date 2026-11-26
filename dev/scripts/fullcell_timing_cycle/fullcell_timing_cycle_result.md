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

`diff -r` of the two `src/` trees is 164 lines (90 changed lines). Seven files differ, and five
exist on one side only:
- **On one side only:** `Cell/params.jl`, `CME/params.jl`, `ODE/params.jl` and `Cell/cell.jl`
  (syn3a-ensemble); `Cell/run.jl` and `SSA/` (standalone). `cell.jl` and `run.jl` are the same
  file apart from the lines below.
- **`Cell/cell.jl` against `Cell/run.jl`:** `LowPrecSSA.Model` in place of `SSA.Model`; the
  `params`, `stop` and `ode_failure` keywords and the `stopped`/`failure` fields of `CellRun`;
  docstrings.
- **`CME/segment.jl`** (22 lines each way): `LowPrecSSA.` in place of `SSA.` and comment
  wording; `MODEL_DIR` points at the syn3a-ensemble repo's `data/model/`.
- **`ODE/ODE.jl`** (5 out, 7 in): `rhs_body(exprs=FLUX_EXPRS)`, so the parameterised path can
  generate its own right-hand side; `include("params.jl")`; a comment.
- **`Cell/Cell.jl`, `Cell/hook.jl`, `CME/CME.jl`, `CME/rates.jl`, `SpeedyCell.jl`** (2 to 5
  lines each): exports and includes for the above, and the module docstring.

On the published path (`params=nothing`), the code each copy runs is the same apart from the
SSA's module name. The tables it reads are byte-identical: every file under `data/model/cme/`
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
  the paper's criterion (Thornburg et al. 2022): division begins when the surface area has
  doubled.
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

"wall" is that simulated minute's wall time. It grows from 0.37 s at the start to 0.75 s at
minute 120, as the cell and its event rate grow. So a cycle's cost is not 105 × the
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
`run_cell`'s own per-minute timers. The 120-minute wall is the outer wall clock around
`run_cell`. dave23 was shared with other jobs.

| per cell, one core | cells | median | range |
|---|---|---|---|
| **to simulated minute 105 (the verdict number)** | 20 | **49.5 s** | 28.3–52.9 s |
| to minute 105, cells that reach t_SA2 | 15 | 50.6 s | 43.1–52.9 s |
| to t_SA2 | 15 | 46.8 s | 40.8–51.0 s |
| to minute 120 (the whole run) | 20 | 58.1 s | 28.4–62.3 s |
| to minute 120, cells that reach t_SA2 | 15 | 59.6 s | 49.4–62.3 s |

The cells below 40 s are those that collapse early (procids 2, 4, 13 at minutes 55 to 63): a
collapsed cell fires few CME events, so its later minutes are cheap. It is not a full cycle.
The rows that reach t_SA2 are the cost of a healthy cycle.

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

**PASS against the 2026-11-26 threshold** (at most 120 s per cell cycle on one core). The
median full cycle to minute 105 takes 49.5 s, and the slowest of the 20 takes 52.9 s. The
verdict does not depend on how the cycle is defined. The slowest cell under each definition
is still well under 120 s: 51.0 s to t_SA2 and 62.3 s to minute 120.

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
the full 120 minutes, with no errors and no ODE failures. 84 reach t_SA2 and 19 collapse.
- **Each cell is slower threaded.** Per cell, the whole run took a median of 90.6 s (range
  35.4 to 113.0 s), or 80.9 s to minute 105. For procids 1 to 20 that is a median of 1.60× the
  one-core time.
- **The cells' output is the same.** All 20 of those procids fire the same number of CME
  events threaded as on one core.
- **GC takes 50.8 s, 27.0% of the wall time.** A collection stops every thread.
- **100 cells do not fill 64 threads twice.** 36 threads ran two cells and 28 ran one. The last
  cell started at 113.1 s, and the threads were busy 71.8% of the node's thread-seconds. With a
  batch large enough to keep every thread busy, 64 threads at 90.6 s a cell would give about
  2,540 cells per node-hour. That figure is arithmetic, not a measurement.

**Context, not measured here.** At the one-core time of section 4 (58.1 s to minute 120), 64
independent one-core cells would make about 3,970 cells per node-hour. This needs one worker
process per core, and assumes they keep their one-core speed. syn3a-ensemble measured the gap on
20-minute cells (README, Phase 6 lane 6a, F-B3). There, one worker per core gave 1.79× the cells
per node-hour of one 64-thread process. GC took 44% of the wall time at 64 threads, and workers
scaled to 85% of linear at 64 cores. So KIM-58 should be sized on workers, not on the threaded
figure above.

For scale only: KIM-58's 10³ draws × 100 cells is 10⁵ cells. That is about 1,600 core-hours at
58.1 s a cell, about 52 node-hours (3,350 CPU-h) at the threaded 1,912 cells per node-hour, or
about 25 node-hours at the workers' rate if it holds.

The per-minute record of the 100 cells stays in the log directory, outside the repo:
`/fred/oz022/tkimpson/infercell-kim183-2026-11-26/ensemble/ensemble_minutes.csv`, sha256
`e6c416ce00217a3f4ba6e5de35438999079012be512dcdda71f125fac06ffdb5`.

## Reproduce

From the worktree root, at a pushed commit. Set the environment up once, on the login node
(this is the state committed here):

```bash
L=/fred/oz022/tkimpson/infercell-kim183-2026-11-26
export JULIA_DEPOT_PATH=$L/depot:/fred/oz022/tkimpson/julia_depot JULIA_PKG_OFFLINE=true JULIA_PKG_PRECOMPILE_AUTO=0
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
```

Each cell's output depends only on (procid, seed, minutes), and is bit-exact only on an
EPYC 7543, as syn3a-ensemble's README notes. The wall times are not reproducible to the digit.
They vary with the node and with what else runs on it.

`results/` holds the run tables these numbers come from: `serial_cells.csv`,
`serial_minutes.csv` and `ensemble_cells.csv`, and each run's `*_meta.txt`. The meta files
record the host, CPU, thread counts, both git commits, `versioninfo()`, and the sha256 of
every file under `syn3a-ensemble/data/model/`. The tables the published path reads are:

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

**Compute used:** about 5 CPU-h. That is 21664418 (0.30), 21664422 (64 × 4 min, 4.3) and the
precompile and smoke jobs (0.5).
