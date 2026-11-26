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
