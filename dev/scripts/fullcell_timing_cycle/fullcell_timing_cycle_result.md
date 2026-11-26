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
