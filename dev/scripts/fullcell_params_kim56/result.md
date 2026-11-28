# KIM-56: the parameter-override hook in the syn3a-ensemble port

**Verdict.** With every multiplier at 1, an override run reproduces the default run bit for bit. The check covered 6 cells over 181 simulated minutes, each run three ways, and found 0 mismatches (Slurm 21853120). The hook exposes 838 named kinetic constants (Slurm 21853119).

All runs were at InferCell c4f62ee, from the detached worktree `.worktrees/kim56-run-c4f62ee`. syn3a-ensemble was at 2075809, with `code/SpeedyCell` and `data/model` clean. The node was dave2 (AMD EPYC 7543). The archive is `/fred/oz022/tkimpson/infercell-kim56-2026-11-28/`, with `SHA256SUMS`, `JOBS` and `COMMIT`.

## 1. The copy, and how it runs

**The copy** is syn3a-ensemble's, `code/SpeedyCell`: src at 72fedc8, repo HEAD 2075809. It already has the hook. Lane 10a (c141e78) added three things there:
- `Params(model; multipliers)`;
- `param_names()` and `param_defaults()`;
- `run_cell(...; params=P)`.

The standalone `SpeedyCell.jl` (0823b53) has none of these. KIM-183 timed this same copy (its note, §1). Nothing in syn3a-ensemble is edited here. That follows the Decision log of 2026-11-28, which moved KIM-56's Touches into InferCell.

**How it runs.** The environment in this directory is a copy of KIM-183's (`dev/scripts/fullcell_timing_cycle/`):
- Julia 1.11.9 at `/fred/oz022/tkimpson/julia-1.11.9`;
- SpeedyCell developed by path;
- KIM-183's depot, stacked on the shared one that holds LowPrecSSA.

`params_check.slurm` refuses to run unless syn3a-ensemble is at the pin and clean.

From precompile job 21853118:

```
Julia Version 1.11.9
Commit 53a02c0720c (2026-02-06 00:27 UTC)
syn3a_ensemble_git 2075809cbb825140cb31a7bec8bc3529770c8e13 (pin 2075809, code/SpeedyCell and data/model clean)
```

InferCell's own suite runs on julia/1.10.5, and SpeedyCell needs 1.11. So the Spec's check is this Slurm script, not part of `Pkg.test`.

## 2. The bitwise check (the Spec's full-cell correctness check)

> An override run of the simulator given the default values reproduces the default run bitwise.

`check.jl` runs each procid three ways, with seed 1, `minutes=181` and `ode_failure=:record`. Those are KIM-180's options of record. The three paths are:
- the published path, `params=nothing`;
- lane 10a's parameterised path, `Params(m)`;
- the vector adapter, `params_from_vector(m, ones(838))`.

Each pair must agree under `same_output`, which compares the procid, seed, minutes, event count, minimum count, stop, failure, the per-minute record, the spans and the hooks per minute. Each pair must also agree under `output_digest`.

As a control, a run with `PGI.kcatF` × 2 must give a different digest. That shows the check can fail.

| procid | events | collapsed (hooks < 59) | digest, all three paths |
|---|---|---|---|
| 1 | 408,009,519 | no | c79dce983bab4e3a |
| 2 | 68,855,337 | yes | ce2e7f91d15a1939 |
| 3 | 341,724,320 | no | f2f09c597e3bf64b |
| 4 | 86,496,110 | yes | ed6aaff8779dc954 |
| 5 | 311,918,951 | no | b166c83febb06dc2 |
| 6 | 352,929,209 | no | 9adc39134a320405 |

- "6 cells x 3 paths, 0 mismatches; control differs", then `CHECK PASS` (job 21853120, 31 min).
- The control gave 9efe6c681eb8b8ed for procid 1, against c79dce983bab4e3a.
- Procid 1's published run was repeated for the control and gave the same digest.
- The full record is in `results/check.csv` and `results/meta.txt`.

**What this adds to lane 10a's checks.** syn3a-ensemble's C-D1 found the parameterised path byte-identical over 20 minutes on 100 cells, at 688b470. This run checks three further things:
- the current commit, which includes the early stop (562396e);
- KIM-180's 181-minute window, on cells that complete and cells that collapse;
- the vector entry point.

This run does not cover multipliers other than 1. Lane 10a's C-D2 and C-D3 check those: every one of the 838 constants acts, at 0.5 and 2, to 1.8e-14 against the reference's own builders.

## 3. A parameter vector per run

`adapter.jl`:

```julia
include("dev/scripts/fullcell_params_kim56/adapter.jl")
P = params_from_vector(m, v)   # v[i] multiplies param_names()[i]; length(v) == 838
run_cell(m, procid; seed, minutes=181, params=P, ode_failure=:record)
```

- The interface is multiplicative. `ones(838)` is the published model.
- A length mismatch is a `DimensionMismatch`.
- `Params` rejects a multiplier that is non-finite or negative.
- Pass `names=param_names()` once to avoid rereading the constants for each draw.

## 4. The parameter count and names

`results/param_names.csv` (job 21853119) has the columns `index,name,module,reaction,key,default`, in the adapter's order: 805 ODE constants, then 33 CME constants.

| | count |
|---|---|
| named constants | **838** (838 unique) |
| ODE, as `<reaction>.<key>` | 805, over 174 of 175 reactions |
| CME, by source (`file:variable`, `file:line:col`, or 2 labelled sites) | 33: 22 in `MinCell_restart.py`, 8 in `rep_restart.py` (replication), 3 in `translation_rate_restart.py` (ribosome) |
| defaults that are zero, negative or non-finite | 0 (range 4.5e-11 to 2e+05) |

The ODE constants by key: KmProd1 159, KmSub1 159, kcatF 158, KmProd2 100, KmSub2 89, kcatR 86, KmProd3 16, Const1 11, KmSub3 8, Const2 5, KmProd4 4, Kcat\_R 3, and one each of KcatR, Km1, KmProd5, KmSub4, KmSub5, Kprod and Ksub.

The count of 838 matches KIM-180's frozen eligible set ("all 838 SpeedyCell constants"). It is the row count KIM-57's inventory has to match.

**Not named, so not reachable by a multiplier.** The export's `data/model/params/README.md` and `classification.csv` list these, by class:
- physical constants and geometry;
- unit conversions;
- odecell's `onoff` switch;
- structural numbers: count floors, rate-law exponents, and sequence-derived counts;
- fixed concentrations in the ODE rate laws;
- enzyme levels and copy numbers;
- kinetic-degenerate literals that only rescale a named constant, such as the ribosome loading of 0.25;
- media.

A key is also left unnamed in a reaction where its value is 0, which is why no default is zero. KIM-57 decides whether any of these classes needs a prior. If one does, the hook needs a new named constant in syn3a-ensemble's export, which would be an amendment.

## Reproduce

From a detached worktree at this commit, at the repo root:

```bash
L=/fred/oz022/tkimpson/infercell-kim56-<date>; X="gina[1-18,301-305],dave148,dave149,gina18"
S=dev/scripts/fullcell_params_kim56/params_check.slurm
j=$(sbatch --parsable --partition=milan --exclude="$X" --cpus-per-task=1 --mem=8G --time=01:00:00 \
      --output=$L/logs/%x_%j.log --export=ALL,MODE=precompile $S)
sbatch --dependency=afterok:$j --partition=milan --exclude="$X" --cpus-per-task=1 --mem=8G --time=00:30:00 \
      --output=$L/logs/%x_%j.log --export=ALL,MODE=table,OUT=$L/table $S
sbatch --dependency=afterok:$j --partition=milan --exclude="$X" --cpus-per-task=1 --mem=2G --time=01:00:00 \
      --output=$L/logs/%x_%j.log --export=ALL,MODE=check,PROCIDS=1-6,SEED=1,MINUTES=181,OUT=$L/check $S
```
