# 16a.9b, V6: the M0 reference chain converges at 1,650 sweeps per chain

Sized by §12 2026-10-08. Runs at `281d548`, from `.worktrees/v6-16a9b`, with
the code of the pilot (`6f8bcbc`) apart from the test extras' compat.

**The jobs.**
- **Stage 2:** chains 18276511 to 18276514 extended stage 1's checkpoints from
  700 to 1,552 sweeps, at 20 threads, in 48.6 to 51.8 h each. Merge 18276515
  (`merge 100`) missed V6 on one value: MCSE(q05)/SD on the GAPD promoter,
  `S_JCVISYN3A_0607`, was 0.051.
- **The top-up:** chains 18420986 to 18420989 continued to **1,650 sweeps**, in
  5.4 to 6.9 h each. Merge 18420990 (`merge 100`) is the table below.
- **Archives:** the sweep-1,552 checkpoints are in `stage2_1552/` and the
  pilot's sweep-200 starts in `pilot_200/`, both untracked under the V6
  worktree's `dev/scripts/gibbs_pilot_16a9a/`.

The merge drops 100 sweeps per chain and keeps 6,200 draws. The script's header
says "Pilot of 16a.9a" because it is the pilot's script.

| Parameter (log) | R̂ (rank) | Bulk ESS | Tail ESS | Bulk ESS per sweep | Tail ESS per sweep | MCSE(q05)/SD | MCSE(q95)/SD | Posterior mean ± SD | Truth | Truth's quantile |
|---|---|---|---|---|---|---|---|---|---|---|
| S_JCVISYN3A_0607 | 1.007 | 1081 | 2323 | 0.174 | 0.375 | 0.049 | 0.044 | 1.545 ± 0.116 | 1.489 | 0.31 |
| S_JCVISYN3A_0779 | 1.001 | 3441 | 4056 | 0.555 | 0.654 | 0.030 | 0.022 | 0.700 ± 0.199 | 1.012 | 0.95 |
| krnadeg | 1.006 | 1088 | 2318 | 0.176 | 0.374 | 0.043 | 0.049 | 2.225 ± 0.076 | 2.243 | 0.60 |
| kcatF_R_ENO | 1.002 | 1531 | 2470 | 0.247 | 0.398 | 0.039 | 0.035 | 4.242 ± 0.042 | 4.182 | 0.07 |
| kcatF_R_FBA | 1.001 | 5750 | 4652 | 0.927 | 0.750 | 0.035 | 0.040 | 4.408 ± 0.184 | 4.207 | 0.13 |
| sigma_b | 1.000 | 4248 | 4520 | 0.685 | 0.729 | 0.030 | 0.030 | -2.625 ± 0.060 | -2.622 | 0.52 |

## Against V6: met

| Criterion | Limit | Worst | On |
|---|---|---|---|
| Rank-normalised R̂ | < 1.01 | 1.007 | the GAPD promoter, `S_JCVISYN3A_0607` |
| Bulk ESS | ≥ 1,000 | 1,081 | the GAPD promoter |
| Tail ESS | ≥ 1,000 | 2,318 | `krnadeg` |
| MCSE of the 5% and 95% quantiles | ≤ 0.05 SD | 0.049 | the GAPD promoter's 5% quantile and `krnadeg`'s 95% |

**The quantile MCSE passes by 0.001, and is the binding criterion.** It was
0.051 at 1,552 sweeps. The merge's own projection, from the slowest ESS rate, asks
for 1,610 sweeps per chain, and the chains ran 1,650.

**The slowest parameter is the GAPD promoter, not ptsG.** The handoffs of
2026-10-08 and 2026-10-11 called `S_JCVISYN3A_0607` the ptsG promoter. It is GAPD
(`M0_GENES`, `src/organisms/coreA/m0.jl`). The ptsG promoter, `S_JCVISYN3A_0779`,
mixes three times faster.

**Against stage 1.** The posterior means agree with stage 1's burn-in-100 merge
to 0.04 SD or better. The ESS rates per sweep hold: the GAPD promoter's is 0.174,
against 0.180.

**Update rates:** over 414 live cell-windows the lowest is 0.37, and none is
below 0.10. The 86 empty cell-windows average 0.15.

**Cost:** 209 s per sweep at 20 threads (B3 133, B1 77), or 1.16 CPU-h.

**The truth's quantiles are not a calibration check.** They come from one
dataset. V7 is the calibration check.

## Compute

From Slurm's CPUTimeRAW (`phase16_ledger.md`):
- stage 2, four chains and the merge: 4,038.7 CPU-h;
- the top-up, four chains and the merge: 462.6 CPU-h;
- job 18421022, which regenerated stage 1's table: 0.07 CPU-h.

**V6 in all, with stage 1 (2,308.3 and 0.07):** about 6,809.7 CPU-h.

## Reproduce

From a worktree at `281d548`. A chain resumes bitwise from its checkpoint, so each
stage is the same command with a larger `SWEEPS`. Repeat it until all four chains
reach that length. Each `.sh` call also queues a merge with no burn-in, which is
not the merge of record, so cancel it.

```
SWEEPS=200  TIME=2-00:00:00 bash dev/scripts/gibbs_pilot_16a9a.sh   # 16a.9a pilot
SWEEPS=700  TIME=6-00:00:00 bash dev/scripts/gibbs_pilot_16a9a.sh   # stage 1
SWEEPS=1552 TIME=6-00:00:00 bash dev/scripts/gibbs_pilot_16a9a.sh   # stage 2
SWEEPS=1650 TIME=14:00:00   bash dev/scripts/gibbs_pilot_16a9a.sh   # top-up
sbatch --time=00:45:00 --cpus-per-task=1 --mem=16G -J v6_merge100 \
  dev/scripts/gibbs_pilot_16a9a.slurm merge 100                      # of record
```

The chains run at 20 threads (the `.slurm` file's `--cpus-per-task`). The
checkpoints at 1,650 sweeps, the archives `stage2_1552/` and `pilot_200/`, and
every chain and merge log are archived, untracked, at
`/fred/oz022/tkimpson/infercell-16a9b-archive-2026-10-17/v6/gibbs_pilot_16a9a/`.
