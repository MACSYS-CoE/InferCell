# 16a.9b, V7: the M0 sampler passes simulation-based calibration on 50 prior-drawn replicates

Sized by §12 2026-10-08. The pass rule and the complete-set guard are §12
2026-10-12, fixed before any rank was computed.

**The jobs.**
- **Replicates:** array 18276794 at `dfe570a`, from `.worktrees/v7-16a9b`. 34 of
  50 completed. The other 16 hit the 6-day limit on 14–15 Oct and resumed from
  their checkpoints as array 18628508 (replicates 2, 4, 10, 15, 17, 18, 20, 21,
  24, 34, 36, 45 to 49), from the same worktree and commit. The last finished at
  07:48 on 17 Oct.
- **Merge:** job 18689602 at `0ff6cff` (pushed), from a worktree at that commit,
  with `V7_DIR` pointing at the v7 worktree's checkpoints. Its run path is
  unchanged from `dfe570a`, and its merge path is the one §12 2026-10-12 fixed.
- **Each replicate:** one chain from a prior-drawn start, 50 cells, PG at N = 20,
  5 threads, 656 sweeps. The merge drops 100 sweeps and thins the rest to
  L = 100 draws.
- **Cost:** 33.9k CPU-h (33.3k in the first array, 0.63k on resumption), against
  25.6k sized and 35.7k projected on 2026-10-09. Phase 16 is at about 42k of the
  50k.

## Against V7: passes

Each parameter's rank-ECDF is against a simultaneous band at
1 − 0.05/6 = 0.9917 (pointwise γ = 0.00076), so V7 as a whole is at 95%.

| Parameter (log) | Outside the band at | Largest |ECDF − z| | Verdict |
|---|---|---|---|
| S_JCVISYN3A_0607 | nowhere | 0.104 | within |
| S_JCVISYN3A_0779 | nowhere | 0.084 | within |
| krnadeg | nowhere | 0.164 | within |
| kcatF_R_ENO | nowhere | 0.094 | within |
| kcatF_R_FBA | nowhere | 0.075 | within |
| sigma_b | nowhere | 0.094 | within |

**V7 passes on all 50 replicates: 6 of 6 parameters are within their bands.**
No parameter leaves its band anywhere.
- **The largest deviation** is `krnadeg`'s, 0.164 at z ≈ 0.36, where the band
  allows up to 0.244.
- **The nearest approach to a band edge** is at the extremes, where the bands are
  a few counts wide: the band at z ≈ 0.01 allows at most 4 ranks of 0, and
  `krnadeg` has 2. The extremes are where autocorrelated draws put extra mass
  (below).

**Sensitivity (not the verdict).** Four replicates have a smallest bulk ESS below
30: 4 (15), 17 (28), 24 (19) and 47 (19). Without them, 6 of 6 are within, on 46
replicates.

**F8:** `dev/notes/figures/v7-16a9b/f8_rank_ecdf.{pdf,png}` (primary) and
`f8_rank_hist.{pdf,png}` (secondary), from `dev/scripts/plot_v7_16a9b.py`.

## What the pass does and does not show

- **The design target of 100 effective draws per replicate was not met by 23 of
  the 50.** The smallest bulk ESS over the kept sweeps ranges from 15 to 325,
  median 109: 4 below 30, 19 from 30 to 100, 27 at 100 or more. §12 2026-10-08
  sized 656 sweeps from the slowest ESS rate measured on V6's one dataset, 0.180
  per sweep. Across prior-drawn θ, some replicates mix more slowly than that.
- **Autocorrelated draws widen the rank distribution.** A thinned set with
  fewer than L effective draws puts extra mass at the extreme ranks and makes the
  ECDF noisier. A large miscalibration would still bend the ECDF out of its band,
  but a small one that a well-mixed R = 50 could detect may be lost in that
  noise, or offset by it.
- **So the claim is:** at R = 50 and the mixing these chains reach, there is no
  evidence of miscalibration in any of the six parameters, and the
  pre-registered rule passes. It does not bound how small a miscalibration it
  would miss.
- **FBA,** flagged as reading high in 16a.5 (R16.4's watch item), is the flattest
  of the six (0.075).

## Reproduce

```
# replicates (run path), from a worktree at dfe570a:
ARRAY=1-50 bash dev/scripts/sbc_v7_16a9b.sh            # then cancel its merge
# merge, from a worktree at 0ff6cff:
V7_DIR=<v7 worktree>/dev/scripts/sbc_v7_16a9b \
  sbatch --export=ALL --time=01:00:00 --cpus-per-task=1 --mem=16G -J v7_merge \
  dev/scripts/sbc_v7_16a9b.slurm merge
# F8:
julia dev/scripts/v7_16a9b_export.jl <v7 worktree>/dev/scripts/sbc_v7_16a9b
python3 dev/scripts/plot_v7_16a9b.py
```

The checkpoints (`replicate*.jls`, `ranks.jls`) are untracked under
`.worktrees/v7-16a9b/dev/scripts/sbc_v7_16a9b/`. The ranks, per-replicate ESS
and bands that the plots read are committed as `v7_16a9b_ranks.tsv` and
`v7_16a9b_bands.tsv`.

## Per replicate (merge 18689602)

Ranks are in the order of the table above.

| Replicate | Seed | Min bulk ESS (kept) | Ranks (of 100) |
|---|---|---|---|
| 1 | 16200001 | 39 | 1, 83, 23, 88, 93, 78 |
| 2 | 16200002 | 37 | 24, 31, 54, 24, 2, 39 |
| 3 | 16200003 | 147 | 67, 38, 31, 14, 0, 18 |
| 4 | 16200004 | 15 | 3, 47, 20, 73, 18, 13 |
| 5 | 16200005 | 97 | 81, 14, 49, 42, 59, 34 |
| 6 | 16200006 | 122 | 9, 74, 14, 39, 15, 26 |
| 7 | 16200007 | 141 | 20, 68, 35, 82, 34, 23 |
| 8 | 16200008 | 211 | 97, 66, 93, 51, 67, 42 |
| 9 | 16200009 | 143 | 44, 18, 43, 66, 65, 95 |
| 10 | 16200010 | 113 | 94, 73, 99, 36, 96, 90 |
| 11 | 16200011 | 71 | 34, 50, 54, 7, 53, 67 |
| 12 | 16200012 | 104 | 1, 32, 27, 26, 83, 35 |
| 13 | 16200013 | 96 | 3, 2, 82, 86, 59, 28 |
| 14 | 16200014 | 157 | 87, 42, 91, 98, 85, 23 |
| 15 | 16200015 | 100 | 52, 89, 56, 79, 12, 19 |
| 16 | 16200016 | 200 | 22, 95, 30, 64, 18, 79 |
| 17 | 16200017 | 28 | 23, 19, 6, 19, 42, 47 |
| 18 | 16200018 | 92 | 81, 92, 93, 26, 83, 3 |
| 19 | 16200019 | 196 | 98, 92, 35, 25, 82, 53 |
| 20 | 16200020 | 93 | 29, 100, 51, 51, 65, 91 |
| 21 | 16200021 | 39 | 17, 9, 54, 70, 26, 33 |
| 22 | 16200022 | 170 | 35, 5, 8, 11, 98, 25 |
| 23 | 16200023 | 187 | 77, 88, 92, 98, 14, 66 |
| 24 | 16200024 | 19 | 41, 40, 78, 8, 69, 17 |
| 25 | 16200025 | 167 | 69, 26, 27, 39, 49, 30 |
| 26 | 16200026 | 57 | 65, 15, 100, 29, 30, 71 |
| 27 | 16200027 | 90 | 93, 31, 34, 39, 81, 31 |
| 28 | 16200028 | 131 | 53, 35, 60, 1, 96, 39 |
| 29 | 16200029 | 70 | 60, 36, 61, 18, 76, 100 |
| 30 | 16200030 | 203 | 37, 66, 12, 67, 52, 12 |
| 31 | 16200031 | 211 | 58, 37, 24, 21, 61, 26 |
| 32 | 16200032 | 55 | 33, 32, 25, 15, 39, 40 |
| 33 | 16200033 | 168 | 5, 27, 53, 95, 52, 53 |
| 34 | 16200034 | 90 | 55, 69, 50, 96, 40, 92 |
| 35 | 16200035 | 181 | 84, 82, 28, 63, 71, 78 |
| 36 | 16200036 | 101 | 66, 66, 28, 58, 33, 8 |
| 37 | 16200037 | 140 | 46, 30, 80, 39, 13, 50 |
| 38 | 16200038 | 54 | 97, 82, 86, 72, 56, 78 |
| 39 | 16200039 | 32 | 27, 18, 10, 48, 22, 86 |
| 40 | 16200040 | 115 | 14, 0, 15, 62, 37, 92 |
| 41 | 16200041 | 240 | 39, 37, 0, 99, 26, 81 |
| 42 | 16200042 | 164 | 36, 36, 17, 39, 84, 7 |
| 43 | 16200043 | 139 | 95, 43, 87, 33, 74, 83 |
| 44 | 16200044 | 209 | 0, 63, 12, 32, 78, 31 |
| 45 | 16200045 | 93 | 88, 69, 96, 39, 82, 63 |
| 46 | 16200046 | 221 | 28, 94, 88, 64, 98, 55 |
| 47 | 16200047 | 19 | 33, 79, 28, 41, 47, 33 |
| 48 | 16200048 | 44 | 27, 84, 8, 27, 1, 34 |
| 49 | 16200049 | 76 | 19, 84, 9, 47, 98, 44 |
| 50 | 16200050 | 325 | 94, 18, 0, 69, 5, 78 |
