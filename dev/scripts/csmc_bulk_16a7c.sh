#!/bin/bash
# Submit 16a.7c's measurement: particle Gibbs under bulk metabolites on M0, at
# N in {5, 10, 20, 50}, then the merge. Run from a worktree at a pushed commit.
set -euo pipefail
mkdir -p dev/scripts/csmc_bulk_16a7c
S=dev/scripts/csmc_bulk_16a7c.slurm
SC=16                     # scans of the five updated cells
jobs=()
for N in 5 10 20 50; do
    t=$(( N * 5 * SC * 3 / 60 + 40 ))     # minutes: about N cycles of ~1.5 s per cell-sweep, 2x margin
    jobs+=($(sbatch --parsable --time=$t -J b7c_N$N $S run $N $SC))
done
dep=$(IFS=:; echo "${jobs[*]}")
m=$(sbatch --parsable --time=00:15:00 -J b7c_merge --dependency=afterany:$dep $S merge)
echo "runs ${jobs[*]}; merge $m"
