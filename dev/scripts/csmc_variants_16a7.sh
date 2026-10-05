#!/bin/bash
# Submit 16a.7's KP measurement (spec/phases/16-recovery.md D16.3 and §8 KP):
# every variant at N = 50 on M0's five cells, and particle Gibbs at N = 50 on the
# full-scale cell, then the merge. Run from a worktree at a pushed commit.
set -euo pipefail
mkdir -p dev/scripts/csmc_variants_16a7
S=dev/scripts/csmc_variants_16a7.slurm
N=50
SW=8                      # sweeps per task
jobs=()
for v in pg l1 l5 full; do
    case $v in pg) t=01:30:00;; l1) t=02:30:00;; l5) t=05:00:00;; full) t=06:00:00;; esac
    jobs+=($(sbatch --parsable --time=$t --array=1-5 -J v_m0_${v} $S run m0 $v $N $SW))
done
jobs+=($(sbatch --parsable --time=08:00:00 --array=1-2 -J v_full_pg $S run full pg $N $((SW / 2))))
dep=$(IFS=:; echo "${jobs[*]}")
m=$(sbatch --parsable --time=00:15:00 -J v_merge --dependency=afterany:$dep $S merge)
echo "runs ${jobs[*]}; merge $m"
