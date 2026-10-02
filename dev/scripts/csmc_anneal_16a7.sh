#!/bin/bash
# Submit 16a.7a's measurement (spec/phases/16-recovery.md §12 2026-10-02): particle
# Gibbs with annealed windows, K ∈ {20, 50, 100}, at N ∈ {2, 5} on M0's five cells,
# then the merge. Run from a worktree at a pushed commit.
set -euo pipefail
mkdir -p dev/scripts/csmc_variants_16a7
S=dev/scripts/csmc_variants_16a7.slurm
SW=8                      # sweeps per cell
jobs=()
for K in 20 50 100; do
    for N in 2 5; do
        t=$(( (N * (K + 1) * SW * 3 / 60) + 30 ))     # minutes: ~1.5 s a cycle, 2× margin
        jobs+=($(sbatch --parsable --time=$t --array=1-5 -J va_m0_a${K}_N${N} $S run m0 a$K $N $SW))
    done
done
dep=$(IFS=:; echo "${jobs[*]}")
m=$(sbatch --parsable --time=00:15:00 -J va_merge --dependency=afterany:$dep $S merge)
echo "runs ${jobs[*]}; merge $m"
