#!/bin/bash
# Submit V5's bulk "on" half (task 16a.7c, spec/phases/16-recovery.md §12
# 2026-10-05) for one case (m0 or toy). The rejection draws are the unannealed
# V5 run's: pass the directory holding them (its dev/scripts/csmc_v5_16a7).
# Run from a worktree at a pushed commit.
set -euo pipefail
case=$1
src=$2
mkdir -p dev/scripts/csmc_v5_16a7/$case
# In a single clone the source can be this tree; then there is nothing to copy.
if [ "$(cd "$src" && pwd -P)" != "$(cd dev/scripts/csmc_v5_16a7 && pwd -P)" ]; then
    cp "$src"/$case/reject_*.jls dev/scripts/csmc_v5_16a7/$case/
fi
S=dev/scripts/csmc_v5_16a7.slurm
p=$(sbatch --parsable --time=01:00:00 -J v5b_${case}_prep $S bulkprep $case)
k=$(sbatch --parsable --time=06:00:00 --array=0-49 -J v5b_${case}_kernel --dependency=afterok:$p $S bulkkernel $case)
m=$(sbatch --parsable --time=00:30:00 -J v5b_${case}_merge --dependency=afterok:$k $S bulkmerge $case)
echo "$case bulk: prep $p, kernel $k, merge $m"
