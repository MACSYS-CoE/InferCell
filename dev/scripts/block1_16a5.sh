#!/bin/bash
# Submit V4's stages in order (spec/phases/16-recovery.md task 16a.5).
# Run from the tree the jobs should read, ideally a worktree at a pushed commit.
set -euo pipefail
mkdir -p dev/scripts/block1_16a5
S=dev/scripts/block1_16a5.slurm
g1=$(sbatch --parsable --time=00:40:00 --array=0-20 -J v4_grid1 $S grid 1)
m1=$(sbatch --parsable --time=00:20:00 -J v4_merge1 --dependency=afterok:$g1 $S merge 1)
g2=$(sbatch --parsable --time=01:00:00 --array=0-40 -J v4_grid2 --dependency=afterok:$m1 $S grid 2)
m2=$(sbatch --parsable --time=00:20:00 -J v4_merge2 --dependency=afterok:$g2 $S merge 2)
ch=$(sbatch --parsable --time=06:00:00 --array=0-19 -J v4_chain --dependency=afterok:$m2 $S chain)
fi=$(sbatch --parsable --time=00:20:00 -J v4_final --dependency=afterok:$ch $S final)
echo "grid1 $g1, merge1 $m1, grid2 $g2, merge2 $m2, chains $ch, final $fi"
