#!/bin/bash
# Submit V4 on block 1's bulk form (task 16a.7c, spec/phases/16-recovery.md §12
# 2026-10-05). As block1_16a5.sh, with V4_BULK=1, and with the FBA-fine grid run
# from the start, since a coarse FBA axis biased the per-cell V4. If merge 2 finds
# mass at a box edge, it writes box3 and the chains fail for want of exact.jls;
# then run grid 3 and merge 3 by hand, as the per-cell V4 did.
# Run from a worktree at a pushed commit.
set -euo pipefail
mkdir -p dev/scripts/block1_16a5_bulk dev/scripts/block1_16a5
S=dev/scripts/block1_16a5.slurm
E="--export=ALL,V4_BULK=1"
g1=$(sbatch --parsable $E --time=00:40:00 --array=0-20 -J v4b_grid1 $S grid 1)
m1=$(sbatch --parsable $E --time=00:20:00 -J v4b_merge1 --dependency=afterok:$g1 $S merge 1)
g2=$(sbatch --parsable $E --time=01:00:00 --array=0-40 -J v4b_grid2 --dependency=afterok:$m1 $S grid 2)
m2=$(sbatch --parsable $E --time=00:20:00 -J v4b_merge2 --dependency=afterok:$g2 $S merge 2)
fb=$(sbatch --parsable $E --time=01:30:00 --array=0-40 -J v4b_fba --dependency=afterok:$m2 $S fba)
fm=$(sbatch --parsable $E --time=00:20:00 -J v4b_fbamerge --dependency=afterok:$fb $S fbamerge)
ch=$(sbatch --parsable $E --time=06:00:00 --array=0-19 -J v4b_chain --dependency=afterok:$m2 $S chain)
ck=$(sbatch --parsable $E --time=00:30:00 -J v4b_check --dependency=afterok:$m2 $S check)
fi=$(sbatch --parsable $E --time=00:20:00 -J v4b_final --dependency=afterok:$ch:$fm $S final)
echo "grid1 $g1, merge1 $m1, grid2 $g2, merge2 $m2, fba $fb, fbamerge $fm, chains $ch, check $ck, final $fi"
