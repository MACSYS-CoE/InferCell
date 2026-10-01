#!/bin/bash
# Submit V5's stages for one case (m0 or toy), in order
# (spec/phases/16-recovery.md task 16a.7). Run from a worktree at a pushed commit.
set -euo pipefail
case=$1
mkdir -p dev/scripts/csmc_v5_16a7
S=dev/scripts/csmc_v5_16a7.slurm
t=$([ "$case" = m0 ] && echo 00:40:00 || echo 01:00:00)
r=$(sbatch --parsable --time=$t --array=0-99 -J v5_${case}_reject $S reject $case)
p=$(sbatch --parsable --time=00:30:00 -J v5_${case}_prep --dependency=afterok:$r $S prep $case)
koff=$(sbatch --parsable --time=02:00:00 --array=0-49 -J v5_${case}_off --dependency=afterok:$p $S kernel $case off)
kon=$(sbatch --parsable --time=02:00:00 --array=0-49 -J v5_${case}_on --dependency=afterok:$p $S kernel $case on)
m=$(sbatch --parsable --time=00:30:00 -J v5_${case}_merge --dependency=afterok:$koff:$kon $S merge $case)
echo "$case: reject $r, prep $p, kernel off $koff, on $kon, merge $m"
