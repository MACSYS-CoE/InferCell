#!/bin/bash
# Submit the fresh V5 reference for one case (m0 or toy) and check earlier
# kernel draws against it (task 16a.7c's watch item). Each further argument is
# an earlier run's csmc_v5_16a7/<case> directory holding kernel draws.
# Run from a worktree at a pushed commit.
set -euo pipefail
case=$1
shift
mkdir -p dev/scripts/csmc_v5_16a7/${case}_fresh
S=dev/scripts/csmc_v5_16a7.slurm
E="--export=ALL,V5_FRESH=1"
t=$([ "$case" = m0 ] && echo 00:40:00 || echo 01:00:00)
r=$(sbatch --parsable $E --time=$t --array=0-99 -J v5f_${case}_reject $S reject $case)
for d in "$@"; do
    c=$(sbatch --parsable $E --time=00:30:00 -J v5f_${case}_check --dependency=afterok:$r $S refcheck $case "$d")
    echo "$case: refcheck $c against $d"
done
echo "$case: reject $r"
