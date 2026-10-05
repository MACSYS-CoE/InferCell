#!/bin/bash
# Submit the fresh V5 reference for one case (m0 or toy) and check earlier
# kernel draws against it (task 16a.7c's watch item). Each further argument is
# an earlier run's csmc_v5_16a7/<case> directory holding kernel draws.
# POOL=<the original run's csmc_v5_16a7/<case>> also submits the pooled checks,
# against the fresh and original rejection sets combined, which are the V5 result
# of record for 16a.7c (jobs 18061130 to 18061135):
#   POOL=<v5 dir>/m0 csmc_v5_fresh_16a7c.sh m0 <v5 dir>/m0 <annealed dir>/m0 <bulk dir>/m0
# The three kernel directories are the per-cell run (csmc_v5_16a7.sh), the
# annealed run (csmc_v5_anneal_16a7.sh) and the bulk run (csmc_v5_bulk_16a7c.sh).
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
    if [ -n "${POOL:-}" ]; then
        p=$(sbatch --parsable --export=ALL,V5_FRESH=1,V5_POOL=$POOL --time=01:00:00 -J v5p_${case}_check \
            --dependency=afterok:$r $S refcheck $case "$d")
        echo "$case: pooled refcheck $p against $d"
    fi
done
echo "$case: reject $r"
