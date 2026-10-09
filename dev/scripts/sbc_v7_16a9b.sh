#!/bin/bash
# Submit 16a.9b's V7: R replicate chains as an array, then the merge. Run from a
# worktree at a pushed commit. A replicate resumes from its checkpoint, so one
# that hits its time limit is resubmitted with the same command (ARRAY=<ids>).
# Smoke run: V7_DIR, V7_NCELLS, V7_HORIZON, V7_SWEEPS, V7_BURN and V7_L are
# passed through to the jobs.
set -euo pipefail
mkdir -p dev/scripts/sbc_v7_16a9b
S=dev/scripts/sbc_v7_16a9b.slurm
R=${R:-50}
ARRAY=${ARRAY:-1-$R}
TIME=${TIME:-6-00:00:00}
a=$(sbatch --parsable --export=ALL --time=$TIME --array=$ARRAY $S)
m=$(sbatch --parsable --export=ALL --time=01:00:00 --cpus-per-task=1 --mem=16G -J v7_merge --dependency=afterany:$a $S merge)
echo "replicates $a (array $ARRAY); merge $m"
