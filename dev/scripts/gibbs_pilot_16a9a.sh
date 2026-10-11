#!/bin/bash
# Submit 16a.9a's pilot: four Gibbs chains on M0, then the merge. Run from a
# worktree at a pushed commit. A chain resumes from its checkpoint, so a job that
# hits its time limit is resubmitted with the same command.
set -euo pipefail
mkdir -p dev/scripts/gibbs_pilot_16a9a
S=dev/scripts/gibbs_pilot_16a9a.slurm
SWEEPS=${SWEEPS:-200}
TIME=${TIME:-2-00:00:00}
jobs=()
for k in 1 2 3 4; do
    jobs+=($(sbatch --parsable --time=$TIME -J p9a_c$k $S run $k $SWEEPS))
done
dep=$(IFS=:; echo "${jobs[*]}")
m=$(sbatch --parsable --time=00:30:00 --cpus-per-task=1 --mem=16G -J p9a_merge --dependency=afterany:$dep $S merge)
echo "chains ${jobs[*]}; merge $m"
