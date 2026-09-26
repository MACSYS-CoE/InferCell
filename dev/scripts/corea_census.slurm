#!/bin/bash
#SBATCH --job-name=infercell_corea_census
#SBATCH --partition=milan
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=1
#SBATCH --mem=8G
#SBATCH --time=04:00:00
#SBATCH --array=0-19
#SBATCH --output=dev/scripts/census/corea_census_%A_%a.log

# Check 7's census (spec §11 task 14.8): 10 published seeds and 201 draws over
# 20 tasks, about 11 cycles each. Merge afterwards with
#   julia --project dev/scripts/corea_census_merge.jl
module load julia/1.10.5
echo "Core A′ census task ${SLURM_ARRAY_TASK_ID} at $(date) on $(hostname), HEAD $(git rev-parse --short HEAD)$(test -z "$(git status --porcelain -- src test dev/scripts/corea_census.jl)" || echo ' (tree dirty)')"
echo "SMOKE=${SMOKE:-0}"
export NTASK=20
julia --project dev/scripts/corea_census.jl "${SLURM_ARRAY_TASK_ID}" 20
status=$?
echo "Exit code: $status  Finished at $(date)"
exit $status
