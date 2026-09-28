#!/bin/bash
#SBATCH --job-name=infercell_census_14c
#SBATCH --partition=milan
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=1
#SBATCH --mem=8G
#SBATCH --time=06:00:00
#SBATCH --array=0-19
#SBATCH --output=dev/scripts/corea_census_14c_%A_%a.log

# Check 7's census rerun for phase 14c (spec §11 tasks 14c.2, 14c.3, 14c.6).
# Usage: sbatch dev/scripts/corea_census_14c.slurm <d11|broad|gtp>
# Merge afterwards with julia --project dev/scripts/corea_census_14c_merge.jl
module load julia/1.10.5
MODE=$1
echo "Core A′ census 14c, mode ${MODE}, task ${SLURM_ARRAY_TASK_ID} at $(date) on $(hostname), HEAD $(git rev-parse --short HEAD)$(test -z "$(git status --porcelain -- src test dev/scripts/corea_census_14c.jl)" || echo ' (tree dirty)')"
julia --project dev/scripts/corea_census_14c.jl "${MODE}" "${SLURM_ARRAY_TASK_ID}" 20
status=$?
echo "Exit code: $status  Finished at $(date)"
exit $status
