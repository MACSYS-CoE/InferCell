#!/bin/bash
#SBATCH --job-name=infercell_v1_m0
#SBATCH --partition=milan
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=1
#SBATCH --mem=4G
#SBATCH --time=01:00:00
#SBATCH --output=dev/scripts/density_16a2_m0/slurm_%A_%a.log

# V1 importance identity on M0 (spec/phases/16-recovery.md task 16a.2).
module load julia/1.10.5
echo "16a.2 V1 on M0, task $SLURM_ARRAY_TASK_ID, at $(date) on $(hostname), HEAD $(git rev-parse --short HEAD)$(test -z "$(git status --porcelain -- src test dev/scripts/path_density_16a2_m0.jl)" || echo ' (tree dirty)')"
julia --project dev/scripts/path_density_16a2_m0.jl "$SLURM_ARRAY_TASK_ID" 100
status=$?
echo "Exit code: $status  Finished at $(date)"
exit $status
