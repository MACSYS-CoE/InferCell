#!/bin/bash
#SBATCH --job-name=infercell_enoinfo
#SBATCH --partition=milan
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=1
#SBATCH --mem=4G
#SBATCH --time=00:40:00
#SBATCH --output=dev/scripts/eno_path_information_16a7/slurm_%x_%A_%a.log

# F6 ENO cell: complete-data information from full paths (spec/phases/16-recovery.md
# task 16a.7, handoff 2026-10-03). Arguments: run (array index = cell) or merge.
module load julia/1.10.5
echo "16a.7 ENO path information $* at $(date) on $(hostname), HEAD $(git rev-parse --short HEAD)$(test -z "$(git status --porcelain -- src dev/scripts/eno_path_information_16a7.jl)" || echo ' (tree dirty)')"
julia --project dev/scripts/eno_path_information_16a7.jl "$1" ${SLURM_ARRAY_TASK_ID:+"$SLURM_ARRAY_TASK_ID"}
status=$?
echo "Exit code: $status  Finished at $(date)"
exit $status
