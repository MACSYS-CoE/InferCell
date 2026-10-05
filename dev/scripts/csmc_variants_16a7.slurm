#!/bin/bash
#SBATCH --job-name=infercell_variants
#SBATCH --partition=milan
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=1
#SBATCH --mem=4G
#SBATCH --output=dev/scripts/csmc_variants_16a7/slurm_%x_%A_%a.log

# One task of 16a.7's variant choice (spec/phases/16-recovery.md D16.3);
# csmc_variants_16a7.sh submits them. Arguments: run <scale> <variant> <N>
# <sweeps>, with the array index as the task, or merge.
module load julia/1.10.5
echo "16a.7 variants $* task ${SLURM_ARRAY_TASK_ID:-none} at $(date) on $(hostname), HEAD $(git rev-parse --short HEAD)$(test -z "$(git status --porcelain -- src test dev/scripts/csmc_variants_16a7.jl)" || echo ' (tree dirty)')"
if [ "$1" = run ]; then
    julia --project dev/scripts/csmc_variants_16a7.jl run "$2" "$3" "$4" "${SLURM_ARRAY_TASK_ID}" "$5"
else
    julia --project dev/scripts/csmc_variants_16a7.jl "$@"
fi
status=$?
echo "Exit code: $status  Finished at $(date)"
exit $status
