#!/bin/bash
#SBATCH --job-name=infercell_corea_14b
#SBATCH --partition=milan
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=1
#SBATCH --mem=8G
#SBATCH --time=04:00:00
#SBATCH --output=dev/scripts/validation_14b/corea_14b_%x_%A_%a.log

# Phase 14b's full-scale checks (spec §11). One section per job; pass it with
# --export, since a bare VAR=1 before sbatch does not reach the job:
#
#   sbatch --export=ALL,SECTION=floor  --array=0 dev/scripts/corea_validation_14b.slurm
#   sbatch --export=ALL,SECTION=check8 --array=0 dev/scripts/corea_validation_14b.slurm
#   sbatch --export=ALL,SECTION=check6 --array=0 dev/scripts/corea_validation_14b.slurm
#   sbatch --export=ALL,SECTION=f5     --array=0 dev/scripts/corea_validation_14b.slurm
#   # after floor:
#   sbatch --export=ALL,SECTION=band,NTASK=20 --array=0-19 dev/scripts/corea_validation_14b.slurm
#   # after band:
#   sbatch --export=ALL,SECTION=merge-band --array=0 dev/scripts/corea_validation_14b.slurm
#   # last:
#   sbatch --export=ALL,SECTION=report --array=0 dev/scripts/corea_validation_14b.slurm
module load julia/1.10.5
echo "Core A′ phase 14b SECTION=${SECTION} task ${SLURM_ARRAY_TASK_ID:-0}/${NTASK:-1} at $(date) on $(hostname), HEAD $(git rev-parse --short HEAD)$(test -z "$(git status --porcelain -- src test dev/scripts/corea_validation_14b.jl)" || echo ' (tree dirty)')"
echo "SMOKE=${SMOKE:-0}"
julia --project dev/scripts/corea_validation_14b.jl "${SLURM_ARRAY_TASK_ID:-0}" "${NTASK:-1}"
status=$?
echo "Exit code: $status  Finished at $(date)"
exit $status
