#!/bin/bash
#SBATCH --job-name=infercell_sens88
#SBATCH --partition=milan
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=1
#SBATCH --mem=8G
#SBATCH --time=01:30:00
#SBATCH --output=dev/scripts/sensitivity_88_%A_%a.log

# Issue #88 A1/A2: one array task of the ensemble sensitivity. Submit through
# dev/scripts/sensitivity_88.sh, which sets MODEL, NSEEDS, NTASKS, OUTDIR and
# the array range.
module load julia/1.10.5
echo "#88 sensitivity $MODEL, task ${SLURM_ARRAY_TASK_ID} of $NTASKS ($NSEEDS seeds) to $OUTDIR at $(date) on $(hostname), HEAD $(git rev-parse --short HEAD)$(test -z "$(git status --porcelain -- src dev/scripts/sensitivity_88.jl)" || echo ' (tree dirty)')"
julia --project dev/scripts/sensitivity_88.jl "$MODEL" "${SLURM_ARRAY_TASK_ID}" "$NTASKS" "$NSEEDS" "$OUTDIR"
status=$?
echo "Exit code: $status  Finished at $(date)"
exit $status
