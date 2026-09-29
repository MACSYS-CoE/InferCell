#!/bin/bash
#SBATCH --job-name=infercell_sensitivity_15
#SBATCH --partition=milan
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=1
#SBATCH --mem=8G
#SBATCH --time=01:30:00
#SBATCH --array=0-99
#SBATCH --output=dev/scripts/corea_sensitivity_15_%A_%a.log

# The ensemble sensitivity of spec §11 tasks 15.6 and 15.8 (§12 2026-09-29):
# NSEEDS paired seeds, one per task, each through 15 configurations. The pilot
# is 100 seeds; the full run is sized from it and submitted with a matching
# --array and NSEEDS, e.g.
#   sbatch --array=0-399 --export=ALL,NSEEDS=400 dev/scripts/corea_sensitivity_15.slurm
# Merge afterwards with
#   julia --project dev/scripts/corea_sensitivity_15_merge.jl
module load julia/1.10.5
N=${NSEEDS:-100}
echo "Core A′ sensitivity 15, task ${SLURM_ARRAY_TASK_ID} of $N at $(date) on $(hostname), HEAD $(git rev-parse --short HEAD)$(test -z "$(git status --porcelain -- src dev/scripts/corea_sensitivity_15.jl)" || echo ' (tree dirty)')"
julia --project dev/scripts/corea_sensitivity_15.jl "${SLURM_ARRAY_TASK_ID}" "$N" "$N"
status=$?
echo "Exit code: $status  Finished at $(date)"
exit $status
