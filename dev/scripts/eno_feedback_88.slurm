#!/bin/bash
#SBATCH --job-name=infercell_enofb88
#SBATCH --partition=milan
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=1
#SBATCH --mem=4G
#SBATCH --time=01:00:00
#SBATCH --output=dev/scripts/eno_feedback_88_%A_%a.log

# Issue #88 A3: one stage of the reverse-feedback measurement. Submit through
# dev/scripts/eno_feedback_88.sh. Arguments: the stage's own, without the array
# index, which is inserted for the ensemble and path stages.
module load julia/1.10.5
echo "#88 A3 $* (array task ${SLURM_ARRAY_TASK_ID:-none}) at $(date) on $(hostname), HEAD $(git rev-parse --short HEAD)$(test -z "$(git status --porcelain -- src dev/scripts/eno_feedback_88.jl)" || echo ' (tree dirty)')"
case "$1" in
    ensemble) julia --project dev/scripts/eno_feedback_88.jl ensemble "$2" "$3" "$SLURM_ARRAY_TASK_ID" "$4" "$5" ;;
    path)     julia --project dev/scripts/eno_feedback_88.jl path "$2" "$3" "$SLURM_ARRAY_TASK_ID" "$4" ;;
    merge)    julia --project dev/scripts/eno_feedback_88.jl merge "$2" ;;
    *) echo "unknown stage $1"; exit 2 ;;
esac
status=$?
echo "Exit code: $status  Finished at $(date)"
exit $status
