#!/bin/bash
#SBATCH --job-name=infercell_jre89
#SBATCH --partition=milan
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=1
#SBATCH --mem=4G
#SBATCH --time=01:00:00
#SBATCH --output=dev/scripts/jensen_rethink_89_%A_%a.log

# Issue #89, after B3: one stage of the rethink. Submit through
# dev/scripts/jensen_rethink_89.sh. Arguments: the stage's own, with the array
# index inserted for the ensemble stage.
module load julia/1.10.5
echo "#89 rethink $* (array task ${SLURM_ARRAY_TASK_ID:-none}) at $(date) on $(hostname), HEAD $(git rev-parse --short HEAD)$(test -z "$(git status --porcelain -- src dev/scripts/jensen_rethink_89.jl)" || echo ' (tree dirty)')"
s=dev/scripts/jensen_rethink_89.jl
case "$1" in
    ensemble)  julia --project $s ensemble "$2" "$3" "$SLURM_ARRAY_TASK_ID" "$4" "$5" "$6" "$7" ;;
    meanfield) julia --project $s meanfield "$2" "$3" ;;
    popmean)   julia --project $s popmean "$2" "$3" ;;
    score)     julia --project $s score "$2" ;;
    *) echo "unknown stage $1"; exit 2 ;;
esac
status=$?
echo "Exit code: $status  Finished at $(date)"
exit $status
