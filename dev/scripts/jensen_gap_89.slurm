#!/bin/bash
#SBATCH --job-name=infercell_jgap89
#SBATCH --partition=milan
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=1
#SBATCH --mem=4G
#SBATCH --time=01:00:00
#SBATCH --output=dev/scripts/jensen_gap_89_%A_%a.log

# Issue #89 B3: one stage of the Jensen-gap measurement. Submit through
# dev/scripts/jensen_gap_89.sh. Arguments: the stage's own, without the array
# index, which is inserted for the ensemble stage.
module load julia/1.10.5
echo "#89 B3 $* (array task ${SLURM_ARRAY_TASK_ID:-none}) at $(date) on $(hostname), HEAD $(git rev-parse --short HEAD)$(test -z "$(git status --porcelain -- src dev/scripts/jensen_gap_89.jl)" || echo ' (tree dirty)')"
case "$1" in
    ensemble)  julia --project dev/scripts/jensen_gap_89.jl ensemble "$2" "$SLURM_ARRAY_TASK_ID" "$3" "$4" ;;
    meanfield) julia --project dev/scripts/jensen_gap_89.jl meanfield "$2" "$3" ;;
    merge)     julia --project dev/scripts/jensen_gap_89.jl merge "$2" ;;
    *) echo "unknown stage $1"; exit 2 ;;
esac
status=$?
echo "Exit code: $status  Finished at $(date)"
exit $status
