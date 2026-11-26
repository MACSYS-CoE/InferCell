#!/bin/bash
#SBATCH --job-name=infercell_kim90
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=1
#SBATCH --mem=16G

# KIM-90: one stage of the m0_v6 refits or the cross-covariance check.
# Submit through dev/scripts/gap_diag_refit.sh. Arguments: <script> <stage args…>.
export EBROOTJULIA=/apps/modules/software/Julia/1.10.5
export PATH="$EBROOTJULIA/bin:$PATH"
export JULIA_PKG_OFFLINE=1
script=$1; shift
echo "KIM-90 $script $* at $(date) on $(hostname), HEAD $(git rev-parse --short HEAD)$(test -z "$(git status --porcelain -- src dev/scripts/gap_diag_*)" || echo ' (tree dirty)')"
julia --project "dev/scripts/$script" "$@"
status=$?
echo "Exit code: $status  Finished at $(date)"
exit $status
