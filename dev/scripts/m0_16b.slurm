#!/bin/bash
#SBATCH --job-name=infercell_16b
#SBATCH --partition=milan
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=5
#SBATCH --mem=8G

# One stage of 16b.1 or 16b.4 at M0 (spec/phases/16-recovery.md §12): run <sweeps>,
# rates, or summary <v6dir> <dir>.... m0_16b.sh submits the chains. RUN_MODE and
# PILOT_DIR come from the environment.
module load julia/1.10.5
echo "16b $RUN_MODE $* at $(date) on $(hostname), HEAD $(git rev-parse --short HEAD)$(test -z "$(git status --porcelain -- src dev/scripts/m0_16b.jl)" || echo ' (tree dirty)')"
julia --project --threads=$SLURM_CPUS_PER_TASK dev/scripts/m0_16b.jl "$@"
status=$?
echo "Exit code: $status  Finished at $(date)"
exit $status
