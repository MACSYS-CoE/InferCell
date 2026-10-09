#!/bin/bash
#SBATCH --job-name=infercell_pilot9a
#SBATCH --partition=milan
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=20
#SBATCH --mem=32G
#SBATCH --output=dev/scripts/gibbs_pilot_16a9a/slurm_%x_%j.log

# One stage of 16a.9a's pilot (spec/phases/16-recovery.md §12 2026-10-06):
# run <chain> <sweeps>, or merge [burn]. gibbs_pilot_16a9a.sh submits them.
# The 20 threads run a cell's 20 particles, and block 1's cells, in parallel.
module load julia/1.10.5
echo "16a.9a pilot $* at $(date) on $(hostname), HEAD $(git rev-parse --short HEAD)$(test -z "$(git status --porcelain -- src dev/scripts/gibbs_pilot_16a9a.jl)" || echo ' (tree dirty)')"
julia --project --threads=$SLURM_CPUS_PER_TASK dev/scripts/gibbs_pilot_16a9a.jl "$@"
status=$?
echo "Exit code: $status  Finished at $(date)"
exit $status
