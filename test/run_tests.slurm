#!/bin/bash
#SBATCH --job-name=infercell_test
#SBATCH --partition=milan
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=1
#SBATCH --mem=16G
#SBATCH --time=00:30:00
#SBATCH --output=test/test_%j.log

module load julia/1.10.5

# Compute nodes have no network. Skips the one Aqua check that needs the
# registry; see the comment in test/runtests.jl.
export INFERCELL_OFFLINE_TESTS=true

echo "Running InferCell tests at $(date)"
echo "Node: $(hostname)"
echo ""

julia --project -e 'using Pkg; Pkg.test()'

echo ""
echo "Finished at $(date)"
