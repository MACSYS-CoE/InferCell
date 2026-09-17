#!/bin/bash
#SBATCH --job-name=infercell_integration
#SBATCH --partition=milan
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=1
#SBATCH --mem=16G
#SBATCH --time=06:00:00
#SBATCH --output=test/integration_%j.log

# The full suite with the Bayesian integration tests switched on. This is the
# home of record for them: they no longer run on push to `main`, because every
# such run hit the hosted runner's 90-minute ceiling without completing. The
# walltime here is deliberately generous — ABC-SMC is serial and the two SSA
# files dominate.

module load julia/1.10.5

# Compute nodes have no network. Skips the one Aqua check that needs the
# registry; see the comment in test/runtests.jl.
export INFERCELL_OFFLINE_TESTS=true
export INFERCELL_INTEGRATION_TESTS=true

echo "Running InferCell integration suite at $(date)"
echo "Node: $(hostname)"
echo "HEAD: $(git rev-parse --short HEAD)$(test -z "$(git status --porcelain -- src test Project.toml Manifest.toml)" || echo ' (tree dirty)')"
echo ""

julia --project -e 'using Pkg; Pkg.test()'

echo ""
echo "Finished at $(date)"
