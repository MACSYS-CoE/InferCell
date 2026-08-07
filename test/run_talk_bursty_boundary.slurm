#!/bin/bash
#SBATCH --job-name=infercell_isab_boundary
#SBATCH --partition=milan
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=1
#SBATCH --mem=16G
#SBATCH --time=08:00:00
#SBATCH --output=test/isab_boundary_%j.log

# Settings come from the script's own defaults (400 particles, 6 populations,
# 50 replicates, 1000 NUTS samples). Override for a smoke test with:
#   sbatch --time=00:30:00 \
#     --export=ALL,INFERCELL_N_PARTICLES=30,INFERCELL_N_POPULATIONS=2,\
# INFERCELL_N_REPLICATES=10,INFERCELL_N_NUTS=100 \
#     test/run_talk_bursty_boundary.slurm

module load julia/1.10.5

echo "ISAB boundary figure run at $(date)"
echo "Node: $(hostname)"
echo ""

# Compute nodes have no route to the internet -- Pkg cannot reach the registry
# from here, so package DOWNLOADS must be done first from the login node:
#
#   export EBROOTJULIA=/apps/modules/software/Julia/1.10.5
#   export PATH="$EBROOTJULIA/bin:$PATH"
#   JULIA_PKG_PRECOMPILE_AUTO=0 julia --project -e 'using Pkg; Pkg.instantiate()'
#
# (JULIA_PKG_PRECOMPILE_AUTO=0 keeps the heavy precompile off the login node;
# it happens here instead, on the first `using`, and needs no network.)
# The depot under ~/.julia was found empty in Aug 2026 -- it was populated for
# the May integration runs and has been cleared since -- so this is not a
# one-off. OFFLINE=1 makes a missing package fail fast and legibly rather than
# hanging on a DNS timeout.
export JULIA_PKG_OFFLINE=1
echo "Verifying project environment (offline)..."
julia --project -e 'using Pkg; Pkg.instantiate()' || {
    echo "Environment incomplete -- run the login-node download above." >&2
    exit 1
}
echo "Environment ready at $(date)"
echo ""

julia --project test/run_talk_bursty_boundary.jl

echo ""
echo "Finished at $(date)"
