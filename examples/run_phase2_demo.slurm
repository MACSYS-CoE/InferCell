#!/bin/bash
#SBATCH --job-name=infercell_p2
#SBATCH --partition=milan
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=1
#SBATCH --mem=16G
#SBATCH --time=02:00:00
#SBATCH --output=examples/phase2_demo_%j.log

module load julia/1.10.5

echo "Starting Phase 2 demo at $(date)"
echo "Node: $(hostname)"
echo ""

julia --project examples/phase2_demo.jl

echo ""
echo "Finished at $(date)"
