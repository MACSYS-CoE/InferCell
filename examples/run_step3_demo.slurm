#!/bin/bash
#SBATCH --job-name=infercell_step3
#SBATCH --partition=milan
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=1
#SBATCH --mem=16G
#SBATCH --time=04:00:00
#SBATCH --output=examples/step3_demo_%j.log

module load julia/1.10.5

echo "Starting Step 3 boundary protocol demo at $(date)"
echo "Node: $(hostname)"
echo ""

julia --project examples/step3_boundary_demo.jl

echo ""
echo "Finished at $(date)"
