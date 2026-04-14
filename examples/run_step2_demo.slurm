#!/bin/bash
#SBATCH --job-name=infercell_step2
#SBATCH --partition=milan
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=1
#SBATCH --mem=16G
#SBATCH --time=02:00:00
#SBATCH --output=examples/step2_demo_%j.log

module load julia/1.10.5

echo "Starting Step 2 joint inference demo at $(date)"
echo "Node: $(hostname)"
echo ""

julia --project examples/step2_joint_demo.jl

echo ""
echo "Finished at $(date)"
