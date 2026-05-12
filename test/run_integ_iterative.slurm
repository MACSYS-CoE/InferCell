#!/bin/bash
#SBATCH --job-name=infercell_integ_iter
#SBATCH --partition=milan
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=1
#SBATCH --mem=16G
#SBATCH --time=04:00:00
#SBATCH --output=test/integ_iter_%j.log

module load julia/1.10.5

echo "Running iterative-boundary integration at $(date)"
echo "Node: $(hostname)"
echo ""

julia --project test/run_integ_iterative.jl

echo ""
echo "Finished at $(date)"
