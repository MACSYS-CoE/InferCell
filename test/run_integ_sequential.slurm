#!/bin/bash
#SBATCH --job-name=infercell_integ_seq
#SBATCH --partition=milan
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=1
#SBATCH --mem=16G
#SBATCH --time=04:00:00
#SBATCH --output=test/test_integ_%j.log

module load julia/1.10.5

echo "Starting sequential inference integration test at $(date)"
echo "Node: $(hostname)"
echo ""

julia --project test/run_integ_sequential.jl

echo ""
echo "Exit code: $?"
echo "Finished at $(date)"
