#!/bin/bash
#SBATCH --job-name=infercell_calib_smoke
#SBATCH --partition=milan
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=1
#SBATCH --mem=8G
#SBATCH --time=00:15:00
#SBATCH --output=examples/calib_smoke_%j.log

module load julia/1.10.5

julia --project examples/calibration_sweep.jl --smoke
