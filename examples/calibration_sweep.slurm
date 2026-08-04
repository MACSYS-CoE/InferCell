#!/bin/bash
#SBATCH --job-name=infercell_calib
#SBATCH --partition=milan
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=1
#SBATCH --mem=4G
#SBATCH --time=04:00:00
#SBATCH --array=0-3
#SBATCH --output=examples/calib_sweep_%A_%a.log

module load julia/1.10.5

julia --project examples/calibration_sweep.jl
