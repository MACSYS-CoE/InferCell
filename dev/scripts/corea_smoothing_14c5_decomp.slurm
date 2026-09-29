#!/bin/bash
#SBATCH --job-name=infercell_decomp_14c5
#SBATCH --partition=milan
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=1
#SBATCH --mem=8G
#SBATCH --time=01:00:00
#SBATCH --array=0-3
#SBATCH --output=dev/scripts/corea_smoothing_14c5_decomp_%A_%a.log

# The 2×2 decomposition of 14c.5's coupled-half gap (spec §12 2026-09-29), on the
# 5,000-seed run's four worst seeds, each to its decoupling. Analyse afterwards
# with
#   julia --project dev/scripts/corea_smoothing_14c5_decomp_merge.jl dev/scripts/decomp_14c5
module load julia/1.10.5
mkdir -p dev/scripts/decomp_14c5
SEEDS=(18227 16184 19958 16391); TS=(4515 6110 5477 3315)
i=${SLURM_ARRAY_TASK_ID}
echo "Core A′ 14c.5 decomposition, seed ${SEEDS[$i]} at $(date) on $(hostname), HEAD $(git rev-parse --short HEAD)"
julia --project dev/scripts/corea_smoothing_14c5_decomp.jl ${SEEDS[$i]} ${TS[$i]} dev/scripts/decomp_14c5/decomp_${SEEDS[$i]}.jls
