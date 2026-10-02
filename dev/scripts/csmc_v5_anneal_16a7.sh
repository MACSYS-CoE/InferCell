#!/bin/bash
# Submit V5 for the annealed kernel (task 16a.7a, spec/phases/16-recovery.md §12
# 2026-10-02), for one case (m0 or toy), at K = 5 stages. The rejection draws and
# prep are the unannealed run's: pass the directory holding them (its
# dev/scripts/csmc_v5_16a7). No kernel change touches the forward simulation they
# come from. Run from a worktree at a pushed commit.
set -euo pipefail
case=$1
src=$2
K=${3:-5}
mkdir -p dev/scripts/csmc_v5_16a7/$case
cp "$src"/$case/reject_*.jls "$src"/$case/prep.jls dev/scripts/csmc_v5_16a7/$case/
S=dev/scripts/csmc_v5_16a7.slurm
koff=$(V5_ANNEAL=$K sbatch --parsable --export=ALL,V5_ANNEAL=$K --time=08:00:00 --array=0-49 -J v5a_${case}_off $S kernel $case off)
kon=$(V5_ANNEAL=$K sbatch --parsable --export=ALL,V5_ANNEAL=$K --time=08:00:00 --array=0-49 -J v5a_${case}_on $S kernel $case on)
m=$(sbatch --parsable --export=ALL,V5_ANNEAL=$K --time=00:30:00 -J v5a_${case}_merge --dependency=afterok:$koff:$kon $S merge $case)
echo "$case (K = $K): kernel off $koff, on $kon, merge $m"
