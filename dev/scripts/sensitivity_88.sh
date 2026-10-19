#!/bin/bash
# Issue #88 A1/A2: submit the ensemble sensitivity as an array of NTASKS tasks.
#   dev/scripts/sensitivity_88.sh m0    1470   49 /fred/oz022/tkimpson/infercell-88-a1-archive-<date>
#   dev/scripts/sensitivity_88.sh corea 1470 1470 /fred/oz022/tkimpson/infercell-88-a2-archive-<date>
# Set SMOKE=1 in the environment for a 120 s smoke run. Extra arguments go to
# the array's sbatch (e.g. --time). The merge is submitted after the array
# (afterany). Prints both job ids.
set -euo pipefail
model=$1; nseeds=$2; ntasks=$3; outdir=$4; shift 4
case "$outdir" in */.worktrees/*|*/.claude/worktrees/*) echo "outdir is inside a worktree" >&2; exit 1;; esac
mkdir -p "$outdir"
a=$(sbatch --parsable --array=0-$((ntasks - 1)) \
    --export=ALL,MODEL="$model",NSEEDS="$nseeds",NTASKS="$ntasks",OUTDIR="$outdir",SMOKE="${SMOKE:-0}" \
    "$@" dev/scripts/sensitivity_88.slurm)
m=$(sbatch --parsable --dependency=afterany:"$a" dev/scripts/sensitivity_88_merge.slurm "$outdir")
echo "array $a, merge $m"
