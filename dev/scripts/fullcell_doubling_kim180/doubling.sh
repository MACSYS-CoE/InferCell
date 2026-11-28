#!/bin/bash
# KIM-180: the default-run check of the doubling-time definition frozen in
# spec/prereg/kim-180.md §1. It runs KIM-183's driver (dev/scripts/fullcell_timing_cycle/
# timing_cycle.jl, unchanged) on 100 default cells (procids 1–100, seed 1) for 180
# simulated minutes, then doubling_reads.py applies the definition. Submit from the
# worktree root at a pushed commit:
#
#   dev/scripts/fullcell_doubling_kim180/doubling.sh <outdir>
#
# Job ids are appended to <outdir>/jobs.tsv, and logs go to <outdir>/logs. The run's
# compile cache goes in <outdir>/depot, stacked on KIM-183's and the shared depot, so
# nothing is written to another issue's archive.
set -euo pipefail
outdir=$(realpath -m "$1")
source "$(dirname "$0")/../launcher_guard.sh"
# KIM-183's node set: the EPYC 7543 nodes, so the gina nodes are excluded as well.
X="gina[1-18,301-305],$excluded"
slurm=dev/scripts/fullcell_doubling_kim180/doubling.slurm
common=(--partition=milan --exclude="$X" --ntasks=1 --output="$outdir/logs/%x_%j.log")
git rev-parse HEAD >> "$outdir/COMMIT"
pre=$(sbatch --parsable "${common[@]}" --cpus-per-task=1 --mem=8G --time=01:00:00 --job-name=kim180-pre \
             --export=ALL,MODE=precompile,DEPOT="$outdir/depot" "$slurm")
record precompile "$pre"
ens=$(sbatch --parsable "${common[@]}" --dependency=afterok:$pre --nodes=1 --cpus-per-task=32 --mem=100G \
             --time=02:00:00 --job-name=kim180-ensemble \
             --export=ALL,MODE=ensemble,PROCIDS=1-100,SEED=1,MINUTES=181,OUT="$outdir/ensemble/ensemble",DEPOT="$outdir/depot" \
             "$slurm")
record ensemble "$ens"
cat "$outdir/jobs.tsv"
