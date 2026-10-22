#!/bin/bash
# Submit 16a.10's production arm for F10 (spec/phases/16-recovery.md §12
# 2026-10-13): one fresh chain on V6's dataset at production length, run by the
# pilot's driver as chain 5 at 5 threads. Run from a worktree at a pushed commit.
# The chain resumes from its checkpoint, so one that hits its time limit is
# resubmitted with the same command. No merge is chained: F10 reads the chain
# once it has all its sweeps.
# Smoke run: PILOT_DIR, PILOT_NCELLS and PILOT_HORIZON are passed through, with
# SWEEPS and TIME.
set -euo pipefail
export PILOT_DIR=${PILOT_DIR:-$PWD/dev/scripts/f10_16a10}   # its own: a V6 merge globs chain*.jls
mkdir -p "$PILOT_DIR"
SWEEPS=${SWEEPS:-656}
TIME=${TIME:-6-00:00:00}
j=$(sbatch --parsable --export=ALL --time=$TIME --cpus-per-task=5 --mem=8G -J f10_c5 \
    --output="$PILOT_DIR/slurm_%x_%j.log" dev/scripts/gibbs_pilot_16a9a.slurm run 5 $SWEEPS)
echo "F10 production chain $j, $SWEEPS sweeps, in $PILOT_DIR"
