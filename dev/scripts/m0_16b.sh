#!/bin/bash
# Submit 16b.1's and 16b.4's production chains at M0 (spec/phases/16-recovery.md
# §12 "16a.11" and "16b at M0"), side by side, one per mode, at F10's settings.
# Run from a worktree at a pushed commit. A chain resumes from its checkpoint, so
# one that hits its time limit is resubmitted with the same command.
# MODES picks a subset. A smoke run passes ROOT, PILOT_NCELLS, PILOT_HORIZON,
# SWEEPS and TIME through.
set -euo pipefail
ROOT=${ROOT:-$PWD/dev/scripts}
SWEEPS=${SWEEPS:-656}
TIME=${TIME:-6-00:00:00}
for mode in ${MODES:-perturbed tx_only met_only}; do
    export RUN_MODE=$mode PILOT_DIR=$ROOT/m0_16b_$mode
    mkdir -p "$PILOT_DIR"
    j=$(sbatch --parsable --export=ALL --time=$TIME -J m0_16b_$mode \
        --output="$PILOT_DIR/slurm_%x_%j.log" dev/scripts/m0_16b.slurm run $SWEEPS)
    echo "16b $mode chain $j, $SWEEPS sweeps, in $PILOT_DIR"
done
