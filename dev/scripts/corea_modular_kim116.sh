#!/bin/bash
# KIM-116: submit dev/scripts/corea_modular_kim116.jl's stages to Slurm, from a clean
# worktree at a pushed commit.
#   dev/scripts/corea_modular_kim116.sh <mode> <outdir>
# Modes:
#   profile   item 2: one job, each stage's density and gradient (profile.md)
# Job ids go to <outdir>/jobs.tsv, logs to <outdir>/logs, and the commit to <outdir>/COMMIT.
# SMOKE=1 in the environment passes through to the job.
set -euo pipefail
mode=${1:?give the mode}
outdir=$(realpath -m "${2:?give the output directory}")
source "$(dirname "$0")/launcher_guard.sh"
script=dev/scripts/corea_modular_kim116.jl
git rev-parse HEAD >> "$outdir/COMMIT"
smoke=${SMOKE:-0}
case "$mode" in
profile)
    [ -e "$outdir/profile.md" ] && { echo "$outdir/profile.md already exists" >&2; exit 1; }
    j=$(sub --cpus-per-task=1 --mem=16G --time=04:00:00 --job-name=kim116_profile \
            --export=ALL,SMOKE="$smoke" --wrap "$run $script profile $outdir")
    record profile "$j" ;;
*)
    echo "unknown mode $mode" >&2; exit 1 ;;
esac
echo "submitted at $(git rev-parse --short HEAD):"
cat "$outdir/jobs.tsv"
