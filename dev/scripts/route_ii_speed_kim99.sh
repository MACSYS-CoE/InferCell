#!/bin/bash
# KIM-99: submit route (ii)'s profile stages (dev/scripts/route_ii_speed_kim99.jl), from
# a worktree at a pushed commit.
#   dev/scripts/route_ii_speed_kim99.sh <outdir> profile   precompile, then split, rhs
#                                                          and variants, one job each
#   dev/scripts/route_ii_speed_kim99.sh <outdir> grad      precompile, then grad
#   dev/scripts/route_ii_speed_kim99.sh <outdir> all       both
# SMOKE=1 in the environment: a plumbing pass at tiny sizes. Job ids go to
# <outdir>/jobs.tsv, logs to <outdir>/logs.
set -euo pipefail
outdir=$(realpath -m "$1")
what=${2:?say what to submit: profile, grad or all}
case $what in profile|grad|all) ;; *) echo "unknown: $what" >&2; exit 1 ;; esac
# The tree this launcher lives in, wherever it is called from.
cd "$(dirname "$(realpath "$0")")/../.."
source dev/scripts/launcher_guard.sh
smoke=${SMOKE:-0}
script=dev/scripts/route_ii_speed_kim99.jl
ex="--export=ALL,SMOKE=$smoke"
if [ "$smoke" = 1 ]; then t=01:00:00; t_grad=01:00:00; else t=12:00:00; t_grad=1-00:00:00; fi

pre=$(sub "$ex" --cpus-per-task=1 --mem=16G --time=01:00:00 --job-name=kim99_pre \
          --wrap "$run -e 'using InferCell'")
record precompile "$pre"
stage() {   # $1 stage, $2 wall
    local id
    id=$(sub "$ex" --dependency=afterok:$pre --cpus-per-task=1 --mem=32G --time=$2 \
             --job-name=kim99_$1 --wrap "$run $script $1 $outdir")
    record "$1" "$id"
}
case $what in
profile) stage split $t; stage rhs $t; stage variants $t ;;
grad) stage grad $t_grad ;;
all) stage split $t; stage rhs $t; stage variants $t; stage grad $t_grad ;;
*) echo "unknown: $what" >&2; exit 1 ;;
esac
echo "submitted $what at $(git rev-parse --short HEAD), SMOKE=$smoke:"
cat "$outdir/jobs.tsv"
