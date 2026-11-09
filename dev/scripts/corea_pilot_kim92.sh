#!/bin/bash
# KIM-92: submit the reduced-model pilot's stages (dev/scripts/corea_pilot_kim92.jl),
# from a worktree at a pushed commit.
#   dev/scripts/corea_pilot_kim92.sh <outdir> profile
# precompiles on milan, draws the dataset, then profiles the density and gradient
# (item 2). With SMOKE=1 in the environment it is a plumbing pass at tiny sizes.
# Job ids go to <outdir>/jobs.tsv, logs to <outdir>/logs.
set -euo pipefail
outdir=$(realpath -m "$1")
what=${2:?say what to submit: profile}
source "$(dirname "$0")/launcher_guard.sh"
smoke=${SMOKE:-0}
script=dev/scripts/corea_pilot_kim92.jl
ex="--export=ALL,SMOKE=$smoke"
if [ "$smoke" = 1 ]; then t_ds=01:00:00; t_prof=02:00:00; else t_ds=12:00:00; t_prof=2-00:00:00; fi

case $what in
profile)
    pre=$(sub "$ex" --cpus-per-task=1 --mem=16G --time=01:00:00 --job-name=kim92_pre \
              --wrap "$run -e 'using InferCell'")
    record precompile "$pre"
    ds=$(sub "$ex" --dependency=afterok:$pre --cpus-per-task=1 --mem=32G --time=$t_ds \
             --job-name=kim92_dataset --wrap "$run $script dataset $outdir")
    record dataset "$ds"
    pr=$(sub "$ex" --dependency=afterok:$ds --cpus-per-task=1 --mem=32G --time=$t_prof \
             --job-name=kim92_profile --wrap "$run $script profile $outdir")
    record profile "$pr"
    ;;
*)
    echo "unknown: $what" >&2; exit 1 ;;
esac
echo "submitted $what at $(git rev-parse --short HEAD), SMOKE=$smoke:"
cat "$outdir/jobs.tsv"
