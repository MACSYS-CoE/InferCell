#!/bin/bash
# KIM-181: submit dev/scripts/corea_eigenspectrum/corea_eigenspectrum.jl, from a worktree
# at a pushed commit.
#
#   corea_eigenspectrum.sh <outdir>             the MAP, −H and the spectrum, one job
#   SMOKE=1 corea_eigenspectrum.sh <outdir>     the same at tiny sizes, never a result
#
# Job ids are appended to <outdir>/jobs.tsv, and logs go to <outdir>/logs.
set -euo pipefail
outdir=$(realpath -m "$1")
source "$(dirname "$0")/../launcher_guard.sh"
smoke=${SMOKE:-0}
script=dev/scripts/corea_eigenspectrum/corea_eigenspectrum.jl
ex="--export=ALL,SMOKE=$smoke"
if [ "$smoke" = 1 ]; then
    t=03:00:00
else
    # KIM-178's 21-coordinate fitted MAP took 5187 s and its −H 1872 s (job 21727400).
    # At 29 coordinates a gradient costs more, and −H has 29 columns: about 5 h expected.
    t=1-00:00:00
fi

pre=$(sub "$ex" --cpus-per-task=1 --mem=16G --time=01:00:00 --job-name=kim181_pre \
          --wrap "$run -e 'using InferCell'")
record precompile "$pre"
git rev-parse HEAD >> "$outdir/COMMIT"
record eigenspectrum "$(sub "$ex" --dependency=afterok:$pre --cpus-per-task=1 --mem=32G --time=$t \
                          --job-name=kim181_eig --wrap "$run $script $outdir")"
cat "$outdir/jobs.tsv"
