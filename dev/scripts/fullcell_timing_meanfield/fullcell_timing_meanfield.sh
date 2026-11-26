#!/bin/bash
# KIM-184: submit the full-cell mean-field timing (fullcell_timing_meanfield.jl), from a
# worktree at a pushed commit.
#   dev/scripts/fullcell_timing_meanfield/fullcell_timing_meanfield.sh <outdir> [what]
# what = all (default): precompile, then check and time (one thread each), then
# threaded (4 threads, after time, whose gradients it compares against).
# what = check | time | threaded: that stage alone, after its own precompile.
# SMOKE=1 in the environment: a plumbing pass (300 s window, 2 repeats).
# Job ids go to <outdir>/jobs.tsv, logs to <outdir>/logs.
set -euo pipefail
outdir=$(realpath -m "$1")
what=${2:-all}
case $what in all|check|time|threaded) ;; *) echo "unknown: $what" >&2; exit 1 ;; esac
# The tree this launcher lives in, wherever it is called from.
cd "$(dirname "$(realpath "$0")")/../../.."
source dev/scripts/launcher_guard.sh
smoke=${SMOKE:-0}
env=dev/scripts/fullcell_timing_meanfield
script=$env/fullcell_timing_meanfield.jl
nthr=4
if [ "$smoke" = 1 ]; then t=01:00:00; else t=12:00:00; fi
jl='export EBROOTJULIA=/apps/modules/software/Julia/1.10.5; export PATH="$EBROOTJULIA/bin:$PATH"; export JULIA_PKG_OFFLINE=1; julia --project='"$env"
ex() { echo "--export=ALL,SMOKE=$smoke,OPENBLAS_NUM_THREADS=1${1:+,$1}"; }

pre=$(sub "$(ex)" --cpus-per-task=1 --mem=16G --time=01:00:00 --job-name=kim184_pre \
          --wrap "$jl -e 'using OrdinaryDiffEq, LinearSolve, SparseArrays, ForwardDiff, SciMLBase'")
record precompile "$pre"
serial() {   # $1 stage, $2 extra dependency
    local id
    id=$(sub "$(ex)" --dependency=afterok:$pre${2:+:$2} --cpus-per-task=1 --mem=16G --time=$t \
             --job-name=kim184_$1 --wrap "$jl $script $1 $outdir")
    record "$1" "$id"
    echo "$id"
}
threaded() {   # $1 extra dependency (afterany)
    local id dep="--dependency=afterok:$pre"
    [ -n "${1:-}" ] && dep="$dep,afterany:$1"
    id=$(sub "$(ex)" "$dep" --cpus-per-task=$nthr --mem=16G --time=$t --job-name=kim184_threaded \
             --wrap "$jl --threads=$nthr --gcthreads=$nthr,1 $script threaded $outdir")
    record threaded "$id"
}
case $what in
all) serial check >/dev/null; tid=$(serial time); threaded "$tid" ;;
check|time) serial "$what" >/dev/null ;;
threaded) threaded ;;
esac
echo "submitted $what at $(git rev-parse --short HEAD), SMOKE=$smoke:"
cat "$outdir/jobs.tsv"
