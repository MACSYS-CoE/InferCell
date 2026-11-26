#!/bin/bash
# KIM-185: submit dev/scripts/fullcell_timing_transcript_kim185.jl's stages to Slurm, from a
# clean worktree at a pushed commit.
#   dev/scripts/fullcell_timing_transcript_kim185.sh <mode> <outdir>
# Modes:
#   data       item 1: the 455-gene inputs and counts (one CPU; also the first precompile)
#   serial     item 2: the serial timing (one CPU), after data.jls exists
#   threaded   item 3: the threaded timing (THREADS CPUs, default 32), after data.jls exists
#   project    the projection, after serial.jls and threaded.jls exist
#   all        data, then serial and threaded side by side, then project
# Job ids go to <outdir>/jobs.tsv, logs to <outdir>/logs, each submission's commit to
# <outdir>/COMMIT. SMOKE=1 in the environment passes through to the jobs.
set -euo pipefail
mode=${1:?give the mode}
outdir=$(realpath -m "${2:?give the output directory}")
source "$(dirname "$0")/launcher_guard.sh"
script=dev/scripts/fullcell_timing_transcript_kim185.jl
git rev-parse HEAD >> "$outdir/COMMIT"
smoke=${SMOKE:-0}
threads=${THREADS:-32}
workers=${WORKERS:-}
env=(--export=ALL,SMOKE="$smoke",WORKERS="$workers",OPENBLAS_NUM_THREADS=1,JULIA_PKG_OFFLINE=1)
jl='export EBROOTJULIA=/apps/modules/software/Julia/1.10.5; export PATH="$EBROOTJULIA/bin:$PATH"; julia --project'
after() { [ -n "$1" ] && echo "--dependency=afterok:$1" || true; }

submit() {  # name cpus mem time dependency julia-flags stage
    local j
    j=$(sub --cpus-per-task="$2" --mem="$3" --time="$4" --job-name="kim185_$1" "${env[@]}" $(after "$5") \
            --wrap "$jl $6 $script $7 $outdir")
    record "$1" "$j"
    echo "$j"
}
data() { submit data 1 8G 02:00:00 "${1:-}" "" data; }
serial() { submit serial 1 16G 06:00:00 "${1:-}" "" serial; }
threaded() { submit threaded "$threads" 32G 04:00:00 "${1:-}" "--threads=$threads --gcthreads=$((threads / 2 > 0 ? threads / 2 : 1)),1" threaded; }
project() { submit project 1 8G 01:00:00 "${1:-}" "" project; }

case "$mode" in
    data) data ;;
    serial) serial ;;
    threaded) threaded ;;
    project) project ;;
    all)
        d=$(data)
        s=$(serial "$d")
        t=$(threaded "$d")
        project "$s:$t" ;;
    *) echo "unknown mode $mode" >&2; exit 1 ;;
esac
