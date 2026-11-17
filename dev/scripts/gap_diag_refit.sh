#!/bin/bash
# KIM-90: submit the m0_v6 refits (items 3, 4a) and the cross-covariance check (4b).
#   dev/scripts/gap_diag_refit.sh <outdir outside the worktrees> <arm>...
# arms: control n1600 nojensen postmean (gap_diag_refit.jl). Run from a worktree
# at a pushed commit. SMOKE=1 runs the small version. Per arm: prep, then two
# chains; the merge waits for every chain. The cross-covariance ensemble runs
# beside them and its analysis after it. Job ids go to <outdir>/jobs.tsv.
set -euo pipefail
outdir=$1; shift
case "$outdir" in */.worktrees/*|*/.claude/worktrees/*) echo "outdir is inside a worktree" >&2; exit 1;; esac
[ $# -gt 0 ] || { echo "name at least one arm" >&2; exit 1; }
mkdir -p "$outdir/logs"
smoke=${SMOKE:-0}
if [ "$smoke" = 1 ]; then prep_time=01:00:00; chain_time=02:00:00; else prep_time=06:00:00; chain_time=${CHAIN_TIME:-1-00:00:00}; fi
sub() { sbatch --parsable --partition=milan --exclude='dave[148-149]' --export=ALL,SMOKE="$smoke" \
               --output="$outdir/logs/%x_%j.log" "$@"; }
job=dev/scripts/gap_diag_refit.slurm
rec() { printf '%s\t%s\n' "$1" "$2" >> "$outdir/jobs.tsv"; }
: > "$outdir/jobs.tsv"
all=()
for a in "$@"; do
    p=$(sub --time=$prep_time --job-name=kim90_prep_$a $job gap_diag_refit.jl prep "$outdir" "$a")
    c1=$(sub --dependency=afterok:$p --time=$chain_time --job-name=kim90_chain_$a $job gap_diag_refit.jl chain "$outdir" "$a" 1)
    c2=$(sub --dependency=afterok:$p --time=$chain_time --job-name=kim90_chain_$a $job gap_diag_refit.jl chain "$outdir" "$a" 2)
    rec "$a" "$p"; rec "$a" "$c1"; rec "$a" "$c2"; all+=("$c1" "$c2")
done
m=$(sub --dependency=afterok:$(IFS=:; echo "${all[*]}") --time=00:30:00 --job-name=kim90_merge $job gap_diag_refit.jl merge "$outdir")
rec merge "$m"
e=$(sub --time=$prep_time --job-name=kim90_crossens $job gap_diag_crosscov.jl ensemble "$outdir")
x=$(sub --dependency=afterok:$e --time=01:00:00 --job-name=kim90_crosscov $job gap_diag_crosscov.jl analyse "$outdir")
rec crossens "$e"; rec crosscov "$x"
cat "$outdir/jobs.tsv"
