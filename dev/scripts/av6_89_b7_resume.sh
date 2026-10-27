#!/bin/bash
# Issue #89 B7: resubmit what an AV6 run is missing, after jobs fail (a prep past its
# wall limit: the bounded MAP's Fminbox runs up to 200 inner iterations per barrier
# round, so a hard replicate's MAP can take far longer than the median's).
#   dev/scripts/av6_89_b7_resume.sh <out>
# Run it from a worktree at the run's own commit, with the run's SMOKE and AV6_REPS,
# once the failed jobs have left the queue. A job of <out>/jobs.tsv counts as live
# only while RUNNING: every PENDING one is taken to wait on a failure, is relabelled
# <fit>~superseded in jobs.tsv (so KC's sacct leaves it out) and is printed, to be
# cancelled. So are the failed ones. Same code, seeds and stages as
# av6_89_b7.sh; only the prep's wall limit (PREP_TIME, default 3 days) is longer.
# DRY=1 prints what it would submit and changes nothing.
set -euo pipefail
outdir=$1
smoke=${SMOKE:-0}
reps=${AV6_REPS:-1:50}
prep_time=${PREP_TIME:-3-00:00:00}
if [ "$smoke" = 1 ]; then chain_time=04:00:00; else chain_time=${CHAIN_TIME:-3-00:00:00}; fi
job=dev/scripts/av6_89_b7.slurm
tsv=$outdir/jobs.tsv
if [ "${DRY:-0}" = 1 ]; then
    cp "$tsv" "${TMPDIR:-/tmp}/av6_resume_dry.tsv"; tsv=${TMPDIR:-/tmp}/av6_resume_dry.tsv
    sbatch() { local a; for a in "$@"; do case $a in --job-name=*|--dependency=*|--time=*) printf '%s ' "$a" >&2;; esac; done
               echo >&2; echo "DRY$RANDOM"; }
fi
sub() { sbatch --parsable --partition=milan --exclude='dave[148-149]' \
               --export=ALL,SMOKE="$smoke",AV6_REPS="$reps",NADAPT="${NADAPT:-}",NDRAW="${NDRAW:-}" \
               --output="$outdir/logs/%x_%j.log" "$@"; }
rec() { printf '%s\t%s\n' "$1" "$2" >> "$tsv"; }

# Live = RUNNING. Everything else not COMPLETED is superseded.
declare -A state
while IFS='|' read -r id st; do state[$id]=$st; done < <(
    sacct -j "$(cut -f2 "$tsv" | paste -sd,)" -X -n -P --format=JobID,State)
live_ids() {   # fit -> its RUNNING job ids, ':'-joined
    awk -F'\t' -v f="$1" '$1 == f {print $2}' "$tsv" |
        while read -r id; do if [ "${state[$id]:-}" = RUNNING ]; then echo "$id"; fi; done | paste -sd:
}
superseded=()
tmp=$(mktemp)
while IFS=$'\t' read -r f id; do
    case "${state[$id]:-}" in
        COMPLETED|RUNNING) printf '%s\t%s\n' "$f" "$id" ;;
        *) printf '%s~superseded\t%s\n' "${f%~superseded}" "$id"
           [[ $f == *~superseded ]] || superseded+=("$id:${state[$id]:-unknown}") ;;
    esac
done < "$tsv" > "$tmp"
mv "$tmp" "$tsv"

have() { [ -f "$outdir/$1/$2" ]; }
ids=$(julia_r=$reps python3 -c '
import os, re
m = re.fullmatch(r"\s*(\d+)\s*:\s*(\d+)\s*", os.environ["julia_r"])
if m: print(" ".join(str(i) for i in range(int(m[1]), int(m[2]) + 1)))
else: print(" ".join(re.findall(r"\d+", os.environ["julia_r"])))')
all_chains=(); base_chains=(); new=0
chains_for() {   # fit, dependency ('' for none): resubmit missing chains, collect ids
    local f=$1 dep=$2 c id run
    run=$(live_ids "$f")
    for c in 1 2; do
        have "$f" "chain_$c.jls" && continue
        if [ -n "$run" ]; then continue; fi
        id=$(sub ${dep:+--dependency=afterok:$dep} --time=$chain_time --job-name=av6_chain_$f $job chain "$outdir" "$f" $c)
        rec "$f" "$id"; all_chains+=("$id"); new=$((new + 1))
        if [ "${3:-}" = base ]; then base_chains+=("$id"); fi
    done
    if [ -n "$run" ]; then
        all_chains+=(${run//:/ })
        if [ "${3:-}" = base ]; then base_chains+=(${run//:/ }); fi
    fi
}
for f in $(for r in $ids; do echo r$r; done) v6; do
    dep=""
    if ! have "$f" prep.jls; then
        if [ -n "$(live_ids "$f")" ]; then chains_for "$f" "" base; continue; fi
        dep=$(sub --time=$prep_time --job-name=av6_prep_$f $job prep "$outdir" "$f")
        rec "$f" "$dep"; new=$((new + 1))
    fi
    chains_for "$f" "$dep" base
done
bdep=$(IFS=:; echo "${base_chains[*]}")
for kf in ks_v6 ks_low1 ks_low2; do
    dep=""
    if ! have "$kf" prep.jls && [ -z "$(live_ids "$kf")" ]; then
        dep=$(sub ${bdep:+--dependency=afterok:$bdep} --time=$prep_time --job-name=av6_ksprep_$kf $job ksprep "$outdir" "$kf")
        rec "$kf" "$dep"; new=$((new + 1))
    fi
    chains_for "$kf" "$dep"
done
for sf in tx_v6 bulk_v6; do chains_for "$sf" ""; done
if [ ! -f "$outdir/av6_summary.md" ] || [ "$new" -gt 0 ]; then
    cdeps=$(IFS=:; echo "${all_chains[*]}")
    m=$(sub --time=02:00:00 ${cdeps:+--dependency=afterok:$cdeps} --job-name=av6_merge $job merge "$outdir")
    rec merge "$m"; new=$((new + 1))
fi
echo "resubmitted $new jobs; superseded (cancel the pending ones): ${superseded[*]:-none}"
