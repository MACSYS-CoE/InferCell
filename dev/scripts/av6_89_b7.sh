#!/bin/bash
# Issue #89 B7: submit every stage of AV6 (spec/phases/approximate-method.md §12
# 2026-10-24, "B7 as planned").
#   dev/scripts/av6_89_b7.sh /fred/oz022/tkimpson/infercell-89-b7-run-<date>
# Run it from a worktree at a pushed commit. SMOKE=1 runs the small version;
# AV6_REPS (a Julia range, default 1:50) picks V7's replicates; NADAPT and NDRAW
# set the chains' length. One job precompiles first, on milan. Per base fit (r*,
# v6): prep, then two chains. KS fits wait for their base's chains; the stream
# arms wait for v6's prep. The merge waits for every chain. Each fit's job ids go
# to <out>/jobs.tsv, which the merge reads for KC's sacct.
set -euo pipefail
outdir=$1
case "$outdir" in */.worktrees/*|*/.claude/worktrees/*) echo "outdir is inside a worktree" >&2; exit 1;; esac
mkdir -p "$outdir/logs"
smoke=${SMOKE:-0}
reps=${AV6_REPS:-1:50}
if [ "$smoke" = 1 ]; then prep_time=02:00:00; chain_time=04:00:00; else prep_time=12:00:00; chain_time=${CHAIN_TIME:-3-00:00:00}; fi
sub() { sbatch --parsable --partition=milan --exclude='dave[148-149]' \
               --export=ALL,SMOKE="$smoke",AV6_REPS="$reps",NADAPT="${NADAPT:-}",NDRAW="${NDRAW:-}" \
               --output="$outdir/logs/%x_%j.log" "$@"; }
job=dev/scripts/av6_89_b7.slurm
rec() { printf '%s\t%s\n' "$1" "$2" >> "$outdir/jobs.tsv"; }
: > "$outdir/jobs.tsv"
pre=$(sub --time=01:00:00 --mem=8G --job-name=av6_pre --wrap 'module load julia/1.10.5; julia --project -e "using InferCell"')
# Julia's range syntax, expanded by Julia's own rules (an empty range is refused).
ids=$(julia_r=$reps python3 -c '
import os, re
m = re.fullmatch(r"\s*(\d+)\s*:\s*(\d+)\s*", os.environ["julia_r"])
if m: print(" ".join(str(i) for i in range(int(m[1]), int(m[2]) + 1)))
else: print(" ".join(re.findall(r"\d+", os.environ["julia_r"])))')
[ -n "$ids" ] || { echo "AV6_REPS=$reps names no replicate" >&2; exit 1; }
declare -A prep chains
all_chains=()
base_fit() {   # fit
    local f=$1 p c1 c2
    p=$(sub --dependency=afterok:$pre --time=$prep_time --job-name=av6_prep_$f $job prep "$outdir" "$f")
    c1=$(sub --dependency=afterok:$p --time=$chain_time --job-name=av6_chain_$f $job chain "$outdir" "$f" 1)
    c2=$(sub --dependency=afterok:$p --time=$chain_time --job-name=av6_chain_$f $job chain "$outdir" "$f" 2)
    rec "$f" "$p"; rec "$f" "$c1"; rec "$f" "$c2"
    prep[$f]=$p; chains[$f]="$c1:$c2"; all_chains+=("$c1" "$c2")
}
for r in $ids; do base_fit "r$r"; done
base_fit v6
# KS: the bases are v6 and the two replicates with the lowest kcatF_R_ENO truths,
# which the driver names (stage lowest_eno) when the job runs.
base_chains=$(IFS=:; echo "${all_chains[*]}")      # the low-ENO bases are known only to Julia
for kf in ks_v6 ks_low1 ks_low2; do
    case $kf in
        ks_v6) dep=${chains[v6]} ;;
        *) dep=$base_chains ;;
    esac
    k=$(sub --dependency=afterok:$dep --time=$prep_time --job-name=av6_ksprep_$kf $job ksprep "$outdir" "$kf")
    c1=$(sub --dependency=afterok:$k --time=$chain_time --job-name=av6_chain_$kf $job chain "$outdir" "$kf" 1)
    c2=$(sub --dependency=afterok:$k --time=$chain_time --job-name=av6_chain_$kf $job chain "$outdir" "$kf" 2)
    rec "$kf" "$k"; rec "$kf" "$c1"; rec "$kf" "$c2"
    all_chains+=("$c1" "$c2")
done
for sf in tx_v6 bulk_v6; do
    c1=$(sub --dependency=afterok:${prep[v6]} --time=$chain_time --job-name=av6_chain_$sf $job chain "$outdir" "$sf" 1)
    c2=$(sub --dependency=afterok:${prep[v6]} --time=$chain_time --job-name=av6_chain_$sf $job chain "$outdir" "$sf" 2)
    rec "$sf" "$c1"; rec "$sf" "$c2"
    all_chains+=("$c1" "$c2")
done
cdeps=$(IFS=:; echo "${all_chains[*]}")
m=$(sub --time=02:00:00 --dependency=afterok:$cdeps --job-name=av6_merge $job merge "$outdir")
rec merge "$m"
echo "precompile $pre; merge $m; $(wc -l < "$outdir/jobs.tsv") jobs in $outdir/jobs.tsv"
