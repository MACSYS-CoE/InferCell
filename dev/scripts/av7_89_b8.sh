#!/bin/bash
# Issue #89 B8: submit every stage of AV7 (spec/phases/approximate-method.md §12
# 2026-10-29, "B8 as planned").
#   dev/scripts/av7_89_b8.sh /fred/oz022/tkimpson/infercell-89-b8-run-<date> [<bound dir>]
# Run it from a worktree at a pushed commit. SMOKE=1 runs the small version;
# AV7_REPS (a Julia range, default 1:200) picks the injections; NADAPT and NDRAW
# set the chains' length; PREP_TIME, ENS_TIME and CHAIN_TIME the limits (set from
# B8.4's pilot, at 3× its measured times). Per base fit: prep, the ensemble as an
# array of NCHUNK chunks, ref, then two chains. The stream arms on r = 1–20 wait
# for their base's ref. KS fits and the matched nulls wait for every base chain
# (their bases are chosen by the driver from the drawn truths). The merge waits
# for every chain, and reads <bound dir> if given. Job ids go to <out>/jobs.tsv,
# which the merge reads for the cost.
set -euo pipefail
outdir=$1; bound=${2:-}
case "$outdir" in */.worktrees/*|*/.claude/worktrees/*) echo "outdir is inside a worktree" >&2; exit 1;; esac
mkdir -p "$outdir/logs"
smoke=${SMOKE:-0}
reps=${AV7_REPS:-1:200}
if [ "$smoke" = 1 ]; then
    nchunk=2; prep_time=02:00:00; ens_time=02:00:00; chain_time=04:00:00
else
    nchunk=16
    prep_time=${PREP_TIME:?set PREP_TIME from B8.4}; ens_time=${ENS_TIME:?set ENS_TIME from B8.4}
    chain_time=${CHAIN_TIME:?set CHAIN_TIME from B8.4}
fi
sub() { sbatch --parsable --partition=milan --exclude='dave[148-149]' \
               --export=ALL,SMOKE="$smoke",AV7_REPS="$reps",NADAPT="${NADAPT:-}",NDRAW="${NDRAW:-}" \
               --output="$outdir/logs/%x_%A_%a.log" "$@"; }
job=dev/scripts/av7_89_b8.slurm
rec() { printf '%s\t%s\n' "$1" "$2" >> "$outdir/jobs.tsv"; }
: > "$outdir/jobs.tsv"
pre=$(sub --time=01:00:00 --mem=8G --job-name=av7_pre --wrap 'module load julia/1.10.5; julia --project -e "using InferCell"')
ids=$(julia_r=$reps python3 -c '
import os, re
m = re.fullmatch(r"\s*(\d+)\s*:\s*(\d+)\s*", os.environ["julia_r"])
if m: print(" ".join(str(i) for i in range(int(m[1]), int(m[2]) + 1)))
else: print(" ".join(re.findall(r"\d+", os.environ["julia_r"])))')
[ -n "$ids" ] || { echo "AV7_REPS=$reps names no injection" >&2; exit 1; }
declare -A ref
all_chains=()
# An ensemble (array), its ref, and two chains, after dependency $2.
ensemble_and_chains() {   # fit dep
    local f=$1 dep=$2 e r c1 c2
    e=$(sub --dependency=afterok:$dep --array=1-$nchunk --time=$ens_time --job-name=av7_ens_$f $job ens "$outdir" "$f" @K)
    r=$(sub --dependency=afterok:$e --time=01:00:00 --mem=32G --job-name=av7_ref_$f $job ref "$outdir" "$f")
    c1=$(sub --dependency=afterok:$r --time=$chain_time --job-name=av7_chain_$f $job chain "$outdir" "$f" 1)
    c2=$(sub --dependency=afterok:$r --time=$chain_time --job-name=av7_chain_$f $job chain "$outdir" "$f" 2)
    rec "$f" "$e"; rec "$f" "$r"; rec "$f" "$c1"; rec "$f" "$c2"
    ref[$f]=$r; all_chains+=("$c1" "$c2")
}
chains_only() {   # fit dep
    local f=$1 dep=$2 c1 c2
    c1=$(sub --dependency=afterok:$dep --time=$chain_time --job-name=av7_chain_$f $job chain "$outdir" "$f" 1)
    c2=$(sub --dependency=afterok:$dep --time=$chain_time --job-name=av7_chain_$f $job chain "$outdir" "$f" 2)
    rec "$f" "$c1"; rec "$f" "$c2"
    arm_chains+=("$c1" "$c2")
}
for r in $ids; do
    p=$(sub --dependency=afterok:$pre --time=$prep_time --job-name=av7_prep_r$r $job prep "$outdir" "r$r")
    rec "r$r" "$p"
    ensemble_and_chains "r$r" "$p"
done
base_chains=$(IFS=:; echo "${all_chains[*]}")
arm_chains=()
for r in $ids; do
    [ "$r" -le 20 ] || continue
    chains_only "tx_r$r" "${ref[r$r]}"
    chains_only "bulk_r$r" "${ref[r$r]}"
done
for kf in ks_ptsg_lo ks_eno_lo ks_ptsg_mid; do ensemble_and_chains "$kf" "$base_chains"; done
for nf in null_ptsg_lo null_eno_lo null_ptsg_mid; do chains_only "$nf" "$base_chains"; done
every=("${all_chains[@]}" "${arm_chains[@]}")
cdeps=$(IFS=:; echo "${every[*]}")
m=$(sub --time=04:00:00 --mem=32G --dependency=afterok:$cdeps --job-name=av7_merge $job merge "$outdir" $bound)
rec merge "$m"
echo "precompile $pre; merge $m; $(wc -l < "$outdir/jobs.tsv") jobs in $outdir/jobs.tsv"
