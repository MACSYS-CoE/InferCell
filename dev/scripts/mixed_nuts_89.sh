#!/bin/bash
# Issue #89 B5: submit every stage of the mixed likelihood's smoke on m0_v6.
#   dev/scripts/mixed_nuts_89.sh /fred/oz022/tkimpson/infercell-89-b5-archive-<date>
# SMOKE=1 runs the small version. NADAPT and NDRAW set the chains' length. One job
# precompiles first, on milan (handoff 2026-10-23). Stages chain by afterok, and
# the merge waits for every chain. Prints the job ids.
set -euo pipefail
outdir=$1
case "$outdir" in */.worktrees/*|*/.claude/worktrees/*) echo "outdir is inside a worktree" >&2; exit 1;; esac
mkdir -p "$outdir/logs"
smoke=${SMOKE:-0}
nchunk=$([ "$smoke" = 1 ] && echo 2 || echo 16)
nchain=4
chain_time=${CHAIN_TIME:-2-00:00:00}
sub() { sbatch --parsable --partition=milan --exclude='dave[148-149]' \
               --export=ALL,SMOKE="$smoke",NADAPT="${NADAPT:-}",NDRAW="${NDRAW:-}" \
               --output="$outdir/logs/%x_%j.log" "$@"; }
job=dev/scripts/mixed_nuts_89.slurm
pre=$(sub --time=01:00:00 --mem=8G --wrap 'module load julia/1.10.5; julia --project -e "using InferCell"')
ref=$(sub --dependency=afterok:$pre --time=08:00:00 $job theta_ref "$outdir")
ens=()
for k in $(seq 1 "$nchunk"); do
    ens+=("$(sub --dependency=afterok:$ref --time=12:00:00 $job ensemble "$outdir" "$k")")
done
edeps=$(IFS=:; echo "${ens[*]}")
chains=()
for c in $(seq 1 "$nchain"); do
    chains+=("$(sub --dependency=afterok:$edeps --time="$chain_time" $job nuts "$outdir" "$c")")
done
cdeps=$(IFS=:; echo "${chains[*]}")
m=$(sub --time=00:30:00 --dependency=afterok:$cdeps $job merge "$outdir")
echo "precompile $pre; theta_ref $ref; ensemble ${ens[*]}; chains ${chains[*]}; merge $m"
