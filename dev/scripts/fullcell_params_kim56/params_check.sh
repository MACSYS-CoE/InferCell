#!/bin/bash
# KIM-56: submits the three jobs of params_check.slurm (precompile, then the name table
# and the bitwise check) on EPYC 7543 nodes. Submit from the worktree root at a pushed
# commit:
#
#   dev/scripts/fullcell_params_kim56/params_check.sh <outdir>
#
# <outdir> must not hold another run's jobs.tsv, COMMIT or results. dev/scripts/launcher_guard.sh refuses an
# outdir inside a worktree, a dirty tree and an unpushed HEAD, and makes <outdir>/logs.
# Job ids go to <outdir>/jobs.tsv, the commit to <outdir>/COMMIT, the name table to
# <outdir>/table, the check to <outdir>/check, and the run's compile cache to
# <outdir>/depot.
set -euo pipefail
outdir=$(realpath -m "$1")
for f in check table jobs.tsv COMMIT; do
    [ ! -e "$outdir/$f" ] || { echo "$outdir already holds $f from another run; use a new outdir" >&2; exit 1; }
done
source "$(dirname "$0")/../launcher_guard.sh"
# KIM-183's node set: the EPYC 7543 nodes, so the gina nodes are excluded as well.
X="gina[1-18,301-305],$excluded"
slurm=dev/scripts/fullcell_params_kim56/params_check.slurm
common=(--partition=milan --exclude="$X" --ntasks=1 --cpus-per-task=1 --mem=8G
        --output="$outdir/logs/%x_%j.log")
git rev-parse HEAD >> "$outdir/COMMIT"
pre=$(sbatch --parsable "${common[@]}" --time=01:00:00 --job-name=kim56-pre \
             --export=ALL,MODE=precompile,DEPOT="$outdir/depot" "$slurm")
record precompile "$pre"
tab=$(sbatch --parsable "${common[@]}" --dependency=afterok:$pre --time=00:30:00 --job-name=kim56-table \
             --export=ALL,MODE=table,OUT="$outdir/table",DEPOT="$outdir/depot" "$slurm")
record table "$tab"
chk=$(sbatch --parsable "${common[@]}" --dependency=afterok:$pre --time=03:00:00 --job-name=kim56-check \
             --export=ALL,MODE=check,PROCIDS=1-6,SEED=1,MINUTES=181,OUT="$outdir/check",DEPOT="$outdir/depot" \
             "$slurm")
record check "$chk"
cat "$outdir/jobs.tsv"
