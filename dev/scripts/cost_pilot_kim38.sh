#!/bin/bash
# KIM-38: cost the evidence package (dev/scripts/cost_pilot_kim38.jl), from a
# worktree at a pushed commit.
#   dev/scripts/cost_pilot_kim38.sh /fred/oz022/tkimpson/infercell-m4-kim38-<date>
#       precompiles on milan, then submits cells, threads (16 and 64) and sizing;
#   dev/scripts/cost_pilot_kim38.sh <outdir> sizing_resolved
#       submits the resolved M0 sizing (added after the PR #123 review);
#   dev/scripts/cost_pilot_kim38.sh <outdir> project <kim35_outdir>
#       submits the projection, once those and KIM-35's arms have finished.
# Job ids go to <outdir>/jobs.tsv, logs to <outdir>/logs.
set -euo pipefail
outdir=$(realpath -m "$1")
case "$outdir" in "$(git rev-parse --show-toplevel)"/*|*/.worktrees/*|*/.claude/worktrees/*)
    echo "outdir $outdir is inside a worktree" >&2; exit 1;; esac
if [ -n "$(git status --porcelain -- src test dev/scripts Project.toml Manifest.toml)" ]; then
    echo "the tree is dirty; run from a clean worktree at a pushed commit" >&2; exit 1
fi
if [ -z "$(git branch -r --contains HEAD)" ]; then
    echo "HEAD $(git rev-parse --short HEAD) is on no remote branch; push it first" >&2; exit 1
fi
mkdir -p "$outdir/logs"
sub() { sbatch --parsable --partition=milan --exclude='dave[148-149]' --ntasks=1 \
               --output="$outdir/logs/%x_%j.log" "$@"; }
run="module load julia/1.10.5; julia --project"
if [ "${2:-}" = sizing_resolved ]; then
    # Needs `cells` from the first run in the same outdir (its per-run times).
    j=$(sub --cpus-per-task=16 --mem=48G --time=02:00:00 --job-name=kim38_sizing_resolved \
            --wrap "$run --threads=16 dev/scripts/cost_pilot_kim38.jl sizing_resolved $outdir")
    printf 'sizing_resolved\t%s\n' "$j" >> "$outdir/jobs.tsv"
    cat "$outdir/jobs.tsv"; exit 0
fi
if [ "${2:-}" = project ]; then
    k35=$(realpath -m "$3")
    j=$(sub --cpus-per-task=1 --mem=8G --time=00:30:00 --job-name=kim38_project \
            --wrap "$run dev/scripts/cost_pilot_kim38.jl project $outdir $k35")
    printf 'project\t%s\n' "$j" >> "$outdir/jobs.tsv"
    cat "$outdir/jobs.tsv"; exit 0
fi
: > "$outdir/jobs.tsv"
pre=$(sub --cpus-per-task=1 --mem=16G --time=01:00:00 --job-name=kim38_pre \
          --wrap "$run -e 'using InferCell'")
printf 'precompile\t%s\n' "$pre" >> "$outdir/jobs.tsv"
j=$(sub --dependency=afterok:$pre --cpus-per-task=1 --mem=16G --time=03:00:00 --job-name=kim38_cells \
        --wrap "$run dev/scripts/cost_pilot_kim38.jl cells $outdir")
printf 'cells\t%s\n' "$j" >> "$outdir/jobs.tsv"
for n in 16 64; do
    j=$(sub --dependency=afterok:$pre --cpus-per-task=$n --mem=$((n * 3))G --time=02:00:00 \
            --job-name=kim38_threads$n \
            --wrap "$run --threads=$n dev/scripts/cost_pilot_kim38.jl threads $outdir")
    printf 'threads%s\t%s\n' "$n" "$j" >> "$outdir/jobs.tsv"
done
j=$(sub --dependency=afterok:$pre --cpus-per-task=8 --mem=32G --time=02:00:00 --job-name=kim38_sizing \
        --wrap "$run --threads=8 dev/scripts/cost_pilot_kim38.jl sizing $outdir")
printf 'sizing\t%s\n' "$j" >> "$outdir/jobs.tsv"
echo "submitted at $(git rev-parse --short HEAD):"
cat "$outdir/jobs.tsv"
