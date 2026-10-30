# Sourced by a dev/scripts launcher once it has set `outdir` (KIM-68):
#   outdir=$(realpath -m "$1"); source "$(dirname "$0")/launcher_guard.sh"
# Before anything is made, it refuses an output directory inside any of the repo's
# worktrees (the main checkout included), a dirty tree, and a HEAD on no remote
# branch. Only then does it make $outdir/logs. It also defines:
#   sub ...           sbatch on milan (dave148–149 excluded), one task, logging to
#                     $outdir/logs, printing the job id;
#   record name id    appends "name<TAB>id" to $outdir/jobs.tsv;
#   run               the command prefix for Julia in a job's --wrap, with BLAS on
#                     one thread (the site's profile sets OMP_NUM_THREADS=1, which
#                     every KIM-68 run inherited; this does not rely on it).
# The M4 launchers carry their own copy of the checks.
[ -n "${outdir:-}" ] || { echo "set outdir before sourcing launcher_guard.sh" >&2; exit 1; }
case "$outdir/" in */.worktrees/*|*/.claude/worktrees/*)
    echo "outdir $outdir is inside a worktree" >&2; exit 1;; esac
while read -r wt; do
    case "$outdir/" in "$wt"/*) echo "outdir $outdir is inside the worktree $wt" >&2; exit 1;; esac
done < <(git worktree list --porcelain | sed -n 's/^worktree //p')
# `:/`: from the tree's top, wherever the launcher is run from.
if [ -n "$(git status --porcelain -- :/src :/test :/dev/scripts :/Project.toml :/Manifest.toml)" ]; then
    echo "the tree is dirty; run from a clean worktree at a pushed commit" >&2; exit 1
fi
if [ -z "$(git branch -r --contains HEAD)" ]; then
    echo "HEAD $(git rev-parse --short HEAD) is on no remote branch; push it first" >&2; exit 1
fi
mkdir -p "$outdir/logs"
sub() { sbatch --parsable --partition=milan --exclude='dave[148-149]' --ntasks=1 \
               --output="$outdir/logs/%x_%j.log" "$@"; }
record() { printf '%s\t%s\n' "$1" "$2" >> "$outdir/jobs.tsv"; }
run="module load julia/1.10.5; OPENBLAS_NUM_THREADS=1 julia --project"
