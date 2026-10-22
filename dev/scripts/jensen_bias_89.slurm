#!/bin/bash
#SBATCH --job-name=infercell_jbias89
#SBATCH --partition=milan
#SBATCH --exclude=dave[148-149]
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=1
#SBATCH --mem=16G
#SBATCH --time=01:00:00

# Issue #89, B3 close-out: the implied-bias report, read from the rethink's archive.
#   sbatch --output=<outdir>/%x_%j.log dev/scripts/jensen_bias_89.slurm <outdir outside the worktrees>
case "$1" in */.worktrees/*|*/.claude/worktrees/*) echo "outdir is inside a worktree" >&2; exit 1;; esac
module load julia/1.10.5
echo "#89 bias $* at $(date) on $(hostname), HEAD $(git rev-parse --short HEAD)$(test -z "$(git status --porcelain -- src dev/scripts)" || echo ' (tree dirty)')"
julia --project dev/scripts/jensen_bias_89.jl "$1"
status=$?
echo "Exit code: $status  Finished at $(date)"
exit $status
