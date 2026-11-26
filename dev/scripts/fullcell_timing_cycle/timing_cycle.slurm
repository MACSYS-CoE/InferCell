#!/bin/bash
# KIM-183: SpeedyCell (syn3a-ensemble copy) full-cycle timing on milan EPYC 7543 nodes.
# Submit from the worktree root at a pushed commit (full commands: the result note's
# Reproduce section), with the node flags of
# dev/scripts/launcher_guard.sh extended to syn3a-ensemble's EPYC 7543 set:
#
#   L=/fred/oz022/tkimpson/infercell-kim183-2026-11-26
#   X="gina[1-18,301-305],$excluded"     # $excluded as dev/scripts/launcher_guard.sh builds it
#   precompile: sbatch --partition=milan --exclude="$X" --cpus-per-task=1 --mem=8G --time=01:00:00 \
#       --output=$L/logs/%x_%j.log --export=ALL,MODE=precompile dev/scripts/fullcell_timing_cycle/timing_cycle.slurm
#   serial:     ... --cpus-per-task=1 --mem=8G --export=ALL,MODE=serial,PROCIDS=1-20,OUT=$L/serial/serial ...
#   ensemble:   ... --nodes=1 --cpus-per-task=64 --mem=200G --export=ALL,MODE=ensemble,PROCIDS=1-100,OUT=$L/ensemble/ensemble ...
#
# Julia 1.11.9 (SpeedyCell needs 1.11); the depot is this issue's own, stacked on the
# shared one that already holds LowPrecSSA (unregistered, SSH-pinned at cb0abb6).
#SBATCH --job-name=kim183-cycle
#SBATCH --ntasks=1
#SBATCH --time=02:00:00

set -euo pipefail
lscpu | grep -q "AMD EPYC 7543" || { echo "not an EPYC 7543 node: $(hostname)" >&2; exit 1; }
echo "host $(hostname)"; lscpu | grep 'Model name'
export PATH=/fred/oz022/tkimpson/julia-1.11.9/bin:$PATH
export JULIA_DEPOT_PATH=/fred/oz022/tkimpson/infercell-kim183-2026-11-26/depot:/fred/oz022/tkimpson/julia_depot
export JULIA_PKG_OFFLINE=true
export OMP_NUM_THREADS=1 OPENBLAS_NUM_THREADS=1 MKL_NUM_THREADS=1
# sbatch copies the script, so $0 is not in the repo; use the submit directory.
ENV_DIR=$SLURM_SUBMIT_DIR/dev/scripts/fullcell_timing_cycle
[ -f "$ENV_DIR/timing_cycle.jl" ] || { echo "submit from the worktree root" >&2; exit 1; }
echo "git $(git -C "$ENV_DIR" rev-parse HEAD)"
NT=${SLURM_CPUS_PER_TASK:-1}
case "${MODE:-}" in
  precompile)
    julia --project="$ENV_DIR" -e 'using Pkg; Pkg.precompile(); using SpeedyCell.Cell, InteractiveUtils; versioninfo()' ;;
  serial)
    julia -t 1 --project="$ENV_DIR" "$ENV_DIR/timing_cycle.jl" ;;
  ensemble)
    julia -t "$NT" --gcthreads="$NT,1" --project="$ENV_DIR" "$ENV_DIR/timing_cycle.jl" ;;
  *) echo "MODE must be precompile, serial or ensemble" >&2; exit 1 ;;
esac
