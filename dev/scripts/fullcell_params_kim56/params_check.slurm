#!/bin/bash
# KIM-56: the parameter-override hook of syn3a-ensemble's SpeedyCell, checked from
# InferCell. syn3a-ensemble is read-only here and pinned: the job refuses to run unless
# it is at SYN3A_PIN with code/SpeedyCell and data/model clean, because this directory's
# environment (copied from dev/scripts/fullcell_timing_cycle/, KIM-183) develops
# /fred/oz022/tkimpson/syn3a-ensemble/code/SpeedyCell by path.
#
# Submit from the worktree root at a pushed commit, on an EPYC 7543 node (byte identity
# is only claimed there), with the node flags of dev/scripts/launcher_guard.sh:
#
#   L=/fred/oz022/tkimpson/infercell-kim56-2026-11-28
#   X="gina[1-18,301-305],$excluded"     # $excluded as dev/scripts/launcher_guard.sh builds it
#   precompile: sbatch --partition=milan --exclude="$X" --cpus-per-task=1 --mem=8G --time=01:00:00 \
#       --output=$L/logs/%x_%j.log --export=ALL,MODE=precompile dev/scripts/fullcell_params_kim56/params_check.slurm
#   table:      ... --export=ALL,MODE=table,OUT=$L/table ...
#   check:      ... --time=03:00:00 --export=ALL,MODE=check,PROCIDS=1-6,SEED=1,MINUTES=181,OUT=$L/check ...
#
# Julia 1.11.9 (SpeedyCell needs 1.11; InferCell's suite is on julia/1.10.5, so this
# check is a Slurm script, not part of Pkg.test). Depot: the run's own (DEPOT), stacked
# on KIM-183's and the shared one that holds LowPrecSSA (unregistered, SSH-pinned at
# cb0abb6), as in KIM-180's doubling.slurm, so no compile cache lands in another
# issue's archive.
#SBATCH --job-name=kim56-params
#SBATCH --ntasks=1
#SBATCH --time=03:00:00

set -euo pipefail
SYN3A=/fred/oz022/tkimpson/syn3a-ensemble
# sbatch copies the script, so $0 is not in the repo; use the submit directory.
ENV_DIR=$SLURM_SUBMIT_DIR/dev/scripts/fullcell_params_kim56
[ -f "$ENV_DIR/params_check.slurm" ] || { echo "submit from the worktree root" >&2; exit 1; }
# The pin's one copy, also read by adapter.jl.
SYN3A_PIN=$(cat "$ENV_DIR/SYN3A_PIN")
lscpu | grep -q "AMD EPYC 7543" || { echo "not an EPYC 7543 node: $(hostname)" >&2; exit 1; }
echo "host $(hostname)"; lscpu | grep 'Model name'
head=$(git -C "$SYN3A" rev-parse HEAD)
[ "$head" = "$(git -C "$SYN3A" rev-parse "$SYN3A_PIN^{commit}")" ] ||
    { echo "syn3a-ensemble is at $head, not $SYN3A_PIN" >&2; exit 1; }
[ -z "$(git -C "$SYN3A" status --porcelain -- code/SpeedyCell data/model)" ] ||
    { echo "syn3a-ensemble code/SpeedyCell or data/model is not clean" >&2; exit 1; }
echo "syn3a_ensemble_git $head (pin $SYN3A_PIN, code/SpeedyCell and data/model clean)"
export PATH=/fred/oz022/tkimpson/julia-1.11.9/bin:$PATH
export JULIA_DEPOT_PATH=${DEPOT:?set DEPOT}:/fred/oz022/tkimpson/infercell-kim183-2026-11-26/depot:/fred/oz022/tkimpson/julia_depot
export JULIA_PKG_OFFLINE=true
export OMP_NUM_THREADS=1 OPENBLAS_NUM_THREADS=1 MKL_NUM_THREADS=1
echo "infercell_git $(git -C "$ENV_DIR" rev-parse HEAD)"
case "${MODE:-}" in
  precompile)
    julia --project="$ENV_DIR" -e 'using Pkg; Pkg.precompile(); using SpeedyCell.Cell, InteractiveUtils; versioninfo()' ;;
  table)
    julia -t 1 --project="$ENV_DIR" "$ENV_DIR/table.jl" ;;
  check)
    julia -t 1 --project="$ENV_DIR" "$ENV_DIR/check.jl" ;;
  *) echo "MODE must be precompile, table or check" >&2; exit 1 ;;
esac
