#!/bin/bash
# KIM-180: KIM-183's timing_cycle.slurm with the depot as a setting (DEPOT), so this
# run's compile cache is its own. Submitted by doubling.sh; see there.
#SBATCH --job-name=kim180-cycle
#SBATCH --ntasks=1
#SBATCH --time=02:00:00

set -euo pipefail
lscpu | grep -q "AMD EPYC 7543" || { echo "not an EPYC 7543 node: $(hostname)" >&2; exit 1; }
echo "host $(hostname)"; lscpu | grep 'Model name'
export PATH=/fred/oz022/tkimpson/julia-1.11.9/bin:$PATH
export JULIA_DEPOT_PATH=${DEPOT:?set DEPOT}:/fred/oz022/tkimpson/infercell-kim183-2026-11-26/depot:/fred/oz022/tkimpson/julia_depot
export JULIA_PKG_OFFLINE=true
export OMP_NUM_THREADS=1 OPENBLAS_NUM_THREADS=1 MKL_NUM_THREADS=1
ENV_DIR=$SLURM_SUBMIT_DIR/dev/scripts/fullcell_timing_cycle
[ -f "$ENV_DIR/timing_cycle.jl" ] || { echo "submit from the worktree root" >&2; exit 1; }
echo "git $(git -C "$ENV_DIR" rev-parse HEAD)"
NT=${SLURM_CPUS_PER_TASK:-1}
case "${MODE:-}" in
  precompile)
    julia --project="$ENV_DIR" -e 'using Pkg; Pkg.precompile(); using SpeedyCell.Cell, InteractiveUtils; versioninfo()' ;;
  ensemble)
    julia -t "$NT" --gcthreads="$NT,1" --project="$ENV_DIR" "$ENV_DIR/timing_cycle.jl" ;;
  *) echo "MODE must be precompile or ensemble" >&2; exit 1 ;;
esac
