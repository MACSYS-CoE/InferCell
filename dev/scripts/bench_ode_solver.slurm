#!/bin/bash
#SBATCH --job-name=infercell_bench_ode_solver
#SBATCH --partition=milan
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=16
#SBATCH --mem=32G
#SBATCH --time=03:00:00
#SBATCH --output=dev/scripts/bench_ode_solver_%j.log

# The serial stages use one thread; THREADS copies run at once only in the
# concurrent stage. Ask for that many cores, so the copies do not share them.
# Pass settings with --export (this cluster does not forward the caller's environment):
#   sbatch --export=ALL,NSTEPS=60,REPS=1,THREADS=2,RESULT=dev/scripts/bench_ode_solver_smoke.md --cpus-per-task=2 ...
module load julia/1.10.5
THREADS=${THREADS:-16}
export THREADS
echo "ODE solver benchmark at $(date) on $(hostname), HEAD $(git rev-parse --short HEAD)$(test -z "$(git status --porcelain -- src test dev/scripts/bench_ode_solver.jl)" || echo ' (tree dirty)')"
echo "MODEL=${MODEL:-both} NSTEPS=${NSTEPS:-full} VARIANTS=${VARIANTS:-all} REPS=${REPS:-3} THREADS=$THREADS"
julia --project -t "$THREADS" dev/scripts/bench_ode_solver.jl
status=$?
echo "Exit code: $status  Finished at $(date)"
exit $status
