# Does V7's commit `dfe570a` regenerate the bulk data the figure overlays? The fit's file does
# not hold them, so `trajectories_predict.jl`'s `data` stage regenerates them at its own commit.
# A difference goes in the result note and the caption, and does not stop the work, since the
# posterior is unchanged. It stands alone because it runs in a worktree at `dfe570a`, which has
# none of the trajectories scripts. Started by `trajectories_predict.sh <outdir> bulk <worktree>`:
#   julia --project <this file> <outdir>
using InferCell
using Serialization

out = ARGS[1]
dest = joinpath(out, "bulk_check.md")
isfile(dest) && error("$dest exists; this run never overwrites an output")
commit = strip(read(`git rev-parse HEAD`, String))
startswith(commit, "dfe570a") || error("this tree is at $commit, not V7's dfe570a")

x = deserialize(joinpath(out, "data.jls"))
ds = m0_dataset(16_200_001; ncells = 50)
same = ds.times == x.result.times && ds.bulk == x.result.bulk
worst = size(ds.bulk) == size(x.result.bulk) ? maximum(abs.(ds.bulk ./ x.result.bulk .- 1)) : NaN
open(dest, "w") do io
    println(io, "# The bulk data at V7's commit against the run's\n")
    println(io, "- at ", commit, ", job ", get(ENV, "SLURM_JOB_ID", "none"), "; the run's data at ", x.commit)
    println(io, "- ", same ? "EQUAL, exactly" : "DIFFERENT", "; largest relative difference ", worst)
end
println("bulk data ", same ? "equal" : "differ", ", largest relative difference ", worst)
