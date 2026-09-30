# Merge path_density_16a2_m0.jl's tasks and report V1's importance identity on M0.
# Usage: julia --project dev/scripts/path_density_16a2_m0_merge.jl

using Serialization
using Statistics
using Printf

const DIR = joinpath(@__DIR__, "density_16a2_m0")
parts = [deserialize(joinpath(DIR, f)) for f in readdir(DIR) if endswith(f, ".jls")]
all(p -> p.commit == parts[1].commit && p.truth == parts[1].truth, parts) ||
    error("the tasks come from more than one commit or truth")
n = sum(length(p.paths) for p in parts)
allunique(reduce(vcat, [p.paths for p in parts])) || error("a path repeats")
f = reduce(vcat, [p.firings for p in parts])
@printf("V1 on M0: %d paths from %d tasks, commit %s, job %s; firings per path median %d (%d to %d)\n",
        n, length(parts), parts[1].commit, parts[1].job, median(f), minimum(f), maximum(f))
for k in (:cme, :ode)
    w = reduce(vcat, [getfield(p.w, k) for p in parts])
    se = std(w) / sqrt(length(w))
    ess = sum(w)^2 / sum(w .^ 2)
    @printf("  %s: mean %.5f, SE %.5f, |mean - 1| / SE = %.2f (pass at 3: %s); weight sd %.4f, max share %.2e, ESS %.0f\n",
            k, mean(w), se, abs(mean(w) - 1) / se, abs(mean(w) - 1) <= 3se, std(w),
            maximum(w) / sum(w), ess)
end
