# Export serialized results to CSV for Python plotting
using InferCell
using Serialization
using DelimitedFiles

results = deserialize("examples/results/step3_3block_results.jls")
ode_chain = results["ode_chain"]
cond_result = results["cond_result"]
uncond_result = results["uncond_result"]

mkpath("examples/results")

# ODE chain samples (8 params x n_samples)
param_names = ["k_tx", "k_tl", "gamma_mRNA", "gamma_protein",
               "k_atp", "k_ntp", "k_aa", "sigma_obs"]
ode_matrix = hcat([vec(ode_chain[p].data) for p in param_names]...)
header = join(param_names, ",")
open("examples/results/ode_chain.csv", "w") do f
    println(f, header)
    writedlm(f, ode_matrix, ',')
end

# Conditioned ABC particles + weights
ssa_names = ["k_tx", "k_tl", "gamma_mRNA", "gamma_protein"]
open("examples/results/cond_particles.csv", "w") do f
    println(f, join(ssa_names, ","))
    writedlm(f, cond_result.particles', ',')
end
writedlm("examples/results/cond_weights.csv", cond_result.weights, ',')

# Unconditioned ABC particles + weights
open("examples/results/uncond_particles.csv", "w") do f
    println(f, join(ssa_names, ","))
    writedlm(f, uncond_result.particles', ',')
end
writedlm("examples/results/uncond_weights.csv", uncond_result.weights, ',')

println("Exported to examples/results/*.csv")
