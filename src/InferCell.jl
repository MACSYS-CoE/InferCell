module InferCell

using OrdinaryDiffEq: Tsit5, solve, ODEProblem
using SciMLSensitivity: ForwardDiffSensitivity
using SciMLBase: ReturnCode, remake, DiscreteProblem
using JumpProcesses: JumpProblem, JumpSet, ConstantRateJump, Direct, SSAStepper
using Turing
using Distributions
using StaticArrays
using Random
using MCMCChains: replacenames
using ForwardDiff
using LinearAlgebra: rank
using Statistics: quantile, mean, std

include("parameters.jl")
# The Core A′ registry is included before the framework files because edges.jl,
# resolver.jl and loader.jl still call its functions (held_value, is_registered,
# species_index, TRANSCRIPTION_ATOL) directly. That dependency runs the wrong way
# for a second organism; parameterising them over a registry is deferred until
# there is one.
include("organisms/coreA/registry.jl")
include("edges.jl")
include("interface.jl")
include("resolver.jl")
include("labels.jl")
include("loader.jl")
include("likelihoods.jl")
include("orchestrator.jl")
include("models/transcription_translation.jl")
include("models/stochastic_gene_expression.jl")
include("models/bursty_gene_expression.jl")
include("models/light_metabolism.jl")
include("models/tier_b_metabolism.jl")
include("summary_statistics.jl")
include("abc_smc.jl")
include("inference.jl")
include("boundary.jl")

export InferParameter, free_params, rate_params, ic_params, obs_params,
       model_free_params, ode_free_params, obs_free_params, unique_params
export ParameterSource, provenance_of, source_file, informedness,
       asserted_prior_params, uninformed_params, provenance_conflicts,
       was_chosen_over_alternative
export AbstractSubModel, SubModelContext
export states, parameters, dynamics, inputs, formalism, inference_mode, reactions,
       coupling, reduction_notes, module_id
# The Core A′ interface contract (registry, edge kinds, resolver, labels,
# loader) exports from its own file, so seven parallel wave-1 branches do not
# all append to one export block here.
export ObservedData, PosteriorPredictive, ABCPosterior
export build_problem
export TranscriptionTranslation, StochasticGeneExpression, BurstyGeneExpression, LightMetabolism, TierBMetabolism
export build_turing_model, infer, observe, posterior_predictive, check_identifiability
export compute_summary_stats, summary_distance
export abc_smc
export KDEPrior, boundary_condition, sequential_infer, iterative_infer, kl_divergence, chain_kl
export Tsit5, solve, NUTS, SSAStepper

end # module InferCell
