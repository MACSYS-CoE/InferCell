module InferCell

using OrdinaryDiffEq: Tsit5, solve, ODEProblem, Rodas5P
using SciMLSensitivity: ForwardDiffSensitivity
using SciMLBase: ReturnCode, remake, DiscreteProblem, init, step!, u_modified!
using JumpProcesses: JumpProblem, JumpSet, ConstantRateJump, Direct, SSAStepper,
                     reset_aggregated_jumps!
using Turing
using Distributions
using StaticArrays
using Random
using MCMCChains: replacenames
using ForwardDiff
using LinearAlgebra: rank
using Statistics: quantile, mean, std

include("parameters.jl")
# Reading order, not a load-time requirement: resolver.jl and loader.jl call
# registry functions only from inside function bodies, which Julia resolves at
# call time, so this order could legally be anything. edges.jl and labels.jl
# name no registry symbol at all. Making resolver and loader take a registry
# rather than reach for the global one is what a second organism costs.
include("organisms/coreA/registry.jl")
include("edges.jl")
include("interface.jl")
include("resolver.jl")
include("labels.jl")
include("loader.jl")
include("likelihoods.jl")
include("orchestrator.jl")
include("handshake.jl")
# The Core A′ modules. Unlike the files above these must come after
# interface.jl: a sub-model's `<: AbstractSubModel` is resolved when the
# struct is defined, not when a method is called. Each exports its own
# names from its own file, so concurrent module branches append one
# include here and touch nothing else.
#
# PtsTransport needs handshake.jl in scope too, not merely interface.jl: it
# defaults a constructor argument to COREA_INITIAL_RADIUS_NM.
include("organisms/coreA/central_glycolysis.jl")
include("organisms/coreA/pts_transport.jl")
# NucleotideRecycling reads COREA_DATA_DIR, which pts_transport.jl defines, so
# it is included after it rather than beside the models below.
include("organisms/coreA/nucleotide_recycling.jl")
# Transcription reads its own extract via `@__DIR__`, so it has no ordering
# constraint against the three above; it needs handshake.jl in scope because
# its rate law converts a polymerase copy number at the registry's volume.
include("organisms/coreA/transcription.jl")
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
export AbstractSubModel, SubModelContext, Reaction
export states, parameters, dynamics, inputs, formalism, inference_mode, reactions,
       coupling, reduction_notes, module_id, contributed_states, contributions,
       written_states, rebuilt_params, rate_constants,
       membrane_protein_states, extracellular_states
# The Core A′ interface contract (registry, edge kinds, resolver, labels,
# loader) exports from its own files, so seven parallel wave-1 branches do not
# all append to one export block here.
export ObservedData, PosteriorPredictive, ABCPosterior
export build_problem
export TranscriptionTranslation, StochasticGeneExpression, BurstyGeneExpression, LightMetabolism, TierBMetabolism
export build_turing_model, infer, observe, posterior_predictive, check_identifiability
export compute_summary_stats, summary_distance
export abc_smc
export KDEPrior, boundary_condition, sequential_infer, iterative_infer, kl_divergence, chain_kl
export Tsit5, solve, NUTS, SSAStepper, Rodas5P

end # module InferCell
