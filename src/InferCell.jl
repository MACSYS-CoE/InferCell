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
include("corea/registry.jl")
include("corea/edges.jl")
include("interface.jl")
include("corea/resolver.jl")
include("corea/labels.jl")
include("corea/loader.jl")
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
export ParameterSource, INFORMEDNESS, provenance_of, source_file, informedness,
       asserted_prior_params, uninformed_params, provenance_conflicts,
       was_chosen_over_alternative
export AbstractSubModel, SubModelContext
export states, parameters, dynamics, inputs, formalism, inference_mode, reactions,
       coupling, reduction_notes, module_id
# Core A′ interface contract: the state registry, the seven edge kinds, the
# resolver, and the provenance-carrying loader.
export SpeciesEntry, COREA_SPECIES, COREA_SPECIES_INDEX,
       species_index, species_entry, is_registered, species_group,
       is_chemostatted, is_dynamic, dynamic_species, chemostat_species,
       species_in_group, held_value, is_informed, uninformed_species,
       n_dynamic_states, n_chemostats
export CouplingEdge, MassEdge, CurrencyEdge, DeferredCounterEdge, CatalyticEdge,
       RateConstantEdge, VolumeEdge, ClampedEdge, EDGE_KINDS, EDGE_DIRECTIONS,
       CLIP_POLICIES, RATE_CADENCES, EDGE_ORIGINS,
       edge_kind, edge_species, edge_direction, edge_peer,
       is_consumer, is_producer, carries_mass, mass_contribution,
       obstructs_gradients, deviates_from_published, deviation_reason
export ResolvedEdge, DeadEnd, CouplingGraph, resolve_coupling,
       dead_end_report, gradient_report, check_gradient_safety
export ReductionLabel, reduction_declarations, reduction_report
export SourceTable, AmbiguousValue, read_source_table, ambiguity_report,
       disagreements, load_parameter, governing_choices
export ObservedData, PosteriorPredictive, ABCPosterior
export build_problem
export TranscriptionTranslation, StochasticGeneExpression, BurstyGeneExpression, LightMetabolism, TierBMetabolism
export build_turing_model, infer, observe, posterior_predictive, check_identifiability
export compute_summary_stats, summary_distance
export abc_smc
export KDEPrior, boundary_condition, sequential_infer, iterative_infer, kl_divergence, chain_kl
export Tsit5, solve, NUTS, SSAStepper

end # module InferCell
