"""
Labelling what is ours.

Core A′ carries treatments the published model does not make: chemostatted CTP,
UTP and amino-acid pools, a lumped tRNA charging step, asserted priors on the
PTS mass-action constants, and — where an author selects them — a smoothed or
unclamped expression-cost drain or a continuous rebuild.

Each is a place where a result could depend on our choice rather than on the
published model, so each has to be able to reach the report that depends on it.
This file is the enumeration that makes that possible.
"""

"""
The label categories, in the order [`reduction_report`](@ref) prints them.
A label outside this vocabulary would be counted by the report's header and
silently dropped from its body, so membership is enforced at construction.
"""
const REDUCTION_CATEGORIES = (:clamp, :smoothed_counter, :unclamped_counter,
                              :continuous_rebuild, :asserted_prior, :lumping)

"""
    ReductionLabel

One declaration that departs from the published model.

- `category` — one of [`REDUCTION_CATEGORIES`](@ref).
- `subject` — the species, parameter or module the label attaches to.
- `description` — one line, written for a human reading a report.
"""
struct ReductionLabel
    category::Symbol
    subject::Symbol
    description::String

    function ReductionLabel(category, subject, description)
        _check_vocab(:ReductionLabel, :category, category, REDUCTION_CATEGORIES)
        return new(category, subject, description)
    end
end

"""
    reduction_declarations(models) -> Vector{ReductionLabel}

Every declaration in a composition that is this reduction's rather than the
published model's.

Collects, in order: clamps introduced by the reduction, deferred counters using
a smoothed or unclamped clip, rate-constant edges refreshing continuously,
parameters whose prior this project asserted, and lumpings a module registered
through [`reduction_notes`](@ref).

A deviation that reaches a result unlabelled is the failure this exists to
prevent, so prefer calling it over remembering what was declared where.
"""
function reduction_declarations(models::Vector{<:AbstractSubModel})
    labels = ReductionLabel[]

    graph = resolve_coupling(models)
    for r in graph.deviations
        reason = deviation_reason(r.edge)   # non-nothing for every deviation
        push!(labels, ReductionLabel(deviation_category(r.edge), r.species, reason))
    end

    # Deduplicate by name first: a parameter shared by two modules is one
    # asserted prior, not two.
    all_params = unique_params(InferParameter[p for m in models for p in parameters(m)])
    for p in asserted_prior_params(all_params)
        push!(labels, ReductionLabel(
            :asserted_prior, p.name,
            "prior on :$(p.name) is asserted by this project — its source " *
            "carries a point value with no quantified uncertainty"))
    end

    for m in models
        for note in reduction_notes(m)
            push!(labels, ReductionLabel(:lumping, module_id(m), note))
        end
    end

    return labels
end

reduction_declarations(model::AbstractSubModel) = reduction_declarations([model])

"""
    reduction_report(models) -> String

[`reduction_declarations`](@ref) written out for a human, grouped by category.
Suitable for printing beside any result whose interpretation depends on a
treatment that is ours rather than the published model's.
"""
function reduction_report(models::Vector{<:AbstractSubModel})
    labels = reduction_declarations(models)
    isempty(labels) && return "Nothing in this composition departs from the published model."

    lines = ["$(length(labels)) declaration(s) that are this reduction's, not the published model's:"]
    for category in REDUCTION_CATEGORIES
        of_category = filter(l -> l.category === category, labels)
        isempty(of_category) && continue
        push!(lines, "  $category:")
        for l in of_category
            push!(lines, "    - $(l.description)")
        end
    end
    return join(lines, "\n")
end

reduction_report(model::AbstractSubModel) = reduction_report([model])

export ReductionLabel, REDUCTION_CATEGORIES, reduction_declarations,
       reduction_report
