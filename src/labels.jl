"""
Labelling what is ours.

Core A′ carries treatments the published model does not make: chemostatted CTP,
UTP and amino-acid pools, a lumped tRNA charging step, asserted priors on the
PTS mass-action constants, and — where an author selects them — a smoothed or
unclamped expression-cost drain or a continuous rebuild.

`:driver_policy` and `:discarded_prior` need a driver too: the first records
the drain granularity and the rounding when they are *not* departures, so both
reach every report, and the second is an asserted prior on a slot the driver
overwrites. Four of the categories are the *driver's* policy rather than any
module's declaration, and so are enumerated by [`driver_declarations`](@ref) rather than
by [`reduction_declarations`](@ref): a drain coarser than the handshake
(`:coarse_drain`), a rounding policy other than fractional carry
(`:rounding_policy`), and the two the growth chain carries — the membrane
footprint being a calibrated constant of the published model rather than a
measurement (`:calibrated_constant`) and the frozen non-ptsG membrane baseline
(`:exogenous_growth`). They share this vocabulary because they share the
report.

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
                              :continuous_rebuild, :coarse_drain, :rounding_policy,
                              :calibrated_constant, :exogenous_growth,
                              :asserted_prior, :discarded_prior, :lumping,
                              :formalism, :model_note,
                              :capped_rate_law_geometry, :driver_policy)

# `:driver_policy` records the driver's policy where it is *not* a departure —
# the published 1 s drain, and fractional carry — so that the drain granularity
# and the rounding reach every report whether or not they depart (spec §11 task
# 13.6). The report prints it apart and leaves it out of its departure count.
const _POLICY_RECORDS = (:driver_policy,)

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
parameters whose prior this project asserted, and what each module registered
through [`reduction_notes`](@ref) — under the category a note names, or
`:model_note`.

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

    # A plain note is a module's own simplification; a `category => text` pair
    # names its category, which is how the lumping is reported as `:lumping`
    # and nothing else is (spec §11 task 13.6).
    for m in models
        for note in reduction_notes(m)
            push!(labels, note isa Pair ?
                  ReductionLabel(first(note), module_id(m), last(note)) :
                  ReductionLabel(:model_note, module_id(m), note))
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

Pass a [`HandshakeDriver`](@ref) as a second argument to include the driver's
own policy — the drain granularity and the rounding — which is a departure the
sub-models do not declare and which this method therefore cannot see. §6 F11's
panel and §6 T2 are generated from the two-argument form; the one-argument form
says so where it finds nothing.
"""
reduction_report(models::Vector{<:AbstractSubModel}) =
    _reduction_report(reduction_declarations(models))

# The rendering, shared with the two-argument form in `handshake.jl`, which adds
# the driver's own policy to the same list.
function _reduction_report(labels::Vector{ReductionLabel}; driver_seen = false)
    isempty(labels) && return driver_seen ?
        "Nothing in this composition's declarations, and nothing in this driver's " *
        "policy, departs from the published model." :
        "Nothing in this composition's declarations departs from the published " *
        "model. A driver's own policy — the drain granularity, the rounding — is " *
        "not a declaration and is enumerated by `driver_declarations`; pass the " *
        "driver to see both."

    departures = count(l -> !(l.category in _POLICY_RECORDS), labels)
    lines = ["$departures declaration(s) that are this reduction's, not the published model's:"]
    for category in REDUCTION_CATEGORIES
        of_category = filter(l -> l.category === category, labels)
        isempty(of_category) && continue
        push!(lines, category in _POLICY_RECORDS ?
              "  $category (recorded whether or not it departs; not counted above):" :
              "  $category:")
        for l in of_category
            push!(lines, "    - $(l.description)")
        end
    end
    return join(lines, "\n")
end

reduction_report(model::AbstractSubModel) = reduction_report([model])

export ReductionLabel, REDUCTION_CATEGORIES, reduction_declarations,
       reduction_report
