# `predict` (KIM-27): one public entry point for probabilistic prediction.
#
# It computes the objects of dev/notes/2026-10-29-cut-predictive-distribution.md,
# and keeps five things apart: the posterior source (a `PosteriorDraws`), the
# forward description (`forward`, `delta`), the completion policy (`completion`,
# `held`, `width`, `viability`), the estimand, and the intervention.
#
# Its layout is one crossed matrix: M posterior draws θ_I^m × K held draws θ_H^k,
# the K drawn once and shared by every row (note §5.2). Each matrix cell (m, k)
# runs its own fresh cells, since the E[σ²]/n correction holds only if no two
# draws share cell seeds (note §2). An intervention's arm reuses its baseline
# cell's seeds, so the pairing is within a draw only (note §4).

"""
    Prediction

What [`predict`](@ref) returns. Rows are posterior draws (m = 1…M), columns held
draws (k = 1…K). Without completion there is one column.

**The draws:**
- `draws`: each row's index into the posterior;
- `theta_I`: each row's posterior values;
- `sigma_b`: each row's σ_b, carried and never simulated, or `nothing`;
- `held` and `theta_H`: the held names, and each column's drawn values;
- `seeds`: the cell seeds, n_cells × K × M. They are empty for the mean field.

**The runs:**
- `cells`: every run, M × K, each a vector of [`CellRun`](@ref);
- `mu`: μ̂_mk, the mean of Q over each matrix cell's finished runs. For the mean
  field this is m(θ), uncorrected;
- `sigma2`: σ̂²_mk, the sample variance of Q over those cells. It is `NaN` for the
  mean field, which has no single-cell distribution, and below two cells;
- `n_ok`: the runs that finished.

**The intervention**, when one is given:
- `cells_intervened`: the paired runs, on the same seeds. Pairs that V rejected
  are not run;
- `effect` and `effect_sigma2`: the paired effect Δ's mean and variance.

**The conditioning:**
- `viability`: the rule;
- `viable`: V's decision per pair, all true without a rule;
- `decisions`: each pair's `(; viable, fraction, counts)`
  ([`assess_viability`](@ref)), so a rejected pair keeps its reason;
- `viability_seeds` and `viability_cells`: V's own cells, n_V × K × M seeds and
  M × K runs, when `n_viability_cells` gives V cells of its own. They are empty
  otherwise, and V reads the cells in `cells`;
- `weights`: the weight sets over the matrix, `uniform`, or `per_theta_I` and
  `joint` under V ([`viability_weights`](@ref)), each with its acceptance and ESS
  in `conditioning`.

**The results:**
- `delta`: the mean field's [`DeltaCorrection`](@ref), or `nothing`;
- `summary` and `effect_summary`: each weight set's population-mean results. See
  [`predict`](@ref);
- `statuses`: the count of every baseline run that gives μ̂ in each
  [`CELL_STATUSES`](@ref) category, read without a rule. This is the unfiltered
  ensemble's validity diagnostic. V's own cells are not counted here; V's
  decisions count them, and their failures are in provenance;
- `statuses_intervened`: the same count over the intervened runs, all zero
  without an intervention;
- `provenance`: what was run, and from what (KIM-33):
  - `fork_seed`: the seed of the stream forked from `rng`, from which every draw,
    held value and cell seed came;
  - `version`: InferCell's;
  - `settings`: the call's arguments, as run;
  - `source`: the posterior's;
  - `reductions`: the composition's declared reductions, the driver's or route
    (ii)'s included ([`reduction_declarations`](@ref));
  - `skipped`: the parameters the completion left out, and why
    ([`held_parameters`](@ref));
  - `pairs`: one record per pair, M × K. Each gives the posterior draw (index,
    chain, iteration), where θ_H came from (`:none`, `:prior` or
    `:viability_conditioned`), its seeds and V's own seeds, V's decision, and every
    failed run (of μ̂'s cells, V's and the intervened ones) with its seed and
    reason;
  - `threads`: the threads the cells ran on, 1 when not threaded. The result
    does not depend on it.
"""
struct Prediction
    forward::Symbol
    completion::Symbol
    estimand::Symbol
    quantity::Symbol
    horizon::Float64
    draws::Vector{Int}
    theta_I::Vector{Vector{Pair{Symbol, Float64}}}
    sigma_b::Union{Nothing, Vector{Float64}}
    held::Vector{Symbol}
    theta_H::Vector{Vector{Pair{Symbol, Float64}}}
    seeds::Array{Int, 3}
    cells::Matrix{Vector{CellRun}}
    mu::Matrix{Float64}
    sigma2::Matrix{Float64}
    n_ok::Matrix{Int}
    intervention::Union{Nothing, Intervention}
    cells_intervened::Matrix{Vector{CellRun}}
    effect::Matrix{Float64}
    effect_sigma2::Matrix{Float64}
    viability::Union{Nothing, ViabilityRule}
    viable::Matrix{Bool}
    decisions::Matrix{NamedTuple}
    viability_seeds::Array{Int, 3}
    viability_cells::Matrix{Vector{CellRun}}
    weights::NamedTuple
    conditioning::NamedTuple
    delta::Union{Nothing, DeltaCorrection}
    summary::NamedTuple
    effect_summary::Union{Nothing, NamedTuple}
    statuses::Dict{Symbol, Int}
    statuses_intervened::Dict{Symbol, Int}
    provenance::NamedTuple
end

function Base.show(io::IO, p::Prediction)
    M, K = size(p.mu)
    print(io, "Prediction($(p.forward), $(p.completion), Q = :$(p.quantity); ",
          "$M posterior × $K held draws",
          p.intervention === nothing ? "" : ", with an intervention",
          p.delta === nothing ? "" : ", δ-corrected", ")")
end

const PREDICT_FORWARDS = (:hybrid, :mean_field)
const PREDICT_COMPLETIONS = (:none, :prior, :viability_conditioned_prior)

# Q over the runs that finished: the mean, and the sample variance from two up.
function _moments(qs)
    q = filter(isfinite, qs)
    n = length(q)
    return n == 0 ? (NaN, NaN, 0) : (mean(q), n >= 2 ? var(q) : NaN, n)
end

# A weighted population-mean summary over the matrix (note §2, §5.1). The pairs
# with no finished run carry no Q, so they leave the weights, and the summary
# counts them.
function _summary(mu, s2, nok, w, shift, forward)
    ok = isfinite.(mu) .& (w .> 0)
    W = sum(w .* ok)
    W > 0 || return (; mean = NaN, var = NaN, mc = NaN, parameter_var = NaN,
                     flagged = false, mc_undefined = false, intrinsic_var = NaN, ess = 0.0,
                     n_pairs = 0, n_no_q = count(w .> 0))
    v = (w .* ok) ./ W
    m = sum(v[i] * mu[i] for i in eachindex(mu) if ok[i])
    # The reliability-weighted variance: the unbiased sample variance when the
    # weights are equal.
    d = 1 - sum(abs2, v)
    raw = d > 0 ? sum(v[i] * (mu[i] - m)^2 for i in eachindex(mu) if ok[i]) / d : NaN
    if forward === :hybrid
        # The noise E[σ̂²/n] adds to the reliability-weighted variance:
        # Σ v(1 − v)τ / (1 − Σv²), which is Σ vτ only for equal weights (KIM-66).
        mc = d > 0 ? sum(v[i] * (1 - v[i]) * s2[i] / nok[i] for i in eachindex(mu) if ok[i]) / d :
                     NaN
        intrinsic = sum(v[i] * s2[i] for i in eachindex(mu) if ok[i])
    else
        mc, intrinsic = 0.0, NaN
    end
    pv = raw - mc
    # Below two finished cells a pair has no σ̂², so E[σ̂²]/n, and the parameter
    # variance with it, is undefined (NaN), and said to be.
    few = forward === :hybrid && any(nok[i] < 2 for i in eachindex(mu) if ok[i])
    return (; mean = m + shift, var = raw, mc, parameter_var = max(pv, 0.0),
            flagged = pv < 0, mc_undefined = few, intrinsic_var = intrinsic, ess = _kish(v),
            n_pairs = count(ok), n_no_q = count((w .> 0) .& .!isfinite.(mu)))
end

# The runs that finished at baseline and failed under the intervention, over the
# pairs `w` weighs. They leave ΔQ, but they are part of the prediction (note §4(4)).
_n_failed_intervened(cells, cells_i, w) =
    sum((count(((b, c),) -> b.failure === :none && c.failure !== :none, zip(cells[i], cells_i[i]))
         for i in eachindex(w) if w[i] > 0); init = 0)

# How many of a matrix's runs fall in each of CELL_STATUSES, read without a rule.
function _status_counts(cells)
    n = Dict(s => 0 for s in CELL_STATUSES)
    for cs in cells, c in cs
        n[cell_status(c)] += 1
    end
    return n
end

function _select_draws(rng, pd::PosteriorDraws, n::Integer)
    1 <= n || throw(ArgumentError("n_parameter_draws must be at least 1, got $n"))
    pd.weights === nothing || return [_weighted_sample(pd.weights, rng) for _ in 1:n]
    n <= length(pd) || throw(ArgumentError(
        "n_parameter_draws = $n, but the posterior holds $(length(pd)) draws"))
    return n == length(pd) ? collect(1:n) : sort(randperm(rng, length(pd))[1:n])
end

# Seeds as positive Ints from the stream, refused if one repeats.
function _cell_seeds(rng, n::Integer)
    s = [Int(x >> 1) for x in rand(rng, UInt64, n)]
    allunique(s) || error("two cell seeds coincide; draw again with another rng")
    return s
end

# The reference point for δ: the posterior mean of log θ_I, held at nominal.
function _theta_ref(pd::PosteriorDraws)
    w = pd.weights === nothing ? fill(1 / length(pd), length(pd)) : pd.weights
    return [n => exp(sum(w .* log.(pd.values[:, j]))) for (j, n) in enumerate(pd.names)]
end

function _check_delta(dc::DeltaCorrection, ms, horizon, quantity, intervention)
    dc.modules == [module_id(m) for m in ms] || throw(ArgumentError(
        "the DeltaCorrection was computed on modules $(dc.modules), not this composition"))
    dc.fingerprint == _fingerprint(ms) || throw(ArgumentError(
        "the DeltaCorrection was computed on these modules with other settings, so it " *
        "is not this composition's"))
    dc.horizon == horizon || throw(ArgumentError(
        "the DeltaCorrection's horizon is $(dc.horizon) s, not $horizon s"))
    dc.quantity === _quantity_label(quantity) || throw(ArgumentError(
        "the DeltaCorrection is for Q = :$(dc.quantity), not :$(_quantity_label(quantity))"))
    isequal(dc.quantity_def, quantity) || throw(ArgumentError(
        "the DeltaCorrection is for another Q with the same label, :$(dc.quantity); " *
        "compute it with this Q"))
    if intervention !== nothing
        dc.intervention == intervention || throw(ArgumentError(
            "the DeltaCorrection does not cover this intervention; compute it with " *
            "intervention = $intervention"))
    end
    return dc
end

"""
    predict(models, posterior::PosteriorDraws; horizon, forward = :hybrid,
            completion = :none, held = nothing, width = 1.0, n_held = 1,
            viability = nothing, quantity = :lactate_out, estimand = :population_mean,
            intervention = nothing, delta = :none, theta_ref = nothing, n_ref = 100,
            rng = Random.default_rng(), n_parameter_draws = length(posterior),
            draws = nothing, n_cells = 1, n_viability_cells = nothing,
            threads = Threads.nthreads() > 1) -> Prediction

The predictive distribution of the quantity Q given the data, for the composition
`models`, as dev/notes/2026-10-29-cut-predictive-distribution.md defines it. The
composition is [`metered`](@ref) first.

**The posterior source** is a [`PosteriorDraws`](@ref) ([`posterior_draws`](@ref)).
`n_parameter_draws` of its draws are used, chosen without replacement (or by
weight, for an ABC posterior), or the indices `draws` exactly. σ_b is carried
beside each draw and never simulated. Two calls with the same `rng` draw the same
held values and the same cell seeds whatever their `draws`. So splitting one
prediction's rows over calls that share an `rng` reruns the same cells on
different rows. Calls with different `rng`s draw different held values, and do
not form one crossed matrix.

**The forward description:**
- `forward = :hybrid` runs `n_cells` stochastic cells for each parameter set;
- `forward = :mean_field` runs one route (ii) solve.

For the mean field, `delta = :ensemble` adds δ from a hybrid ensemble of `n_ref`
cells at `theta_ref` ([`delta_correction`](@ref)). The default `theta_ref` is the
posterior mean of log θ_I, with every held parameter at its nominal value. A
[`DeltaCorrection`](@ref) computed once can be passed as `delta` instead, and is
checked against this prediction. The uncorrected mean field is kept beside the
corrected one, since the gap between them is part of the result (spec §12
2026-10-29, O4).

**The completion policy:**
- `completion = :none` runs the posterior as it stands, as on M0, where no
  parameter is held;
- `completion = :prior` completes the `held` parameters (by default
  [`held_parameters`](@ref)) from their priors, `width` times as wide.
  `n_held` = K draws are crossed with every posterior draw;
- `completion = :viability_conditioned_prior` does the same, and then reads the
  rule `viability` (a [`ViabilityRule`](@ref)) on each pair's baseline runs.
  Both conditioned variants are returned as weights ([`viability_weights`](@ref)),
  with each pair's decision and reason. Under V, `n_viability_cells` must be
  given, so that the choice is never a default:
  - `n_viability_cells = 0`: V reads the same cells that give μ̂. Then, where Q
    differs between viable and unviable cells and a θ's viable fraction is near
    the rule's threshold, the conditioned μ̂ is biased;
  - `n_viability_cells` = n_V > 0 gives V cells of its own: n_V more cells per
    pair, on seeds drawn after every other, that V reads and μ̂ does not, which
    removes the bias (KIM-36). M4's freeze uses 4.

  The mean field has one deterministic run per θ, which V reads, so it takes 0.

**The estimand** is `:population_mean`, the mean of Q over cells, μ(θ). For each
weight set, `summary` gives:
- `mean`: the predicted population mean, δ-corrected when δ is given;
- `var`: the spread of μ̂ over the matrix's pairs;
- `mc`: the Monte Carlo part of that spread, which the hybrid's finite cells
  add: Σ v(1 − v)τ / (1 − Σv²) for normalised weights v and τ = σ̂²/n, which is
  E[σ̂²]/n when the weights are equal;
- `parameter_var`: `var − mc`. A negative difference is reported as 0, and
  `flagged`. Under two finished cells per pair (the default `n_cells = 1` among
  them), σ̂² does not exist, so `mc` and `parameter_var` are NaN and
  `mc_undefined` is set;
- `intrinsic_var`: E[σ̂²], the cell-to-cell spread (hybrid only). It is a separate
  output, not part of the estimand's uncertainty;
- `ess`, `n_pairs`, and `n_no_q`, the pairs with no finished run.

Without completion, the rows are independent posterior draws, and
`parameter_var` estimates the parameter variance. With completion it does not:
every row shares the same K held draws, so with equal weights the held
parameters' main effect enters `var` only in proportion M(K−1)/(MK−1), about
(K−1)/K, and not at all at `n_held = 1`. Then `var` is a spread over the matrix, not the parameter variance.
The total and its θ_I/θ_H split on the crossed matrix are
[`variance_split`](@ref)'s (note §5.2, KIM-66).

**The intervention** (an [`Intervention`](@ref), such as a [`Knockdown`](@ref))
is applied after the parameter draw. Its arm runs on the same seeds, and
`effect_summary` gives the paired effect ΔQ. V is read at baseline only. ΔQ is
taken over the cells that finish in both arms. A cell that fails only under the
intervention leaves ΔQ, but it is part of the prediction (note §4(4)), so it is
counted: in `statuses_intervened` with every intervened run, in each weight set's
`effect_summary` as `n_failed_intervened` (the cells that finished at baseline
and failed under the intervention, over the pairs the set weighs), and in its
pair's `failed_intervened` (provenance). A knockdown of a reverse constant, or of
a name the composition lacks, is refused before any cell runs.

The cells reseed the global generator. `predict` draws everything from a stream
forked from `rng` and puts the global generator back afterwards, as
[`infer`](@ref) does. Two calls with equal `rng` states give the same result.
The fork seed is recorded, but `predict` cannot be called with one. To replay a
call, keep the `rng`'s seed and the commit.

**Parallel runs.** With `threads` (the default when Julia has more than one
thread), the cells run threaded: every baseline cell or solve of the matrix, μ̂'s
and V's own alike, as one pool, then the intervened ones of the pairs V keeps, and
δ's ensemble likewise. Each thread takes the next cell as it frees up. An error
stops the others taking more and is thrown as itself; a cell whose metabolic block
stalls is not an error but a failed run, kept with its reason.
Each hybrid cell reseeds the generator of the task it runs on before its build
and its run, so the result is the serial one. A crossed matrix with completion
is run by one call, threaded on one node: calls cannot share it, as above.
Without completion (M0), the rows can also be split over calls, each with its
own `rng` and disjoint `draws`, and a corrected mean field given one shared
[`DeltaCorrection`](@ref).
A reduced-model cell at 6,300 s allocates about 25 GiB, so at many threads garbage
collection takes most of the time. Start Julia with one concurrent sweep thread,
`--gcthreads=N,1`; the number of mark threads N did not matter (8 and 32 at 64
threads). On 128 such cells it cut the wall from 464–477 s to 150–166 s at 64
threads, in three runs on two nodes, and from 551 s to 394 s at 16 (KIM-68,
`dev/scripts/threading_kim68_result.md`).
"""
function predict(models::AbstractVector{<:AbstractSubModel}, pd::PosteriorDraws;
                 horizon::Real, forward::Symbol = :hybrid, completion::Symbol = :none,
                 held::Union{Nothing, AbstractVector{Symbol}} = nothing, width::Real = 1.0,
                 n_held::Integer = 1, viability::Union{Nothing, ViabilityRule} = nothing,
                 quantity = :lactate_out, estimand::Symbol = :population_mean,
                 intervention::Union{Nothing, Intervention} = nothing,
                 delta::Union{Symbol, DeltaCorrection} = :none, theta_ref = nothing,
                 n_ref::Integer = 100, rng::AbstractRNG = Random.default_rng(),
                 n_parameter_draws::Integer = length(pd),
                 draws::Union{Nothing, AbstractVector{<:Integer}} = nothing,
                 n_cells::Integer = 1,
                 n_viability_cells::Union{Nothing, Integer} = nothing,
                 threads::Bool = Threads.nthreads() > 1)
    # Refuse before any work.
    forward in PREDICT_FORWARDS || throw(ArgumentError(
        "forward must be one of $PREDICT_FORWARDS, not :$forward"))
    completion in PREDICT_COMPLETIONS || throw(ArgumentError(
        "completion must be one of $PREDICT_COMPLETIONS, not :$completion"))
    estimand === :population_mean || throw(ArgumentError(
        "estimand must be :population_mean, not :$estimand"))
    _check_quantity(quantity)
    conditioned = completion === :viability_conditioned_prior
    conditioned == (viability !== nothing) || throw(ArgumentError(conditioned ?
        "completion = :viability_conditioned_prior needs a ViabilityRule as `viability`" :
        "a viability rule is read only under completion = :viability_conditioned_prior"))
    if completion === :none
        held === nothing && n_held == 1 && width == 1 || throw(ArgumentError(
            "completion = :none draws no held parameters, so `held`, `n_held` and " *
            "`width` do not apply"))
    end
    n_held >= 1 || throw(ArgumentError("n_held must be at least 1, got $n_held"))
    forward === :hybrid && n_cells < 1 && throw(ArgumentError(
        "n_cells must be at least 1, got $n_cells"))
    conditioned && n_viability_cells === nothing && throw(ArgumentError(forward === :hybrid ?
        "under a viability rule, give n_viability_cells: 0 reads V on the cells that " *
        "give μ̂, which biases the conditioned μ̂ where acceptance is partial; n_V > 0 " *
        "gives V cells of its own (M4's freeze uses 4)" :
        "under a viability rule, give n_viability_cells = 0: the mean field has one " *
        "deterministic run per θ, which V reads"))
    n_viability_cells = something(n_viability_cells, 0)
    n_viability_cells >= 0 || throw(ArgumentError(
        "n_viability_cells must be at least 0, got $n_viability_cells"))
    n_viability_cells == 0 || conditioned && forward === :hybrid || throw(ArgumentError(
        "n_viability_cells gives V cells of its own, so it needs " *
        "completion = :viability_conditioned_prior and the hybrid; the mean field has " *
        "one deterministic run per θ, which V reads"))
    delta isa Symbol && !(delta in (:none, :ensemble)) && throw(ArgumentError(
        "delta must be :none, :ensemble or a DeltaCorrection, not :$delta"))
    delta === :none || forward === :mean_field || throw(ArgumentError(
        "δ corrects the mean field; the hybrid needs none"))
    delta === :ensemble && n_ref < 2 && throw(ArgumentError("n_ref must be at least 2"))
    theta_ref === nothing || delta === :ensemble || throw(ArgumentError(
        "theta_ref is read only under delta = :ensemble"))

    ms = metered(models)
    # A knockdown of a reverse constant, or of a name the composition lacks, is
    # refused here, not after the baseline cells have run.
    intervention === nothing || intervene(intervention, ms, Pair{Symbol, Float64}[])
    # Built for its slot lists and labels; a build may read the global generator.
    driver = _restoring_global_rng(() -> build_problem(ms; tspan = (0.0, 1.0), complete = true))
    hp = held_parameters(ms, pd; driver)
    unused = [s for s in hp.skipped if s.name in pd.names && s.reason !== :inferred]
    isempty(unused) || throw(ArgumentError(
        "the posterior names $unused: the driver writes these or they are derived, so " *
        "a posterior value for them would never run"))
    intervention isa Knockdown &&
        any(s -> s.name === intervention.param && s.reason === :driver_written, hp.skipped) &&
        throw(ArgumentError("the driver writes :$(intervention.param) at every handshake, " *
                            "so a knockdown of it would never run"))
    inferred = Set(s.name for s in hp.skipped if s.reason === :inferred)
    bad = setdiff(pd.names, inferred)
    isempty(bad) || throw(ArgumentError(
        "the posterior's $bad are not free parameters this composition reads"))
    hnames = completion === :none ? Symbol[] : held === nothing ? hp.held : collect(Symbol, held)
    badh = setdiff(hnames, hp.held)
    isempty(badh) || throw(ArgumentError(
        "$badh cannot be held: each must be a free parameter that is not inferred, not " *
        "written by the driver, and not a derived reverse constant ($(hp.skipped))"))
    delta isa DeltaCorrection && _check_delta(delta, ms, Float64(horizon), quantity, intervention)

    fork = rand(rng, UInt64)
    stream = Xoshiro(fork)
    return _restoring_global_rng() do
        rows = draws === nothing ? _select_draws(stream, pd, n_parameter_draws) :
               collect(Int, draws)
        all(in(1:length(pd)), rows) || throw(ArgumentError(
            "draws must index the posterior's $(length(pd)) draws"))
        M, K = length(rows), completion === :none ? 1 : n_held
        θH = completion === :none || isempty(hnames) ?
             [Pair{Symbol, Float64}[] for _ in 1:K] :
             draw_held(stream, ms, hnames, K; width)
        nc = forward === :hybrid ? n_cells : 0
        seeds = reshape(_cell_seeds(stream, nc * K * M), nc, K, M)

        dc = if delta === :ensemble
            ref = theta_ref === nothing ? _theta_ref(pd) : collect(theta_ref)
            rs = _cell_seeds(stream, n_ref)
            isdisjoint(rs, seeds) || error("δ's seeds meet the prediction's; draw again")
            delta_correction(ms, haldane_values(ms, ref), rs; horizon, quantity, intervention,
                             threads)
        else
            delta === :none ? nothing : delta
        end

        # V's own cells, drawn after every other seed, so that without them the
        # draws are as before.
        nv = n_viability_cells
        vseeds = reshape(_cell_seeds(stream, nv * K * M), nv, K, M)
        isdisjoint(vseeds, seeds) && (dc === nothing || isdisjoint(vseeds, dc.seeds)) ||
            error("V's seeds meet the prediction's; draw again")

        f = _forward_setup(ms, forward; horizon)
        θI = [draw_values(pd, r) for r in rows]
        plain = [vcat(θI[m], θH[k]) for m in 1:M, k in 1:K]
        base = map(x -> haldane_values(ms, x), plain)
        # Every baseline run as one pool (μ̂'s cells, and V's own after them on
        # each pair's seeds), then V on each pair, then the intervened runs of the
        # pairs V keeps. Each phase runs its cells threaded or not, with the same
        # result.
        both = _run_matrix(f, base, cat(seeds, vseeds; dims = 1), quantity; threads)
        cells, vcells = nv == 0 ? (both, fill(CellRun[], M, K)) :
                        (map(c -> c[1:nc], both), map(c -> c[(nc + 1):end], both))
        decisions = Matrix{NamedTuple}(undef, M, K)
        for i in eachindex(cells)
            decisions[i] = conditioned ?
                           assess_viability(nv == 0 ? cells[i] : vcells[i], viability) :
                           (; viable = true, fraction = 1.0, counts = Dict{Symbol, Int}())
        end
        cells_i = if intervention === nothing
            fill(CellRun[], M, K)
        else
            vals_i = [decisions[i].viable ? intervene(intervention, ms, plain[i]) : nothing
                      for i in CartesianIndices(plain)]
            _run_matrix(f, vals_i, seeds, quantity; threads)
        end
        labels = forward === :hybrid ? reduction_declarations(ms, driver) :
                                       reduction_declarations(ms, f.route)
        settings = (; forward, completion, estimand, quantity = _quantity_label(quantity),
                    horizon = f.horizon, n_parameter_draws = M, n_held = K, width,
                    n_cells = nc, n_viability_cells = nv, intervention,
                    delta = delta isa DeltaCorrection ? :cached : delta,
                    n_ref = delta === :ensemble ? n_ref : 0,
                    viability = viability === nothing ? nothing : viability.label)
        _assemble(ms, pd, f, quantity, rows, θI, hnames, θH, seeds, cells, cells_i, decisions,
                  vseeds, vcells, completion, viability, intervention, dc, hp, fork, settings,
                  labels, threads ? Threads.nthreads() : 1)
    end
end

# What each pair ran, from what, and what became of it (KIM-33).
function _pair_provenance(pd, rows, seeds, cells, cells_i, decisions, vseeds, vcells,
                          completion, m, k)
    failed(cs) = [(seed = c.seed, failure = c.failure, message = c.message)
                  for c in cs if c.failure !== :none]
    r = rows[m]
    return (; draw = r, chain = pd.chain[r], iteration = pd.iteration[r],
            theta_H_source = completion === :none ? :none :
                             completion === :prior ? :prior : :viability_conditioned,
            held_draw = k, seeds = seeds[:, k, m], viability_seeds = vseeds[:, k, m],
            viable = decisions[m, k].viable, failed = failed(cells[m, k]),
            failed_viability = failed(vcells[m, k]), failed_intervened = failed(cells_i[m, k]))
end

function _assemble(ms, pd, f, quantity, rows, θI, hnames, θH, seeds, cells, cells_i,
                   decisions, vseeds, vcells, completion, viability, intervention, dc, hp,
                   fork, settings, labels, nthreads)
    M, K = size(cells)
    mom = [_moments([c.q for c in cells[i]]) for i in CartesianIndices(cells)]
    mu, s2, nok = first.(mom), getindex.(mom, 2), last.(mom)
    viable = [d.viable for d in decisions]
    if viability === nothing
        u = fill(1 / (M * K), M, K)
        weights = (; uniform = u)
        conditioning = (; uniform = (; acceptance = 1.0, ess = _kish(u)))
    else
        vw = viability_weights(viable)
        weights = (; per_theta_I = vw.per_theta_I, joint = vw.joint)
        conditioning = (; per_theta_I = (; acceptance = vw.acceptance, ess = vw.ess_per_theta_I,
                                         failed_rows = vw.failed_rows,
                                         row_acceptance = vw.row_acceptance),
                        joint = (; acceptance = vw.acceptance, ess = vw.ess_joint))
    end
    shift = dc === nothing ? 0.0 : dc.delta
    summary = map(w -> _summary(mu, s2, nok, w, shift, f.forward), weights)

    effect = fill(NaN, M, K)
    es2 = fill(NaN, M, K)
    esummary = nothing
    if intervention !== nothing
        enok = zeros(Int, M, K)
        for i in CartesianIndices(cells)
            isempty(cells_i[i]) && continue
            effect[i], es2[i], enok[i] = _moments(_paired_delta(cells[i], cells_i[i]))
        end
        eshift = dc === nothing ? 0.0 : dc.delta_intervened - dc.delta
        esummary = map(weights) do w
            (; _summary(effect, es2, enok, w, eshift, f.forward)...,
             n_failed_intervened = _n_failed_intervened(cells, cells_i, w))
        end
    end
    statuses, statuses_intervened = _status_counts(cells), _status_counts(cells_i)
    pairs = [_pair_provenance(pd, rows, seeds, cells, cells_i, decisions, vseeds, vcells,
                              completion, m, k) for m in 1:M, k in 1:K]
    provenance = (; fork_seed = fork, version = pkgversion(@__MODULE__), settings,
                  source = pd.source, reductions = labels, skipped = hp.skipped, pairs,
                  threads = nthreads)
    return Prediction(f.forward, completion, :population_mean, _quantity_label(quantity),
                      f.horizon, rows, θI,
                      pd.sigma_b === nothing ? nothing : pd.sigma_b[rows], hnames, θH,
                      seeds, cells, mu, s2, nok, intervention, cells_i, effect, es2,
                      viability, viable, decisions, vseeds, vcells, weights, conditioning, dc,
                      summary, esummary, statuses, statuses_intervened, provenance)
end

export Prediction, predict
