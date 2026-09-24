"""
Core A′ transcript decay: seventeen first-order reactions, one per transcript.
Spec §11 phase 12.

The third Core A′ module in the stochastic block. It owns no transcript: the
seventeen it decrements belong to [`CoreATranscription`](@ref), and it writes
them through [`written_states`](@ref) rather than through an edge, because a
transcript is not a registry species and no edge can name one (§12, amendment
of 2026-09-04).

**One rate parameter for all seventeen genes.** The published law is one global
catalytic constant over transcript length, so every half-life is a function of
length alone (`MinCell_CMEODE.py:444-463`):

```python
kcat = (18/452)*88 #1/s # INSTEAD OF 18 or 20
k_deg = kcat / n_tot
```

The trailing comment is the authors' own record of hand-tuning. It is inherited
here as a published value, not as a measured one.

**Each firing returns the transcript's four monomers and accrues its energy
cost**, as the published hook does (`in_out.py:178-186`): one `ATP_mRNAdeg` per
nucleotide, and one NMP per base. The NMP counters credit their pools at the
handshake. AMP and GMP are live pools that nucleotide recycling closes, through
ADK1 and GK1. CMP and UMP have no Core A′ species, so their counters credit the
CTP and UTP chemostats, which absorb them. That is ours, and `reduction_notes`
says so.
"""

"""
The published decay constant's numerator, `(18/452)·88` nt/s
(`MinCell_CMEODE.py:457`). The upstream line carries the comment `# INSTEAD OF
18 or 20`, which records hand-tuning rather than a measurement.

**Naming trap.** The upstream file also defines a variable called `krnadeg`,
as `0.00578/2` at line 272, and never uses it. Spec §4 D11's `krnadeg` target
means *this* constant, the one that runs. An equally dead `rnaDegRate` of the
same value sits at line 295, and that is the constant the archived
transcription design's predicted band was computed from (§12, 2026-09-10 B).
"""
const RNADEG_KCAT = (18 / 452) * 88

"""
The five decay counters, the registry pool each one touches, and the
direction.

`ATP_mRNAdeg` is the decay energy, one ATP per nucleotide, debited from ATP. Its
products, ADP and phosphate, are credited one each per ATP actually paid (spec
§11 task 13.10).

The four NMP counters are productions, with no consumer. AMP and GMP credit the
pools nucleotide recycling owns. CMP and UMP credit the CTP and UTP chemostats,
because the registry has no CMP or UMP species; the chemostat absorbs them with
no owner behind it (task 13.9), and `chemostat_census` counts what it took.
"""
const DECAY_COUNTERS = (
    (counter = :ATP_mRNAdeg, species = :M_atp_c, direction = :in,
     base = nothing, produces = (:M_adp_c, :M_pi_c)),
    (counter = :AMP_mRNAdeg, species = :M_amp_c, direction = :out,
     base = :A, produces = ()),
    (counter = :GMP_mRNAdeg, species = :M_gmp_c, direction = :out,
     base = :G, produces = ()),
    (counter = :CMP_mRNAdeg, species = :M_ctp_c, direction = :out,
     base = :C, produces = ()),
    (counter = :UMP_mRNAdeg, species = :M_utp_c, direction = :out,
     base = :U, produces = ()),
)

"""
    CoreATranscriptDecay(; genes, counters, clip, edges, written)

Seventeen decay jumps, `mRNA_g → n_g·ATP_mRNAdeg + #A·AMP + #G·GMP + #C·CMP +
#U·UMP`, at propensity `krnadeg / n_g · mRNA_g`.

`genes` defaults to the phase 10 extract, so the two modules cannot disagree on
a base count. `written` is overridable so a test can remove the declaration and
watch the refusal it guards; `edges` and `clip` mirror transcription's
constructor.

`counters` may be any subset of [`DECAY_COUNTERS`](@ref), and each entry must
match it exactly: a counter's increment is fixed by its name, so the pool it
touches and the direction are too.
"""
struct CoreATranscriptDecay <: AbstractSubModel
    genes::Vector{TranscriptionGene}
    params::Vector{InferParameter}
    counters::Vector{NamedTuple}
    edges::Vector{CouplingEdge}
    written::Vector{Symbol}
end

function CoreATranscriptDecay(; genes = read_transcription_genes(),
                              counters = DECAY_COUNTERS,
                              clip = :clamped_deficit_carried,
                              edges = nothing,
                              written = nothing)
    genes = collect(TranscriptionGene, genes)
    counters = collect(NamedTuple, counters)
    # Mirrors transcription's check (§12, 2026-09-23 B review): a GMP counter
    # aimed at AMP would otherwise credit guanine into the adenylate pool.
    allunique(c.counter for c in counters) || throw(ArgumentError(
        "Decay counters must be named once each, got " *
        join((string(":", c.counter) for c in counters), ", ")))
    for c in counters
        i = findfirst(d -> d.counter === c.counter, DECAY_COUNTERS)
        i === nothing && throw(ArgumentError(
            "Counter :$(c.counter) is not one decay accrues. The per-event " *
            "amounts are defined for " *
            join((string(":", d.counter) for d in DECAY_COUNTERS), ", ") * " only"))
        c == DECAY_COUNTERS[i] || throw(ArgumentError(
            "Counter :$(c.counter) must be declared as $(DECAY_COUNTERS[i]), got " *
            "$c. Its increment is fixed by its name, so its pool and direction are too"))
    end

    # The one rate parameter, and free: it is a §4 D11 target, identified by
    # the transcript autocorrelation independently of promoter amplitude. The
    # prior width is ours, because the published value carries no uncertainty.
    params = InferParameter[InferParameter(
        RNADEG_KCAT, LogNormal(log(RNADEG_KCAT), log(2.0)), false,
        :krnadeg, :CoreATranscriptDecay, :rate,
        ParameterSource("MinCell_CMEODE.py"; identifier = "DegradationRate",
                        informedness = :asserted))]
    for c in counters
        push!(params, InferParameter(0.0, Poisson(1.0), true,
                                     Symbol(c.counter, "0"),
                                     :CoreATranscriptDecay, :initial_condition))
    end

    default = CouplingEdge[]
    for c in counters
        push!(default, DeferredCounterEdge(species = c.species, direction = c.direction,
                                           counter = c.counter, clip = clip))
        for p in c.produces
            push!(default, DeferredCounterEdge(species = p, direction = :out,
                                               counter = c.counter, clip = clip))
        end
    end

    names = [transcript_state(g.locus) for g in genes]
    return CoreATranscriptDecay(
        genes, params, counters,
        collect(CouplingEdge, edges === nothing ? default : edges),
        collect(Symbol, written === nothing ? names : written))
end

# --- the protocol -----------------------------------------------------------

# The counters follow the constructor's `counters`, not `DECAY_COUNTERS`, so a
# custom set gets its own states — the defect phase 10b fixed in transcription
# (§12, 2026-09-23 H).
states(m::CoreATranscriptDecay) = [c.counter for c in m.counters]
parameters(m::CoreATranscriptDecay) = m.params
formalism(::CoreATranscriptDecay) = :jump
inference_mode(::CoreATranscriptDecay) = :simulation
coupling(m::CoreATranscriptDecay) = m.edges
inputs(m::CoreATranscriptDecay) = [transcript_state(g.locus) for g in m.genes]
written_states(m::CoreATranscriptDecay) = m.written

"""
    reactions(m::CoreATranscriptDecay) -> Vector{Reaction}

One first-order jump per gene, reading and decrementing the transcript
transcription owns. The propensity is zero at zero copies, so decay cannot
drive a count negative. Every reaction reads the same parameter slot, the one
free parameter `krnadeg`.
"""
function reactions(m::CoreATranscriptDecay)
    idx = Dict(c.counter => k for (k, c) in enumerate(m.counters))
    atp = get(idx, :ATP_mRNAdeg, 0)
    bases = [(k, c.base) for (k, c) in enumerate(m.counters) if c.base !== nothing]
    rxns = Reaction[]
    for (i, g) in enumerate(m.genes)
        len = g.length
        incs = [(k, g.counts[b]) for (k, b) in bases]
        push!(rxns, Reaction(
            (u, p, t, w) -> p[1] / len * w[i],
            (u, w) -> begin
                w[i] -= 1
                atp == 0 || (u[atp] += len)
                for (k, n) in incs
                    u[k] += n
                end
            end))
    end
    return rxns
end

function reduction_notes(::CoreATranscriptDecay)
    return [
        "Decay's CMP and UMP are credited into the CTP and UTP chemostats. " *
        "Core A′ has no CMP or UMP species, and the published model's " *
        "M_cmp_c and M_ump_c sit in modules this reduction cuts, so the " *
        "chemostat absorbs them. A monophosphate credited to a triphosphate " *
        "pool is not a reaction; it is an exemption, and it stops being " *
        "defensible the moment either pool is made live.",
    ]
end

# --- queries ----------------------------------------------------------------

"""
    counter_drains(m::CoreATranscriptDecay) -> Vector{NamedTuple}

Each decay counter, the registry species it touches, its direction, and what a
debit produces. Only `ATP_mRNAdeg` is a debit; the other four are productions.
"""
counter_drains(m::CoreATranscriptDecay) = copy(m.counters)

"""
    transcript_decay_constants(m) -> Vector{Pair{Symbol,Float64}}

Each gene's first-order constant, `krnadeg / n_g`, keyed by locus.
"""
transcript_decay_constants(m::CoreATranscriptDecay) =
    [g.locus => m.params[1].value / g.length for g in m.genes]

"""
    transcript_half_lives(m) -> Vector{Pair{Symbol,Float64}}

Each transcript's half-life in seconds, `ln 2 · n_g / krnadeg`. A function of
length alone.
"""
transcript_half_lives(m::CoreATranscriptDecay) =
    [g.locus => log(2) * g.length / m.params[1].value for g in m.genes]

"""
    monomer_closure(genes, n0, n_end, polymerised, returned) -> NamedTuple

The per-moiety residual of decay's own check (task 12.6):

```
polymerised_x + Σ_g n0_g·x_g − returned_x − Σ_g n_end_g·x_g
```

for each base `x`. Transcripts present at the start were never polymerised in
the run, and transcripts present at the end have not yet been returned, so the
stocks enter on both sides. Every term is an integer, so a closed composition
gives exactly zero.

`genes` is the extract, the transcription side's view of each transcript.
`polymerised` and `returned` are `(A, C, G, U)` named tuples of totals, read
from transcription's monomer counters and decay's NMP counters.

**Jump-only runs.** Those counters are run totals only when nothing clears
them. In a hybrid composition the hook clears every counter at each drain, so
their end values are one interval's accrual and the residual here is
meaningless.
"""
function monomer_closure(genes, n0, n_end, polymerised, returned)
    stock(n, b) = sum(n[i] * genes[i].counts[b] for i in eachindex(genes))
    return NamedTuple{(:A, :C, :G, :U)}(Tuple(
        polymerised[b] + stock(n0, b) - returned[b] - stock(n_end, b)
        for b in (:A, :C, :G, :U)))
end

"""
    assert_monomer_closure(genes, n0, n_end, polymerised, returned)

Throw, naming every moiety whose residual is not exactly zero, or return the
residuals.
"""
function assert_monomer_closure(genes, n0, n_end, polymerised, returned)
    r = monomer_closure(genes, n0, n_end, polymerised, returned)
    bad = [b for b in keys(r) if r[b] != 0]
    isempty(bad) || throw(ErrorException(
        "Monomer closure fails for moiety $(join(bad, ", ")): " *
        join(("$b residual $(r[b])" for b in bad), ", ") * ". What decay " *
        "returns does not equal what transcription polymerised, net of the " *
        "transcripts still alive"))
    return r
end

# Exported from this file, as transcription does, so parallel module branches
# do not all append to one list in `InferCell.jl`.
export CoreATranscriptDecay, DECAY_COUNTERS, RNADEG_KCAT
export transcript_decay_constants, transcript_half_lives, monomer_closure,
       assert_monomer_closure
