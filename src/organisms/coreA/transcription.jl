"""
Core A′ transcription: seventeen genes, one constitutive reaction each.

The first Core A′ module in the stochastic block, and the only one that reads
the ODE block's pools — through the 60 s rate-constant rebuild rather than
through `inputs`, because the two blocks hold separate state vectors. Spec §11
phase 10.

Three things about this module are not derivable from the published source, and
each is registered rather than absorbed:

  * **The base-to-nucleotide mapping is corrected** (§4 D4). `TranscriptRate`
    assigns three of four base counts to the wrong pool, because its local
    names run A, U, C, G while the dictionary keys run A, C, G, U. Correcting
    it is a departure, so `:corrected` is what `reduction_notes` labels while
    the published permutation — available behind `base_mapping = :published` —
    does not appear there. This is the unusual case where the default is the
    labelled one.
  * **CTP and UTP are chemostatted by this reduction**, not by the published
    model, because the modules that make them are outside Core A′.
  * **Promoter strength is a proxy** — protein copy number over 180 — inherited
    verbatim, and conditioning on proteomics would close a loop through it
    (§4 D5). Registered twice, once as the proxy and once as the circularity.

**The transcripts are not registry species and must not become them.** The
registry is frozen at 32 dynamic states plus 5 chemostats; a transcript is a
count, not a concentration, and `_resolve_edges` refuses an edge naming a
non-registry species. That is why decay, in phase 12, writes a transcript
through [`written_states`](@ref) rather than through an edge (§12, amendment of
2026-09-04).

**Genes are a fixed quantity, not a state.** Core A′ cuts replication, so
nothing can change a gene's copy number, and seventeen constant states would be
seventeen wasted dimensions in every stored trajectory. The published model
carries them as species because it has `rep_start.py`; adding replication later
must promote them, which `reduction_notes` records.
"""

# ---------------------------------------------------------------------------
# The published rate law's constants (`MinCell_CMEODE.py:259-292`).
# ---------------------------------------------------------------------------

"Polymerase turnover, `0.155·187/493·20` nt/s."
const RNAPOL_KCAT = 0.155 * 187 / 493 * 20

"`rnaPolK0`, mM."
const RNAPOL_K0 = 1e-4

"`rnaPolKd`, mM."
const RNAPOL_KD = 0.1

"RNA polymerase copies per cell, converted to mM at the registry's volume."
const RNAP_COPIES = 187

"""
The promoter proxy's divisor. Not a fitted quantity — protein copy number over
this is the published back-calculation of promoter strength (§4 D5).
"""
const PROMOTER_DIVISOR = 180

"""
The turnover ceiling, `2*90` nt/s. It never binds for these seventeen genes —
the largest promoter-scaled turnover is GAPD's 8.85 — and the module asserts
that rather than assuming it, via [`turnover_headroom`](@ref).
"""
const TURNOVER_CEILING = 2 * 90

"""
The five deferred cost counters, and the registry pool each debits.

`ATP_trsc` is the polymerisation energy, one per nucleotide; the four monomer
counters are the incorporated bases. **ATP is therefore debited twice per
transcript** — 1284 + 521 = 1805 for PGI. That is the published model's
accounting, not double-counting, and `reduction_notes` says so because the
alternative is someone "fixing" it.
"""
const TRANSCRIPTION_COUNTERS = (
    (counter = :ATP_trsc, species = :M_atp_c, produces = (:M_adp_c, :M_pi_c)),
    (counter = :ATP_mRNA, species = :M_atp_c, produces = (:M_ppi_c,)),
    (counter = :GTP_mRNA, species = :M_gtp_c, produces = (:M_ppi_c,)),
    (counter = :CTP_mRNA, species = :M_ctp_c, produces = (:M_ppi_c,)),
    (counter = :UTP_mRNA, species = :M_utp_c, produces = (:M_ppi_c,)),
)

"""
The four nucleotide pools in the order the rate-constant edges declare them,
which is the order [`rate_constants`](@ref) receives `pools` in.
"""
const TRANSCRIPTION_POOLS = (:M_atp_c, :M_ctp_c, :M_gtp_c, :M_utp_c)

"""
Which base count is charged against which nucleotide.

`:corrected` charges each base against the nucleotide it is actually
polymerised from. `:published` reproduces `MinCell_CMEODE.py:377-383`, where
the UTP term takes the cytosine count, the CTP term the guanine count and the
GTP term the uracil count. The permutation moves each rate constant by about
one percent, but it weights the GTP term by the uracil count rather than the
guanine count, which makes each constant more sensitive to the live GTP pool
than it should be — by that gene's U-to-G count ratio, a median 1.85× over the
seventeen and spanning 1.53 to 2.45. GTP is one of the four reverse channels
into this module; ATP, CTP and UTP enter the same rate law.
"""
const BASE_MAPPINGS = (
    corrected = (A = :A, C = :C, G = :G, U = :U),
    published = (A = :A, C = :G, G = :U, U = :C),
)

"""
The first-minute rate-constant inputs of `MinCell_CMEODE.py:275-278`, in
(ATP, CTP, GTP, UTP) order. **Never the concentrations this module runs on**
(§4 D2) — the 60 s rebuild replaces them with live pools, and the nucleotide
file's `Compound` table repeats two of them, which makes the wrong choice look
corroborated. Held here so a test can assert the module took none of them.
"""
const SETUP_CONSTANTS = (M_atp_c = 1.04, M_ctp_c = 0.34, M_gtp_c = 0.68, M_utp_c = 0.68)

"""
    TranscriptionGene

One row of `src/organisms/coreA/data/transcription_genes.tsv`, plus the
metabolic reaction the gene encodes. The reaction is carried so task 10.7's
cross-check against the metabolic modules has something to name; the join key
is the locus.
"""
struct TranscriptionGene
    locus::Symbol
    reaction::Symbol
    length::Int
    counts::NamedTuple{(:A, :C, :G, :U), NTuple{4, Int}}
    first_two::NTuple{2, Symbol}
    ptn_count::Int
    mean_mrna::Float64
end

"""
The metabolic reaction each locus encodes: ten glycolytic, four
phosphotransferase, three nucleotide-recycling. Not upstream data — this is
Core A′'s own scope, and it is what lets task 10.7 report a mismatch against
the metabolic modules by reaction name rather than by bare locus.
"""
const TRANSCRIPTION_LOCI = (
    (:JCVISYN3A_0445, :PGI), (:JCVISYN3A_0220, :PFK), (:JCVISYN3A_0131, :FBA),
    (:JCVISYN3A_0727, :TPI), (:JCVISYN3A_0607, :GAPD), (:JCVISYN3A_0606, :PGK),
    (:JCVISYN3A_0729, :PGM), (:JCVISYN3A_0213, :ENO), (:JCVISYN3A_0221, :PYK),
    (:JCVISYN3A_0475, :LDH_L), (:JCVISYN3A_0233, :ptsI), (:JCVISYN3A_0234, :Crr),
    (:JCVISYN3A_0694, :ptsH), (:JCVISYN3A_0779, :ptsG), (:JCVISYN3A_0651, :ADK1),
    (:JCVISYN3A_0344, :PPA), (:JCVISYN3A_0203, :GK1),
)

"""
    TRANSCRIPTION_EXTRACT

Path to this module's vendored gene extract. Regenerated, never fetched at run
time — compute nodes have no network. See
`src/organisms/coreA/data/README.md`.
"""
const TRANSCRIPTION_EXTRACT =
    joinpath(@__DIR__, "data", "transcription_genes.tsv")

# The extract's nine read columns, in the order `read_transcription_genes`
# indexes them. Asserted against the file's own header so a reorder fails
# loudly instead of loading one base count as another.
const TRANSCRIPTION_COLUMNS = ("ID", "Length", "A", "C", "G", "U", "First2",
                               "PtnCount", "MeanMRNA")

"""
    read_transcription_genes([path]) -> Vector{TranscriptionGene}

The vendored per-gene extract, in [`TRANSCRIPTION_LOCI`](@ref) order.

Read directly rather than through [`load_parameter`](@ref): each column has
exactly one upstream source, so there is no cross-file ambiguity for the
governing machinery to arbitrate (design D8). Every failure names the locus,
because a silently absent gene would surface only as a model with sixteen
transcripts.
"""
function read_transcription_genes(path::AbstractString = TRANSCRIPTION_EXTRACT)
    isfile(path) || throw(ArgumentError(
        "No transcription extract at $path. Regenerate it with " *
        "`julia --project=dev/scripts dev/scripts/extract_transcription_genes.jl " *
        "<Minimal_Cell checkout> $path`; see src/organisms/coreA/data/README.md"))

    rows = Dict{Symbol, NTuple{9, String}}()
    header = nothing
    for line in eachline(path)
        s = strip(line)
        (isempty(s) || startswith(s, "!!") || startswith(s, "%")) && continue
        f = split(s, '\t')
        if startswith(s, "!")
            header = [strip(c, '!') for c in f]
            # The columns are read positionally below, so a reordered extract
            # would load silently wrong data: swapping !A and !C leaves the
            # base-count sum invariant intact and shifts every rate constant.
            # Checking the names is what makes the positions safe.
            header[1:length(TRANSCRIPTION_COLUMNS)] == collect(TRANSCRIPTION_COLUMNS) ||
                throw(ArgumentError(
                    "$path: the first $(length(TRANSCRIPTION_COLUMNS)) columns " *
                    "are $(header[1:min(end, length(TRANSCRIPTION_COLUMNS))]), " *
                    "not $(collect(TRANSCRIPTION_COLUMNS)). The reader takes " *
                    "them by position, so a reorder would load one base count " *
                    "as another"))
            continue
        end
        header === nothing && throw(ArgumentError(
            "$path: a data row appears before the `!`-prefixed header"))
        length(f) == length(header) || throw(ArgumentError(
            "$path: row '$(f[1])' has $(length(f)) fields against the header's " *
            "$(length(header)); a truncated row would otherwise load as a gene " *
            "with a missing base count"))
        rows[Symbol(f[1])] = (f[1], f[2], f[3], f[4], f[5], f[6], f[7], f[8], f[9])
    end

    genes = TranscriptionGene[]
    for (locus, reaction) in TRANSCRIPTION_LOCI
        haskey(rows, locus) || throw(ArgumentError(
            "$path holds no row for $locus ($reaction). Core A′ carries " *
            "$(length(TRANSCRIPTION_LOCI)) genes and this extract would build " *
            "a model with $(length(rows) < length(TRANSCRIPTION_LOCI) ? "fewer" : "different") ones"))
        r = rows[locus]
        n = parse(Int, r[2])
        counts = (A = parse(Int, r[3]), C = parse(Int, r[4]),
                  G = parse(Int, r[5]), U = parse(Int, r[6]))
        total = counts.A + counts.C + counts.G + counts.U
        total == n || throw(ArgumentError(
            "$locus ($reaction): the four base counts sum to $total but the " *
            "transcript length is $n"))
        f2 = r[7]
        length(f2) == 2 || throw(ArgumentError(
            "$locus ($reaction): First2 is '$f2', which is not two bases"))
        push!(genes, TranscriptionGene(locus, reaction, n, counts,
                                       (Symbol(f2[1]), Symbol(f2[2])),
                                       parse(Int, r[8]), parse(Float64, r[9])))
    end
    return genes
end

"""
    transcript_state(locus) -> Symbol

The state name carrying one gene's transcript count. Prefixed so it cannot
collide with a registry species name, and keyed on the locus so phases 11 and
12 name the same state without sharing code with this module.
"""
transcript_state(locus::Symbol) = Symbol("mRNA_", locus)

"""
    promoter_param(locus) / rate_param(locus) -> Symbol

Names of a gene's two free parameters: the promoter strength, which is the
inference target, and the rate constant, which the 60 s rebuild fills.
"""
promoter_param(locus::Symbol) = Symbol("S_", locus)
rate_param(locus::Symbol) = Symbol("k_tx_", locus)
@doc (@doc promoter_param) rate_param

"""
    CoreATranscription(; base_mapping = :corrected, seed = nothing, ...)

Seventeen constitutive transcription reactions.

`base_mapping` selects the corrected mapping (the default, and the labelled
one) or the published permutation. `seed`, when given, draws each transcript's
initial count from `Poisson(measured mean)` as `MinCell_CMEODE.py:630` does;
without it every transcript starts at the rounded measured mean, so a test that
is not about initialisation is deterministic.

Every keyword that names a declaration is overridable so a test can build a
mis-declared variant and witness the refusal.
"""
struct CoreATranscription <: AbstractSubModel
    genes::Vector{TranscriptionGene}
    params::Vector{InferParameter}
    edges::Vector{CouplingEdge}
    rebuilt::Vector{Symbol}
    base_mapping::Symbol
    conc::NamedTuple{(:M_atp_c, :M_ctp_c, :M_gtp_c, :M_utp_c), NTuple{4, Float64}}
    rnap_conc::Float64
    counters::Tuple   # the configured cost counters, in state order
end

"""
    TRANSCRIPTION_NTP_SOURCES

The four balanced nucleotide concentrations of spec §4 D3, each with its
geometric standard deviation and the file that governs it. ATP is the central
file's; the other three are the nucleotide file's, which governs every
nucleotide-module species (§4 D2). Reading the wrong file understates the
guanylate pool by an order of magnitude.
"""
const TRANSCRIPTION_NTP_SOURCES = (
    (species = :M_atp_c, value = 3.6529, gstd = 1.2825, file = "central_balanced"),
    (species = :M_ctp_c, value = 0.6874, gstd = 2.0574, file = "nucleotide_balanced"),
    (species = :M_gtp_c, value = 1.6627, gstd = 1.5684, file = "nucleotide_balanced"),
    (species = :M_utp_c, value = 2.7681, gstd = 1.3664, file = "nucleotide_balanced"),
)

function CoreATranscription(; genes = read_transcription_genes(),
                            base_mapping::Symbol = :corrected,
                            seed = nothing,
                            interval = 60.0,
                            cadence = :piecewise_constant,
                            pools = TRANSCRIPTION_POOLS,
                            counters = TRANSCRIPTION_COUNTERS,
                            clip = :clamped_deficit_carried,
                            edges = nothing,
                            rebuilt = nothing)
    haskey(BASE_MAPPINGS, base_mapping) || throw(ArgumentError(
        "base_mapping must be :corrected or :published, got :$base_mapping. " *
        ":corrected charges each base against the nucleotide it is polymerised " *
        "from; :published reproduces the permutation of MinCell_CMEODE.py:377-383"))

    counters = Tuple(counters)
    for c in counters
        haskey(COUNTER_INCREMENTS, c.counter) || throw(ArgumentError(
            "Counter :$(c.counter) is not one transcription can charge. Each " *
            "event adds a fixed per-gene amount to each counter, and those " *
            "amounts are defined for $(join(string.(":", keys(COUNTER_INCREMENTS)), ", ")) only"))
    end

    genes = collect(TranscriptionGene, genes)
    conc = NamedTuple{(:M_atp_c, :M_ctp_c, :M_gtp_c, :M_utp_c)}(
        Tuple(s.value for s in TRANSCRIPTION_NTP_SOURCES))
    rnap_conc = RNAP_COPIES / corea_particles_per_mM()

    params = InferParameter[]

    # The promoter strengths: seventeen free parameters, and the seventeen the
    # inverse problem is actually about. Marked `:asserted` — the proxy carries
    # a point value with no quantified uncertainty behind it, so the prior
    # width is ours (§4 D5).
    for g in genes
        s = g.ptn_count / PROMOTER_DIVISOR
        push!(params, InferParameter(
            s, LogNormal(log(s), log(2.0)), false,
            promoter_param(g.locus), :CoreATranscription, :rate,
            ParameterSource("proteomics"; table = "Proteomics",
                            identifier = string(g.locus),
                            informedness = :asserted)))
    end

    # The seventeen rate constants, at this module's own concentrations. A
    # composed run overwrites them at every rebuild; a standalone run keeps
    # them, which is what makes the module self-consistent with the pool owner
    # before the driver exists.
    for g in genes
        k = transcription_rate_constant(g, conc; base_mapping = base_mapping,
                                        rnap_conc = rnap_conc)
        push!(params, InferParameter(
            k, LogNormal(log(k), log(2.0)), false,
            rate_param(g.locus), :CoreATranscription, :rate))
    end

    # The four nucleotide concentrations, with provenance. Fixed: they are the
    # pool owner's states, not this module's unknowns, and a composed run reads
    # them live. They are held on the struct as well, because a fixed parameter
    # does not reach the composed parameter vector and so cannot be read from
    # `p` inside the rate law.
    for s in TRANSCRIPTION_NTP_SOURCES
        push!(params, InferParameter(
            s.value, LogNormal(log(s.value), log(s.gstd)), true,
            Symbol("tx_conc_", s.species), :CoreATranscription, :rate,
            ParameterSource(s.file; table = "Quantity",
                            identifier = "conc_$(s.species)",
                            informedness = :balanced)))
    end

    # The polymerase globals. Fixed, and D13 is why: the published stochastic
    # block is a star in which ~19 scalars set all 455 rate constants, and this
    # module inherits its hub. A free `rnaPolKcat` would be a ridge through all
    # seventeen constants at once.
    for (nm, v) in ((:tx_rnapol_kcat, RNAPOL_KCAT), (:tx_rnapol_k0, RNAPOL_K0),
                    (:tx_rnapol_kd, RNAPOL_KD), (:tx_rnap_conc, rnap_conc))
        push!(params, InferParameter(v, LogNormal(log(v), log(2.0)), true,
                                     nm, :CoreATranscription, :rate,
                                     ParameterSource("MinCell_CMEODE.py";
                                                     identifier = string(nm),
                                                     informedness = :asserted)))
    end

    # Initial conditions, fixed. Transcripts start at the measured mean, or at
    # a seeded Poisson draw from it as the published model does.
    rng = seed === nothing ? nothing : MersenneTwister(seed)
    for g in genes
        n0 = rng === nothing ? round(Int, g.mean_mrna) :
             rand(rng, Poisson(g.mean_mrna))
        push!(params, InferParameter(
            float(n0), Poisson(g.mean_mrna), true,
            Symbol(transcript_state(g.locus), "0"), :CoreATranscription,
            :initial_condition))
    end
    for c in counters
        push!(params, InferParameter(0.0, Poisson(1.0), true,
                                     Symbol(c.counter, "0"),
                                     :CoreATranscription, :initial_condition))
    end

    default = CouplingEdge[]
    for s in pools
        push!(default, cadence === :continuous ?
              RateConstantEdge(species = s, direction = :in, cadence = :continuous) :
              RateConstantEdge(species = s, direction = :in, cadence = cadence,
                               interval = interval))
    end
    for c in counters
        push!(default, DeferredCounterEdge(species = c.species, direction = :in,
                                           counter = c.counter, clip = clip))
    end
    # CTP and UTP only. ATP and GTP are live pools the recycling module owns,
    # and clamping them would be a claim this module has no right to make.
    for s in TRANSCRIPTION_NTP_SOURCES
        s.species in (:M_ctp_c, :M_utp_c) || continue
        push!(default, ClampedEdge(species = s.species, direction = :in,
                                   origin = :ours, held_value = s.value))
    end

    return CoreATranscription(
        genes, params,
        collect(CouplingEdge, edges === nothing ? default : edges),
        collect(Symbol, rebuilt === nothing ?
                [rate_param(g.locus) for g in genes] : rebuilt),
        base_mapping, conc, rnap_conc, counters)
end

# What one transcription event of gene `g` adds to each counter: the transcript
# length as polymerisation energy, and each base count as incorporated monomer.
const COUNTER_INCREMENTS = Dict{Symbol, Function}(
    :ATP_trsc => g -> g.length,
    :ATP_mRNA => g -> g.counts.A,
    :GTP_mRNA => g -> g.counts.G,
    :CTP_mRNA => g -> g.counts.C,
    :UTP_mRNA => g -> g.counts.U,
)

"""
    transcription_rate_constant(g, conc; base_mapping, rnap_conc) -> Float64

The published rate law of `MinCell_CMEODE.py:351-391`, equation 3:

```
kcat_mod = min(rnaPolKcat · S_g , 180)              S_g = ptnCount / 180
NMonoSum = K_d · Σ_x  N_x / [NTP_x]
k_g      = kcat_mod / ( (1 + K₀/[RNAP])·K_d²/(C₁C₂) + NMonoSum + n_g − 1 )
```

`C₁` and `C₂` are the concentrations of the nucleotides the transcript's first
two bases are polymerised from. `N_x` is the count of the base charged against
NTP `x`, which is where [`BASE_MAPPINGS`](@ref) enters and where the published
code permutes three of four.

Written once, as a free function, so every test asserts the module's own
arithmetic rather than a second copy of it.
"""
function transcription_rate_constant(g::TranscriptionGene, conc;
                                     base_mapping::Symbol = :corrected,
                                     rnap_conc::Float64 = RNAP_COPIES /
                                                          corea_particles_per_mM())
    s = g.ptn_count / PROMOTER_DIVISOR
    kcat_mod = min(RNAPOL_KCAT * s, TURNOVER_CEILING)
    m = BASE_MAPPINGS[base_mapping]
    nmono = 0.0
    for (ntp, base) in ((:M_atp_c, m.A), (:M_ctp_c, m.C),
                        (:M_gtp_c, m.G), (:M_utp_c, m.U))
        nmono += g.counts[base] * RNAPOL_KD / conc[ntp]
    end
    c1 = conc[_ntp_of(g.first_two[1])]
    c2 = conc[_ntp_of(g.first_two[2])]
    denom = (1 + RNAPOL_K0 / rnap_conc) * RNAPOL_KD^2 / (c1 * c2) +
            nmono + g.length - 1
    return kcat_mod / denom
end

_ntp_of(base::Symbol) = base === :A ? :M_atp_c :
                        base === :C ? :M_ctp_c :
                        base === :G ? :M_gtp_c :
                        base === :U ? :M_utp_c :
                        throw(ArgumentError("no nucleotide for base :$base"))

# --- the protocol -----------------------------------------------------------

states(m::CoreATranscription) =
    vcat([transcript_state(g.locus) for g in m.genes],
         [c.counter for c in m.counters])

parameters(m::CoreATranscription) = m.params
formalism(::CoreATranscription) = :jump
inference_mode(::CoreATranscription) = :simulation
coupling(m::CoreATranscription) = m.edges
rebuilt_params(m::CoreATranscription) = m.rebuilt

"""
    rate_constants(p, t, m::CoreATranscription, pools) -> SVector{17}

The seventeen constants the 60 s rebuild writes, recomputed from the live
pools. `pools` arrive in [`TRANSCRIPTION_POOLS`](@ref) order, in mM.

Reads the promoter strengths from `p` and nothing else from it. The slots this
fills are read back from `p` on the next call, so a law that touched them would
compound its own previous value.

**This module never calls this.** The propensities read their parameter slot;
only the driver writes there. That is the whole content of "the module provides
the means and the driver does the scheduling".
"""
function rate_constants(p, t, m::CoreATranscription, pools)
    conc = NamedTuple{TRANSCRIPTION_POOLS}(Tuple(pools))
    n = length(m.genes)
    return SVector{n, Float64}(ntuple(n) do i
        g = m.genes[i]
        s = p[i]
        kcat_mod = min(RNAPOL_KCAT * s, TURNOVER_CEILING)
        mp = BASE_MAPPINGS[m.base_mapping]
        nmono = 0.0
        for (ntp, base) in ((:M_atp_c, mp.A), (:M_ctp_c, mp.C),
                            (:M_gtp_c, mp.G), (:M_utp_c, mp.U))
            nmono += g.counts[base] * RNAPOL_KD / conc[ntp]
        end
        c1 = conc[_ntp_of(g.first_two[1])]
        c2 = conc[_ntp_of(g.first_two[2])]
        kcat_mod / ((1 + RNAPOL_K0 / m.rnap_conc) * RNAPOL_KD^2 / (c1 * c2) +
                    nmono + g.length - 1)
    end)
end

"""
    reactions(m::CoreATranscription) -> Vector{Reaction}

One constant-rate jump per gene:

```
gene -> gene + mRNA + n·ATP_trsc + #A·ATP_mRNA + #C·CTP_mRNA
                    + #G·GTP_mRNA + #U·UTP_mRNA
```

The gene is reactant and product, so its count is invariant and the propensity
is the rate constant alone — transcription is zeroth order, which is part of
what makes the stochastic block's conditional likelihood closed form (§3).
"""
function reactions(m::CoreATranscription)
    n = length(m.genes)
    # The configured counters sit immediately after the transcripts, in the
    # order `states` lists them.
    nc = length(m.counters)
    rxns = Reaction[]
    for (i, g) in enumerate(m.genes)
        incs = ntuple(k -> COUNTER_INCREMENTS[m.counters[k].counter](g), nc)
        push!(rxns, Reaction(
            (u, p, t, _) -> p[n + i],
            (u, _) -> begin
                u[i] += 1
                for k in 1:nc
                    u[n + k] += incs[k]
                end
            end))
    end
    return rxns
end

function reduction_notes(m::CoreATranscription)
    loci = join(string.(rate_param(g.locus) for g in m.genes), ", ")
    mapping = m.base_mapping === :corrected ?
        "The base-to-nucleotide mapping is corrected: each base is charged " *
        "against the nucleotide it is polymerised from (:corrected), rather " *
        "than against the permutation of MinCell_CMEODE.py:377-383 " *
        "(:published), where the UTP term takes the cytosine count, the CTP " *
        "term the guanine count and the GTP term the uracil count. The rate " *
        "constants move about one percent; the sensitivity to the live GTP " *
        "pool moves by a median 1.85× across the seventeen, spanning 1.53 to " *
        "2.45 — it is the per-gene U-to-G count ratio, not one factor. GTP is " *
        "one of four reverse channels, not the only one: ATP, CTP and UTP " *
        "enter the same rate law, and the measured ATP elasticity matches " *
        "GTP's. Both mappings are computable." :
        "The base-to-nucleotide mapping is the published permutation of " *
        "MinCell_CMEODE.py:377-383 (:published), reproduced deliberately: the " *
        "UTP term takes the cytosine count, the CTP term the guanine count and " *
        "the GTP term the uracil count, so bases are not charged against the " *
        "nucleotides they are polymerised from. Relative to the corrected " *
        "mapping (:corrected, spec §4 D4) the rate constants differ by about " *
        "one percent, and the sensitivity to the live GTP pool is inflated by " *
        "the per-gene U-to-G count ratio, a median 1.85× spanning 1.53 to 2.45. " *
        "Both mappings are computable."
    counter_names = Set(c.counter for c in m.counters)
    double_debit = (:ATP_trsc in counter_names && :ATP_mRNA in counter_names) ?
        ["ATP is debited twice per transcript — once as polymerisation energy " *
         "at one per nucleotide (ATP_trsc), once as an incorporated monomer " *
         "(ATP_mRNA). This is the published model's accounting, not an error, " *
         "and it is recorded so that it is not \"fixed\"."] : String[]
    return [
        mapping,

        "CTP and UTP are held constant at 0.6874 and 2.7681 mM. This " *
        "reduction chemostats them, not the published model: the reactions " *
        "that make them are in modules Core A′ cuts. The assumption stops " *
        "being defensible the moment either pool is made live.",

        "Promoter strength is the published proteomics proxy, protein copy " *
        "number over $PROMOTER_DIVISOR, inherited verbatim. It is not a " *
        "measured parameter — it is back-calculated from the steady-state " *
        "protein abundance the model is meant to predict.",

        "Conditioning on proteomics would be circular. The same protein " *
        "counts that set these seventeen promoter strengths also set every " *
        "enzyme concentration in the metabolic modules, so a protein-level " *
        "observable closes a loop through them. Affected: $loci.",

        double_debit...,

        "Genes are carried as a fixed quantity rather than as states, because " *
        "Core A′ cuts replication and nothing can change them. The published " *
        "model carries them as species; adding replication must promote them.",

        "Transcription is constitutive: one rate constant per gene, no " *
        "promoter switching. The two-state identifiability result does not " *
        "transfer, because there is no on/off rate pair to be degenerate.",
    ]
end

# --- queries ----------------------------------------------------------------

"The genes this module carries, in canonical order."
transcription_genes(m::CoreATranscription) = m.genes

"""
    promoter_strengths(m) -> Vector{Pair{Symbol,Float64}}

Each gene's promoter strength, keyed by locus.
"""
promoter_strengths(m::CoreATranscription) =
    [g.locus => g.ptn_count / PROMOTER_DIVISOR for g in m.genes]

"""
    protein_copy_numbers(m) -> Dict{Symbol,Int}

The protein copy number behind each gene's promoter strength, keyed by locus.
One number with three consumers: this proxy and the enzyme concentrations of
the metabolic modules. Task 10.7's cross-check reads this.
"""
protein_copy_numbers(m::CoreATranscription) =
    Dict(g.locus => g.ptn_count for g in m.genes)

"""
    transcription_rate_constants(m[, conc]) -> Vector{Pair{Symbol,Float64}}

The seventeen rate constants keyed by locus, at this module's own
concentrations or at any four supplied.
"""
transcription_rate_constants(m::CoreATranscription, conc = m.conc) =
    [g.locus => transcription_rate_constant(g, conc;
                                            base_mapping = m.base_mapping,
                                            rnap_conc = m.rnap_conc)
     for g in m.genes]

"""
    turnover_headroom(m) -> Vector{Pair{Symbol,Float64}}

Each gene's promoter-scaled turnover, `rnaPolKcat · S_g`, against the ceiling
of $TURNOVER_CEILING nt/s. Reported rather than assumed, so "the cap never
binds" is a measurement.
"""
turnover_headroom(m::CoreATranscription) =
    [g.locus => RNAPOL_KCAT * g.ptn_count / PROMOTER_DIVISOR for g in m.genes]

"""
    counter_drains(m) -> Vector{NamedTuple}

Each cost counter, the registry species it debits, and what that drain
produces. `ATP_trsc` yields ADP and phosphate; the four monomer counters yield
pyrophosphate — a source the project's earlier accounting attributed to
amino-acid charging alone, and one phase 8's phosphate closure must see.
"""
counter_drains(m::CoreATranscription) = collect(m.counters)

"""
    transcript_decay_constant(g) -> Float64

The published first-order decay constant, `(18/452)·88 / n_g`
(`MinCell_CMEODE.py:444-460`) — one global catalytic constant over transcript
length.

**Phase 12 owns decay; this is here because phase 10's calibration check needs
it.** Without a consumer the transcripts only accumulate and every gene's
time-averaged count is an order of magnitude above its measured mean, so the
check would be vacuous. Recorded as a §12 amendment dated 2026-09-10; phase 12
supersedes this for composed runs and keeps it only for standalone ones.
"""
transcript_decay_constant(g::TranscriptionGene) = (18 / 452) * 88 / g.length

# Exported from this file rather than from `InferCell.jl`'s block, so that
# parallel Core A′ module branches do not all append to one list and conflict.
export CoreATranscription, TranscriptionGene
export read_transcription_genes, transcription_genes, transcription_rate_constant,
       transcription_rate_constants, promoter_strengths, protein_copy_numbers,
       turnover_headroom, counter_drains, transcript_decay_constant
export transcript_state, promoter_param, rate_param
export TRANSCRIPTION_LOCI, TRANSCRIPTION_COUNTERS, TRANSCRIPTION_POOLS,
       TRANSCRIPTION_NTP_SOURCES, TRANSCRIPTION_EXTRACT, BASE_MAPPINGS,
       SETUP_CONSTANTS
export RNAPOL_KCAT, RNAPOL_K0, RNAPOL_KD, RNAP_COPIES, PROMOTER_DIVISOR,
       TURNOVER_CEILING
