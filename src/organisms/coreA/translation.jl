"""
Core A′ translation: seventeen reactions, one per transcript, plus the
translocation of ptsG into the membrane. Spec §11 phase 11.

The second Core A′ module in the stochastic block. It reads the transcripts
[`CoreATranscription`](@ref) owns and publishes the protein counts the rest of
the model reads: the thirteen enzyme counts fill the metabolic modules' enzyme
slots through catalytic edges those modules declare (phase 11a), and the four
phosphotransferase carriers are credited to [`PtsTransport`](@ref)'s states
through producer counters.
"""

"""
    read_translation_residues([path]) -> Dict{Symbol,Int}

Each locus's residue count, from the `!Residues` column of the gene extract
phase 10 vendors (`src/organisms/coreA/data/transcription_genes.tsv`).

A residue is an amino acid translation charges for: the transcript translated
under NCBI table 4, less its one stop codon. The extract's generator counts the
translation and asserts it equals `length/3 − 1`, so this reader asserts the
same rather than trusting either (spec §12, 2026-09-23 D). The column is read
by name, not position, so phase 10's positional columns are untouched.
"""
function read_translation_residues(path::AbstractString = TRANSCRIPTION_EXTRACT)
    header = nothing
    out = Dict{Symbol, Int}()
    for line in eachline(path)
        s = strip(line)
        (isempty(s) || startswith(s, "!!") || startswith(s, "%")) && continue
        f = split(s, '\t')
        if startswith(s, "!")
            header = [strip(c, '!') for c in f]
            continue
        end
        col = header === nothing ? nothing : findfirst(==("Residues"), header)
        col === nothing && throw(ArgumentError(
            "$path has no !Residues column. Regenerate it with " *
            "dev/scripts/extract_transcription_genes.jl; see " *
            "src/organisms/coreA/data/README.md"))
        out[Symbol(f[1])] = parse(Int, f[col])
    end
    for g in read_transcription_genes(path)
        haskey(out, g.locus) || throw(ArgumentError(
            "$path holds no residue count for $(g.locus) ($(g.reaction))"))
        out[g.locus] == g.length ÷ 3 - 1 && g.length % 3 == 0 || throw(ArgumentError(
            "$(g.locus) ($(g.reaction)): $(out[g.locus]) residues, but a " *
            "$(g.length)-nucleotide transcript with one terminal stop has " *
            "$(g.length ÷ 3 - 1)"))
    end
    return out
end

# ---------------------------------------------------------------------------
# The published rate law's constants: `translation_rate_restart.py`, the law
# the 60 s rebuild runs for every minute of the cycle after the first.
# ---------------------------------------------------------------------------

"""
Ribosome turnover, `riboKcat`, in residues per second.

12, as both translation-rate files that run set it
(`translation_rate_start.py:18`, `translation_rate_restart.py:19`). The 10 at
`MinCell_CMEODE.py:332` is a module variable no rate law reads: `TranslatRate`
is imported from the rate files, which define their own.
"""
const RIBO_KCAT = 12.0

"`riboK0`, mM (`translation_rate_restart.py:20`)."
const RIBO_K0 = 4 * 25e-6

"""
`riboKd`, mM: 1e-3, the restart file's value (`translation_rate_restart.py:21`).

The start file has 1e-4, and it governs only the first minute. The published
model runs `MinCell_CMEODE.py` for one minute and then `MinCell_restart.py`
for every minute after, and the restart file's law is the one the 60 s
rebuild evaluates. This module's rebuild is that channel, so it takes the
restart law whole: this `riboKd` and the restart form of `kcat_mod` in
[`translation_kcat`](@ref). Chosen 2026-09-24 and recorded in spec §12.
"""
const RIBO_KD = 1e-3

"Ribosomes per cell, converted to mM at the registry's initial volume."
const RIBOSOME_COPIES = 503

"""
The number of per-amino-acid charged-tRNA pools the lumped pool stands in for.

The published law sums `n_aa · riboKd / [aa-tRNA]` over twenty pools and adds
`riboKd² / [fMet-tRNA]²`. Core A′ has one lumped charged pool, and each of
those twenty-one concentrations reads its per-amino-acid share,
`[M_trna_chg_c] / 20`. That keeps the published law's structure: at the
nominal 0.2 mM charged pool the share is 0.01 mM against upstream's 150
copies, 0.0074 mM. Substituting the whole pool for each per-amino-acid pool
would move each term's concentration 27-fold, and with it the pool's
elasticity, from about 0.09 to about 0.005, for no reason but the lumping.
**Ours**, chosen 2026-09-24 (spec §12), and `reduction_notes` says so.
"""
const TL_AA_TYPES = 20

"Ribosomes per transcript never exceed this (`translation_rate_restart.py:103`)."
const POLYSOME_CAP = 15

"""
    ribosomes_per_transcript(length) -> Int

The published polysome size, `min(15, max(1, round(length/125 − 1)))`.
`round` ties to even in both Python 3 and Julia, and no transcript length
lands on a tie, since that would need `length/125` to end in .5.
"""
ribosomes_per_transcript(n::Integer) = min(POLYSOME_CAP, max(1, round(Int, n / 125 - 1)))

"""
    translation_kcat(length) -> Float64

The restart law's `kcat_mod`: `(0.25·n + 0.2)·riboKcat` for a polysome of
`n > 1` ribosomes, and `0.45·riboKcat` for one. The start file's `+0.25` would
apply for the first minute only (see [`RIBO_KD`](@ref)).
"""
function translation_kcat(n::Integer)
    r = ribosomes_per_transcript(n)
    return r > 1 ? (0.25 * r + 0.2) * RIBO_KCAT : 0.45 * RIBO_KCAT
end

"""
    translation_rate_constant(length, residues, chg_mM; ribo_conc) -> Float64

The published restart law (`translation_rate_restart.py:46-115`), equation 3,
with the lumped pool's per-amino-acid share `c = chg_mM / 20` standing in for
each of the twenty-one charged-tRNA concentrations:

```
k = kcat_mod / ( (1 + K₀/[ribosome])·K_d²/c² + residues·K_d/c + residues )
```

`residues` is the published `n_tot − 1`, since `n_tot` counts the stop codon,
and it is also `Σ n_aa`, since the stop is no amino acid. So the per-gene
composition drops out of the lumped law, which needs only the length and the
residue count.

The share is floored at one particle, at the registry's volume, as upstream
floors each pool with `max(1, count)`. Without it an exhausted charged pool
would divide by zero. `ribo_conc` is frozen at the registry's volume, as
upstream computes `ribosomeConc` once at the initial radius.

Written once, as a free function, so tests assert the module's own arithmetic.
"""
function translation_rate_constant(n::Integer, residues::Integer, chg_mM::Real;
                                   ribo_conc::Real = RIBOSOME_COPIES /
                                                     corea_particles_per_mM())
    c = max(chg_mM / TL_AA_TYPES, 1 / corea_particles_per_mM())
    denom = (1 + RIBO_K0 / ribo_conc) * RIBO_KD^2 / c^2 +
            residues * RIBO_KD / c + residues
    return translation_kcat(n) / denom
end

"""
    translation_elasticity(length, residues, chg_mM; ribo_conc) -> Float64

`∂ ln k / ∂ ln [M_trna_chg_c]`, analytically: the charged-pool terms' share of
the denominator, with the fMet term counted twice because it is squared. This
is the charged-tRNA reverse channel's gain (spec §9, R1, task 11.10).

Below the one-particle floor of [`translation_rate_constant`](@ref) the
constant no longer depends on the pool, so the elasticity there is zero. It is
not a small-signal quantity near the floor: at 0.02 mM it is about 0.5, not
the nominal pool's 0.09.
"""
function translation_elasticity(n::Integer, residues::Integer, chg_mM::Real;
                                ribo_conc::Real = RIBOSOME_COPIES /
                                                  corea_particles_per_mM())
    c = chg_mM / TL_AA_TYPES
    c < 1 / corea_particles_per_mM() && return 0.0
    a = (1 + RIBO_K0 / ribo_conc) * RIBO_KD^2 / c^2
    b = residues * RIBO_KD / c
    return (2a + b) / (a + b + residues)
end

# ---------------------------------------------------------------------------
# The module.
# ---------------------------------------------------------------------------

"ptsG, the one membrane protein, and so the one gene with a translocation step."
const TL_PTSG = :JCVISYN3A_0779

"""
The upstream constant of `TranslocRate`, `50 / ptnLen` per second
(`MinCell_CMEODE.py:574`), a typical secY translocation rate. `ptnLen` is
`len(aasequence)`, which counts the stop codon.
"""
const TRANSLOC_KCAT = 50.0

"""
The seven deferred cost counters translation accrues, what each touches, and
what fills it.

- `GTP_translat` is the elongation energy, two GTP per residue
  (`MinCell_CMEODE.py:1008-1010`; upstream names it `ATP_translat`, and its
  hook drains GTP, `in_out.py:199-218`). Its GDP and phosphate are recorded in
  `produces` and not yet credited: one accrual cannot credit two products until
  task 13.10 (spec §12, 2026-09-23). Two per residue, and residues exclude the
  stop, so this is two GTP per protein fewer than upstream's
  `2·len(aasequence)`. That is ours.
- `tRNA_translat` is the charged tRNA consumed, one per residue. It **debits
  the charged pool and credits the uncharged one**, one counter with two
  edges. The hook credits what the debit actually paid, so a clipped debit
  cannot create tRNA (spec §11 task 11.6).
- `ATP_transloc` is translocation's energy, `int(len(aasequence)/10)` ATP per
  ptsG inserted (`MinCell_CMEODE.py:1012-1013`), debited from ATP, products
  pending 13.10 like GTP's.
- The four carrier counters credit one new carrier each to the
  unphosphorylated state `PtsTransport` owns (spec §12, 2026-09-23 C). The
  three cytosolic carriers are credited as they are translated; ptsG as it is
  translocated, since only a membrane ptsG is a carrier.

`per` says what one event adds: `:twice_residues`, `:residues`, `:transloc_atp`
or `:one`. `event` says which event: `:translation`, or `:translocation`.
"""
const TRANSLATION_COUNTERS = (
    (counter = :GTP_translat, debits = :M_gtp_c, credits = nothing,
     produces = (:M_gdp_c, :M_pi_c), event = :translation, locus = nothing,
     per = :twice_residues),
    (counter = :tRNA_translat, debits = :M_trna_chg_c, credits = :M_trna_c,
     produces = (), event = :translation, locus = nothing, per = :residues),
    (counter = :ATP_transloc, debits = :M_atp_c, credits = nothing,
     produces = (:M_adp_c, :M_pi_c), event = :translocation, locus = TL_PTSG,
     per = :transloc_atp),
    (counter = :ptsI_translat, debits = nothing, credits = :M_ptsi_c,
     produces = (), event = :translation, locus = :JCVISYN3A_0233, per = :one),
    (counter = :ptsH_translat, debits = nothing, credits = :M_ptsh_c,
     produces = (), event = :translation, locus = :JCVISYN3A_0694, per = :one),
    (counter = :Crr_translat, debits = nothing, credits = :M_crr_c,
     produces = (), event = :translation, locus = :JCVISYN3A_0234, per = :one),
    (counter = :ptsG_transloc, debits = nothing, credits = :M_ptsg_c,
     produces = (), event = :translocation, locus = TL_PTSG, per = :one),
)

"""
    protein_state(locus) -> Symbol

The state carrying one gene's protein count, `P_<locus>`: the name phase 11a's
catalytic edges read (`default_protein_sources`, `recycling_enzymes`). For ptsG
it is the *membrane* count, the one translocation increments.
"""
protein_state(locus::Symbol) = Symbol("P_", locus)

"The cytosolic ptsG count, made by translation and consumed by translocation."
const TL_PTSG_CYTO = Symbol("Pcyto_", TL_PTSG)

"""
    translation_rate_param(locus) -> Symbol

The name of a gene's translation rate constant, which the 60 s rebuild fills.
"""
translation_rate_param(locus::Symbol) = Symbol("k_tl_", locus)

"""
    CoreATranslation(; genes, residues, counters, clip, interval, chg_mM, edges, rebuilt)

Seventeen translation reactions and one translocation.

```
mRNA_g → mRNA_g + P_g + 2r_g·GTP_translat + r_g·tRNA_translat [+ carrier credit]
Pcyto_ptsG → P_ptsG + int((r+1)/10)·ATP_transloc + ptsG_transloc
```

Translation is first order in the transcript, which it reads and leaves
unchanged: transcripts are [`CoreATranscription`](@ref)'s states, read through
`inputs`. Translocation is first order in the cytosolic ptsG count.

`residues` defaults to the extract's `!Residues` column and is overridable so
task 11.8's mutation can halve one gene's. `counters` may be any subset of
[`TRANSLATION_COUNTERS`](@ref), each entry matching it exactly, and the states
and edges follow the subset. `chg_mM` is the charged pool the constants start
at, the nominal 0.2 mM that `TrnaCharging` initialises; a composed run
overwrites them at every rebuild. `edges` and `rebuilt` are overridable so a
test can build a mis-declared variant and watch the refusal.
"""
struct CoreATranslation <: AbstractSubModel
    genes::Vector{TranscriptionGene}
    residues::Vector{Int}
    params::Vector{InferParameter}
    counters::Vector{NamedTuple}
    edges::Vector{CouplingEdge}
    rebuilt::Vector{Symbol}
    ribo_conc::Float64
    transloc_k::Float64
end

function CoreATranslation(; genes = read_transcription_genes(),
                          residues = nothing,
                          counters = TRANSLATION_COUNTERS,
                          clip = :clamped_deficit_carried,
                          interval = 60.0,
                          chg_mM = CHARGING_POOL_DEFAULTS.charged_fraction *
                                   CHARGING_POOL_DEFAULTS.pool_mM,
                          edges = nothing,
                          rebuilt = nothing)
    genes = collect(TranscriptionGene, genes)
    res = residues === nothing ?
          (r = read_translation_residues(); [r[g.locus] for g in genes]) :
          collect(Int, residues)
    length(res) == length(genes) || throw(ArgumentError(
        "$(length(res)) residue counts for $(length(genes)) genes"))
    counters = collect(NamedTuple, counters)
    allunique(c.counter for c in counters) || throw(ArgumentError(
        "Translation counters must be named once each, got " *
        join((string(":", c.counter) for c in counters), ", ")))
    for c in counters
        i = findfirst(d -> d.counter === c.counter, TRANSLATION_COUNTERS)
        i === nothing && throw(ArgumentError(
            "Counter :$(c.counter) is not one translation accrues. The per-event " *
            "amounts are defined for " *
            join((string(":", d.counter) for d in TRANSLATION_COUNTERS), ", ") * " only"))
        c == TRANSLATION_COUNTERS[i] || throw(ArgumentError(
            "Counter :$(c.counter) must be declared as $(TRANSLATION_COUNTERS[i]), " *
            "got $c. Its increment is fixed by its name, so what it debits and " *
            "credits is too"))
    end
    iptsg = findfirst(g -> g.locus === TL_PTSG, genes)
    iptsg === nothing && throw(ArgumentError(
        "Translation needs ptsG ($TL_PTSG) among its genes: it is the gene the " *
        "translocation reaction belongs to"))

    ribo_conc = RIBOSOME_COPIES / corea_particles_per_mM()
    transloc_k = TRANSLOC_KCAT / (res[iptsg] + 1)

    params = InferParameter[]
    # The seventeen rate constants, free because the rebuild fills them, at
    # the nominal charged pool. Not inference targets (§4 D11): a rebuilt slot
    # is overwritten at every refresh, and its prior is ours.
    for (g, r) in zip(genes, res)
        k = translation_rate_constant(g.length, r, chg_mM; ribo_conc = ribo_conc)
        push!(params, InferParameter(
            k, LogNormal(log(k), log(2.0)), false,
            translation_rate_param(g.locus), :CoreATranslation, :rate))
    end

    # The ribosome globals. Fixed, for D13's reason: they are the hub of the
    # published star, and a free one would be a ridge through all seventeen
    # constants. The ribosome concentration is computed once at the initial
    # volume, as upstream computes it, and is one of the genuinely frozen
    # quantities; K₀ and K_d are dissociation constants in mM, so volume does
    # not enter them. They are held on the struct and in constants, and listed
    # here so they are enumerable with their provenance.
    src(id) = ParameterSource("translation_rate_restart.py"; identifier = id,
                              informedness = :asserted)
    for (nm, v, id) in ((:tl_ribo_kcat, RIBO_KCAT, "riboKcat"),
                        (:tl_ribo_k0, RIBO_K0, "riboK0"),
                        (:tl_ribo_kd, RIBO_KD, "riboKd"),
                        (:tl_ribo_conc, ribo_conc, "ribosomeConc"))
        push!(params, InferParameter(v, LogNormal(log(v), log(2.0)), true,
                                     nm, :CoreATranslation, :rate, src(id)))
    end
    push!(params, InferParameter(
        transloc_k, LogNormal(log(transloc_k), log(2.0)), true,
        :tl_transloc_k, :CoreATranslation, :rate,
        ParameterSource("MinCell_CMEODE.py"; identifier = "TranslocRate",
                        informedness = :asserted)))

    # Initial conditions, fixed: the proteomics counts (§4 D8 holds protein
    # initial conditions fixed), no cytosolic ptsG, and empty counters.
    for g in genes
        push!(params, InferParameter(
            float(g.ptn_count), Poisson(g.ptn_count), true,
            Symbol(protein_state(g.locus), "0"), :CoreATranslation,
            :initial_condition))
    end
    push!(params, InferParameter(0.0, Poisson(1.0), true, Symbol(TL_PTSG_CYTO, "0"),
                                 :CoreATranslation, :initial_condition))
    for c in counters
        push!(params, InferParameter(0.0, Poisson(1.0), true, Symbol(c.counter, "0"),
                                     :CoreATranslation, :initial_condition))
    end

    default = CouplingEdge[
        RateConstantEdge(species = :M_trna_chg_c, direction = :in,
                         cadence = :piecewise_constant, interval = interval)]
    for c in counters
        c.debits === nothing ||
            push!(default, DeferredCounterEdge(species = c.debits, direction = :in,
                                               counter = c.counter, clip = clip))
        c.credits === nothing ||
            push!(default, DeferredCounterEdge(species = c.credits, direction = :out,
                                               counter = c.counter, clip = clip))
    end

    return CoreATranslation(
        genes, res, params, counters,
        collect(CouplingEdge, edges === nothing ? default : edges),
        collect(Symbol, rebuilt === nothing ?
                [translation_rate_param(g.locus) for g in genes] : rebuilt),
        ribo_conc, transloc_k)
end

# --- the protocol -----------------------------------------------------------

states(m::CoreATranslation) =
    vcat([protein_state(g.locus) for g in m.genes], [TL_PTSG_CYTO],
         [c.counter for c in m.counters])
parameters(m::CoreATranslation) = m.params
formalism(::CoreATranslation) = :jump
inference_mode(::CoreATranslation) = :simulation
coupling(m::CoreATranslation) = m.edges
inputs(m::CoreATranslation) = [transcript_state(g.locus) for g in m.genes]
rebuilt_params(m::CoreATranslation) = m.rebuilt

"""
    rate_constants(p, t, m::CoreATranslation, pools) -> SVector{17}

The seventeen constants the 60 s rebuild writes, from the live charged pool,
`pools[1]` in mM. Reads nothing from `p`, so it cannot compound its own
previous value. **This module never calls this**; only the driver does.
"""
function rate_constants(p, t, m::CoreATranslation, pools)
    chg = pools[1]
    n = length(m.genes)
    return SVector{n, Float64}(ntuple(n) do i
        translation_rate_constant(m.genes[i].length, m.residues[i], chg;
                                  ribo_conc = m.ribo_conc)
    end)
end

# What one event adds to counter `c`: gene `i`'s translation, or ptsG's
# translocation (`i = 0`).
function _tl_increment(c, m::CoreATranslation, event::Symbol, i::Int)
    c.event === event || return 0
    if event === :translation
        c.locus === nothing || c.locus === m.genes[i].locus || return 0
        r = m.residues[i]
        return c.per === :twice_residues ? 2r : c.per === :residues ? r : 1
    else
        r = m.residues[findfirst(g -> g.locus === TL_PTSG, m.genes)]
        return c.per === :transloc_atp ? (r + 1) ÷ 10 : 1
    end
end

"""
    reactions(m::CoreATranslation) -> Vector{Reaction}

Seventeen translation jumps at `k_g · mRNA_g`, reading the transcript through
`inputs` and leaving it unchanged, then one translocation jump at
`tl_transloc_k · Pcyto_ptsG`. Translating ptsG adds to the cytosolic count;
every other gene adds to its own `P_<locus>`.
"""
function reactions(m::CoreATranslation)
    n = length(m.genes)
    cyto = n + 1
    nc = length(m.counters)
    iptsg = findfirst(g -> g.locus === TL_PTSG, m.genes)
    rxns = Reaction[]
    for i in 1:n
        target = i == iptsg ? cyto : i
        incs = [(cyto + k, _tl_increment(m.counters[k], m, :translation, i)) for k in 1:nc]
        incs = filter(x -> x[2] != 0, incs)
        push!(rxns, Reaction(
            (u, p, t, w) -> p[i] * w[i],
            (u, w) -> begin
                u[target] += 1
                for (k, a) in incs
                    u[k] += a
                end
            end))
    end
    kt = m.transloc_k
    tincs = filter(x -> x[2] != 0,
                   [(cyto + k, _tl_increment(m.counters[k], m, :translocation, 0)) for k in 1:nc])
    push!(rxns, Reaction(
        (u, p, t, w) -> kt * u[cyto],
        (u, w) -> begin
            u[cyto] -= 1
            u[iptsg] += 1
            for (k, a) in tincs
                u[k] += a
            end
        end))
    return rxns
end

function reduction_notes(::CoreATranslation)
    return [
        "Translation's rate law reads the lumped charged-tRNA pool through its " *
        "per-amino-acid share, [M_trna_chg_c]/$TL_AA_TYPES, in each of the " *
        "twenty-one concentrations the published law reads per amino acid " *
        "(translation_rate_restart.py). This keeps the law's structure: 0.01 mM " *
        "at the nominal pool against upstream's 150 copies, 0.0074 mM. Reading " *
        "the whole pool instead would cut the charged-tRNA elasticity about " *
        "twenty-fold. The lumping itself is TrnaCharging's declaration, not this one.",

        "Translation runs the published restart law (riboKd 1e-3, kcat_mod " *
        "(0.25n + 0.2)·riboKcat) from the start of the cycle. The published " *
        "model runs the start law (riboKd 1e-4, +0.25) for its first minute only.",

        "A residue excludes the stop codon, the convention k_chg was calibrated " *
        "on. Upstream charges GTP on len(aasequence), which counts the stop, so " *
        "this module charges two GTP per protein fewer than the published model.",
    ]
end

# --- queries ----------------------------------------------------------------

"The genes this module carries, in canonical order."
translation_genes(m::CoreATranslation) = m.genes

"""
    residue_counts(m) -> Vector{Pair{Symbol,Int}}

Each gene's residue count, keyed by locus: what one translation event charges
for.
"""
residue_counts(m::CoreATranslation) = [g.locus => r for (g, r) in zip(m.genes, m.residues)]

"""
    translation_rate_constants(m[, chg_mM]) -> Vector{Pair{Symbol,Float64}}

The seventeen rate constants keyed by locus, at the nominal charged pool or at
any supplied.
"""
translation_rate_constants(m::CoreATranslation,
                           chg_mM = CHARGING_POOL_DEFAULTS.charged_fraction *
                                    CHARGING_POOL_DEFAULTS.pool_mM) =
    [g.locus => translation_rate_constant(g.length, r, chg_mM; ribo_conc = m.ribo_conc)
     for (g, r) in zip(m.genes, m.residues)]

"""
    counter_drains(m::CoreATranslation) -> Vector{NamedTuple}

Each configured counter, the pool it debits and the pool it credits (either
may be `nothing`), and what a debit produces but does not yet credit.
"""
counter_drains(m::CoreATranslation) = copy(m.counters)

"""
    proteins_made(m, u0, u) -> Vector{Int}

Each gene's proteins translated between two states of this module's own slice:
the rise in its `P_<locus>` count, plus, for ptsG, the rise in the cytosolic
count, since a translated ptsG is cytosolic until it is translocated.

**Jump-only runs.** Nothing lowers a protein count but translocation, which
moves ptsG from one of its two states to the other, so this is exact there. In
a hybrid run it is too, since the hook reads and clears the counters but never
the protein counts.
"""
function proteins_made(m::CoreATranslation, u0, u)
    n = length(m.genes)
    made = [Int(u[i] - u0[i]) for i in 1:n]
    iptsg = findfirst(g -> g.locus === TL_PTSG, m.genes)
    made[iptsg] += Int(u[n + 1] - u0[n + 1])
    return made
end

"""
    residue_energy_closure(m, made, gtp; residues) -> Int

Task 11.8's residual: the GTP counter's accrual less twice the residues
translated,

```
GTP_translat − 2 · Σ_g made_g · r_g
```

where `r_g` comes from the **extract**, not from the module, so a module that
charges the wrong residue count cannot close its own check. Every term is an
integer, so a closed module gives exactly zero. `made` is
[`proteins_made`](@ref) over the same interval `gtp` accrued over; in a
jump-only run that is the whole run, since nothing clears the counter.
"""
function residue_energy_closure(m::CoreATranslation, made, gtp;
                                residues = read_translation_residues())
    ref = [residues[g.locus] for g in m.genes]
    return Int(gtp) - 2 * sum(made[i] * ref[i] for i in eachindex(ref))
end

"""
    assert_residue_energy_closure(m, made, gtp; residues) -> Int

Throw if [`residue_energy_closure`](@ref) is not exactly zero, naming every
gene whose single firing charges other than twice its extract residue count,
or return the zero residual. The counter holds one aggregate, so the gene is
found by firing each reaction once on an empty state rather than inferred from
the sum.
"""
function assert_residue_energy_closure(m::CoreATranslation, made, gtp;
                                       residues = read_translation_residues())
    r = residue_energy_closure(m, made, gtp; residues = residues)
    r == 0 && return r
    k = findfirst(c -> c.counter === :GTP_translat, m.counters)
    k === nothing && throw(ArgumentError(
        "This CoreATranslation carries no GTP_translat counter, so there is " *
        "no energy accrual to close against"))
    idx = length(m.genes) + 1 + k
    bad = String[]
    for (i, g) in enumerate(m.genes)
        u = zeros(Int, length(states(m)))
        reactions(m)[i].affect!(u, zeros(Int, length(m.genes)))
        u[idx] == 2 * residues[g.locus] ||
            push!(bad, "$(g.locus) ($(g.reaction)): $(u[idx]) GTP per protein " *
                       "against 2 × $(residues[g.locus])")
    end
    throw(ErrorException(
        "Residue-to-energy closure fails by $r GTP: the counter does not equal " *
        "twice the residues translated. " *
        (isempty(bad) ? "Every gene's single firing charges twice its residues, " *
                        "so the proteins-made tally disagrees with the run" :
                        "Genes charging the wrong amount: " * join(bad, "; "))))
end

export proteins_made, residue_energy_closure, assert_residue_energy_closure
export CoreATranslation, TRANSLATION_COUNTERS, TL_PTSG, TL_PTSG_CYTO, TRANSLOC_KCAT
export protein_state, translation_rate_param, translation_genes, residue_counts,
       translation_rate_constants
export read_translation_residues, ribosomes_per_transcript, translation_kcat,
       translation_rate_constant, translation_elasticity
export RIBO_KCAT, RIBO_K0, RIBO_KD, RIBOSOME_COPIES, TL_AA_TYPES, POLYSOME_CAP
