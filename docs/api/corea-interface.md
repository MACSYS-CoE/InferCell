# Core A′ interface contract

Source: [`src/corea/registry.jl`](https://github.com/MACSYS-CoE/InferCell/blob/main/src/corea/registry.jl), [`src/corea/edges.jl`](https://github.com/MACSYS-CoE/InferCell/blob/main/src/corea/edges.jl), [`src/corea/resolver.jl`](https://github.com/MACSYS-CoE/InferCell/blob/main/src/corea/resolver.jl), [`src/corea/loader.jl`](https://github.com/MACSYS-CoE/InferCell/blob/main/src/corea/loader.jl), [`src/corea/labels.jl`](https://github.com/MACSYS-CoE/InferCell/blob/main/src/corea/labels.jl).

Core A′ is the reduced JCVI-syn3A model scoped in `dev/notes/reduced-syn3a-scoping.md`: glycolysis through lactate, the PTS glucose cascade, lactate export, a live ATP/GTP energy interface, and three moiety-recycling reactions. It is built by seven independent modules that touch each other only through a coupling boundary.

This page documents the contract that boundary is written in. It exists so that a module author working alone can declare how their block touches the others, and have that declaration checked before anything runs — rather than discovering at the join that two modules meant different things by the same word.

## The state registry

Every Core A′ species is named once and ordered once. Modules index into the registry; none defines its own species names.

```julia
COREA_SPECIES              # Vector{SpeciesEntry}, 32 dynamic states then 5 chemostats
species_index(:M_atp_c)    # canonical position
species_entry(:M_gtp_c)    # the registry row
species_group(:M_amp_c)    # :adenylate — the conserved moiety
is_chemostatted(:M_ctp_c)  # true
held_value(:M_glc__D_e)    # 40.0 mM
```

Names are the published model's BiGG-style identifiers (`M_g6p_c`, `M_13dpg_c`). They map one-to-one onto the source model, and unlike the display names of the scoping note's tables (`13DPG`, `3PG`) every one is a valid Julia identifier.

`SpeciesEntry` fields:

- `name::Symbol` — the canonical identifier.
- `group::Symbol` — `:glycolytic`, `:adenylate`, `:guanylate`, `:redox`, `:other`, `:trna`, `:pts`, `:chemostat`. This is the conserved moiety a dead-end report names.
- `treatment::Symbol` — `:dynamic` or `:chemostatted`. A module that declares dynamics for a chemostat is rejected at composition.
- `regime::Symbol` — `:mrna`, `:protein` or `:metabolite`. Core A′ spans nearly five orders of magnitude of copy number, which is why it needs several formalisms at once.
- `initial_value`, `gstd`, `source_file`, `informedness` — the imported value and where it came from.

!!! warning "Read initial conditions from the file that owns the species"
    GTP, GDP, AMP and GMP are nucleotide-module species. The central balanced file shows them at the 0.1 mM prior default because they are out of *its* module's scope; the nucleotide file has real balanced values (GTP at 1.6627 mM). Reading the wrong file understates the guanylate pool by an order of magnitude. The registry records the source file per entry for exactly this reason.

Three states are marked `:prior_default` — G3P, PPi and cytosolic lactate sit at the prior median at prior width, and nothing informed them. `is_informed` and `uninformed_species` distinguish them from measurements.

## The seven edge kinds

The legend of [`fig1r_state_graph_reduced.pdf`](https://github.com/MACSYS-CoE/InferCell/blob/main/dev/notes/figures/reduced-syn3a-coupling/fig1r_state_graph_reduced.pdf) enumerates seven genuinely different ways state crosses an interface. The figure is the source of truth for this list: if the figure and these types disagree, one of them is a bug.

| Kind | What crosses | Beyond species and direction |
|---|---|---|
| `MassEdge` | shared state, continuous; gradients cross inside the ODE block | — |
| `CurrencyEdge` | routed via a shared pool node | `pool` |
| `DeferredCounterEdge` | a cost accrued in one block, debited a step later | `counter`, `clip` |
| `CatalyticEdge` | counts entering a rate law as parameters; no mass flows | `param_slot` |
| `RateConstantEdge` | pools re-entering the stochastic block as rate constants | `cadence`, `interval` |
| `VolumeEdge` | counts set surface area, hence volume, hence every concentration | — |
| `ClampedEdge` | a real dependence replaced by a constant | `held_value`, `origin` |

Every edge names a registry species and a direction — `:in` where the declaring module consumes, `:out` where it produces — and optionally a `peer` module. A `peer` of `nothing` resolves against whichever module owns the species, which is what a single module under test should use.

```julia
coupling(m::MyModule) = [
    MassEdge(species=:M_atp_c, direction=:out),
    DeferredCounterEdge(species=:M_atp_c, direction=:in, counter=:ATP_trsc),
]
```

Omitting a field a kind requires fails at construction, naming both the field and the kind — not at composition, and not at run time.

### The two boundary hazards are fields, not hard-coded behaviour

The scoping note leaves two questions open, and both are properties of the boundary rather than of any module. Each becomes a declared field whose default is what the published model does, so departing from the published model is something an author has to write down.

**The clipped drain.** The published expression-cost drain floors the pool at zero and carries the shortfall forward — a `max(0, ·)` on the interface that will obstruct gradients if the ODE block is sampled with NUTS.

```julia
DeferredCounterEdge(species=:M_atp_c, direction=:in, counter=:ATP_trsc)
# clip defaults to :clamped_deficit_carried — the published model
# alternatives: :unclamped, :smoothed
```

`obstructs_gradients` reports whether an edge is non-differentiable, and `check_gradient_safety` warns when a composition containing a differentiable sub-model carries one. It warns rather than throws: a clamped interface is a legitimate model, just not a differentiable one.

**The rebuild cadence.** The published ODE→stochastic channel is piecewise-constant — rate constants recomputed once a minute, not propensities reading pools continuously.

```julia
RateConstantEdge(species=:M_gtp_c, direction=:out)
# cadence defaults to :piecewise_constant at interval = 60.0 — the published model
# alternative: cadence=:continuous, which carries no interval and is marked a deviation
```

## Resolving a composition

`resolve_coupling` validates a set of modules and returns a `CouplingGraph`. It is standalone: it needs no `ODEProblem`, works on a single module, and produces its reports without running a simulation.

```julia
graph = resolve_coupling([central, transport, nucleotide])
graph.edges                  # every declared edge, resolved
graph.dead_ends              # consumed species with no declared producer
graph.unowned_states         # registry states no module integrates
graph.chemostat_exemptions   # costs exempted because the registry holds the pool
graph.gradient_obstructions
graph.deviations
```

It **throws** when the boundary is inconsistent: an edge naming an unregistered species, an edge whose named peer is absent from the composition, two modules declaring the same species and direction with different kinds, a module integrating a chemostat, a state integrated twice, a declared cost with no state that can pay it, or an `inputs` list that has drifted from the module's own inbound edges.

It **reports** what is merely incomplete: unowned states and dead ends, which a partial composition is expected to have.

!!! note "Every declared cost needs a state that can pay it and a reaction that returns it"
    This is the mechanical form of a rule the scoping note derived from two errors of record. The lumped charging step converts ATP to AMP, and without ADK1 nothing returns it — 3.48M charging events against an adenylate pool of ~79,800 particles exhausts it after 2.3% of the cell cycle. Transcription buries GTP as GMP and, without GK1, that GMP is stranded. `dead_end_report` names the conserved moiety of each stranded species, because those two are one class of error, not two.

    The check is structural: it validates declarations, not stoichiometry. A module that declares a producer and does not implement one still passes. That is why each wave-1 module carries its own conservation check as part of its own work.

## Parameter provenance

Every imported value records the file it came from. This is not bookkeeping neatness: the central-versus-nucleotide cross-file trap has already produced two errors of record, once in kinetics and once in initial concentrations.

```julia
central    = read_source_table("path/to/central.tsv";    file="central_balanced")
nucleotide = read_source_table("path/to/nucleotide.tsv"; file="nucleotide_balanced")

disagreements([central, nucleotide])
# AmbiguousValue("kcat_fwd_PGK3", ["central_balanced" => 319.5,
#                                  "nucleotide_balanced" => 140.8], false)

load_parameter([central, nucleotide], "kcat_fwd_PGK3";
               name=:kcat_PGK3, module_id=:Nucleotide,
               prior=LogNormal(0.0, 1.0),
               governing="nucleotide_balanced")
```

An identifier present in more than one file **must** declare which file governs. Without it the load fails, naming the identifier, both files and both values — an unresolved ambiguity is an error at load time, never a silent choice. With it, the parameter's provenance records both the file chosen and the ones rejected, so the decision is visible in the assembled model rather than only in the module's source. `governing_choices` enumerates them.

Parsing is SBtab-shaped tab-separated text, handled with Base only: the `!!SBtab` declaration and `%` comments are skipped, and the header's leading `!` is stripped. The value column defaults to `Mode`, because that is the column the published simulator reads — an earlier version of the scoping note used a balancing-distribution column instead and produced a fictitious capacity bottleneck.

How well a value is known is recorded as `:balanced`, `:prior_default`, `:asserted` or `:not_imported`. `:asserted` means the source carries a point value with no quantified uncertainty, so any prior on it is this project's rather than the model's — the ten PTS mass-action constants are the case that matters.

## Labelling what is ours

Several Core A′ treatments are ours rather than the published model's, and each is a place where a result could depend on our choice. `reduction_declarations` enumerates them so the label can reach any report that depends on one.

```julia
reduction_declarations(models)   # Vector{ReductionLabel}
reduction_report(models)         # the same, written for a human
```

It collects reduction-introduced clamps, smoothed deferred counters, continuous rate-constant edges, parameters whose prior this project asserted, and any lumping a module registers through `reduction_notes`. A composition that follows the published model throughout reports that nothing departs from it.
