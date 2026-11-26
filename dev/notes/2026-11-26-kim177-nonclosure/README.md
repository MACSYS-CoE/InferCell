# What the full syn3A cell has that the mean field cannot close (KIM-177)

2026-11-26. Linear: [KIM-177](https://linear.app/kimpson/issue/KIM-177/non-closure-inventory-of-the-full-cell),
milestone M11. It feeds the design of M13's exact check (F3 panel c), the Spec's
assumptions and the route decision (KIM-187).

The counts here are printed by `nonclosure_inventory.py` in this directory. Its output
is committed beside it as `nonclosure_inventory.out.txt`. The few numbers it does not
print were checked by hand; "Checked by hand" at the end lists each with its source. To
rerun it on the login node:

```bash
/fred/oz022/tkimpson/miniconda3/bin/python3 -I nonclosure_inventory.py > nonclosure_inventory.out.txt
```

It uses the standard library only. It asserts the issue's checks: 363 second-order
reactions, 493 replication reactions, 1,434 count-dependent rate rules, and every
reaction in exactly one family. It also asserts the exact-closure counts of section 1.

Revised on 2026-11-26 after the PR #204 review. The review found that 177 of the DnaA
reactions, and the first-round initiations, close exactly although `MeanField` rejects
them. This version keeps two terms for closure apart throughout (see Terms).

## Summary

- The full cell has 363 second-order reactions: 279 DnaA at the origins, 62 tRNA
  charging, 20 translation cost payment and 2 replication initiations. M0 and the
  reduced model have none. 361 have constant rule nodes.
- None of the 363 closes syntactically: `MeanField` rejects every one. Not all of them
  fail to close exactly.
  - 177 DnaA reactions are site × site within one origin's cascade. Each origin's
    cascade holds a single 0/1 token, so every such product equals one of its factors.
    They close exactly.
  - 1 (reaction 1720) can never fire.
  - The 2 initiations close exactly in the first round only.
  - The rest do not close exactly: 99 DnaA bindings to a 0/1 site, 2 DnaA bindings to
    a multi-copy site, 62 charging and 20 cost payment. That is 183, plus the
    initiations from the second round on. For the 82 charging and payment reactions
    this was not enumerated.
- The 493 replication reactions are 2 of those initiations and 491 first-order fork
  steps. The steps' form closes; their constants are rebuilt from the dNTP pools and a
  rounded, capped volume.
- All 1,434 count-dependent rate rules reach a `max(1, ·)` floor. 29 transcription
  promoters read their own protein through `min(·, 180)`. 493 replication rules reach
  `round`, and 948 rules reach the volume cap. 18 ODE reactions reach a `min` AND rule.
  None of the 29 protein-count genes is in the reduced model.
- Gene copy number multiplies transcription and nothing else. It keeps the mean field
  affine and breaks the transcript factor: each gene's birth rate steps at a latent
  time, and all 486 genes share the same two fork paths.
- In the dependency graph the script builds, every kind reaches an observable.
  - In 1 hop: tRNA charging and cost payment (bulk metabolites), transcription
    (transcripts) and the ODE AND rules (bulk metabolites).
  - In 2 hops: the replication reactions, through a gene and transcription, or
    through the `*_DNArep` counters the hook pays from the ODE pools.
  - In 3 hops: the DnaA bindings.
  - Doubling time is 3 to 6 hops from every kind.
- The proposed M13 model, M0 with one origin, contains one second-order reaction
  (`D + O → F₁`) and gene copy number for both M0 genes. It adds 6 species and 4
  reactions to M0 and changes 2. Its R2 stands for the 99 live DnaA bindings to a 0/1
  site, and for the 2 bindings to a multi-copy site. It does not stand for the 177
  reactions that close exactly, and it does not stand for charging or cost payment.

## Sources

The exported reference model in `syn3a-ensemble` (read only), at HEAD
`2075809cbb825140cb31a7bec8bc3529770c8e13`. The export itself was made at
syn3a-ensemble `514fa79` from Minimal_Cell `db048ac` (`export_meta.json`).

| Table | sha256 |
|---|---|
| `data/model/cme/reactions.csv` | `63fbeca180fe8bdebec646df0ddbffb0149958fc2ff9fc3e6ec8803cd24c6aaf` |
| `data/model/cme/rate_rules.csv` | `a95619a05e47b32e7e97e8f65128e31a25403c386e80b5ccfb298d905aa3af0e` |
| `data/model/cme/species.csv` | `c1acb80878e8a500a86e0e1d9e195a491b46a01c337788c0a6554df19cbe6ec9` |
| `data/model/cme/initial_counts.csv` | `1149da1f0ac53d53d1eb8d63c302e4d638e5e15666bddbb66cfd1f08276c67b5` |
| `data/model/nodes.csv` | `1bd4b35e1c33f3547fd5683213b5b4c32d7211dd9046d8a4b7ee628ee1f39b89` |
| `data/model/ode/species.csv` | `77d003243ec52770287ba76658f4deeb00cb3e8ad2e51172d0d3a7517b02bd55` |
| `data/model/ode/reactions.csv` | `4da4941402a0e176ad6adad2c35b6ad0f67f992c023a306ea43f5edb3b8b9a27` |
| `data/model/ode/stoichiometry.csv` | `be4293c536cac2ef1a73cabf0144e5a64b696a7bff5f07ab15c871c8a710bf2f` |
| `data/model/export_meta.json` | `e3c21a138aa318eb39f5e21a2e1f5e31476d2d10581ecec8e552b1128ce66251` |
| `data/model/inventory.md` | `2cbf4eae3b624a6b310c8eb8e67d30ea01ebd8b6a99c1a7fb12e13187399cca6` |
| `code/SpeedyCell/src/Cell/writeback.jl` | `ae6d3124ca31fddaad5b0a9c08f150ab601b1ae52727fce2ba0c4ec8802e0a41` |
| `code/SpeedyCell/src/Cell/hook.jl` | `1cec291491c64202495cd6820d439657e333f2228001c8cefad0b285d9788720` |

The four CME tables (`reactions`, `rate_rules`, `species`, `initial_counts`) are
byte-identical to SpeedyCell.jl's `data/model/cme/` (`/fred/oz022/tkimpson/SpeedyCell.jl`,
HEAD `0823b53`; the script compares them). SpeedyCell.jl's own `inventory.md` has
sha256 `23d6b162…`. It differs from syn3a-ensemble's in wording only: it drops the phase
and spec references, and the numbers and rules are the same. That was checked by hand
with `diff`. This note cites syn3a-ensemble's copy.

**Access.** syn3a-ensemble and SpeedyCell.jl are private repositories. On the cluster
the script finds them at its default paths. Elsewhere, check out syn3a-ensemble at
`2075809` and pass `--model <checkout>/data/model` and
`--speedycell <checkout>/code/SpeedyCell/src`. `--speedycell-jl` is optional: without
it the comparison is skipped. The sha256 table pins every input.

**Reproducibility.** At syn3a-ensemble `2075809` the output reproduces byte for byte
across runs and hash seeds. Its line 2 is the checkout's HEAD, which is provenance. It
changes after a pull even if the tables do not, so compare the sha256 lines, not the
whole file.

## Terms and assumptions

**Two meanings of "closes".** This note uses two terms and keeps them apart.
- **Closes syntactically:** `MeanField`'s check accepts the reaction as the export
  writes it. `MeanField` (`src/mean_field.jl`) builds mean equations only for a network
  whose every propensity is affine in the state and whose every state change is
  constant. It throws on anything else ("A propensity … is not affine in the state",
  `src/mean_field.jl:105`). Every type-2 reaction fails this check.
- **Closes exactly:** the mean equations are exact given the rate constants. That holds
  when the propensity is an affine function of the state over the states the cell can
  reach. A product of two counts can collapse there. For example, if a conserved 0/1
  token means y = 1 whenever x = 1, then x·y = x. Such a reaction closes exactly
  although `MeanField` rejects it as written. A re-encoding as first order would
  remove it.

A rate constant that is rebuilt from the counts at a minute boundary is outside both
checks. The mean-field route evaluates it at the mean pools
(`src/mean_field_route.jl:15`), so its error is a Jensen gap in what reads the mean,
not a failure of the mean equations.

**The transcript factor is exact** when each gene's transcript count is a birth–death
chain whose birth rate is known and piecewise constant: constant over each window
except the switch 1 s before its end (`src/transcript_factor.jl`, and the cSMC's
bridges, `src/csmc.jl:140` and `src/transcript_bridge.jl`). These are two different
failures, and this note keeps them apart.

**Observables** are the paper's three.
- *Doubling time* is taken as `t_SA2`: the first minute at which `CellSA` reaches twice
  its minute-0 value (`code/cell_cycle/events.py`, the reference paper's cycle end). It
  is itself a threshold crossing of a stochastic path. `CellSA` is computed in code from
  9 lipid counts and 93 membrane-protein ids, 103 species in all (`writeback.jl`,
  `lipid_sa!`).
- *Transcripts* are each cell's `M_RNA_*` counts.
- *Bulk metabolites* are taken as population means of the 146 ODE species. This is an
  assumption: the full-cell panel is not chosen yet.

**Candidate blocks** are assumed, since M12's variance decomposition picks the block.
- The *gene-expression block* is the constants of transcription (486 reactions),
  translation (455) and mRNA degradation (455).
- The *metabolic block* is the kinetic constants in the 175 ODE rate laws.
- The constants of replication, DnaA binding and tRNA charging are in neither. They
  would join a block only if M12 picks them.

**Between a block and an observable.** The cell is one connected network with
feedback, so a yes/no "is it on a path" says nothing. The script reports a hop count
instead. It builds a dependency graph on the species. An edge a → b means that one CME
reaction, ODE reaction or hook rule reads a and changes b. "Reads" covers reactants,
the counts a rate rule reads and the species an ODE rate law reads. The hook rules are
in `writeback.jl`: the cost counters paid from their pools (lines 120, 128 and 139),
the surface area from lipids and membrane proteins, and the volume from the surface
area.

A kind's hop count to an observable is 1 for the kind's own reaction, plus the shortest
path from what it changes. The script prints the route. The graph counts every read,
however weak, so a short route says that a coupling exists, not that it is strong.

The upstream side, from a block to the kind, is short for every kind, and the note
states it in prose. Every kind reads a block's output directly or one hop away: the
ODE pools, the proteins translation makes, or `M_DnaA_c`.

## 1. The 363 second-order reactions

**Order.** `cme/reactions.csv` gives each reaction an LM type, 0 to 3, with
propensities k, k·x, k·x·y and k·x(x−1)/2 (`data/model/README.md`). The table has
2,002 reactions of type 1 and 363 of type 2. It has none of type 0 or 3. Every type-2
reaction has two distinct reactants. So "second order" here means type 2, a product of
two different counts, and the count is 363, as the issue says.

**Families.** Each reaction is classified by an explicit rule on its reactant and
product names (`second_order_family` in the script). Nothing is left over.

| Family | n | Subfamily | n |
|---|---:|---|---:|
| DnaA at the origins | 279 | DnaA leaves the filament (unbinding cascade, k = 600) | 90 |
| | | DnaA extends the ssDNA filament | 87 |
| | | DnaA leaves the filament tip (k = 0.55) | 87 |
| | | DnaA binds a low-affinity site at oriC | 6 |
| | | DnaA binds a high-affinity site off oriC | 3 |
| | | DnaA binds the high-affinity site at oriC | 3 |
| | | DnaA binds the last ssDNA site (makes `Initiator_C` and `Initiator_CC`) | 3 |
| tRNA charging | 62 | synthetase–ATP–amino acid binds its tRNA | 20 |
| | | synthetase binds ATP | 19 |
| | | synthetase–ATP binds its amino acid | 19 |
| | | Gln-tRNA(Gln) amidotransferase (GatCAB, `JCVISYN3A_0689`) steps | 4 |
| Translation cost payment | 20 | residue counter + charged tRNA → uncharged tRNA | 20 |
| Replication initiation | 2 | `Initiator_C + chromosome_C`, `Initiator_CC + chromosome_CC` | 2 |
| **Total** | **363** | | |

The counts sum to 363. With the first-order families they reproduce `inventory.md`'s
table: its "Other" row (383) is the 279 DnaA, 62 charging and 20 cost-payment reactions
plus 22 first-order charging release steps, and the whole table sums to 2,365.

The DnaA family has three copies of the same cascade: the starting origin and the two
daughter origins (`High_Affinity_Site_oriC`, `…_oriC2`, `…_oriC3`). The daughter oriC
sites start at zero and are made by replication reaction 2364 (section 3), along with 16
copies of the off-oriC site `High_Affinity_Site2`.

Reaction 2364 makes no `High_Affinity_Site3`, and nothing else makes it either. So
reaction 1720 (`M_DnaA_c + High_Affinity_Site3`) can never fire. The script finds it as
the only reaction with a reactant that starts at 0 and that nothing makes. This is
reported as found.

**The 363 and the 493 overlap.** The two replication initiations are second order and
are also two of the 493 replication reactions. The full cell has 363 + 491 = 854
reactions in these two sets, not 856.

**Rate constants.** 361 of the 363 point at a `const` node, so their rule does not
change after minute 1. One of them changes once, at minute 2. Reaction 1494
(`M_glutrnagln_c` + GatCAB) has a first-minute rate of 0.01, and its rule from minute 2
on is the constant 0.1. That is a tenfold step, and the script finds it as the only
reaction whose first-minute rate differs from its constant rule. The two count-dependent
exceptions are the replication initiations (reactions 1873 and 2098), whose constants
are rebuilt each minute from the dNTP pools and `CellSA`. The constant second-order
rates do not scale with the cell volume, although the cell doubles its volume.

**Copy-number scales**, at procid 1's initial counts (`cme/species.csv`):

| Family | Smaller reactant | Larger reactant |
|---|---|---|
| DnaA at the origins | 0 to 16 (sites, filament states) | 0 to 148 (`M_DnaA_c` is 148) |
| tRNA charging | 1 to 736 (synthetases and their intermediates) | 40 to 73,715 (tRNA, amino acids, ATP) |
| Translation cost payment | 0 (the residue counters) | 160 to 480 (charged tRNA) |
| Replication initiation | 0 (`Initiator_*`) | 1 (`chromosome_*`) |

**A finding in the table.** All 20 cost-payment reactions credit `FMET_cost_paid`. None
credits its own `<AA>_cost_paid` counter. The table has 21 `*_cost_paid` species, and
no reaction writes the other 20. This is reported as found. The script checks it; this
note does not say whether it matters, since nothing in the tables reads a
`*_cost_paid` count.

**Closure.** None of the 363 closes syntactically. Each propensity is a product of two
counts as written, so `MeanField` rejects it. Whether each closes exactly, the script
decides for the DnaA cascade by enumerating every state one origin can reach. The
origin starts with its oriC site at 1 and the rest of its cascade at 0. DnaA is treated
as always available, and the off-oriC multi-copy sites are left out of the origin's
state.
- Each origin reaches 64 states.
- In every one of them, each of the origin's 59 site × site reactions has a product
  equal to one of its two factors. So these reactions close exactly: 3 × 59 = 177. The
  cascade carries a single 0/1 token. One final unbinding per origin removes it
  (reaction 1812 for the starting origin).
- Each origin's 33 DnaA bindings are to a site that is 0 or 1 in every reachable state:
  3 × 33 = 99. Their propensity is k·D·s with D fluctuating, so they do not close
  exactly.
- For the daughter origins this holds on the condition that reaction 2364 fires once,
  so that each daughter oriC is made once.

| DnaA reactions (279) | n | Closes syntactically? | Closes exactly? |
|---|---:|---|---|
| site × site within an origin's cascade (filament step-back 87, unbinding 90) | 177 | no | **yes**, by the 0/1 token |
| DnaA × a 0/1 site (oriC 3, low-affinity 6, ssDNA 90) | 99 | no | no |
| DnaA × a multi-copy off-oriC site (1594 at 16 copies; 1657 at 16 after reaction 2364) | 2 | no | no |
| DnaA × `High_Affinity_Site3` (1720) | 1 | no | trivially: it never fires |

For the other 84 second-order reactions:
- **The 2 initiations** (`Initiator × chromosome`) close exactly in the first round. The
  script enumerates a hand-written abstraction of the tables: one origin makes one
  `Initiator_C` (reaction 1656), and one `chromosome_C` is consumed (reaction 1873). It
  has 3 reachable states, (O, I, C) = (1, 0, 1), (0, 1, 1), (0, 0, 0), and I·C = I in
  each. From the second round on they do not close exactly. After reaction 2364,
  `chromosome_C` is 2 and two daughter origins each make an `Initiator_C`. That
  abstraction has 8 states, and the script finds I·C is not affine over them.
- **The 62 charging and 20 cost-payment reactions** have multi-copy reactants on both
  sides. Charging is an enzyme state × ATP, an amino acid or a tRNA pool; payment is a
  residue counter × a charged tRNA pool. No conservation law in their family collapses
  the product. They were not enumerated, and they are counted here as not closing
  exactly.

In total, 177 close exactly and 1 never fires. The 2 initiations close exactly in the
first round only. 183 do not close exactly: 99 + 2 + 62 + 20.

M0 and the reduced model have no second-order reaction. The reduced model lumps tRNA
charging into one deterministic ODE reaction
(`src/organisms/coreA/trna_charging.jl:1-19`), and it pays translation's residue cost
through deferred counters drained at the 1 s handshake, not through a reaction.

**Between a block and an observable.** The table gives hops to the nearest observable
of each type, from the script ("Hops to the observables" in its output). The arrows name
the shortest route found.

| Kind | Upstream: reads | Transcripts | Bulk metabolites | Doubling time |
|---|---|---:|---:|---:|
| tRNA charging (62) | ATP and the amino acids (ODE species), synthetases (translation's proteins), uncharged tRNA. Its type-2 steps do not touch the charged pools; the 22 first-order releases make them | 2: ATP → a transcription rule | 1: it consumes ATP and the amino acids | 4 |
| Translation cost payment (20) | the residue counters translation writes, the charged tRNA pools | 3: charged tRNA → translation → an own-protein promoter | 1: one of the 20 consumes `M_mettrna_c`, an ODE species | 4 |
| DnaA × 0/1 site (99) | `M_DnaA_c`, made by translation reaction 9 | 3: → `Initiator_C` → 1873 → `JCVISYN3A_0001_gene` → transcription | 3: → `Initiator_C` → 1873 → `ATP_DNArep` → paid from ATP at the hook | 5 |
| DnaA × multi-copy site (2) | as above | 4 | 4 | 6 |
| Cascade site × site (177) | the cascade's own sites | 4: they release DnaA | 4 | 6 |
| Replication initiation (2) | `Initiator_*`, `chromosome_*`; the dNTP pools and `CellSA` through the rule | 2: → a gene → transcription | 2: → the `*_DNArep` counters → paid from the ODE pools (`writeback.jl:120, 128`) | 4 |

No second-order family writes a `CellSA` input directly. Doubling time is 4 to 6 hops
away. The shortest routes run through an ODE pool, a lipid reaction and the hook's
surface-area rule.

## 2. Rate rules with `min`, `max`, `round` or a threshold

**Where the rules live.** A CME reaction's rate constant from minute 2 on is the node
`rate_rules.csv` names, evaluated on the counts at each minute boundary
(`data/model/README.md`). An ODE reaction's rate refers to count-dependent nodes as
`N__<id>` (enzyme concentrations), and each ODE species' starting value at each 1 s
step is node `init_node`. All are nodes of one DAG, `nodes.csv`. The enzyme rules
`inventory.md` describes (AND rules by `min` over subunits, OR rules by sum) are among
those `N__` nodes.

**Counting unit.** `export_meta.json`'s op histogram counts nodes over the whole DAG.
The script splits it by what reaches each node. A rule "uses" an op if the op is
reachable from its root.

| Op | Nodes in the DAG | Reached from CME rate rules | From ODE rates | From ODE `init_node` |
|---|---:|---:|---:|---:|
| `min` | 38 | 29 | 9 | 0 |
| `max` | 29 | 29 | 0 | 0 |
| `round` | 1 | 1 | 0 | 0 |
| `select_gt` (the only threshold op) | 2 | 2 | 1 | 1 |

`floor`, `trunc`, `floordiv`, `mod` and `abs` do not occur. The 70 nodes fall into five
kinds. The "at t = 0" column evaluates the DAG at procid 1's initial counts with the
README's Python semantics. It says which branch is taken at that one state only.

| Kind | Nodes | Reached by | At t = 0 |
|---|---:|---|---|
| `max(1, pool)`: a floor on a pool before division. The pools are the 4 NTPs, the 4 dNTPs and the 21 charged tRNAs | 29 | CME rules | no floor active |
| `min(c · P/180, 180)`: a promoter strength from a live protein count P, capped at 180 | 29 | 29 transcription rules | no ceiling binding |
| `min` over subunit terms: a complex's AND rule (2 to 4 arguments) | 9 | 18 ODE reactions | — |
| `select_gt(V, 6.7e-17, …)`: the volume cap | 2 | CME rules, ODE rates, ODE init | below the cap |
| `round(V · 10¹⁹)`: `CellV` recomputed from `CellSA` | 1 | replication rules | `CellV` = 335 |

**Every count-dependent rate rule, classified.** 1,434 CME rules depend on the
counts, as `export_meta.json` says; the other 931 are constants. They are exactly four
families, and every one of the 1,434 reaches at least one of the four ops:

| Family | Rules | Reads | `min` | `max` | `round` | `select_gt` |
|---|---:|---|---:|---:|---:|---:|
| Transcription | 486 | 4 NTP pools and the `CellV` count; 29 also read their gene's own protein | 29 | 486 | 0 | 0 |
| Translation | 455 | all 21 charged tRNA pools (each rule reads all 21), volume from `CellSA` | 0 | 455 | 0 | 455 |
| Replication fork steps | 491 | 4 dNTP pools, `CellSA` | 0 | 491 | 491 | 491 |
| Replication initiation | 2 | 4 dNTP pools, `CellSA` | 0 | 2 | 2 | 2 |
| **Total** | **1,434** | 60 count species in all | 29 | 1,434 | 493 | 948 |

The 931 constant rules are mRNA degradation (455), translocation (93), the 361
constant second-order reactions of section 1 and the 22 charging release steps.

Transcription reaches neither `round` nor `select_gt` in the DAG because it reads the
`CellV` count directly. That count is rounded and capped in code at every hook
(below). So all 486 transcription rules read a rounded, capped volume too.

On the ODE side, 158 of the 175 reactions read a count-dependent node. All 158 reach
the volume cap. 18 reach a `min` AND rule. None reaches `max` or `round`.

**The 29 protein-count promoters are new.** Each of the 29 reads the count of its own
gene's protein, so each is a feedback from a gene's protein onto its own transcription.
The genes are `JCVISYN3A_` 0027 0077 0094 0154 0214 0230 0235 0259 0263 0379 0398 0401
0503 0525 0593 0622 0623 0643 0648 0649 0795 0851 0852 0873 0913 0918 0930 0931 0932.
None is one of the reduced model's 17 genes.

**Rules that are code, not data.** These are in SpeedyCell's hook, not in the tables.
The script lists their lines but does not count them as rate rules.

| Rule | Where | What it does |
|---|---|---|
| Volume cap and `CellV` | `Cell/writeback.jl:98-99` | `V > 6.70e-17` returns the cap and `CellV` = 670; else `CellV = round(V · 10¹⁹)` |
| ODE → counts | `Cell/writeback.jl:114` | each of the 146 ODE species is rounded half to even to a count, every 1 s |
| Cost counters | `Cell/writeback.jl:120, 128, 139` | 22 counters (5 ATP/GTP hydrolysis, 12 NTP, 4 dNTP, 1 `FMET_cost`) are paid with `min(counter, pool)`, the remainder carried; the 4 NMP counters are credited whole |
| Surface area | `Cell/writeback.jl:154, 159` | `CellSA_Lip` is rounded, `CellSA_Prot` truncated, `CellSA` their sum |
| ODE span | `Cell/hook.jl:26` | the span is `max(round(Δt, 2), 0.1)`, stepped in 0.1 s; a function of time, not of counts |

**What M0 and the reduced model already have.**

| Rule | Full cell | M0 and the reduced model |
|---|---|---|
| `max(1, pool)` floor | 29 pools, in all 1,434 rules | translation floors its pools (`src/organisms/coreA/translation.jl:152`); transcription has no floor (`transcription.jl:522-540`) |
| Promoter ceiling `min(·, 180)` | 29 rules read a live protein count | the ceiling is on a parameter, and it never binds (`transcription.jl:61-65`, `turnover_headroom`) |
| Volume cap | reached mid-cycle: the reference paper reports the volume doubling at 56 to 72 min, mean 64 (quoted in `code/cell_cycle/events.py`, `t_Vcap`) | applied (`src/handshake.jl:1732`), but the reduced cell grows to about 1.07× and never nears it (`handshake.jl:1721-1722`) |
| Rounded `CellV` | read by all 486 transcription and 493 replication rules | absent: the reduced model has no `CellV` count |
| Complex AND rules | 9 nodes, 18 ODE reactions | absent: each enzyme is one protein count |
| ODE → counts rounding | round half to even, 146 species | `:fractional_carry` by default (`handshake.jl:194-226`) |
| Counter clamp | 22 counters, `min(counter, pool)` | the hybrid clamps and carries the deficit (`handshake.jl:1878-1879`); the mean-field route does not clamp (`mean_field_route.jl:12-14`) |

**Closure.** None of these rules makes a reaction fail to close syntactically, since a
rebuilt rate constant is outside `MeanField`'s check. Each one adds a Jensen gap where the route evaluates the
rule at the mean counts. A floor, ceiling or cap adds no gap of its own while it is
far from its kink, but the rule beneath it is still nonlinear in the counts. `round`
and `CellV` step at every count, and the volume cap is reached mid-cycle. This note
gives structure and scale, not the size of any gap.

**Between a block and an observable.** These are hops from the script.
- **Transcription and translation rules.** These are the gene-expression block's own
  rate constants, and they read pools the metabolic block writes. Over all 1,434
  count-dependent rules, the nearest transcript and the nearest bulk metabolite are
  1 hop away. Doubling time is 3 hops away: translation makes a cytosolic membrane
  protein, translocation moves it to the membrane, and the hook's surface-area rule
  reads it.
- **The 29 own-protein promoters.** Transcripts are 1 hop away, bulk metabolites 2
  (through the transcription counters paid from the NTP pools), doubling time 4.
- **The replication rules.** They read the dNTP pools and step each gene's copy number
  (2 hops to transcripts and to bulk metabolites, section 3).
- **The 18 ODE reactions with AND rules.** These are the metabolic block's rate laws.
  They change bulk metabolites directly (1 hop), and doubling time is 3 hops away.

## 3. The 493 replication reactions and gene copy number

### How replication is represented

The chromosome is two chains, one per fork. Each of the 493 `*_gene` species starts at
1. Each of the 493 replication reactions adds +1 to exactly one gene, every gene is
incremented by exactly one of them, and no reaction removes a gene. A fork's position
is a chain of `*_inter` species: step i turns `<gene_i>_inter` into `<gene_i+1>_inter`
and copies gene i+1.

| Fork | Starts with | First gene | First-order steps | Ends at |
|---|---|---|---:|---|
| C | reaction 1873, `Initiator_C + chromosome_C` (second order) | `JCVISYN3A_0001` (*dnaA*) | 224 | reaction 2097 |
| CC | reaction 2098, `Initiator_CC + chromosome_CC` (second order) | `JCVISYN3A_0910` | 267 | reaction 2365 |

So the 493 are 2 second-order initiations and 491 first-order fork steps. Reaction
2364 (fork CC, copying `JCVISYN3A_0421`) also makes the daughter origins:
`High_Affinity_Site_oriC2` and `…_oriC3` (+1 each), `High_Affinity_Site2` (+16), and
`chromosome_C` and `chromosome_CC` (+2 each). The daughter origins then run their own
DnaA cascades (section 1), which make `Initiator_C` and `Initiator_CC` again. A second
round reuses the same `*_inter` chain, so a chain species can then hold more than one
fork.

The initiators come from DnaA. `M_DnaA_c` is made by translation reaction 9, from
`M_RNA_JCVISYN3A_0001_c`. It is consumed and released by the 279 DnaA reactions. The
last ssDNA binding at an origin (reactions 1656, 1719 and 1782) makes one
`Initiator_C` and one `Initiator_CC`. Initial counts: `M_DnaA_c` 148,
`High_Affinity_Site_oriC` 1, `chromosome_C` 1, `chromosome_CC` 1.

All 493 rate rules depend on the counts. Each reads the four dNTP pools, floored by
`max(1, ·)`, and the volume recomputed from `CellSA`, through `round` and the
`select_gt` cap. They are rebuilt each minute, like every other count-dependent rule.

### Does the mean field close on them?

- **The 491 fork steps:** the form closes. Each is type 1 with a constant state change,
  so its mean equation is exact given its rate constant. The rate constant is not
  given: it is rebuilt from the dNTP pools and the rounded volume, so the route
  evaluates it at the means (section 2).
- **The 2 initiations:** they do not close syntactically. In the first round, both
  reactants are 0 or 1 and I·C = I in every reachable state, so they close exactly
  (section 1). From the second round on, `chromosome_*` is 2 and there are two
  initiator sources, so they do not.
- **What the means describe.** In the first round, a cell's fork is at one site. Its
  mean is an occupation probability spread over 225 or 268 sites. Each gene's copy
  number steps from 1 to 2 at a random time, and again after a daughter-origin
  initiation; its mean rises smoothly. The mean field gets `E[G_g(t)]` right given
  exact rate constants and the right initiation-time distribution. That distribution
  comes from the DnaA bindings, which do not close exactly. The mean field cannot give
  the step.

### How gene copy number enters transcription

- All 486 transcription reactions are type 1 with one `*_gene` reactant, which the
  reaction leaves unchanged. So the propensity is `k_g(t) · G_g(t)`. The 486 make 455
  mRNAs, 29 tRNAs and 2 rRNA operons (the 2 operon reactions make the 5S, 23S and 16S
  rRNAs each).
- 7 genes have no transcription reaction: `JCVISYN3A_` 0049, 0068, 0069, 0158, 0356,
  0533, 0534. Genes 0067 and 0532 each have one rRNA reaction that makes all three
  rRNAs (reactions 1592 and 1593); 0068, 0069, 0533 and 0534 sit beside them in the
  species table and have none of their own.
- `G_g` enters nothing else. No other reaction has a gene as a reactant. No rate rule
  reads a gene count. The rate constant `k_g` reads the four NTP pools and `CellV`, and
  for 29 genes their own protein count (section 2).
- So `G_g` multiplies the propensity and nothing more. The mean equation
  `dE[m_g]/dt = E[k_g G_g] − μ_g E[m_g]` stays affine in the state given `k_g`, and
  `MeanField`'s check passes. Gene copy number does not break the mean field. It breaks
  the transcript factor.

### What a transcript-only dataset leaves latent

The transcript factor treats each gene's transcripts as a birth–death chain with a
known birth rate that changes only at the 60 s rebuild (`src/transcript_factor.jl:1-22`;
the cSMC reads the birth rate at the window's start state, `src/csmc.jl:140-147`). In
the full cell, the birth rate is `k_g · G_g(t)`, and a transcript-only dataset does not
observe:

- **each gene's copy-number path** `G_g(t)`: when it steps from 1 to 2, and, after a
  daughter-origin initiation, again. The step falls anywhere in a window, not at the
  rebuild.
- **what sets those times:** the two fork positions (the `*_inter` species), the two
  initiation times (reactions 1873 and 2098), and before them the DnaA cascade's state
  (279 reactions over the site species) and `M_DnaA_c`. DnaA's own transcript,
  `M_RNA_JCVISYN3A_0001_c`, is observed; its protein is not.
- **the rate constants' inputs in that cell:** the four NTP pools and `CellV` for every
  gene, and the gene's own protein count for the 29 protein-count promoters.

Every gene's copy number is a function of the same two fork paths. So given the
parameters, the 486 genes' transcript factors are no longer independent: they share
the fork paths. In the reduced model and M0, genes are a fixed quantity, not a state
(`src/organisms/coreA/transcription.jl:32-37`), and this never arose.

M0's two genes are both on fork CC. ptsG (`JCVISYN3A_0779`) is copied at step 68
(reaction 2166). GAPD (`JCVISYN3A_0607`) is copied at step 172 (reaction 2270).

### Between a block and an observable

Replication reads pools the metabolic block writes (the dNTPs). From the script, its
491 fork steps are 2 hops from an observable by two routes:
- through a gene and then transcription, to transcripts;
- through the five `*_DNArep` counters every replication reaction writes. The hook pays
  them from `M_atp_c` and the dNTP pools (`writeback.jl:120, 128`), which are bulk
  metabolites.

A gene copy number itself is 1 hop from a transcript and 2 from a bulk metabolite.
Doubling time is 4 hops from both. Replication does not write a `CellSA` input itself.

## 4. The inventory, and the smallest model containing each kind

"Syntactic" and "exact" are the two meanings of "closes" defined in Terms. "Factor"
asks whether the item breaks the exact transcript factor's assumption. "Hops" gives the
script's hop count to the nearest transcript, bulk metabolite and doubling time, in that
order, and names the shortest route.

| # | Kind | In the full cell | Syntactic? | Exact? | Factor? | Hops (tx / bulk / dt) and route | In M0 / reduced | Smallest model containing it |
|---|---|---|---|---|---|---|---|---|
| 1a | Second order: DnaA × a 0/1 site | 99 bindings | no | no | through G (row 5) | 3 / 3 / 5: → `Initiator_*` → initiation → gene or `ATP_DNArep` | no / no | **the M13 model below** (R2) |
| 1b | Second order: DnaA × a multi-copy site | 2 bindings (1594, 1657); 1720 never fires | no | no | through `M_DnaA_c` | 4 / 4 / 6 | no / no | the M13 model with a second, multi-copy site species that also binds `D` |
| 1c | Second order: site × site in one origin's cascade | 177 | no | **yes**, by the 0/1 token | through `M_DnaA_c` | 4 / 4 / 6 | no / no | none needed: re-encoded as first order in the origin's state, it has no gap |
| 2 | Second order: 0/1 × 0/1, `Initiator × chromosome` | 2 initiations | no | first round yes; second round no | through G | 2 / 2 / 4: → gene → transcription; → `*_DNArep` → ODE pools | no / no | see below the table |
| 3 | Second order: enzyme state × metabolite or tRNA | 62 tRNA charging | no | no (not enumerated) | no | 2 / 1 / 4: consumes ATP and amino acids | no / no: charging is one lumped ODE reaction | M0 with its lumped charging reaction replaced by one synthetase chain as jump reactions (3 second-order bindings, 1 release) |
| 4 | Second order: residue counter × charged tRNA | 20 cost payments | no | no (not enumerated) | no | 3 / 1 / 4: one consumes `M_mettrna_c` | no / no: deferred counters | M0 with one residue counter made a species, paid by one second-order reaction |
| 5 | Gene copy number, a state multiplying transcription | 493 genes, 486 transcription reactions | yes | yes | **yes** | 1 / 2 / 4, from the gene species | no / no: genes fixed | **the M13 model below** (`G_0779`, `G_0607`) |
| 6 | First-order fork steps with count-dependent constants | 491 | yes | form yes; constants at the mean | through G | 2 / 2 / 4: → gene; → `*_DNArep` → ODE pools | no / no | the M13 model with its fork-step constants rebuilt from a dNTP pool; `src/organisms/coreA/` has no dNTP species, so one must be added |
| 7 | `max(1, pool)` floor | 29 nodes, in all 1,434 rules | rebuild | rebuild | no | 1 / 1 / 3, as the rules it sits in | translation only | M0 (already) |
| 8 | Promoter from a live own-protein count, `min(c·P/180, 180)` | 29 transcription rules | rebuild | rebuild | as now, the rebuilt birth rate reads latent state; new is that the state is a stochastic protein count, not an ODE pool | 1 / 2 / 4 | no / no | M0 with one gene's promoter read from its own protein at the rebuild |
| 9 | Rounded `CellV`; the volume cap reached | 1 + 2 nodes; 486 + 493 rules read `CellV` or its `round`, 948 the cap | rebuild; steps | rebuild; steps | as row 8, through the rebuilt k | 1 / 1 / 3, as the rules it sits in | no `CellV`; cap present but never reached | a model that grows to the cap: the full cell, or the reduced model over a full cycle with growth that reaches it |
| 10 | Complex AND rule, `min` over subunits | 9 nodes, 18 ODE reactions | in the ODE, at the mean counts | in the ODE, at the mean counts | no | 2 / 1 / 3 | no / no | M0 with one enzyme made a two-subunit complex |
| 11 | Hook rules: counter clamp, ODE → count rounding, `CellSA` round and trunc | code (section 2) | route-dependent | route-dependent | no | they write the pools and `CellSA` themselves (by construction, not from the script) | clamp yes; rounding by another policy | M0 (already), up to the rounding policy |
| 12 | Doubling time as a threshold crossing of `CellSA` | `t_SA2` | an observable, not a reaction | — | — | it is the observable | no (600 s) / growth, no cycle end | the full cell, or a reduced-model cycle |

**Row 2's smallest model.** Splitting R2 into `D + O → I` and then `I + C → F₁` adds two
species, `I` and `C`, and one reaction. It does not give a 0/1 × 0/1 non-closure. With
O = 1, I = 0 and C = 1 at the start, the model reaches three states, and I·C = I in all
of them (the script's "first round" line). The full cell's first-round initiations
collapse the same way, so the first-round model is degenerate.

The smallest non-degenerate model needs `C` at 2 and two sources of `I`, as after
reaction 2364:
- species `O_a`, `O_b`, `I` and `C`, with C = 2 at the start;
- reactions `D + O_a → I`, `D + O_b → I` and `I + C → F`.

The script's "second round" abstraction of this has 8 states, and I·C is not affine
over them. In the full cell this kind appears only from the second round. The
reference paper puts daughter initiations at 58 to 110 min (quoted in `events.py`),
long after M0's 600 s. It is not proposed for M13.

"Rebuild" means the rule changes a rebuilt rate constant: the mean equations still
close, and the route evaluates the rule at the mean counts (section 2).

## 5. The proposed model for M13: M0 with one origin

The smallest model with one second-order reaction and gene copy number, buildable on
M0's machinery. It is M0 (`src/organisms/coreA/m0.jl`: GAPD and ptsG, 600 s, the four
ODE modules, `HeldEnzymeCounts`) with one new jump module. It mirrors the full cell's
causal chain: DnaA meets the origin, a fork starts, and the fork copies ptsG and then
GAPD, in the full cell's order on fork CC (steps 68 and 172).

**Species added** (all counts):

| Species | Initial | Stands for |
|---|---|---|
| `D` | D₀ | `M_DnaA_c` |
| `O` | 1 | the free origin, `High_Affinity_Site_oriC` |
| `F₁` | 0 | a fork that has not reached ptsG (`*_inter` before step 68) |
| `F₂` | 0 | a fork between ptsG and GAPD (`*_inter` steps 68 to 171) |
| `G_0779` | 1 | `JCVISYN3A_0779_gene`, promoted from a fixed quantity to a state |
| `G_0607` | 1 | `JCVISYN3A_0607_gene`, likewise |

**Reactions added:**

| | Reaction | Propensity | Order | Stands for |
|---|---|---|---|---|
| R1 | ∅ → `D` | k_D | 0 | DnaA's expression (translation reaction 9), lumped |
| R2 | `D` + `O` → `F₁` | k_b · D · O | **2** | the DnaA cascade at one origin and the initiation, lumped (origin 1's cascade is reactions 1595 to 1656, then 2098) |
| R3 | `F₁` → `F₂` + `G_0779` | k_f1 · F₁ | 1 | fork CC's steps 1 to 68 |
| R4 | `F₂` → `G_0607` | k_f2 · F₂ | 1 | fork CC's steps 69 to 172 |

**Reactions changed:** M0's two transcription reactions take the propensity
`k_g · G_g` in place of `k_g`, for g = 0779 and 0607. Their rate constants are rebuilt
from the pools exactly as now.

Everything else in M0 is unchanged: decay, translation, the ODE block, the handshake,
the data and M0's five targets (`M0_TARGETS`). k_D, k_b, k_f1 and k_f2 are fixed at
declared values. Making k_b a sixth target is the first extension.

**Why R1 is needed.** Without a source of `D`, `D` is `D₀ − (1 − O)`, a function of
`O`. Since O² = O, `D · O = D₀ · O` in every state. So E[D·O] is affine in E[O], and the
mean equations close exactly. The check would then test nothing. E[D·O] − E[D]·E[O] is
Var(O), which is not zero, but that does not matter: the mean equation needs E[D·O],
and that is D₀·E[O]. With R1, `D` fluctuates apart from `O`. Then E[D·O] needs a second
moment, and R2 does not close exactly. The script enumerates both cases: without R1
(D₀ = 3) D·O is affine over the reachable states; with three R1 births it is not.

**Lumping the DnaA walk changes the kind of dependence on `D`.** R2 replaces the
origin's ssDNA walk with one irreversible step whose hazard is linear in `D`.
- The walk has a steep threshold in `D`. It binds at about 0.005·D per site, the tip
  unbinds at 0.55/s, it runs over 30 sites, and each binding uses up one DnaA.
- The PR #204 review estimated the walk's first-passage time for one origin, with no
  DnaA production.
  - The lead reviewer's chain gives 235 s at D = 148 and 693 s at D = 129 (CV 0.79),
    with d log T / d log D ≈ −11 at D = 129.
  - The refuter's model gives 225 s, 638 s and −10.
  - The lumped R2 has −1.
- These are the review's illustrative numbers, not this script's. Both assume one
  origin and no production.

So R2 tests the second-order product. It does not test the threshold. Whether M13
keeps the threshold is a choice for KIM-188, and it is Tom's. One way is to keep some
of the walk's steps as reactions. This note makes no recommendation.

**Other design choices left to M13 (KIM-188).**
- The scale of `D`. The check should say which copy-number regime it tests. The full
  cell's is 148 copies with a steep walk (above).
- The rates. k_D, k_b and the fork rates must make the copy number step inside 600 s
  in a stated share of cells, or the factor's new assumption is not exercised.
- Lumping 68 and 104 fork steps into one exponential step each changes the delay's
  shape, not the kind of non-closure.

**What it reuses from `src/`, and what must be built.** This is from reading the code,
not from running it.

| Piece | Status |
|---|---|
| The ODE modules, `CoreATranscriptDecay`, `CoreATranslation`, `HeldEnzymeCounts`, `m0_models` | reused unchanged |
| A jump module for R1 to R4 | new; each is a `Reaction(rate, affect!)` (`src/interface.jl:130-150`), and the hybrid driver takes any propensity |
| `CoreATranscription` | extended: read `G_g` through `inputs` and multiply the propensity by it, as `transcription.jl:32-37` anticipates ("adding replication later must promote them") |
| The exact route (cSMC) | `propose_events` already simulates every non-transcript reaction exactly with the composition's own propensities (`src/csmc.jl:151-200`), so R1 to R4 need no change there. The transcript bridge does: `_birth_death` reads the birth rate once at the window's start (`csmc.jl:140-147`) and bridges at that constant rate (`src/transcript_bridge.jl`). With `G_g` stepping mid-window, the replication path must be proposed first (it does not depend on transcripts here), and each transcript bridged under a birth rate piecewise constant in time, with a break at the step. The weight's transition probabilities become products over the pieces, as the factor already does for its 59 s + 1 s split |
| The path density (`path_logdensity`, `_accrue_exposure!`, `src/path_replay.jl:116-135`) | reused unchanged: it integrates every propensity at the current state between firings, so `k_g · G_g` and `k_b · D · O` enter as they are |
| Gibbs block 2, the promoter and `krnadeg` updates (`rate_update`, `src/conditional_updates.jl:52-77`, which covers `rate_log_conditional` and `rate_update`; `src/organisms/coreA/gibbs.jl:234-238`) | reused unchanged: the firings and exposures come from the replay, so the exposure becomes ∫k_g·G_g dt with no code change. It needs the promoter to stay linear in the propensity, which holds while the turnover ceiling does not bind, as now. A k_b target would take the same form, with exposure ∫D·O dt |
| The mixed route | `MeanField` throws on R2 (`src/mean_field.jl:105`). How R2 is treated, and whether `TranscriptFactor` takes a `G_g` path or its marginal, is the route decision (KIM-187). M13 then measures that choice against the cSMC |
| SBC and the M0 dataset pipeline (`m0_dataset`) | reused, with the new model list |

**What the M13 model does not contain:** rows 1b, 2, 3, 4, 6, 8, 9, 10 and 12 above.
Rows 7 and 11 are already in M0. Row 1c needs no model, since it closes exactly. Of
the 363 second-order reactions, R2 stands for the 99 DnaA bindings to a 0/1 site, and
loosely for the 2 live bindings to a multi-copy site. It does not stand for:
- the 177 cascade reactions that close exactly;
- the initiations, whose 0/1 × 0/1 non-closure appears only from the second round
  (row 2);
- tRNA charging (62) or cost payment (20). These sit on the translation path and are
  second order in an enzyme state or a counter, not in a 0/1 site.

## Checked by hand

The script does not print these. Each was checked against its source.

| Number | Source |
|---|---|
| Origin 1's cascade is reactions 1595 to 1656 | `cme/reactions.csv`, read directly |
| The volume cap 6.70 × 10⁻¹⁷ L and `CellV` = 670 at the cap | `SpeedyCell/src/Cell/writeback.jl:9-10` |
| 225 and 268 sites per fork | 224 + 1 and 267 + 1, from the printed fork steps |
| Volume doubling at 56 to 72 min (mean 64); daughter initiations at 58 to 110 min | Thornburg et al. 2022, as quoted in `code/cell_cycle/events.py:8-25` |
| The reduced cell grows to about 1.07× | `src/handshake.jl:1721` |
| SpeedyCell.jl's `inventory.md` differs in wording only | `diff` of the two files |
| The DnaA walk's first-passage times and elasticities | the PR #204 review, not recomputed |
| Every `src/…:line` citation | reading the file at `e877f29` |

