# What the full syn3A cell has that the mean field cannot close (KIM-177)

2026-11-26. Linear: [KIM-177](https://linear.app/kimpson/issue/KIM-177/non-closure-inventory-of-the-full-cell),
milestone M11. It feeds the design of M13's exact check (F3 panel c), the Spec's
assumptions and the route decision (KIM-187).

Every number here is printed by `nonclosure_inventory.py` in this directory. Its output
is committed beside it as `nonclosure_inventory.out.txt`. To rerun it on the login node:

```bash
/fred/oz022/tkimpson/miniconda3/bin/python3 -I nonclosure_inventory.py > nonclosure_inventory.out.txt
```

It uses the standard library only. It asserts the issue's checks: 363 second-order
reactions, 493 replication reactions, 1,434 count-dependent rate rules, and every
reaction in exactly one family.

## Sources

The exported reference model in `syn3a-ensemble` (read only), at HEAD
`2075809cbb825140cb31a7bec8bc3529770c8e13`. The export itself was made at
syn3a-ensemble `514fa79` from Minimal_Cell `db048ac` (`export_meta.json`).

| Table | sha256 |
|---|---|
| `data/model/cme/reactions.csv` | `63fbeca180fe8bdebec646df0ddbffb0149958fc2ff9fc3e6ec8803cd24c6aaf` |
| `data/model/cme/rate_rules.csv` | `a95619a05e47b32e7e97e8f65128e31a25403c386e80b5ccfb298d905aa3af0e` |
| `data/model/cme/species.csv` | `c1acb80878e8a500a86e0e1d9e195a491b46a01c337788c0a6554df19cbe6ec9` |
| `data/model/nodes.csv` | `1bd4b35e1c33f3547fd5683213b5b4c32d7211dd9046d8a4b7ee628ee1f39b89` |
| `data/model/ode/species.csv` | `77d003243ec52770287ba76658f4deeb00cb3e8ad2e51172d0d3a7517b02bd55` |
| `data/model/ode/reactions.csv` | `4da4941402a0e176ad6adad2c35b6ad0f67f992c023a306ea43f5edb3b8b9a27` |
| `data/model/export_meta.json` | `e3c21a138aa318eb39f5e21a2e1f5e31476d2d10581ecec8e552b1128ce66251` |
| `data/model/inventory.md` | `2cbf4eae3b624a6b310c8eb8e67d30ea01ebd8b6a99c1a7fb12e13187399cca6` |
| `code/SpeedyCell/src/Cell/writeback.jl` | `ae6d3124ca31fddaad5b0a9c08f150ab601b1ae52727fce2ba0c4ec8802e0a41` |
| `code/SpeedyCell/src/Cell/hook.jl` | `1cec291491c64202495cd6820d439657e333f2228001c8cefad0b285d9788720` |

The four CME tables are byte-identical to SpeedyCell.jl's `data/model/cme/`
(`/fred/oz022/tkimpson/SpeedyCell.jl`, HEAD `0823b53`). SpeedyCell.jl's `inventory.md`
differs from syn3a-ensemble's (sha256 `23d6b162…`) in wording only: it drops the
phase and spec references. The numbers and rules are the same. This note cites
syn3a-ensemble's copy.

## Terms and assumptions

**The mean field closes** when its mean equations are exact given the rate constants.
InferCell's `MeanField` (`src/mean_field.jl`) builds them only for a network whose
every propensity is affine in the state and whose every state change is constant. It
throws on anything else ("A propensity … is not affine in the state",
`src/mean_field.jl:105`). A rate constant that is rebuilt from the counts at a minute
boundary is outside that check. The mean-field route evaluates it at the mean pools
(`src/mean_field_route.jl:14`), so its error is a Jensen gap in what reads the mean,
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
feedback, so every item is on some path from every block to every observable. That
answer says nothing. This note reports direct adjacency instead, computed by the
script: does the item read or write a species that a candidate block's reactions
write or read (ODE species, proteins, genes, the pools the rate rules read), and does
it write a species an observable is computed from? The path is then named in prose.

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
daughter origins (`High_Affinity_Site_oriC`, `…_oriC2`, `…_oriC3`). The daughter sites
start at zero and are made by replication reaction 2364 (section 2).

**The 363 and the 493 overlap.** The two replication initiations are second order and
are also two of the 493 replication reactions. The full cell has 363 + 491 = 854
reactions in these two sets, not 856.

**Rate constants.** 361 of the 363 point at a `const` node, so their constants never
change. The two exceptions are the replication initiations (reactions 1873 and 2098),
whose constants are rebuilt each minute from the dNTP pools and `CellSA`. The 361
constant second-order rates do not scale with the cell volume, although the cell
doubles its volume.

**Copy-number scales**, at procid 1's initial counts (`cme/species.csv`):

| Family | Smaller reactant | Larger reactant |
|---|---|---|
| DnaA at the origins | 0 to 16 (sites, filament states) | 0 to 148 (`M_DnaA_c` is 148) |
| tRNA charging | 1 to 736 (enzyme intermediates) | 40 to 73,715 (tRNA, amino acids, ATP) |
| Translation cost payment | 0 (the residue counters) | 160 to 480 (charged tRNA) |
| Replication initiation | 0 (`Initiator_*`) | 1 (`chromosome_*`) |

**A finding in the table.** All 20 cost-payment reactions credit `FMET_cost_paid`. None
credits its own `<AA>_cost_paid` counter. The table has 21 `*_cost_paid` species, and
no reaction writes the other 20. This is reported as found. The script checks it; this
note does not say whether it matters, since nothing in the tables reads a
`*_cost_paid` count.

**Closure.** None of the 363 closes. Each propensity is a product of two counts, so
`MeanField` refuses it. M0 and the reduced model have no second-order reaction. The
reduced model lumps tRNA charging into one deterministic ODE reaction
(`src/organisms/coreA/trna_charging.jl:1-19`), and it pays translation's residue cost
through deferred counters drained at the 1 s handshake, not through a reaction.

**Between a block and an observable.** Adjacency, from the script:

| Family | Reads or writes | Writes an observable's input? | Sits between |
|---|---|---|---|
| tRNA charging | 21 ODE species (ATP and 20 amino acids), 20 synthetase proteins, the charged tRNA pools that all 455 translation rules read | bulk metabolites: ATP and the amino acids are ODE species | yes: translation (gene-expression block) and metabolism (metabolic block) on one side, bulk metabolites and, through translation, doubling time on the other |
| Translation cost payment | the 20 charged tRNA pools the translation rules read; `M_mettrna_c` is an ODE species | indirectly: it sets the charged pools that translation's rates read | yes: translation's own demand feeds back into translation's rates |
| DnaA at the origins | `M_DnaA_c`, made by translation reaction 9 from `M_RNA_JCVISYN3A_0001_c` | no direct write; it makes `Initiator_C` and `Initiator_CC` | yes: from the gene-expression block (DnaA's own expression) to transcripts, through initiation timing and gene copy number (section 4) |
| Replication initiation | `JCVISYN3A_0001_gene`, `JCVISYN3A_0910_gene`; its rate reads the dNTP pools and `CellSA` | genes, which every transcription propensity multiplies | yes: as for DnaA |

None of the four families writes a `CellSA` input directly. They reach doubling time
only through translation (charging, cost payment) or through gene copy number and
then transcription (DnaA, initiation).

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
| `max(1, pool)` floor | 29 pools, in all 1,434 rules | translation floors its pools (`src/organisms/coreA/translation.jl:142-143`); transcription has no floor (`transcription.jl:522-540`) |
| Promoter ceiling `min(·, 180)` | 29 rules read a live protein count | the ceiling is on a parameter, and it never binds (`transcription.jl:61-65`, `turnover_headroom`) |
| Volume cap | reached mid-cycle: the reference paper reports the volume doubling at 56 to 72 min, mean 64 (quoted in `code/cell_cycle/events.py`, `t_Vcap`) | applied (`src/handshake.jl:1732`), but the reduced cell grows to about 1.07× and never nears it (`handshake.jl:1721-1722`) |
| Rounded `CellV` | read by all 486 transcription and 493 replication rules | absent: the reduced model has no `CellV` count |
| Complex AND rules | 9 nodes, 18 ODE reactions | absent: each enzyme is one protein count |
| ODE → counts rounding | round half to even, 146 species | `:fractional_carry` by default (`handshake.jl:194-226`) |
| Counter clamp | 22 counters, `min(counter, pool)` | the hybrid clamps and carries the deficit (`handshake.jl:1878-1879`); the mean-field route does not clamp (`mean_field_route.jl:12-14`) |

**Closure.** None of these rules breaks `MeanField`'s check by itself, since a rebuilt
rate constant is outside it. Each one adds a Jensen gap where the route evaluates the
rule at the mean counts. The gap is zero while a floor, ceiling or cap is far from its
kink. It is not zero for `round` and `CellV`, which step at every count, or for the
volume cap, which a cell reaches mid-cycle. This note gives structure and scale, not the
size of any gap.

**Between a block and an observable.** All four count-dependent families sit between
them, by direct adjacency. Transcription and translation rules are the
gene-expression block's own rate constants. They read pools that the metabolic block
writes (NTPs, charged tRNA through charging) and they write transcripts and proteins,
which feed doubling time through the membrane proteins. The replication rules read
pools the metabolic block writes (dNTPs) and set when each gene's copy number steps,
which multiplies transcription (section 3). The 18 ODE reactions with AND rules are the
metabolic block's rate laws and write bulk metabolites directly.

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
- **The 2 initiations:** they do not close. Both are second order, and both reactants
  are 0 or 1 in a single cell.
- **What the means describe.** A cell's fork is at one site. Its mean is an occupation
  probability spread over 225 or 268 sites. Each gene's copy number steps once, from 1
  to 2, at a random time; its mean rises smoothly. The mean field gets
  `E[G_g(t)]` right given exact rate constants. It cannot give the step.

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

Yes, directly. Replication reads pools the metabolic block writes (the dNTPs). Its
output, `G_g`, multiplies every transcription propensity, which is the gene-expression
block's rate constant, and so writes transcripts, an observable. Through translation
of the membrane proteins it reaches doubling time. Replication does not write a
`CellSA` input itself.
