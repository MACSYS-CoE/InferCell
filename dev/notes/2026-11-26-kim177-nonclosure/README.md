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
