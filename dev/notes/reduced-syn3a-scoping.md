# Reduced syn3A: Scoping the Step 1 Model

Step 1 of the three-step ladder toward running InferCell on JCVI-syn3A. Step 2
is the published well-stirred model of Thornburg et al. 2022; Step 3 is the
spatially resolved 4DWCM. This note decides what Step 1 keeps and what it
throws away, and commits to a target.

Read alongside [well-stirred-minimal-cell.md](well-stirred-minimal-cell.md),
which maps the 2022 model this reduction is carved out of. All numbers here were
derived from that model's own data files at commit `db048ac`.

**Decision: the Step 1 target is Core A', the glycolysis/energy core with a live
ATP/GTP interface, lactate export and closed adenylate/guanylate moieties.**
21 metabolic reactions, 17 genes.
Rationale and the rejected alternative are recorded below, because the reasons
matter more than the choice.

Revised after external review, which caught three substantive errors in the
first specification. Corrections of record are kept inline rather than deleted,
so the same mistakes are not made twice:

| Was | Actually |
|---|---|
| The CME cannot see metabolism; channel 4 is a defect | The production loop rebuilds the CME every 60 s from live NTP and charged-tRNA pools. Freezing them would have made the model unidirectional. |
| Transport reduces to the PTS cascade alone | Lactate export is required; without it lactate reaches 166-333 mM |
| FBA is the capacity bottleneck at 9,104 glucose/s | The simulator reads the SBtab `Mode` column. ENO is tightest at 31,029 glucose/s; there is no meaningful bottleneck |
| Initial ATP is 1.04 mM | 3.6529 mM (~73,700 particles); 1.04 mM is a first-minute rate-constant input |
| 24 ODE states | The declared costs need GTP, GDP, AMP and PPi too — ~31 states, and ~32 once the second review round adds GMP |
| ISAB's promoter identifiability result applies | The published motif is constitutive; there is no two-state promoter to be degenerate |

A second review pass caught three more, the first of which is fatal as
specified:

| Was | Actually |
|---|---|
| AMP and PPi are states; nothing consumes them | The lumped charging step is a **dead end**. 3.48M charging events versus an adenylate pool of ~79,800 particles exhausts it after **2.3% of the cell cycle**, and PPi reaches 173 mM. ADK1 and PPA are required |
| mRNA decay returns NMPs, destination unstated | GMP has no route back to GDP/GTP. ~24,000 GMP over a cycle against a ~39,800-particle guanylate pool — 60%, a real leak. GK1 closes it |
| ~533 glucose/s funds the doubling | That counts translation GTP only. Charging costs another 2 ~P per residue, so demand is **~1,106 glucose/s** — a 2x undercount |
| GTP, GDP, AMP sit at 0.1 mM with gstd 10, uninformed | True of the *central* file. The **nucleotide** file — which governs, since these are nucleotide-module species — has GTP 1.6627 (gstd 1.57), GDP 0.2981 (2.65), AMP 0.0832 (3.76), GMP 0.0117 (5.82). Only PPi is genuinely at prior width |

A third pass, from the module-level design work (2026-09-03), caught five more.
These are recorded in full in `spec/spec.md` §4; the corrections of record are:

| Was | Actually |
|---|---|
| "~35 `K_M` from the balanced file (with priors)" — i.e. the balanced column is what runs | The balanced `Parameter` table carries the priors; the **`Quantity` table is what the published simulator reads**, and the two disagree on **8 of the glycolytic module's 32** constants — PGI's G6P constant by 82x, FBA's G3P constant by 227x. Core A′ takes the balanced column for its priors and provenance, which makes those eight **ours, not the published model's** |
| The cross-file trap affects PGK3 and PYK3 | **All five** nucleotide-module reactions are parameterised in both files. GK1 differs by **54x** and PPA by **902x**; the central file's gstd values of 1e33 to 1e63 identify them mechanically as balancing artefacts with no data behind them. Taking the central PPA constant would run the enzyme 902x too fast — an error that hides because it makes a conservation check pass |
| Transcription's NTP demand is charged against the right pools | `TranscriptRate` permutes three of four base counts against the wrong NTP: U gets #C, C gets #G, G gets #U. Present on both code paths including the 60 s rebuild. The rate constant moves ~1%, but the **GTP term is weighted by #U rather than #G, making every transcription rate constant 1.9x more sensitive to the live GTP pool than it should be** — and that is the only ODE→CME channel. Core A′ corrects it and labels the correction |
| The PTS phospho-state split is unspecified (state list leaves it open) | The model's own data files give it: `protein_metabolites_frac.csv` 0.05/0.95 for ptsi, ptsh and crr; `membrane_protein_metabolites.csv` 0.15/0.85 for ptsg, in a `proteomics_fraction` column read live at startup. Applied to the copy numbers these reproduce the totals exactly (353, 314, 290, 831), so both the totals and the split are the published model's. The transport SBtab's uniform 0.024 mM for all eight is a placeholder contradicting the proteomics for three of four carriers. Note also that its 42.77 mM external glucose is superseded by `setICs_two.py:279`'s 40 mM |
| PPi comes from the lumped charging step, which is why PPA is required | Charging dominates, so the argument survives — but **transcription also produces PPi**, one per incorporated NTP monomer across four counters, and the note's phosphate accounting omits it. The published model charges transcription twice over: `n` ATP hydrolysed per nucleotide polymerised *and* one NTP monomer incorporated per base. Phosphate closure must see this flux |

| The lumped charging step's block is unstated, and the wave plan read it as CME | **This note already puts it on the ODE side** and always did: the CME-side count above is "17 genes x 3 reactions + 1 TranslocRate = 52", with no charging reaction; `trna` and `trna_chg` are in the state list of ~32 dynamic states; and the charged pool sits in the energy-interface table beside ATP, ADP, Pi, AMP, PPi, GTP, GDP and GMP. `dev/plans/reduced-syn3a-wave-plan.md` assigned it to the CME and `spec/spec.md` inherited that. Corrected in the spec's §12 amendment 1. So the ODE block carries **21 metabolic reactions plus one lumped step**, and integrating that step deterministically is ours even though the placement is this note's |

One measurement rather than a correction, from the same pass: the NTP-dependent
term contributes only **4.4-5.1%** of each transcription rate constant's
denominator, the rest being transcript length. So **channel 4's elasticity is
about 0.045** — live and bidirectional, but weak: doubling every NTP pool moves a
transcription rate constant by ~3%. That satisfies design criterion 2 and it
bounds how much information can cross. `spec/spec.md` §8 K2 turns this into a
kill criterion on the posterior rather than on the gain.


## What Step 1 is for

**To demonstrate the architecture, not to reproduce syn3A.**

Core A exists to show that an inference-first design works on a real
whole-cell-derived model: heterogeneous formalisms in one graph, a genuine
ODE/stochastic boundary, parameters shared across it, and uncertainty
propagating through rather than being discarded at it.

It is *not* an attempt to reproduce syn3A's phenotype. We have deliberately
removed most of the cell, so Core A will not reproduce the 105 min doubling time
and is not expected to. That is a consequence of the reduction, not a defect in
it.

**The work is the inverse problem.** In almost every case we will be running the
model backwards, not forwards. The forward model's job is to serve as a
likelihood — cheap to evaluate many times, structurally faithful at the
boundary — not to make a prediction that matches an experiment. This changes what
counts as success:

| Not the bar | The bar |
|---|---|
| Reproduces the 105 min doubling time | Parameters recover from synthetic data with nominal coverage |
| Matches syn3A metabolite concentrations | Posterior is calibrated (SBC passes per module) |
| Predicts the observed proteome | Identifiability is characterised: we know what the data constrains |
| Forward trajectories look biological | Uncertainty crosses the boundary correctly |

Forward realism matters only insofar as a degenerate model would make the
inference problem uninteresting. The energy-budget and capacity checks below
exist to rule that out, not to claim biological fidelity.

See dev/notes/modular-bayesian-inference-heterogeneous-modules.md for discussion re inference. 

## The design criterion

A reduction is only honest if it preserves the features that make the inference
problem what it is, and removes only multiplicity. Four features must survive:

1. A genuine ODE/stochastic split with parameters shared across it.
2. Bidirectional coupling. If information flows only one way the boundary is
   trivial and the whole exercise is two loosely-joined toys.
3. Low-copy discrete species, so the stochastic block is not just an ODE with
   noise bolted on.
4. Autocatalytic closure. The cell must make the enzymes that make the energy
   that makes the enzymes. Without the loop there is nothing for a posterior to
   be constrained by.


## Two findings that set the cut

**Every gene is already in the low-copy regime.** From `mRNA_counts.csv`, over
all 455 mRNA-coding loci:

| statistic | mRNA copies per cell |
|---|---|
| mean | 0.31 |
| median | 0.18 |
| max | 2.29 |
| loci below 1 copy | 444 of 455 |

Transcripts exist in zero-or-one copies. This is why the gene-expression block
must be stochastic, and it is the licence to cut gene count hard: each locus is
the same three-reaction motif with different parameters, so going from 455 genes
to 14 removes dimensions without removing phenomena. Criterion 3 is satisfied
automatically at any gene count. It is also the same regime as the ISAB bursty
transcription work, so the machinery already exists.

**The coupling is sparse.** Across the five named biosynthetic module lists —
central, nucleotide, lipid, cofactor, amino acid — only **69 distinct loci** are
enzymes:

| Module | Reactions | With GPR | Distinct genes |
|---|---|---|---|
| Central | 25 | 22 | 20 |
| Nucleotide | 47 | 47 | 23 |
| Lipid | 17 | 17 | 19 |
| Cofactor | 11 | 11 | 9 |
| Amino acid | 1 | 1 | 1 |

This counts the module lists only, not the 87 transport reactions or the 61
hand-coded rate forms, so it is a lower bound on the model's full enzyme
complement (190 active reactions in total — see the companion note). It is still
enough to make the point: the great majority of the 455 loci are expression
machinery (49 ribosomal proteins, 5 RNAP subunits, 30 tRNA synthetases, 95
membrane proteins), replication, and unmodelled functions. They set the capacity
of the CME rather than coupling into the ODE, so most of the genome can go
without touching the boundary at all.


## The four coupling channels

"Preserve the boundary" has to mean something specific. In the 2022 model the
two blocks touch in four places, on **two different timescales**.

Every 1 s (the hook, `hook.py:107`):

1. **Protein counts to enzyme concentrations.** `createGeneExpression`
   (`defMetRxns.py:1360`) resolves each reaction's GPR rule to protein species
   and feeds `pmap` counts into its rate law, overwriting the balanced enzyme
   concentration entirely. Reactions with no GPR get
   `defaultPtnConcentration = 0.001 mM` and are *decoupled* — cutting one costs
   no boundary.
2. **Expression cost to metabolite drain.** The CME accumulates counters
   (`ATP_trsc`, `ATP_translat`, `ATP_mRNAdeg`, per-amino-acid costs,
   NTP/dNTP/NMP counters) which are drained against the pools in
   `in_out.py:178`.
3. **Metabolite pools back into the particle map**, setting the budget available
   at the next interval.

Every 60 s (the CME rebuild, `MinCell_restart.py`):

4. **Metabolite pools into CME rate constants.** The production driver loops once
   per minute and rebuilds the CME, recomputing transcription rates from live
   NTP counts (`:442-445`) and translation rates from live charged-tRNA counts
   (`translation_rate_restart.py`). This is the ODE→CME influence path.

**Channel 4 is a real coupling, not a defect.** An earlier version of this note
had it backwards — read off `MinCell_CMEODE.py`, the first-minute setup, where
the NTP concentrations genuinely are frozen constants. In the production loop
they are not. The consequence for Core A is direct: **freezing NTP or
charged-tRNA pools would delete the only ODE→CME channel and leave the model
unidirectional**, with no "energy makes enzymes" loop and nothing to justify
criterion 2. Core A needs a live energy interface from the start.

**A differentiability hazard in channel 2.** The drain is a clipped subtraction
(`in_out.py`, ATP branch):

```python
if costCnt > atpCnt:
    pmap[costID] = costCnt - atpCnt   # carry the shortfall forward
    pmap['M_atp_c'] = 0
```

A non-smooth `max(0, .)` on the interface. It will obstruct gradients if the ODE
block is sampled with NUTS, and is a natural candidate for smoothing in the port.

## Free cuts

| Cut | Size | Why it costs nothing |
|---|---|---|
| Gene multiplicity | 455 to 14 | Finding 1: no new physics per locus |
| Lipid module | 17 rxns, 19 genes | Membrane lipid supply is made an exogenous boundary condition instead — see Growth below |
| Cofactor module | 11 rxns, 9 genes | NAD/FAD/folate supply; supporting, not core |
| Nucleotide module | 46 rxns to 5 | Keeps PGK3/PYK3 (GTP production) and ADK1/PPA/GK1 (moiety closure). The salvage network, the NMP kinases beyond GK1, and the whole dNTP branch go — the last of these only because replication is cut, so nothing demands dNTPs |
| DNA replication | `rep_start.py`, 624 lines | Cell-cycle machinery, off the boundary |
| Transport | 87 rxns to 6 | PTS glucose cascade (5) + lactate export (1). External glucose stays a fixed 40 mM constant |
| tRNA charging | 20 x 5 rxns to 1 | Lumped into one effective charging step; the published chain's rates are hard-coded inline (0.01/0.1/30/100) with no priors to inherit |
| Amino acid module | 1 rxn | 22 of its 23 reactions are already commented out upstream |


## Core A': the target

Glycolysis through lactate, the PTS glucose cascade, lactate export, a live
ATP/GTP energy interface, and the three reactions that close the adenylate,
guanylate and phosphate moieties. **Twenty-one metabolic reactions, seventeen
genes.**

This is a revision of an earlier "Core A" spec after two rounds of external
review. The second round added the three moiety-recycling reactions below; the
first round's three substantive changes were a live energy interface (an earlier
version froze the NTP pools, which would have deleted the only ODE→CME channel);
lactate export
(without it carbon balance is open and lactate reaches 166-333 mM); and the
`Mode` column for kinetics (earlier capacity tables used balancing-distribution
columns and produced a fictitious bottleneck).

### Metabolic reactions (18)

**Glycolysis through lactate (10).** Every GPR here is a single gene, so none of
`createGeneExpression`'s AND/OR effective-enzyme logic is needed — the
gene-to-enzyme-to-reaction map is one-to-one.

| Reaction | Gene (MM code) | Locus (syn3A) | Copies | S = cnt/180 |
|---|---|---|---|---|
| PGI | MMSYN1_0445 | JCVISYN3A_0445 | 266 | 1.48 |
| PFK | MMSYN1_0220 | JCVISYN3A_0220 | 458 | 2.54 |
| FBA | MMSYN1_0131 | JCVISYN3A_0131 | 775 | 4.31 |
| TPI | MMSYN1_0727 | JCVISYN3A_0727 | 410 | 2.28 |
| GAPD | MMSYN1_0607 | JCVISYN3A_0607 | 1355 | 7.53 |
| PGK | MMSYN1_0606 | JCVISYN3A_0606 | 411 | 2.28 |
| PGM | MMSYN1_0729 | JCVISYN3A_0729 | 322 | 1.79 |
| ENO | MMSYN1_0213 | JCVISYN3A_0213 | 998 | 5.54 |
| PYK | MMSYN1_0221 | JCVISYN3A_0221 | 551 | 3.06 |
| LDH_L | MMSYN1_0475 | JCVISYN3A_0475 | 1100 | 6.11 |

**NOX is dropped.** LDH already regenerates NAD+, so NOX is redundant for redox
balance; it had no GPR rule (hence no boundary), and its k_cat was
prior-dominated. Removing it costs nothing and deletes a parameter that was never
an inference target.

**Lactate export (1).** `L_LACt2r`, passive diffusion, from
`transport_NoH2O_Zane-TB-DB.tsv`:

```
R_L_LACt2r   M_lac__L_c <=> M_lac__L_e
             P_R_L_LACt2r*(M_lac__L_c - M_lac__L_e)*3/r_cell
P_R_L_LACt2r = 5.00e-09 m/s
```

**This is not optional.** Glycolysis produces two lactate per glucose. At Core
A's corrected ~1,106 glucose/s that is ~13.9M lactate over 105 min, or **~691 mM
at initial volume, ~345 mM even at doubled volume**. (The figures here were ~333
and ~166 mM while the glucose demand omitted tRNA charging; the correction makes
the case for the exporter twice as strong, not weaker.) Without export the core
is not a viable sustained pathway. The earlier spec dropped it; that was an
error, and it is the kind of error a balance check catches immediately, which is
why the validation list below leads with balances.

**GTP regeneration (2).** `PGK3` and `PYK3`, the GDP-phosphorylating analogues of
PGK and PYK:

```
PGK3   GDP + 13DPG <=> GTP + 3PG
PYK3   PEP + GDP + H+ --> GTP + pyruvate
```

Both are catalysed by genes already in the core — MMSYN1_0606 and MMSYN1_0221 —
so **they add reactions but no genes**. This is the published mechanism rather
than a lumped invention, and it makes GTP a live product of glycolysis, which is
what the translation cost needs. It also restores a little of the enzyme
competition that the reduction otherwise removes.

Note a cross-file inconsistency to resolve on implementation: PGK3 and PYK3 are
parameterised in *both* balanced files with different modes (PGK3: 319.5 central
vs 140.8 nucleotide; PYK3: 1874 vs 672.8). They are nucleotide-module reactions,
so the nucleotide file's values are the ones the model uses.

**Moiety recycling (3).** Every conserved moiety needs a closed cycle, and
without these three the core is a set of dead ends rather than a sustained
pathway:

```
ADK1   AMP + ATP <=> 2 ADP          MMSYN1_0651 (adk)
PPA    PPi + H2O  -> 2 Pi           MMSYN1_0344 (ppa)
GK1    GMP + ATP <=> GDP + ADP      MMSYN1_0203 (gmk)
```

As written in the reconstruction. Core A' carries PPA as `PPi -> 2 Pi`, since
this is the `NoH2O` model and water is not a state.

| Reaction | Gene | Locus | Copies | S = cnt/180 | kcat_fwd (`Mode`) |
|---|---|---|---|---|---|
| ADK1 | MMSYN1_0651 | JCVISYN3A_0651 | 213 | 1.18 | 319.16 |
| PPA | MMSYN1_0344 | JCVISYN3A_0344 | 190 | 1.06 | 646.73 |
| GK1 | MMSYN1_0203 | JCVISYN3A_0203 | 186 | 1.03 | 410.23 |

All three are active nucleotide-module reactions (lines 6, 52 and 19 of
`nucleo_rxns_list.txt`), so their kinetics and priors are **inherited from the
balanced file rather than asserted by us** — the reason for preferring them over
a re-lumping of the charging step. Each has a single-gene GPR, so the
one-to-one gene-to-enzyme map holds across the whole core.

Why each is required, with the magnitude that makes it required:

- **ADK1** is the fatal one. The lumped charging step converts ATP to AMP, and
  nothing else in the core touches AMP. Doubling the proteome takes 3,484,518
  charging events against a total adenylate pool of ~79,800 particles (3.95 mM):
  **44 pool-equivalents**, so the cell runs out of adenylate after 144 s of a
  6,300 s cycle. ADK1 returns AMP to ADP, which PGK and PYK then rephosphorylate.
- **PPA** is the same argument in the other product. One PPi per charging event,
  unhydrolysed, is **173 mM** — lactate-scale accumulation against a 0.1 mM
  starting pool. It also removes a spurious product-inhibition term if charging
  is written reversibly.
- **GK1** closes the smaller leak. Transcription buries GTP as GMP in mRNA and
  decay releases it free; with PGK3/PYK3 rephosphorylating only GDP, that GMP is
  stranded. Roughly 24,000 GMP over a cycle against a ~39,800-particle guanylate
  pool — about **60%**, not fatal within one cycle but enough to bias the
  transcription rate constants that channel 4 rebuilds from live NTP pools, which
  is precisely the coupling Core A' exists to exercise. The same argument applies
  to transcription's AMP return, which ADK1 already handles.

The general rule this encodes, and the one the earlier spec violated twice:
**every declared cost needs both a state that can pay it and a reaction that
returns it.** Validation checks 4 and 4b below are the mechanical version.

**PTS glucose uptake (5).** How glucose actually gets in:

```
GLCpts0   ptsi + PEP         <=> ptsi_P + pyruvate
GLCpts1   ptsh + ptsi_P      <=> ptsh_P + ptsi
GLCpts2   crr  + ptsh_P      <=> crr_P  + ptsh
GLCpts3   ptsg + crr_P       <=> ptsg_P + crr
GLCpts4   glucose_e + ptsg_P <=> G6P + ptsg
```

| Gene | Locus | Copies | Length (aa) | Note |
|---|---|---|---|---|
| ptsI | JCVISYN3A_0233 | 353 | 574 | |
| Crr | JCVISYN3A_0234 | 314 | 155 | |
| ptsH | JCVISYN3A_0694 | 290 | 90 | |
| **ptsG** | JCVISYN3A_0779 | 831 | 746 | **membrane protein** |

These are mass-action forward/reverse constants and **carry no SBtab priors** —
unlike the balanced glycolytic parameters. They are point values with no
quantified uncertainty, which makes them poor early inference targets and a place
where a prior has to be asserted rather than inherited.

The PTS consumes PEP rather than ATP, which is why the net yield is 2 ATP per
glucose: of the two PEP produced, one imports the next glucose and one reaches
PYK.

### The energy interface

The minimum needed to keep channel 4 alive, with everything else declared
chemostatted:

| Quantity | Treatment | Why |
|---|---|---|
| ATP / ADP / Pi | **dynamic** | produced by PGK, PYK; drained by transcription and charging |
| AMP | **dynamic** | charging's other product; returned to ADP by ADK1 |
| PPi | **dynamic** | charging's other product; hydrolysed to 2 Pi by PPA |
| GTP / GDP | **dynamic** | produced by PGK3, PYK3; drained by translation |
| GMP | **dynamic** | returned by mRNA decay; phosphorylated to GDP by GK1 |
| Charged tRNA (lumped) | **dynamic**, one effective species | translation propensity reads it; charging consumes ATP |
| CTP, UTP | chemostatted | sources are in the nucleotide module, out of scope |
| Amino acids | chemostatted | supplied in the medium at 0.1 mM |
| RNAP (187), ribosomes (503) | chemostatted | external expression capacity |
| External glucose (40 mM), O2 | chemostatted | medium |

The lumped charging step, replacing the 20 per-amino-acid chains:

```
M_trna_c + ATP  ->  M_trna_chg_c + AMP + PPi
```

so that translation propensity depends on a live pool that costs ATP to
maintain. This is a documented lumping, not a published reaction — flag it as
such wherever results depend on it. Its AMP and PPi are what ADK1 and PPA exist
to clear; the alternative of re-lumping it as `trna + 2 ATP -> trna_chg + 2 ADP
+ 2 Pi` closes the same moieties by fiat, but at the cost of putting two more
entries on the "ours, not the model's" list. We take the published reactions
instead.

**The loop, stated explicitly:** glucose → PEP/13DPG → ATP and GTP → transcription
and translation of the seventeen genes → enzymes → glucose, with ADK1/PPA/GK1
returning the spent adenylate, phosphate and guanylate to the top of the loop.
Both directions are live: protein counts set enzyme concentrations every 1 s, and
metabolite pools set CME rate constants every 60 s.

### State list

Every declared cost must have a state that can pay it. The earlier spec failed
this: it debited GTP for translation and recycled NMPs without carrying GTP, GDP,
AMP or PPi at all.

| Group | Species |
|---|---|
| Glycolytic intermediates (11) | G6P, F6P, FDP, DHAP, G3P, 13DPG, 3PG, 2PG, PEP, PYR, LAC_c |
| Adenylate (4) | ATP, ADP, AMP, Pi |
| Guanylate (3) | GTP, GDP, GMP |
| Redox (2) | NAD+, NADH |
| Other (2) | PPi, LAC_e |
| tRNA (2) | trna, trna_chg |
| PTS phospho-states (8) | ptsi/ptsi_P, ptsh/ptsh_P, crr/crr_P, ptsg/ptsg_P |
| Chemostats (5) | glucose_e, CTP, UTP, amino-acid pool, O2 |

**~32 dynamic states + 5 chemostats.** H2O and H+ are dropped (`NoH2O` model,
hydrogen-ion accounting removed), which is also why PPA is written as
`PPi -> 2 Pi`.

CME side: 17 genes x 3 reactions + 1 `TranslocRate` for ptsG = **52 reactions**.

### Initial conditions, from the SBtab `Mode` column

| Species | mM | gstd | | Species | mM | gstd |
|---|---|---|---|---|---|---|
| ATP | 3.6529 | 1.28 | | Pi | 17.8185 | 1.06 |
| ADP | 0.2178 | 2.90 | | NAD+ | 2.1844 | 1.45 |
| G6P | 3.7076 | 1.28 | | NADH | 0.0253 | 4.96 |
| F6P | 0.8538 | 1.92 | | PYR | 3.3660 | 1.31 |
| FDP | 7.6037 | 1.14 | | 3PG | 1.1015 | 1.78 |
| DHAP | 0.6445 | 2.10 | | 2PG | 0.0272 | 4.88 |
| 13DPG | 0.0098 | 6.02 | | PEP | 0.0409 | 4.46 |

Note the ATP pool is **3.6529 mM (~73,700 particles)**. The 1.04 mM that appeared
in an earlier version of this note is not the metabolic pool — it is an input to
the first-minute transcription rate constant.

The nucleotide species are **not** uninformed, contrary to an earlier version of
this note. That claim came from reading the central file, where they sit at the
0.1 mM / gstd 10 prior default because they are out of that module's scope. They
are nucleotide-module species, so the nucleotide file governs — and it has real
balanced values:

| Species | mM | gstd | | Species | mM | gstd |
|---|---|---|---|---|---|---|
| GTP | 1.6627 | 1.57 | | AMP | 0.0832 | 3.76 |
| GDP | 0.2981 | 2.65 | | GMP | 0.0117 | 5.82 |

This is the same cross-file trap as the PGK3/PYK3 kinetics, in the concentration
table rather than the parameter table, and it cuts the other way: the guanylate
pool is ~39,800 particles rather than the ~4,000 that 0.1 mM would imply. **Read
initial conditions from the file that owns the species**, and record which file
each value came from in the port.

**Two states genuinely have no measured initial condition:** PPi and lactate
(and G3P) sit at 0.1 mM with gstd 10 in both files — the prior median at prior
width. Nothing informed them. For PPi it barely matters once PPA is present,
since the pool is turned over ~1,700 times per cycle and any starting value is
forgotten within seconds; for lactate it matters even less, since it is an
export product. Worth knowing before treating either as a trusted starting
point.

### Three copy-number regimes in one model

With a 200 nm cell radius, 1 mM = 20,180 particles:

| Species class | Copies per cell | Natural formalism |
|---|---|---|
| mRNA | 0 to 2 | discrete, SSA required |
| Enzymes (protein) | 266 to 1355 | SDE / chemical Langevin sits naturally here |
| Metabolites | ~2,000 (0.1 mM) to ~74,000 (ATP) | ODE justified |

Nearly five orders of magnitude of copy number inside twenty-one reactions and
seventeen genes — the clearest available demonstration of why a whole-cell model
needs several formalisms at once, and the reason Core A' is a real test of the
architecture rather than a toy.

### Growth: retained, as fractional growth

The 2022 model has **no division event**: volume growth is capped at 6.70e-17 L
(exactly 2x initial) and stops. Growth is surface-area accounting on a sphere
(`in_out.py:107`):

```
CellSA_Lip  = 0.513 x SUM_lipids(count x headgroup area)
CellSA_Prot = 28.0 nm^2 x SUM(membrane protein counts)
r = sqrt(CellSA / 4pi);  V = (4/3) pi r^3
```

ptsG is Core A's only membrane protein: 831 copies x 28 nm^2 = 23,268 nm^2, so
doubling it moves surface area from 502,831 to 526,099 nm^2 — r from 200.0 to
204.6 nm, volume to ~1.07x.

**Report this as fractional growth or time-to-threshold, never as a syn3A
doubling-time prediction.** Core A' omits ~92% of the cell; it is not built to
reproduce the 105 min doubling time and any comparison to it would be
meaningless. What makes the growth law worth keeping is unrelated to phenotype:
`Rxns.partTomM` calls `calcCellVolume(pmap)`, so volume feeds the
particle-to-mM conversion and **growth dilutes every concentration in the ODE
block** — a second coupling channel from the stochastic block into the
deterministic one, for ten lines of code.

### Sanity check: the model is not degenerate

Forward capacity from the `Mode` column, expressed as glucose per second, with
turnovers-per-glucose accounted for (note LDH needs **two**, not one):

| Enzyme | kcat_fwd (`Mode`) | Copies | Per glucose | Max glucose/s |
|---|---|---|---|---|
| PYK | 3204 | 551 | 1 | 1,765,404 |
| TPI | 759.6 | 410 | 1 | 311,440 |
| PGI | 804.3 | 266 | 1 | 213,954 |
| LDH_L | 388.5 | 1100 | 2 | 213,657 |
| GAPD | 119.6 | 1355 | 2 | 81,024 |
| PGM | 434 | 322 | 2 | 69,874 |
| PFK | 111 | 458 | 1 | 50,838 |
| FBA | 59.7 | 775 | 1 | 46,268 |
| PGK | 220 | 411 | 2 | 45,210 |
| **ENO** | **62.18** | 998 | 2 | **31,029** |

Doubling the seventeen proteins needs **3,484,518 amino acids**. The energy
budget has two terms of equal size, and an earlier version of this note counted
only the first:

| Cost | Per residue | Over the cycle |
|---|---|---|
| Translation elongation GTP | 2 ~P | 6,969,036 ~P |
| tRNA charging (ATP → AMP) | 2 ~P | 6,969,036 ~P |
| **Total** | **4 ~P** | **13,938,072 ~P** |

At a net 2 ~P per glucose that is **~1,106 glucose/s**, not the ~533 reported
before — the old figure omitted charging entirely. **The tightest constraint is
still ENO at 31,029 glucose/s, so Core A' runs at ~3.6% of capacity.** The
conclusion survives the correction with a factor of 28 to spare; there is no
meaningful bottleneck.

The three added enzymes are nowhere near limiting either, which is what makes
them cheap to include:

| Enzyme | Capacity (kcat x copies) | Demand | Utilisation |
|---|---|---|---|
| ADK1 | 67,981 /s | 553 AMP/s | 0.8% |
| PPA | 122,879 /s | 553 PPi/s | 0.5% |
| GK1 | 76,303 /s | ~4 GMP/s | 0.005% |

An earlier version of this note claimed FBA was the bottleneck at 9,104
glucose/s and built an argument on it — that FBA's loose prior therefore governed
the predicted doubling time. That came from the balancing-distribution columns
rather than `Mode`, and it does not survive. **No single kinetic parameter is a
natural inference target on capacity grounds.**

The purpose of this check is narrow: it confirms the autocatalytic loop can fund
itself with large margin, so the coupling is live rather than pinned at a bound.
It is not a fidelity claim.

### Parameters, and which to actually infer

Available: ~10 glycolytic k_cat_fwd/k_cat_rev pairs and ~35 K_M from the balanced
file (with priors); 2 GTP-branch and 3 moiety-recycling reactions (with priors);
10 PTS mass-action constants (**no priors**); 17 promoter strengths; ~6
gene-expression globals; one lumped charging rate.

The five nucleotide-file reactions all carry forward-kcat gstd 1.051 — the
tightest band in the core — so adding them widens the parameter vector without
loosening it.

**Infer 3-6 of them initially, not all of them.** Recommended first set:

1. Two or three promoter strengths — the `/180` proxy is the thing the surrogate
   work aims to replace, so this is the target that matters.
2. One shared gene-expression global (`rnaPolKcat` or `riboKcat`).
3. One glycolytic k_cat with a genuinely informative prior, as a control: it
   should recover tightly, and if it does not, the machinery is wrong.

Everything else fixed at published values for the first runs, then freed
selectively.

Note a semantic correction from review: `rnaPolKcat`, `riboKcat` and `krnadeg` are
**CME-local** parameters whose effects propagate into the ODE through the
coupling. They do not literally appear in both blocks' rate laws, and calling
them "shared" in the ISAB sense overstates it.

### Validation before any inference

In this order, because each catches errors the next would mask:

1. **Non-negativity** of all states over a nominal trajectory.
2. **Carbon balance** — glucose in versus lactate out plus biomass. This is the
   check that would have caught the missing lactate exporter immediately.
3. **Redox balance** — NAD+/NADH pool conserved.
4. **Adenylate and guanylate conservation** — ATP+ADP+AMP and GTP+GDP+GMP,
   allowing for what expression consumes and what ADK1/GK1 return. Run this over
   a *full* cycle, not a short burst: the failure mode it catches is a slow
   one-way drain, and 144 s of a 6,300 s trajectory looked fine before the
   adenylate dead end was found.
4b. **Phosphate closure** — Pi + all phosphorylated species, with PPi hydrolysis
   counted as 1 PPi to 2 Pi. This is the check that catches a missing PPA.
5. **PTS protein conservation** — each of ptsi/ptsh/crr/ptsg conserved across its
   phosphorylated and unphosphorylated forms.
6. **A nominal trajectory** that looks sane, at published parameters.

Only then SBC.


## The inverse problem: what we condition on

The forward model is now specified in detail and the observables are not. Under
an inference-first framing that is backwards, so this section records the design
space and one trap.

### What the model can emit

The 2022 model writes at `write = 60.0` s intervals: per-species counts (every
gene's mRNA and protein, every metabolite), metabolic fluxes for every reaction
(`solver.calcFlux`), and the `CellSA`/`CellV` counters. Core A inherits all of
this. So the candidate observables are time series of counts, fluxes, and
volume, at 60 s resolution, over as many replicate cells as we care to simulate.

### The circularity trap

**Do not condition on proteomics.** Per-gene promoter strength in this model is
`S = ptnCount/180` — a deterministic function of the proteomics counts. Those
counts have already been used to set the parameters. Conditioning on proteomics
would therefore use the same data twice: once implicitly through the prior on
promoter strength, once explicitly in the likelihood. The posterior would come
out impressively tight for entirely spurious reasons, and the tightness would
not announce itself as an artefact.

This is the same prior-double-counting hazard that shows up in modular
inference, in a concrete and easily-missed form. It also applies to initial
protein counts, which come from the same table.

### The options

| Data | Ground truth? | Circular? | Notes |
|---|---|---|---|
| **Synthetic, from known theta** | Yes | No | The only option with ground truth. Coverage and SBC are computable. |
| Metabolite time courses | No | No | Not used to set any Core A parameter. Genuinely independent. |
| Single-cell mRNA distributions | No | No | Matches the low-copy regime; distributional, not mean-based. |
| Proteomics counts | No | **Yes** | Sets the promoter proxy. Avoid. |
| Doubling time / growth rate | No | No | But Core A is not built to reproduce it, so uninformative here. |

**For Step 1a: synthetic-data recovery.** Simulate Core A at known parameter
values, infer, and check that the posterior covers truth at nominal rates. This
tests the machinery, which is what Step 1a is for, and it is the only route that
admits a ground truth to check against. Real data enters at Step 2, where the
model is large enough for the comparison to mean something.

### Which parameters are the targets

Not all ~60. The interesting subsets are the ones that cross or define the
boundary:

- The 17 promoter strengths — currently pinned by the `/180` proxy, and the thing
  the surrogate work aims to replace.
- The shared gene-expression globals (`rnaPolKcat`, `riboKcat`, `krnadeg`) —
  these are the analogue of ISAB's `k_tl`, `gamma_mRNA`, `gamma_prot`, the
  parameters that appear in both blocks.
- FBA's `kcat_fwd` — the loosest prior in the core (gstd 1.747, against 1.05-1.17
  for everything else). Note this is a prior-width argument only: an earlier
  version also called FBA the capacity bottleneck, which the `Mode`-column
  correction retracted. ENO is tightest on capacity, and at ~3.6% utilisation
  neither constrains anything.

The remaining metabolic k_cat and K_M values are well determined (gstd
1.05-1.17) and can reasonably be fixed in early runs, then freed to see whether
the inference notices.

### Summary statistics

**One imported result does not apply here.** The ISAB work found that promoter
rates are identifiable only through the product
`k_tx_burst * k_on/(k_on+k_off)` — but that result is about a *two-state*
promoter. The 2022 model's transcription is **constitutive**: one reaction per
gene, firing on gene presence, with no `k_on`/`k_off` and no bursting. So the
degeneracy that motivated distributional summaries in ISAB is not present in Core
A' as specified.

Distributional summaries are still the right choice, for a different reason:
every gene sits at 0-2 mRNA copies, so the count distribution across replicates
carries information that any mean-based summary destroys. The justification is
the low-copy regime, not the bursting degeneracy.

Worth noting as a modelling option rather than an inherited fact: *adding* a
two-state promoter to Core A' would make the ISAB result directly relevant and
would be a defensible model-comparison question (constitutive versus bursty, by
Bayes factor). That is a deliberate extension, not a reproduction of the
published model.


## An important caveat about Core A'

**Core A' is the best-measured region of the network, not a representative sample
of it.** The glycolytic reactions carry genuinely informative priors, and that is
not an accident: the reactions that made it into the published model are the
better-measured subset of the reconstruction.

Note that the companion note's prior-dominance figures were corrected downward
substantially by the derived audit in `figures/minimal-cell-coupling/`: restricted
to reactions the model actually uses, Michaelis constants sit at prior width for
12-14% and KV for 3-13%, not the 77% and 73% computed over all 244 reactions in
the balanced file. So "most parameters are unidentified" was never a safe claim
at model scope, and Step 1a should not be framed as demonstrating it.

Two consequences:

- **Good:** Step 1a is a well-posed problem with real informative priors. If
  parameter recovery fails here, the machinery is broken rather than the problem
  being hard. That is what you want from a first rung.
- **Constraint:** Step 1a cannot demonstrate an identifiability gap, because
  there is not much of one in this subnetwork. What it demonstrates is that the
  architecture works.

**Two places where Core A' is *easier* than the same subnetwork inside the full
model**, and both should be stated wherever Step 1a's results are reported:

1. **Enzyme competition is partly removed.** MMSYN1_0221 catalyses nine reactions
   in the published model and MMSYN1_0606 four. Core A' keeps PYK/PYK3 and
   PGK/PGK3 — two of nine and two of four — so the GTP branch restores some
   competition but not all of it. The moiety-recycling genes are the same story
   in miniature: MMSYN1_0651 catalyses ADK1 and DADK, and MMSYN1_0203 catalyses
   GK1 and DGK1, of which Core A' keeps one each. Both dropped partners are dNTP
   reactions, so they return with the RNDR branch at Step 1b — which makes the
   competition Core A' removes here exactly the competition Step 1b restores.
2. **The PTS parameters carry no priors at all.** Ten mass-action constants with
   point values and no quantified uncertainty. Any prior on them is asserted by
   us, not inherited from the model.


## Core B: rejected for Step 1a, adopted for Step 1b

The alternative was to build the reduced cell around a low-copy enzyme, so that
noise as well as parameter uncertainty crosses the boundary. Every low-copy
metabolic enzyme sits outside central metabolism:

| Reaction(s) | Gene | Copies/cell |
|---|---|---|
| DUTPDP | MMSYN1_0447 | 10 (at the `defaultPtnCount` floor) |
| CLPNS | MMSYN1_0147 | 15 |
| DADNK, DGSNK, DCYTK | MMSYN1_0382 | 23 |
| 5FTHFPGS | MMSYN1_0823 | 35 |
| PGSA | MMSYN1_0875 | 36 |
| FTHFCL | MMSYN1_0443 | 38 |
| DASYN | MMSYN1_0304 | 40 |
| **RNDR1-4** | **MMSYN1_0772** | **50** |

Median copies by module: central 301 (min 148), nucleotide 306, lipid 130,
cofactor 72.

Ribonucleotide reductase is the interesting one: ~50 copies, catalysing four
reactions, converting NTPs to dNTPs — and the 2026 paper reports replication
rate responding to dNTP pools. Roughly 14% Poisson noise feeding a pool that
gates a downstream discrete process.

**Why not first.** Core A's enzymes are 266 to 1355 copies, so relative Poisson
noise at the boundary is 3 to 6%. The CME-to-ODE channel is therefore nearly
deterministic, and Core A honestly cannot support a claim that stochasticity
crosses the boundary. It supports a claim about *parameter uncertainty*
crossing, which is a different and sufficient first result. Building Core B
first would mean a harder model, a bigger reaction set (dNTPs need an NTP
source) and a less well-conditioned inference problem, all before the machinery
has been shown to work at all.

**Step 1b is Core A plus the RNDR/dNTP branch.** That rung adds the genuinely
stochastic boundary, makes the NTP pools live (fixing channel 4), and is where
both surrogate proposals from the companion note earn their keep.


## What the reduction costs

**What survives — what the architecture demonstration needs.** A genuinely
bidirectional ODE/stochastic boundary: protein counts setting enzyme
concentrations every 1 s, metabolite pools setting CME rate constants every 60 s.
Published parameter values carrying real priors, for the glycolytic core.
Three copy-number regimes spanning nearly five orders of magnitude, so multiple
formalisms are required rather than decorative. A closed autocatalytic loop with
large margin. Volume-mediated dilution as a second coupling channel. The `/180`
promoter proxy the surrogate work aims to replace. And a model small enough to
run thousands of times, which is the property the inverse problem depends on.

**What is lost, and is fine to lose.** No DNA replication, no septum or
partitioning geometry, no cell-wide phenotype, no claim on syn3A as an organism.
Core A' will not reproduce the 105 min doubling time. Under a forward-modelling
framing that is the headline weakness; under an inference-first framing it is out
of scope. Growth is reported as fractional growth or time-to-threshold.

**What is lost and needs watching.** Four things could make the inference
problem *misleading* rather than merely small:

- **Enzyme competition is partly removed** — see the caveat above. Core A' will
  look better-conditioned than the same subnetwork does inside the full model.
- **Core A' is the best-measured region** of the network, so recovery succeeding
  is evidence the machinery works, not that inference on the full model is
  well-posed.
- **Several treatments are ours, not the model's**: the lumped tRNA charging
  step, chemostatted CTP/UTP/amino acids, and exogenous non-ptsG membrane growth.
  Each is a documented simplification, and each is a place where a result could
  depend on our choice rather than on the published model. Label them wherever
  they matter.
- **The nucleotide module is cut to five reactions chosen by what the rest of the
  core happens to demand**, not by a biological criterion. Core A' keeps ADK1,
  PPA and GK1 because the charging lump and transcription would otherwise leave
  dead ends; syn3A's actual nucleotide economy is the salvage network, which
  synthesises nothing de novo and is cut entirely. So the core's NTP pools are
  maintained by a mechanism the real cell does not primarily use. This is
  defensible while CTP and UTP are chemostatted, and stops being defensible the
  moment Step 1b makes the NTP pools live.

None of these is a reason to build something bigger first. All are reasons to be
careful about what Step 1a's success licenses.


## Open questions

- **Wall-clock cost of one Core A' trajectory.** Now the binding practical
  number, because the inverse problem needs thousands. The 2022 full model is
  reported at "a few hours" per cell cycle; Core A' should be orders of magnitude
  cheaper, but this needs measuring before any inference budget can be sized.
  **First measurement, 2026-09-05 (spec phase 3, Slurm job 16212707).** On a
  two-module toy — one jump gene-expression module, one ODE metabolite module,
  exchanging every simulated second — the 1 s handshake costs 1.50–1.51e-06 s of
  wall-clock per simulated second with the pool frozen — steady across 60, 300
  and 600 s horizons — and 1.57e-06 s with the rate law live at 600 s, rising to
  2.16e-06 s at the shortest horizon. The worst row extrapolates to **0.014 s
  per 6,300 s trajectory** (0.010 s at 600 s), three orders of magnitude inside
  the 10 s budget the spec's kill criterion K1 is stated against. Read it as a *floor*
  rather than an estimate: the toy carries one gene and two metabolite states
  against Core A''s seventeen and thirty-two, so the real composition can only
  cost more. What it does establish is that the exchange itself carries no
  fatal per-handshake overhead. The number at full scale is spec task 13.7, and
  that is where K1 is decided. Details in
  `dev/scripts/bench_handshake_result.md`.
- **Does the 60 s CME rebuild need reproducing faithfully?** The published
  ODE→CME coupling is piecewise-constant: rate constants recomputed once a
  minute, not propensities reading pools continuously. A continuous version would
  be a different (arguably better) model. Decide deliberately, and if we deviate,
  measure what it changes.
- **PGK3/PYK3 cross-file inconsistency.** Both appear in the central and
  nucleotide balanced files with different modes (319.5 vs 140.8; 1874 vs 672.8).
  The nucleotide file governs, but the discrepancy is worth understanding — see
  `mu0_cross_file.csv` in the audit for the analogous problem in chemical
  potentials. The same trap has now bitten twice, once in kinetics and once in
  initial concentrations, so the port should carry the source file for every
  imported value rather than resolving it once by hand.
- **Does the lumped charging step's stoichiometry matter?** Real synthetases
  produce AMP + PPi, which is why ADK1 and PPA are needed; a `2 ATP -> 2 ADP`
  lumping would close the moieties without them. The two differ in how much
  traffic they put through ADK1 (553/s versus zero), and so in how much the ADK1
  equilibrium constant can influence the ATP/ADP ratio that the 60 s CME rebuild
  reads. Worth one comparison run, since it is cheap and the answer determines
  whether ADK1 belongs in the inference target set.
  **Measured 2026-09-23 (spec task 9.9, job 17023622,
  `dev/scripts/trna_charging_diagnostics_result.md`).** Matched on steady flux,
  not event count, with both forms as mass-action laws: ADK1 carries 548/s
  under AMP + PPi and 0 under 2 ATP, but the ATP/ADP ratio the rebuild reads
  differs by only 1.02% (8.097 vs 8.180). Pyrophosphate settles at 0.367 against
  0.024 mM. So the stoichiometry decides ADK1's *traffic* entirely and moves its
  lever on ATP/ADP by about a percent, at least against the phase 8 glycolytic
  double. Whether ADK1 belongs in the target set now turns on whether that
  percent is visible in a posterior (spec R8), not on the traffic.
- **Which observable is most informative about the boundary-crossing
  parameters?** Answer by simulation before committing to a data model: simulate
  Core A', then measure how much each candidate observable moves each target
  parameter.
- **Do we inherit the `/180` promoter proxy verbatim?** For synthetic recovery it
  barely matters — it only sets the truth we recover. It matters at Step 2.
  Leaning toward inheriting it as the baseline the surrogate is measured against.
- **The clipped ATP drain** (`max(0, .)` in channel 2) will obstruct gradients if
  the ODE block is sampled with NUTS. Decide whether to smooth it in the port.
- **The lumped charging step is ours.** Whether one effective charged-tRNA pool
  is adequate, or whether translation needs per-amino-acid resolution to behave
  correctly at low copy number, is untested.
