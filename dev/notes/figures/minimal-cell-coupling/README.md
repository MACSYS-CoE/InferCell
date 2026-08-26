# Module coupling graphs — well-stirred minimal cell

Two derived figures for [`dev/notes/well-stirred-minimal-cell.md`](../../well-stirred-minimal-cell.md).
Everything here is computed from the model files, not drawn by hand: node
positions are hand-set, every edge and every number is derived.

**Source:** `github.com/Luthey-Schulten-Lab/Minimal_Cell` at `db048ac`, directory
`CME_ODE/`. Clone it next to these scripts as `mc/` and `mkdir out` to re-run them.

## Figures

| File | What it shows |
|---|---|
| `fig1_state_graph.{pdf,png}` | **State coupling** — which module's species reach which other module, and by what protocol |
| `fig2_param_graph.{pdf,png}` | **Parameter coupling** — which parameters constrain which, within and between modules |

## Derived tables

| File | Contents |
|---|---|
| `state_interfaces.csv` | 177 rows: `(source, target, species, kind, protocol, differentiable, provenance)`. The canonical artefact; fig 1 is a view of it. |
| `param_prior_dominance.csv` | Per module, per quantity type: n, dependence (basic/derived), how many sit at prior width |
| `mu0_cross_file.csv` | Standard chemical potentials for compounds active in more than one separately-balanced file |

## Pipeline

```
sbtab.py       minimal multi-table SBtab reader + ReactionFormula parser
modules.py     module membership, transcribed from defMetRxns.py with line refs
handcoded.py   extracts the 61 hand-written model.addReaction blocks
build_state.py  -> state_interfaces.csv
build_param2.py -> param_prior_dominance.csv, mu0_cross_file.csv
draw_state.py   -> fig1   (neato -n, positions hand-set)
draw_param.py   -> fig2
```

## The three things fig 1 shows that a two-block picture does not

Fig 1 is not "CME on the left, ODE on the right, one arrow between". Three
mechanisms sit between the blocks, on three different footings.

**The hook, every 1.0 s** (`hook.py:107`). Carries three kinds of traffic that
are not interchangeable: *mass* (counts ↔ mM), *deferred counters* (the CME
accrues a cost, the hook debits the ODE pool one step later, zero-clamped with
the deficit carried forward — a `max(0,·)` sitting on the interface), and
*catalytic* (protein counts overwrite every ODE rate law's enzyme concentration
via the GPR rule; no mass moves).

**The CME rebuild, every 60 s** (`MinCell_restart.py`). The production driver
loops once per biological minute and reconstructs the whole CME, recomputing
rate constants from live pools — NTP into `k_transcription`, charged tRNA into
`k_translation`, dNTP into replication elongation. This is the ODE→CME direction
and it passes *parameters, not state* — a piecewise-constant coupling.

**Growth.** Surface-area accounting that is in neither block: lipid headgroup
areas plus membrane-protein counts give `CellSA`, hence `r`, hence `V`, and `V`
rescales every count↔concentration conversion in both directions. Volume is
clamped at 2× initial and there is no division event.

## Corrections of record

Kept inline rather than deleted, so the same mistakes are not made twice.

- **190 ODE reactions**, not ~158. An earlier count summed the module tables
  before the code's own turn-offs and missed the 61 hand-coded rate forms
  (ATPase, PUNP5, GHMT2, NADHK, 19 ABC transporters, 40 amino-acid transporters).
  Active per module: Transport 87, Nucleotide 46, Central 26, Lipid 17,
  Cofactor 13, Amino acid 1.
- **Prior dominance is much weaker than the note reports.** Restricted to the
  reactions actually in the model, Michaelis constants sit at prior width for
  12–14% of central and nucleotide reactions (not 77%), and KV for 3–13% (not
  73%). The note's figures were computed over all 244 reactions in the
  reconstruction file, most of which the model never uses. The reactions that
  made it into the model are the better-measured subset.
- **Enzyme concentration is 100% at prior width in every module** — and it does
  not matter, because the ODE never uses the balanced value: it is overwritten
  every hook call from CME protein counts via the GPR rule.
- **The NTP → transcription edge is not severed.** Earlier versions of fig 1
  drew it as a red tee, and annotated Transport→tRNA charging as "the only live
  ODE→CME edge". Both were read off `MinCell_CMEODE.py`, which is the
  *first-minute setup script*. In the production driver the hard-coded scalars
  are commented out and replaced by live pool reads
  (`MinCell_restart.py:442-445`). Same for charged tRNAs: `translation_rate_
  start.py:88,132` uses the frozen `ctRNAconc = 150`, but
  `translation_rate_restart.py:89,113` reads `pmap[trnaID]` and
  `pmap['M_fmettrna_c']` live, and `ctRNAconc` is declared and never used.
  This matters: without this channel the model is unidirectional.
- **What *is* still frozen through the rebuild:** `RnaPconc = 187`
  (`MinCell_restart.py:462`) and `ribosomeConc = 503` (`:506`). Both are
  computed from a hard-coded `r_cell = 2.0e-7`, i.e. the initial volume.
- **`CellSA` does not drive replication.** Earlier versions drew
  `Lipid → DNA replication` labelled "CellSA → division". `CellSA` has no
  consumer outside volume calculation and file I/O, and there is no division
  event anywhere in `CME_ODE/`. The real edge into replication is
  `translation → DNA replication`: `M_DnaA_c` (= MMSYN1_0001,
  `model_data/protein_metabolites_frac.csv:7`) is a translation product, and
  DnaA filament assembly on oriC is what initiates replication
  (`rep_restart.py:55-95`).
- **`memProtList` holds 93 loci**, not 96 (`in_out.py:31`) — and 93 *loci*, not
  93 `M_PTN_*` species: `JCVISYN3A_0779` expands to `ptsg` and `ptsg_P`, so 94
  terms are summed.

## Second review pass

A consistency check against `db048ac` found ten further problems. All were
confirmed against the source and are fixed in the current generator.

Wrong or imprecise in the drawing:

- **ATP was attached to the wrong charging edge.** The table has always carried
  `Transport → tRNA charging` (20 aa pools) and `Central → tRNA charging`
  (`M_atp_c`) as separate rows; the drawing omitted the Central edge and folded
  ATP into the Transport label. Every edge into charging is now drawn from the
  table.
- **The hook's counter traffic was collapsed onto one edge.** Only
  `HOOK → Nucleotide` was drawn, labelled "debit NTP/dNTP pools", hiding the ATP
  hydrolysis debit to Central, the fMet-tRNA debit to Amino acid, and the NMP
  *credit* from mRNA decay (`in_out.py:178-317`). All four are now drawn, with
  the credit styled distinctly.
- **"Currency pool … shared by all six" was false.** No `currency-mass` edge
  involves the Amino-acid module, and this is a fact about the model, not the
  derivation: its one reaction, `R_FMETTRS`
  (`10fthfglu3 + mettrna <=> fmettrna + thfglu3`), uses no currency metabolite.
  The label now says five, and names the exception.

Couplings the table was missing entirely:

- **`Nucleotide → rebuild → DNA replication`.** Replication elongation constants
  are recomputed from live `M_datp_c`/`M_dttp_c`/`M_dctp_c`/`M_dgtp_c` every
  rebuild, with cell volume entering through `getConc` — the same construction as
  the NTPs, in `rep_restart.py:239-243`.
- **`transcription → tRNA charging`.** tRNA transcription produces the uncharged
  tRNA that charging consumes (`MinCell_restart.py:1533,1566`).
- **`translation → tRNA charging`.** The aminoacyl-tRNA synthetases are explicit
  reactants and products of the charging chain (`:1502-1509,1519-1525`).
- **`DNApol3 = 35` is frozen too** (`rep_restart.py:244`), alongside
  `RnaPconc = 187` and `ribosomeConc = 503`. Three polymerase capacities, not two.

Claims that needed qualification (fixed in the figure and in the note):

- The volume clamp is **not** a termination condition. `in_out.py:118-120` stops
  volume *growth* at 6.70e-17 L; the simulation runs for the `runTime` the user
  passed. "~6,300 hook calls per cell cycle" therefore assumes a 105-minute run.
- **54% protein is the calibrated initial area split**, not a statement about
  which term dominates the growth increment. Establishing that needs
  d(CellSA_Prot)/dt against d(CellSA_Lip)/dt on a trajectory.
- **"Differentiable"** is qualified on the ODE-block label: differentiable in
  principle within the ODE subsystem, but the shipped SciPy/LSODA implementation
  is not AD-enabled and the model carries `min`, `max`, clamps and integer
  rounding.
