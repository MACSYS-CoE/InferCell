# The Well-Stirred Minimal Cell (Thornburg et al. 2022)

A map of the CME-ODE whole-cell model of JCVI-syn3A, read from the source. The
companion note [4D-minimal-cell.md](4D-minimal-cell.md) covers the 2026 4DWCM;
this one is deliberately confined to the 2022 well-stirred model, which is the
realistic near-term inference target because it runs in hours on a CPU rather
than days on two A100s.

**Provenance.** Everything below was read from
`github.com/Luthey-Schulten-Lab/Minimal_Cell` at commit `db048ac` (14 Dec 2021),
directory `CME_ODE/`. Line references are to that commit. The whole well-stirred
model is 6 MB and ~10,400 lines of Python — small enough to read end to end,
which is why this note exists.

Paper: Thornburg et al., *Fundamental behaviors emerge from simulations of a
living minimal cell*, Cell 185(2):345-360.e28 (2022),
[doi:10.1016/j.cell.2021.12.025](https://doi.org/10.1016/j.cell.2021.12.025).

The reduction carved out of this model for Step 1 is scoped in
[reduced-syn3a-scoping.md](reduced-syn3a-scoping.md): Core A', 21 metabolic
reactions and 17 genes, against the 190 reactions and 455 loci mapped here.

A note on scope: the same repository also contains `RDME_gCME_ODE/`, the
spatially resolved model. That one covers only the first 20 minutes of the cell
cycle with a static cell structure and non-diffusing ribosomes. The well-stirred
model in `CME_ODE/` is the one that simulates a complete cell cycle.


## Two solvers and a one-second handshake

The model is a stochastic gene-expression block and a deterministic metabolic
block, exchanging state once per second of biological time.

```
   CME (Gillespie, Lattice Microbes)          ODE (odecell -> LSODA)
   gene expression, replication               190 metabolic reactions
   tRNA charging                              transport/nucleotide/central/
                                              lipid/cofactor
            |                                          ^
            |  counts -> initial conditions            |
            +------------- every delt = 1.0 s ---------+
            |        integrate 1 s, write counts back  |
            |                                          |
            +<--- CME REBUILT every 60 s: rate --------+
            |     constants recomputed from live NTP,   |
            |     dNTP and charged-tRNA pools           |
            |     (MinCell_restart.py, rep_restart.py)  |
            |                                          |
            v                                          v
       +----------------------------------------------------+
       |  GROWTH  (in neither block, recomputed every hook)  |
       |  membrane protein counts + lipid counts -> CellSA   |
       |  -> r -> V -> rescales every count<->mM conversion  |
       +----------------------------------------------------+
```

Note the two distinct timescales. The 1 s hook passes *state* both ways. The 60 s
CME rebuild passes *rate constants* from the ODE's pools into the stochastic
block's propensities — a piecewise-constant coupling, and the mechanism by which
metabolism modulates gene expression at all.

The mechanism is a subclass of the Lattice Microbes Gillespie solver that
overrides `hookSimulation` (`program/hook.py:107`). Timesteps are set at
`program/MinCell_CMEODE.py:1663`:

```python
delt    = 1.0        # s   -- CME/ODE communication interval
odestep = delt/10.0  # s   -- max step given to the adaptive integrator
write   = 60.0       # s   -- trajectory write interval
```

So a 105-minute run is ~6,300 hook calls, each one rebuilding the ODE model,
pulling current counts, integrating one second, and writing back. The metabolic
integration is the dominant cost. Note that 105 minutes is a *choice*: `runTime`
is a command-line argument to `MinCell_restart.py` and nothing in the code ends
the simulation at a cell cycle. See Growth and division below.

The important architectural point for InferCell: **this decomposition is already
ours.** A smooth ODE block and a non-differentiable stochastic block, sharing
parameters and exchanging state at a defined interface, is the ISAB three-module
setup at scale. "Differentiable" here means *in principle, within the ODE
subsystem*: the shipped SciPy/LSODA implementation is not AD-enabled, and even
the ODE side carries `min`, `max`, zero-clamps and integer rounding at its
edges. Making it actually differentiable is work the port has to do, not a
property we inherit. We are not proposing a new architecture, we are
proposing inference over an existing one.


## The stochastic block

Per gene locus, three reactions (`MinCell_CMEODE.py:728`):

```python
sim.addReaction(geneMetID, tuple(TrscProd),     TranscriptRate(...))  # transcription
sim.addReaction(rnaMetID,  tuple(mRNAdegProd),  DegradationRate(...))  # mRNA decay
sim.addReaction(rnaMetID,  tuple(TranslatProd), TranslatRate(...))     # translation
```

Membrane proteins get a fourth, `TranslocRate` (`:1064`). There are 455
mRNA-coding loci (`model_data/mRNA_counts.csv`). rRNA and tRNA operons are
handled by separate code paths with their own rate functions
(`riboTranscriptRate`, `:1238`).

tRNA charging is a five-reaction chain per amino acid with rate constants
written directly into the call (`:1307`). Note what its reactants are: the
synthetase protein `M_PTN_*` comes from **translation**, the uncharged tRNA
comes from **tRNA transcription** (29 tRNA loci mapping to 20 uncharged species,
`model_data/trna_metabolites_synthase.csv`), and the amino acid and ATP come
from the **ODE**. Charging is the one place where all three sources meet:

```python
sim.addReaction(('M_atp_c',synthase_ptn),      (synthase_atp),          0.01)
sim.addReaction((aminoAcid,synthase_atp),      (synthase_atp_aa),       0.01)
sim.addReaction((synthase_atp_aa,unchargedMetID),(synthase_atp_aa_trna),0.1)
sim.addReaction((synthase_atp_aa_trna), (chargedMetID,'M_amp_c',...),   30)
```

DNA replication lives in `program/rep_start.py` (624 lines).


## The deterministic block

`defMetRxns.py` (2,588 lines) builds the metabolic ODE system from explicit
per-module reaction-ID lists. This is the model's most useful property for us:
the modularity is already declared in code, not implicit.

| Module | Active reactions | Where |
|---|---|---|
| Transport | 87 | `transport_NoH2O_Zane-TB-DB.tsv` + hand-coded |
| Nucleotide | 46 | `model_data/nucleo_rxns_list.txt` |
| Central | 26 | `cntrMetRxn`, `defMetRxns.py:133` |
| Lipid | 17 | `lipMetRxn`, `:744` |
| Cofactor | 13 | `cofactMetRxn`, `:916` |
| Amino acid | **1** (`FMETTRS`) | `aaMetRxn`, `:288` |

**190 active ODE reactions**, not the ~158 an earlier version of this note
reported. Two things caused that undercount: summing the module lists before the
code applies its own `lipTurnOff`/`nuclTurnOff` exclusions, and missing the 61
hand-written `model.addReaction` blocks (ATPase, PUNP5, GHMT2, NADHK, 19 ABC
transporters, 40 amino-acid transporters). Counts above are from the derived
audit in `figures/minimal-cell-coupling/`, which extracts them from the code
rather than from the module tables.

Drawn from a 340-reaction / 308-species SBML reconstruction
(`model_data/FBA/iMB155_*.xml`, derived from Breuer et al. 2019).

The amino-acid module is worth dwelling on: 22 of its 23 tRNA-synthetase
reactions are commented out at `:288`, because charging was moved into the CME.
So the ODE/CME split is not a clean biological boundary — it is partly a record
of where the authors found it easier to put things.

Rate laws are generated by `odecell` from the SBtab `KineticLaw` column.
Transport uses irreversible Michaelis-Menten; intracellular reactions use
reversible convenience-kinetics forms, e.g. from
`transport_NoH2O_Zane-TB-DB.tsv`:

```
R_GLCpts0  ptsi+M_pep_c <=> ptsi_P+M_pyr_c
           (KF_0_R_GLCpts0*ptsi*M_pep_c)-(KR_0_R_GLCpts0*ptsi_P*M_pyr_c)
```


## Growth and division

Smaller than expected, and worth knowing before designing any reduction: **the
well-stirred model has no division event at all.** Volume *growth* stops once
volume has doubled — the simulation itself carries on for whatever duration the
user asked for. From `in_out.py:107`:

```python
cellRadius = ((SurfaceArea/4/np.pi)**(1/2))*1e-9
cellVolume = ((4/3)*np.pi*(cellRadius)**3)*(1000)
if (cellVolume > 6.70e-17):     # exactly 2x the initial 3.35e-17 L
    cellVolume = 6.70e-17       # growth stops here
```

There is no septum, no partitioning, no daughter geometry. Those arrive only in
the 4DWCM. Here, growth is surface-area accounting on a sphere — and it is worth
treating as a **third block**, since it lives in neither solver, takes input from
both, and feeds both:

```
CellSA_Lip  = 0.513 x SUM_lipids(count x headgroup area)   # 9 species, 0.35-0.6 nm^2
CellSA_Prot = 28.0 nm^2 x SUM(membrane protein counts)     # 93 loci in memProtList
CellSA      = CellSA_Lip + CellSA_Prot
r = sqrt(CellSA / 4pi);  V = (4/3) pi r^3
```

`0.513` is the outer-leaflet fraction of a 5 nm bilayer. `28.0` nm^2 is an
average membrane-protein footprint, chosen to reproduce 54% protein coverage for
the observed ~9.6k membrane proteins — a calibrated constant, not a measurement.

Initial conditions (`setICs_two.py:342`) are 231,875 nm^2 lipid and 270,956 nm^2
protein, totalling 502,831 nm^2, which returns r = 200.0 nm exactly. Note that
`setICs_two.py` carries four commented-out alternatives for this split (40%,
47%, 53%, 54% protein), so the lipid/protein partition was tuned.

**Membrane protein is 54% of surface area** — but note carefully what that is
and is not. It is the calibrated *initial* area split (`setICs_two.py:342`), not
a measured growth rate. Which term dominates the *increment* over a trajectory
does not follow from it, and would need d(CellSA_Prot)/dt compared against
d(CellSA_Lip)/dt on an actual run. What the code does license is structural:
the growth block takes input from translation and from the lipid module both,
so it is not a lipid-module appendage.

Volume is recomputed every hook call and both conversion directions read it
live: `Rxns.partTomM` (counts to mM, on the way into the ODE) and
`in_out.mMtoPart` (mM back to counts) each call `calcCellVolume(pmap)`. So
growth dilutes every concentration in the ODE block. Note that the module-level
`countToMiliMol` constants in `Rxns.py`, `Simp.py` and `defMetRxns.py` are *not*
this quantity — they are computed from a fixed `r_cell = 200 nm` and are used
only for default concentrations. `MinCell_restart.py`'s `getConc` reads
`specDict['CellV']`, so the live volume reaches the 60 s CME rebuild too.

**`CellSA` does not drive replication.** It has no consumer outside the volume
calculation and file I/O, and there is no division event in `CME_ODE/` at all.
Replication is initiated by DnaA filament assembly on oriC (`rep_start.py:32-72`),
and `M_DnaA_c` is a translation product (`model_data/protein_metabolites_frac.csv`
maps it to MMSYN1_0001). The edge into replication is from translation, not from
the membrane.


## Metabolic parameters: three columns, and only one of them runs

The metabolic parameter files are **SBtab**, the format used by Liebermeister's
parameter-balancing method. Each carries several tables, and the distinction
between them is the single most important thing to get right when reading these
files:

| Table / column | What it is | Used by |
|---|---|---|
| `Parameter prior` | log-normal priors per quantity type: `PriorMedian`, `PriorGeometricStd`, bounds | parameter balancing |
| `Parameter` → **`Mode`** | the balanced point estimate | **the simulator** |
| `Parameter` → `UnconstrainedGeometricMean` / `GeometricStd` | the balancing posterior's location and spread | nothing at runtime |
| `Quantity` | per-reaction local parameter values | model construction |

**`Mode` is what runs.** From `defMetRxns.py:191`:

```python
kcatF = ParDF.loc[ParDF["Reaction:SBML:reaction:id"] == rxnID, "Mode"].values[3]
kcatR = ParDF.loc[ParDF["Reaction:SBML:reaction:id"] == rxnID, "Mode"].values[4]
```

Rows 3 and 4 of each reaction's block are the substrate and product catalytic
constants. Do **not** use `UnconstrainedGeometricMean` for flux, and do not use
`catalytic rate constant geometric mean` at all — the latter is
sqrt(kcat_fwd x kcat_rev), which for a near-irreversible reaction is orders of
magnitude below the forward rate. Worked examples of how far wrong each is:

| Reaction | kV(geo) | Unconstrained geo. mean (substrate kcat) | **`Mode` (what runs)** |
|---|---|---|---|
| PYK | 0.905 | 386.6 | **3204** |
| PGK | 895.1 | 1796 | **220** |
| PGM | 90.46 | 112.5 | **434** |
| FBA | 3.662 | 11.75 | **59.7** |

Reading the wrong column produces a model that either silently starves or has
imaginary headroom. This note got it wrong twice before landing on `Mode`.

### The uncertainty is computed and then discarded

That much stands. Every metabolic parameter ships with a prior, a balanced mode,
and a posterior spread, and the simulator consumes only the mode. For the
central + amino-acid module: 244 reactions, 930 Michaelis constants, 244 k_cat,
244 K_eq, 304 concentrations — each with a geometric standard deviation that
never reaches the integrator.

That is a real gap and it is InferCell's opening. But the *size* of the gap is
much smaller than an earlier version of this note claimed.

### How much is actually prior-dominated: far less than first reported

An earlier version reported 77% of Michaelis constants and 73% of k_cat sitting
at prior width. **Those figures were computed over all 244 reactions in the
balanced file, most of which the model never uses.** Restricted to the reactions
actually active in the model, from the derived audit in
`figures/minimal-cell-coupling/param_prior_dominance.csv`:

| Quantity | Earlier claim (all 244 rxns) | Reactions the model uses |
|---|---|---|
| Michaelis constants | 77% | **12-14%** (central and nucleotide) |
| KV | 73% | **3-13%** |

The reactions that made it into the model are the better-measured subset — which
is unsurprising in hindsight and worth remembering as a general caution about
reading these files at reconstruction scope rather than model scope.

Enzyme concentration remains 100% at prior width in every module, and **it does
not matter**: the ODE never reads the balanced value, because
`createGeneExpression` overwrites it from CME protein counts via the GPR rule at
every hook call.

So the honest statement is: *the parameter uncertainty exists, is quantified in
the shipped files, and is thrown away at simulation time* — not that most
parameters are unidentified. Characterising what is actually identified is our
job, and this note should not pre-empt the answer.

## Gene-expression parameters

Metabolism went through parameter balancing. Gene expression did not — its
constants are written into the source. But the coupling story is subtler than it
first looks, and getting it wrong changes the architecture.

### The constants, and which of them are actually constant

`MinCell_CMEODE.py:259-343` — the **first-minute setup** — reads:

```python
rnaPolKcat  = 0.155*187/493*20   # nt/s
rnaPolKd    = 0.1                # mM
krnadeg     = 0.00578/2          # 1/s
ptnDegRate  = 7.70e-06           # 1/s
ATPconc = 1.04; UTPconc = 0.68; CTPconc = 0.34; GTPconc = 0.68   # mM
RnaPconc     = 187*countToMiliMol
ribosomeConc = 503*countToMiliMol
riboKcat = 10; riboKd = 0.0001
```

**Those NTP concentrations are first-minute initial values, not the model's
metabolic pools, and not constants over the run.** Two corrections follow.

First, they are not the metabolic ATP pool: the balanced initial concentration of
`M_atp_c` is **3.6529 mM** (~73,700 particles), from the SBtab `Mode` column.
The 1.04 mM above is an input to the transcription rate constant for minute one.

Second, and more important: **the CME does see metabolism.** The production
driver is `MinCell_restart.py`, which loops one iteration per minute of biological
time (`for sTime in np.arange(1,runTime+1,1)`) and rebuilds the CME each time.
At `:442-445` the hard-coded values are commented out and replaced with live
pool reads:

```python
ATPconc = getConc(max(1,specDict['M_atp_c']),specDict)  #1.04 #mM
UTPconc = getConc(max(1,specDict['M_utp_c']),specDict)  #0.68 #mM
CTPconc = getConc(max(1,specDict['M_ctp_c']),specDict)  #0.34 #mM
GTPconc = getConc(max(1,specDict['M_gtp_c']),specDict)  #0.68 #mM
```

Translation likewise. `translation_rate_restart.py` builds its monomer term from
live charged-tRNA counts rather than a constant:

```python
NMonoSum = NMonoSum + aaCntPtn*riboKd/partTomM(max(1,pmap[trnaID]),pmap)
...
k_translation = kcat_mod / ((1+riboK0/ribosomeConc)*(riboKd**2)
                / (partTomM(max(1,pmap['M_fmettrna_c']),pmap)**2) + NMonoSum + n_tot - 1)
```

Replication elongation is rebuilt the same way. `rep_restart.py:239-243` reads
`M_datp_c`, `M_dttp_c`, `M_dctp_c` and `M_dgtp_c` live, with the hard-coded
values commented out beside them exactly as for the NTPs — so the rebuild
channel carries dNTP pools as well, and the note's earlier "NTP and charged-tRNA
pools" was incomplete.

What does *not* become live at the rebuild: `RnaPconc = 187`
(`MinCell_restart.py:462`), `ribosomeConc = 503` (`:506`) and `DNApol3 = 35`
(`rep_restart.py:244`). All three are computed from a hard-coded
`r_cell = 2.0e-7`, i.e. the initial volume. The three polymerase capacities are
the genuinely frozen quantities in the stochastic block.

So the ODE→CME direction is real. **Its mechanism is a full CME rebuild at 60 s
boundaries, not a propensity that reads pools continuously.** That is the detail
a port has to reproduce: metabolite pools enter the stochastic block through
recomputed rate *constants* once a minute, giving a piecewise-constant coupling
rather than a smooth one. An earlier version of this note claimed the CME
"cannot see the metabolism it is coupled to" — that was read off the setup
script alone and is wrong.

### What is genuinely hard-coded

**The promoter-strength proxy.** In `TranscriptRate` (`:367`):

```python
kcat_mod = min(rnaPolKcat*(ptnCount/(180)), 2*90)
```

Per-gene transcription strength is a deterministic function of that gene's
initial proteomics count over 180 (the mean count), capped at 180 nt/s, scaling
binding affinity and elongation speed with one prefactor. The identical `/180`
construction survives into the 2026 4DWCM, where the authors call the model
*sensitive* to RNAP-promoter binding and name proteomics-assigned promoter
strengths as a known inaccuracy. Load-bearing, acknowledged, and carrying no
uncertainty.

**mRNA degradation has no per-gene parameter.** In `DegradationRate` (`:457`):

```python
kcat  = (18/452)*88   # 1/s # INSTEAD OF 18 or 20
k_deg = kcat / n_tot
```

One global catalytic constant over transcript length, so every gene's mRNA
half-life is a function of its length alone. The trailing comment is the
authors' own record of hand-tuning.

**Transcription is constitutive.** One transcription reaction per gene, firing on
gene presence. There is no two-state promoter, no `k_on`/`k_off`, and therefore
no transcriptional bursting in this model. Identifiability results derived for
bursty two-state promoters do not transfer to it.

## An apparent index transposition in TranscriptRate

Flagging this because it is in the published model and worth verifying properly
before we either reproduce or correct it.

`baseMap` is defined correctly at `:298`:

```python
baseMap = OrderedDict({ "A":ATPconc, "U":UTPconc, "G":GTPconc, "C":CTPconc })
```

But the per-base monomer counts at `:378` do not follow it:

```python
NMono_A = baseCount["A"]
NMono_U = baseCount["C"]
NMono_C = baseCount["G"]
NMono_G = baseCount["U"]

NMonoSum = NMono_A*rnaPolKd/ATPconc + NMono_C*rnaPolKd/CTPconc \
         + NMono_U*rnaPolKd/UTPconc + NMono_G*rnaPolKd/GTPconc
```

The variable names are permuted against the counts in a 3-cycle
(C -> U -> G -> C), and `baseMap` does not compensate: it is used only to
validate the base alphabet and to pick out the first two bases (`CMono1`,
`CMono2`). So the NTP-demand term pairs the count of C bases with UTP's
concentration, G with CTP's, and U with GTP's.

The numerical consequence is smaller than it looks, and worth stating precisely:

- U and G swap harmlessly, because `UTPconc == GTPconc == 0.68` mM. Inert.
- The real effect is C against G. The C-base term is divided by 0.68 instead of
  0.34, halving it; the G-base term is divided by 0.34 instead of 0.68,
  doubling it.

Net effect per transcript scales with the imbalance between its G and C counts,
so there is partial cancellation but not exact. It perturbs a denominator in the
transcription rate, so the direction of bias varies gene by gene.

Status: this is a code reading, not a tested claim. Before asserting a bug we
should confirm the sequence convention feeding `rnasequence` and check whether
the 2026 code carries the same lines.


## What counts as validation in this model

Worth recording, because it defines the bar we would be raising.

The headline metric is whether protein counts double over one cell cycle. In the
2026 paper this is reported as a distribution with median below 2 and a tail out
to ~6x, with the underproduced set identified as long genes (>3 kb). It is a
diagnosis made by inspecting a histogram.

The only sensitivity analysis in the repository is
`simulations/growthMedia-transportSensitivty/transportFluxes-SensitivityAnalysis.ipynb`:
a one-at-a-time scan of external medium concentration against Michaelis-Menten
transport flux, compared to FBA fluxes. Useful, deterministic, local, and not
UQ. Nothing in the repository computes a posterior, an identifiability
diagnostic, or a global sensitivity index.


## Where an ML surrogate would earn its place

Three candidates. They are not alternatives to each other so much as answers to
different questions, and the first two pair naturally.

**1. Replace the promoter-strength proxy with a learned prior over latent
strengths.** The strongest *scientific* claim available. Today 455 promoter
strengths are pinned by `ptnCount/180` with no uncertainty, in a mechanism the
authors themselves flag as sensitive. Instead treat them as latent variables,
and use a sequence-to-strength model — trained on the genome, which ships in the
repo as `model_data/syn3A.gb` — as an informative prior with calibrated
uncertainty. This replaces a modelling choice with a different description of
the same mechanics, which is exactly the substitution we want to demonstrate,
and the uncertainty has a real destination: the downstream stochastic
translation block. It also answers the 2026 paper's own stated next step
(constrain expression from transcriptomics rather than proteomics) with a
principled version of the same move.

**2. Replace the ODE metabolism with a differentiable learned surrogate.** The
strongest *computational* claim. The metabolic integrator runs ~6,300 times per
cell cycle and dominates runtime; a surrogate mapping (enzyme counts, metabolite
state) to (state after 1 s, fluxes) would make the whole graph differentiable
and fast enough for the sampling budgets inference needs. The risk is the
demonstration's whole weight rests on the surrogate's epistemic uncertainty
being *calibrated* — an overconfident surrogate would pass overconfidence
downstream exactly as the smoke-setting NUTS runs did in the ISAB work, and the
result would look better than it is. Calibration of the surrogate has to be
shown, not assumed.

**3. Learned likelihood (NLE) for the stochastic block.** Already the stated
direction in `docs/positioning.md`. Worth distinguishing from the two above: it
is an inference method, not a replacement modelling choice. It makes the
stochastic block joinable to the differentiable graph; it does not make a claim
about the biology. Useful, and orthogonal.

The pairing to aim for: (1) makes a scientific point the original model cannot
make, (2) makes the computation tractable enough to make it at the scale of the
published model.


## Open questions

- **Port or wrap?** `docs/positioning.md` rules out cross-language
  orchestration, since gradients do not cross language boundaries — which
  implies reimplementation in Julia. The mitigating fact is that the durable
  assets are language-neutral: SBML for structure, SBtab for parameters. A port
  means writing rate laws, not re-deriving parameters. This makes issue #15
  (Catalyst.jl as source of truth) load-bearing rather than housekeeping, since
  SBMLToolkit.jl into Catalyst is the natural ingest path.
- Do the balanced parameter files for the nucleotide and lipid modules show the
  same prior-dominated fraction as central+AA? Both carry the same SBtab tables;
  only central+AA was counted here.
- Is the `NMono_*` transposition present in the 2026 code as well?
- What is the actual wall-clock cost of one well-stirred cell cycle on one core?
  The 2026 paper says well-stirred runs take "as little as a few hours"; we need
  a measured number to size any inference budget.
