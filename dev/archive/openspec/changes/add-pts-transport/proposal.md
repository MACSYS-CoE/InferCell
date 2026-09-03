## Why

Nothing in Core A′ has any carbon until glucose gets in, and in syn3A glucose
gets in exactly one way: the phosphotransferase cascade, which phosphorylates the
sugar with PEP rather than ATP as it crosses the membrane. That detail is not
cosmetic — it is why the net yield is 2 ATP per glucose rather than 3, and why
PEP is a contested intermediate rather than a free product.

The same change carries lactate out. Both halves are the boundary between Core A′
and its medium, they share four genes' worth of protein bookkeeping, and they are
the only two places in the model where a species leaves the cell.

**Wave:** 1, ODE block, per `dev/plans/reduced-syn3a-wave-plan.md`. This is the
`add-pts-transport` row of that file's wave-1 table.

## Dependencies

- **`establish-corea-interface`** (wave 0). Required and already satisfied:
  merged, archived, and its three delta specs promoted to
  `openspec/specs/corea-interface/`. This change consumes the species registry,
  the seven edge kinds and the provenance-carrying loader, and modifies none of
  them.

No other change is required. `add-central-glycolysis` and
`add-nucleotide-recycling` are parallel siblings in wave 1, not prerequisites:
the species shared with them are named through the registry and declared as typed
edges with peers left unnamed, so this module resolves and integrates alone.

## What Changes

- **The five-step PTS glucose cascade.** `GLCpts0` through `GLCpts4`, as
  mass-action forward/reverse pairs, from `transport_NoH2O_Zane-TB-DB.tsv`:

  ```
  GLCpts0   ptsi + M_pep_c        ⇌  ptsi_P + M_pyr_c
  GLCpts1   ptsh + ptsi_P         ⇌  ptsh_P + ptsi
  GLCpts2   crr  + ptsh_P         ⇌  crr_P  + ptsh
  GLCpts3   ptsg + crr_P          ⇌  ptsg_P + crr
  GLCpts4   M_glc__D_e + ptsg_P   ⇌  M_g6p_c + ptsg
  ```

- **`L_LACt2r`, lactate export**, as passive diffusion across the membrane.

- **Four genes**: ptsI `JCVISYN3A_0233`, Crr `JCVISYN3A_0234`, ptsH
  `JCVISYN3A_0694`, ptsG `JCVISYN3A_0779`.

- **`ptsG` is the only membrane protein in Core A′**, and the module flags it as
  such. This is not a label for its own sake: `in_out.py:51` carries a special
  case, `otherNamesDict = {'JCVISYN3A_0779': ['ptsg','ptsg_P']}`, precisely so
  that `getProtSA` can sum both phospho-forms into the membrane-protein surface
  area at 28.0 nm² each. Surface area sets radius, radius sets volume, and volume
  rescales every concentration in the ODE block. So `add-growth-volume-coupling`
  in wave 2 reads this module's `M_ptsg_c` and `M_ptsg_P_c` directly, and the flag
  is the handle it needs.

- **Nine dynamic states owned**: the registry's eight `:pts` phospho-states
  (`M_ptsi_c`/`M_ptsi_P_c`, `M_ptsh_c`/`M_ptsh_P_c`, `M_crr_c`/`M_crr_P_c`,
  `M_ptsg_c`/`M_ptsg_P_c`) and `M_lac__L_e`.

- **External glucose is declared as a clamped edge at 40 mM with origin
  `:published`.** Not this reduction's simplification: `defMetRxns.py:1338` adds
  every `_e` species to the ODE as a constant parameter, so the published model
  clamps it too. It therefore does *not* appear in the "what is ours"
  enumeration, alongside the chemostatted CTP, UTP and amino-acid pool, which do.

- **This module's own conservation check: PTS protein conservation.** Each of the
  four proteins is conserved across its phosphorylated and unphosphorylated forms,
  because every one of the five cascade steps moves a phosphate between adjacent
  carriers and creates no carrier. Four independent invariants, checked here on a
  standalone trajectory, not deferred to wave 3.

### Lactate export is required, and is not a one-liner to prune

`L_LACt2r` is one line of rate law, which makes it exactly the kind of thing a
later reader deletes as redundant. It is not redundant, and the earlier Core A
spec did delete it — that was an error, and it is recorded here so it is not
repeated.

Glycolysis produces two lactate per glucose. At Core A′'s corrected ~1,106
glucose/s, that is ~13.9 M lactate over a 105-minute cycle: **~691 mM at initial
volume, and ~345 mM even at doubled volume**. Cytosolic lactate would exceed the
entire measured phosphate pool by a factor of forty. Without export the core is
not a sustained pathway at all — it is a pathway that poisons itself in the first
few minutes — and carbon balance is open, so wave 3's second validation check
cannot even be posed.

The figures were ~333 and ~166 mM in an earlier version of the scoping note,
while the glucose demand still omitted tRNA charging. Correcting that demand made
the case for the exporter **twice as strong**, not weaker.

### Two decisions of record

**1. The PTS phospho-state split is the published model's, from proteomics.**
The registry marks all eight states `not_imported`, and the source model offers
two candidate initial conditions that disagree. The transport SBtab's `Compound`
table gives a uniform 0.024 mM to all eight — a placeholder that puts 969 copies
behind every protein, contradicting the proteomics counts for three of the four.
The model's own data files give the real answer:

| File | Rows |
|---|---|
| `protein_metabolites_frac.csv` | `ptsi`/`ptsi_P` 0.05/0.95, `ptsh`/`ptsh_P` 0.05/0.95, `crr`/`crr_P` 0.05/0.95 |
| `membrane_protein_metabolites.csv` | `ptsg`/`ptsg_P` 0.15/0.85 |

The column is `proteomics_fraction`, and it is live model data read at startup by
`MinCell_CMEODE.py:199-203`, not a disabled setting. Applying those fractions to
the scoping note's copy numbers at the model's initial volume reproduces the
totals exactly — ptsI 353, Crr 314, ptsH 290, ptsG 831 — so both the totals and
the split are the published model's, and both are recorded as such rather than
asserted by us.

**2. `M_lac__L_e` is dynamic here and constant in the published model, and the
gap is closed with a volume ratio that is ours.** `defMetRxns.py:1338` clamps
every external species, so published lactate export drains into a pool pinned at
zero and never saturates. The Core A′ registry instead lists `M_lac__L_e` among
its 32 dynamic states. Integrating it in a single shared volume would make export
saturate as the two pools equilibrate — the ~691 mM would not leave, it would
merely move outside — which is a different model, not a faithful one.

So the module integrates `M_lac__L_e` as the registry requires, and scales its
accumulation by a medium-to-cell volume ratio, reproducing the published
unsaturated efflux through a state that is genuinely dynamic. **The ratio is
ours**, registered through `reduction_notes` so it reaches any result that
depends on it.

## Capabilities

### New Capabilities

- `corea-pts-transport/pts-cascade`: the five phosphotransferase steps, the eight
  phospho-states they own, the four protein-conservation invariants they imply,
  and the membrane-protein flag on ptsG that wave 2's growth coupling reads.
- `corea-pts-transport/lactate-export`: the passive-diffusion export law, the
  external lactate state, the volume ratio that keeps efflux unsaturated, and the
  requirement that the reaction not be silently removable.
- `corea-pts-transport/transport-kinetics`: the ten mass-action constants and the
  membrane permeability, all carrying priors this project asserts rather than
  inherits, and the initial conditions with their provenance.
- `corea-pts-transport/transport-coupling`: what this module declares at its
  boundary — the clamped glucose edge, the four mass edges on states it does not
  own, the inputs those require — and that it resolves standalone.

### Modified Capabilities

None. `corea-interface` is consumed unchanged.

## Impact

**New source**

| Path | Contents |
|---|---|
| `src/organisms/coreA/pts_transport.jl` | `PtsTransport <: AbstractSubModel` |
| `src/organisms/coreA/data/pts_transport.tsv` | 11 rate constants, no uncertainty column, so every row loads as `:asserted` |
| `src/organisms/coreA/data/pts_initial_conditions.tsv` | 9 initial conditions with their derivation recorded |
| `dev/scripts/extract_pts_transport.jl` | regenerates both from a source-model checkout |

**Changed source**

- `src/InferCell.jl` — one `include` after the registry.

**New tests**

- `test/test_corea_pts_transport.jl`, included from `test/runtests.jl`.

**Read-only to this branch**

`src/edges.jl`, `src/resolver.jl`, `src/loader.jl`, `src/labels.jl`,
`src/interface.jl`, `src/orchestrator.jl`, `src/organisms/coreA/registry.jl`.

**Asserted priors, which this module contributes and no other wave-1 module
does.** The ten PTS mass-action constants carry no balancing distribution
anywhere in the source model — the transport file has no `Parameter` table at all,
only `Compartment`, `Compound`, `Reaction` and `Quantity`. They are point values.
So is the membrane permeability. Every prior on the eleven is **asserted by this
project**, recorded with informedness `:asserted` rather than `:balanced`, and
enumerable from the composed model. The wave-0 contract already names the ten as
the canonical example of this class; the permeability makes eleven.

**The same orchestrator limitation `add-central-glycolysis` records.**
`_build_rhs` (`src/orchestrator.jl:203-215`) lets a module write only its own
state slice, so this module's production of `M_g6p_c` and `M_pyr_c` is declared as
outbound mass edges and left unexecuted until `add-hook-1s-coupling`. Standalone,
`M_pep_c` and `M_lac__L_c` arrive as held inputs.

**Documentation**

`dev/notes/reduced-syn3a-scoping.md` gains the PTS split and its source files,
which the note's state list currently leaves open, and a note that the transport
SBtab's 42.77 mM external glucose is superseded by `setICs_two.py:279`'s 40 mM.
