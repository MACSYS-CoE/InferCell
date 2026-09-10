# Core A′ vendored parameter extracts

Derived extracts, not upstream files. Each is produced by a script under
`dev/scripts/` from a checkout of
[`Luthey-Schulten-Lab/Minimal_Cell`](https://github.com/Luthey-Schulten-Lab/Minimal_Cell)
at commit **`db048ac`**, and **nothing here is fetched at run time** — compute
nodes have no network, and the test suite must be able to read these files with
Base alone.

The raw upstream tables are not vendored. `read_source_table` assumes one
identifier column per row, and the balanced `Parameter` table keys each row on
the triple `(QuantityType, Reaction id, Compound id)`; `src/loader.jl` is
framework and is not amended on a module branch (spec §10 R15). So the reshape
happens here, before the loader sees the file.

**Every extract must round-trip.** Re-running the command below over an
unchanged checkout reproduces the file byte for byte. That is what keeps the
regeneration script honest: if it did not, the vendored file and the script
would be two independent copies of the same numbers, which is the trap spec §4
D2 records this project hitting three times.

## `central_glycolysis.tsv`

| | |
|---|---|
| Upstream file | `CME_ODE/model_data/Central_AA_Zane_Balanced_direction_fixed_nounqATP.tsv` |
| Upstream commit | `db048aca5fe85438e0129819bbf0314b037dd931` |
| Upstream table | `TableType='Quantity' TableName='Parameter'` (the balanced table, `Document='MinCell_params_Central_Zane_direction_fixed.tsv'`) |
| Script | `dev/scripts/extract_central_glycolysis.jl` |
| Consumed by | `src/organisms/coreA/central_glycolysis.jl` (spec §11 phase 6) |

### Regenerate

```
julia dev/scripts/extract_central_glycolysis.jl <path to Minimal_Cell checkout>
```

Base only, no project — it runs on a login node with an empty depot in about two
seconds. The output path defaults to this directory.

The script **refuses a checkout at the wrong revision** rather than regenerating
against it: the row-class counts catch a table that changed shape, but nothing
catches one that changed a value, and a silent regeneration against a different
upstream would round-trip cleanly against the wrong source. A checkout that is
not a git repository is warned about rather than refused, since a tarball is a
legitimate way to have the file.

It prints the SHA-256 of what it wrote. For the committed file that is

```
c5278483a9d8bc9e6c91bb6da7b5065581f050d9eeea9460c215c02a4a2d9322
```

**Nothing in the test suite can check the round-trip**, because CI and the
compute nodes have no `Minimal_Cell` checkout. What the suite does assert is
that the committed file holds exactly the 65 identifiers the model needs, in the
right classes, with four cells matching upstream byte for byte, and that all
thirteen concentrations agree with `registry.jl`. A hand-edit to a *concentration*
therefore fails the loader's registry-agreement check; a hand-edit to a catalytic
or Michaelis constant would not, and the hash above is what catches it.

### The 65 rows, and where each class came from

| Rows | Identifier | Upstream `!QuantityType` |
|---|---|---|
| 10 | `kcatF_R_*` | `substrate catalytic rate constant` |
| 10 | `kcatR_R_*` | `product catalytic rate constant` |
| 32 | `km_R_*_M_*` | `Michaelis constant` |
| 13 | `conc_M_*` | `concentration` |

All from the one `Parameter` table. The reshape **renames identifiers and
selects rows; it changes no value.** Every `Mode` and `GeometricStd` cell is
copied as the string the upstream file holds rather than parsed and reprinted,
so `kcatR_R_PGI` stays `650` and not `650.0`. `!UpstreamRow` carries the
original triple, so a renamed identifier stays traceable to the line it came
from.

The 32 Michaelis constants are 2 (PGI) + 4 (PFK) + 3 (FBA) + 2 (TPI) + 5 (GAPD)
+ 4 (PGK) + 2 (PGM) + 2 (ENO) + 4 (PYK) + 4 (LDH_L), which is the
substrate-plus-product term count of each reaction's rate law.

### The width column keeps its upstream name, and that is not cosmetic

The header reads `!UnconstrainedGeometricStd`, not a bare `!GeometricStd`,
because **there is no balanced width in the source.** The upstream `Parameter`
table carries exactly one spread column and it is the unconstrained one, so
every prior built from this extract pairs a *balanced* median (`Mode`) with an
*unconstrained* width.

For most rows that costs nothing, and the file says which: **all 13
concentrations and all 32 Michaelis constants have `Mode` equal to
`UnconstrainedGeometricMean` exactly**, so the width belongs to the median it is
given. **Ten of the twenty catalytic constants do not** — those are the rows the
thermodynamic balancing actually moved — and there the width is inherited from a
distribution not centred on the prior's median:

| identifier | `Mode` (the prior's median) | `UnconstrainedGeometricMean` | width used |
|---|---|---|---|
| `kcatR_R_TPI` | 4 | 65341.6923 | 1.0513 |
| `kcatR_R_PGK` | 0.128899064884157 | 446.1665 | 1.0513 |
| `kcatF_R_PGK` | 220 | 1795.6486 | 1.0513 |
| `kcatF_R_PYK` | 3204 | 386.6434 | 1.0513 |
| `kcatF_R_PGM` | 434 | 112.5355 | 1.0921 |
| `kcatR_R_PGM` | 14 | 72.7139 | 1.1457 |
| `kcatF_R_FBA` | 59.7 | 11.7469 | 1.7466 |
| `kcatR_R_FBA` | 0.56 | 1.1418 | 8.3563 |
| `kcatR_R_PGI` | 650 | 1001.2634 | 1.0513 |
| `kcatF_R_PFK` | 111 | 111.7687 | 1.093 |

Spec §4 D1 carries this as a labelled departure, and §12's 2026-09-10 entry
records how it was found. Emitting the column as `!GeometricStd` would have
hidden exactly the mode-versus-mean distinction D1 exists to police, which is
why the rename was undone.

### Three columns not read, each for a different reason

- **`!UnconstrainedGeometricMean`** — the unconstrained estimate, not the
  balanced one. Spec §4 D1 records reading the wrong column twice before landing
  on `Mode`. It is read here only to *report* the ten divergences above, never
  to supply a value.
- **The `Quantity` table**, earlier in the same document — the hand-patch layer
  the published simulator actually runs. Eight of the 32 Michaelis constants
  differ from the balanced column, two of them by 82× and 227×. Taking the
  balanced column is a deliberate departure (spec §4 D3): it gives every constant
  a prior and one provenance, where the `Quantity` column carries no uncertainty
  at all. The module registers it through `reduction_notes`, so
  `reduction_declarations` enumerates it.
- **`!KineticLaw`**, in the `Reaction` table — inert for these reactions. The
  published simulator builds every central-metabolism rate from
  `Rxns.Enzymatic(substrates, products)`, not from that string (spec §4 D1).

### The thirteen concentrations are the thirteen *owned* states

Eleven glycolytic intermediates and the redox pair — not `M_atp_c`, `M_adp_c` or
`M_pi_c`, which this module reads but does not integrate. The `conc_<species>`
prefix is what fires the loader's registry-agreement check, so each of the
thirteen is checked against `src/organisms/coreA/registry.jl` on every load and
the two copies of every initial condition check each other.

The archived reference for this phase
(`dev/archive/openspec/changes/add-central-glycolysis/tasks.md`, task 1.1) names
`conc_M_atp_c` among its four spot values while also fixing the count at
thirteen. The two cannot both hold, and the count is what spec §11 task 6.1
carries; `conc_M_g6p_c = 3.7076` (gstd 1.2785) is asserted in its place. Whoever
vendors the adenylate species — phase 8, which needs the central file's rival
values in any case (spec §5) — picks up `M_atp_c` there.
