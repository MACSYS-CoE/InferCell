## 1. The vendored extracts

- [ ] 1.1 Write `dev/scripts/extract_pts_transport.jl`, taking a path to a
  `Minimal_Cell` checkout, and run it to produce
  `src/organisms/coreA/data/pts_transport.tsv`. Verify the output has exactly 11
  data rows, **no `GeometricStd` column**, and the values of design.md D2 —
  `KF_1_R_GLCpts1 = 200000`, `KR_4_R_GLCpts4 = 1.00e-05`, `P_R_L_LACt2r = 5.00e-09`
  among them — each carrying the upstream row it came from.

- [ ] 1.2 Extend the same script to produce
  `src/organisms/coreA/data/pts_initial_conditions.tsv`: the eight phospho-state
  concentrations derived as copy number × `proteomics_fraction` at 20,180
  particles per mM, plus `conc_M_lac__L_e = 0.0`. Verify the nine rows match
  design.md D3 to seven decimals (`conc_M_ptsi_c = 0.0008746`,
  `conc_M_ptsg_P_c = 0.0350025`, …), that each records both the copy number and
  the fraction it was derived from, and that the script fails loudly if either
  upstream input has moved.

- [ ] 1.3 Write `src/organisms/coreA/data/README.md` covering both extracts:
  upstream filenames (`transport_NoH2O_Zane-TB-DB.tsv`,
  `protein_metabolites_frac.csv`, `membrane_protein_metabolites.csv`,
  `setICs_two.py`), commit `db048ac`, the derivation for the phospho-states, and
  three recorded discrepancies — the transport table's uniform 0.024 mM
  placeholder, its 42.77 mM external glucose superseded by `setICs_two.py:279`,
  and `getProtSA`'s docstring saying 35 nm² where its code uses 28.0. Verify by
  re-running the regeneration command from the README and getting unchanged files.
  Append rather than overwrite if `add-central-glycolysis` has already created
  this file.

## 2. The sub-model, its rate laws and its parameters

- [ ] 2.1 Create `src/organisms/coreA/pts_transport.jl` with
  `PtsTransport <: AbstractSubModel`, `formalism = :ode`, and `states` returning
  the registry's eight `:pts` species plus `M_lac__L_e`. Wire one `include` into
  `src/InferCell.jl` after the registry and one into `test/runtests.jl` for
  `test/test_corea_pts_transport.jl`. Verify `using InferCell` loads clean and a
  test asserting `length(states(m)) == 9`, that the set equals
  `species_in_group(:pts) ∪ [:M_lac__L_e]`, and that `species_index.(states(m))`
  is strictly increasing, passes under `sbatch test/run_tests.slurm`.

- [ ] 2.2 Implement the five cascade steps as reversible second-order mass action
  per design.md D1. Verify unit tests: each step's net rate is zero when forward
  and reverse terms balance; doubling both substrates of a step quadruples its
  forward term (confirming no saturation); and no step's stoichiometry names
  `M_atp_c`, `M_adp_c`, `M_amp_c` or `M_pi_c`.

- [ ] 2.3 Implement `L_LACt2r` as `P · (lac_c − lac_e) · 3 / r_cell` with `P` and
  `r_cell` carried separately per design.md D6. Verify the rate is zero when the
  two lactate pools are equal, negative when external exceeds cytosolic, that it
  equals `0.075 · (lac_c − lac_e)` at `r_cell = 200 nm`, and that halving
  `r_cell` doubles the rate constant.

- [ ] 2.4 Import all 20 values through `load_parameter` against the two extracts
  as two `SourceTable`s, with an `asserted_gstd` keyword setting every rate
  constant's prior width and `fixed = true` by default plus a `free` keyword.
  Verify every one of the 11 rate constants loads with informedness `:asserted`;
  that the composed model's asserted-prior enumeration returns all 11 and no
  balanced parameter; that a rate constant reports `"transport"` as its source
  while an initial condition reports `"model_ics"`; and that
  `ambiguity_report` over the two tables is empty.

- [ ] 2.5 Confirm the initial conditions reproduce the published copy numbers.
  Verify a test summing each carrier's two states, converting at 20,180 particles
  per mM, and asserting 353, 314, 290 and 831 exactly; and that each
  phospho-state's provenance reports both the copy number and the
  `proteomics_fraction` behind it.

## 3. Boundary declarations and the membrane flag

- [ ] 3.1 Implement `coupling` returning the clamped `M_glc__D_e` edge at 40 mM
  with `origin = :published` and the four `MassEdge`s of design.md D8, every one
  with `peer` unnamed. Verify a test asserting exactly five edges with the
  tabulated species, kinds and directions; that no currency, deferred-counter,
  rate-constant or volume edge is present; and that constructing the clamp at the
  transport table's 42.77 mM throws against the registry's 40.

- [ ] 3.2 Implement `inputs` returning exactly `M_pep_c` and `M_lac__L_c`. Verify
  `resolve_coupling` on the module alone raises no inputs/edges disagreement in
  either direction; that adding `M_glc__D_e` to `inputs` makes it throw naming
  the chemostat; and that removing `M_pep_c` makes it throw naming that species.

- [ ] 3.3 Implement `membrane_protein_states` as a plain function over
  `PtsTransport` per design.md D7, returning `M_ptsg_c` and `M_ptsg_P_c` with the
  28.0 nm² footprint and a record that the footprint is the published model's
  calibrated constant rather than a measurement. Verify a test asserting exactly
  those two states, that no other Core A′ species is returned, and that the
  footprint carries its calibrated-not-measured label.

- [ ] 3.4 Implement `reduction_notes` carrying two entries: the medium-to-cell
  volume ratio, naming the published constant-external-pool behaviour it
  reproduces and the dynamic state it replaces; and the held metabolites in any
  composition that does not execute the coupling. Verify `reduction_declarations`
  returns both, and — the negative case that matters — that the external glucose
  clamp does **not** appear, because its origin is `:published`.

- [ ] 3.5 Verify standalone resolution reports rather than fails: a test composing
  the module alone asserts `resolve_coupling` succeeds, lists `M_pep_c` and
  `M_lac__L_c` among the unowned states, and distinguishes this deliberate
  partial composition from an incomplete one.

## 4. Dynamics and standalone integration

- [ ] 4.1 Implement `dynamics`, assembling `du` for the nine owned states from the
  six reaction rates, reading `M_pep_c` and `M_lac__L_c` from `u_inputs` and
  `M_glc__D_e` from its fixed held value. Verify a hand-checked unit test at one
  state vector: each species' derivative equals the signed sum of the rates of the
  reactions that touch it, and no derivative is returned for `M_g6p_c`,
  `M_pyr_c`, `M_pep_c` or `M_lac__L_c`.

- [ ] 4.2 Implement external lactate's accumulation as `+v_export / R` with the
  medium-to-cell volume ratio `R` exposed, marked asserted, and defaulting to 1e5
  per design.md D5. Verify a test that `R = 1` conserves total lactate across the
  two pools and drives the export rate toward zero, while the default keeps
  external lactate below one percent of cytosolic over a full-cycle integration.

- [ ] 4.3 Add a `HeldMetabolites <: AbstractSubModel` double to
  `test/test_corea_pts_transport.jl` owning `M_pep_c` and `M_lac__L_c` with zero
  derivatives, integrate the two-module composition from the module's initial
  conditions at published parameters with a stiff solver and explicit
  `abstol`/`reltol`, and verify a trajectory is produced with all nine owned
  states non-negative at every saved point.

## 5. Conservation, and the record

- [ ] 5.1 **PTS protein conservation — this module's own conservation check.**
  Assert each of `M_ptsi_c + M_ptsi_P_c`, `M_ptsh_c + M_ptsh_P_c`,
  `M_crr_c + M_crr_P_c` and `M_ptsg_c + M_ptsg_P_c` equals its initial value at
  every saved point of the 4.3 trajectory, as four independent bounds derived from
  the solver tolerances. Verify all four pass and that a failure message names
  which carrier drifted.

- [ ] 5.2 Show the conservation check can fail, and fails locally. Verify a
  mutation test that perturbs one cascade step so a carrier is created rather than
  transferred, asserts that carrier's check fails, and asserts the other three
  still pass.

- [ ] 5.3 Record why lactate export is not removable, in the code rather than only
  in the planning artifacts: the ~691 mM at initial volume and ~345 mM at doubled
  volume, and that an earlier specification deleted the reaction in error. Verify
  the composed model can be queried for this rationale and returns both magnitudes.

- [ ] 5.4 Add the PTS phospho-state split and its two source files to the state
  list in `dev/notes/reduced-syn3a-scoping.md`, which currently records that the
  split is unknown, together with the 42.77-versus-40 mM external glucose
  discrepancy. Verify the note and `design.md` D3 state the same fractions,
  totals and copy numbers.

- [ ] 5.5 Run the full suite and record the count. Verify `sbatch
  test/run_tests.slurm` passes with no failures, then update `docs/handoff.md`
  with what changed, the job id, and the two things wave 2 inherits: that
  `membrane_protein_states` is a plain function awaiting promotion to the
  protocol, and that the four conservation invariants must be restated as
  "conserved up to what translation adds" once translation feeds the carriers.
