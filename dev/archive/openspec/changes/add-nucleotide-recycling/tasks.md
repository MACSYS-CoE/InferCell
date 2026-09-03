## 1. The two vendored extracts

- [ ] 1.1 Write `dev/scripts/extract_nucleotide_recycling.jl`, taking a path to a
  `Minimal_Cell` checkout, and run it to produce
  `src/organisms/coreA/data/nucleotide_recycling.tsv` from
  `Nucleotide_Kinetic_Parameters.tsv`. Verify 35 data rows — 10 kcat, 17 K_M, 8
  `conc_` — and spot values `kcatF_R_PPA = 646.727` (gstd 1.0513),
  `kcatF_R_GK1 = 410.227`, `conc_M_gtp_c = 1.6627` (gstd 1.5684),
  `conc_M_ppi_c = 0.1` (gstd 10).

- [ ] 1.2 Extend the script to produce
  `src/organisms/coreA/data/nucleotide_recycling_central.tsv` — the central
  balanced file's rows for the same 35 identifiers. Verify it holds all 35, that
  `kcatF_R_PPA = 583611.6071`, `kcatF_R_GK1 = 7.5991`, and that every one of the
  17 K_M rows reads 0.1 at gstd 10.

- [ ] 1.3 Verify the ambiguity is real: a test loading both extracts asserts that
  all 35 identifiers have two holders; that `ambiguity_report` returns 35 entries;
  that `conc_M_atp_c`, `conc_M_adp_c`, `conc_M_pi_c` and `conc_M_ppi_c` are marked
  agreeing and the other four `conc_` rows disagreeing; and that the five forward
  catalytic constants are all marked disagreeing.

- [ ] 1.4 Append to `src/organisms/coreA/data/README.md` (creating it if
  `add-central-glycolysis` has not) the upstream files, commit `db048ac`, why the
  rival extract exists, and the regeneration command. Verify re-running that
  command from the README leaves both files unchanged.

## 2. The sub-model, its rate law and its imports

- [ ] 2.1 Create `src/organisms/coreA/nucleotide_recycling.jl` with
  `NucleotideRecycling <: AbstractSubModel`, `formalism = :ode`, and `states`
  returning the registry's four `:adenylate` and three `:guanylate` species plus
  `M_ppi_c`. Wire one `include` into `src/InferCell.jl` after the registry and one
  into `test/runtests.jl`. Verify `using InferCell` loads clean and a test
  asserting `length(states(m)) == 8`, that the set equals
  `species_in_group(:adenylate) ∪ species_in_group(:guanylate) ∪ [:M_ppi_c]`, and
  that `species_index.(states(m))` is strictly increasing, passes under `sbatch
  test/run_tests.slurm`.

- [ ] 2.2 Implement the five reactions with the modular rate law of design.md D1,
  expanding stoichiometric coefficients above one into repeated terms. Verify unit
  tests: ADK1's product term scales as the square of the ADP saturation ratio and
  PPA's as the square of the phosphate ratio; both use one Michaelis constant per
  species, giving 17 across the five reactions; each reaction is zero at balance
  and negative under product excess; and no rate law names water or a hydrogen ion.

- [ ] 2.3 Import all 35 values through `load_parameter` against both extracts, each
  with an explicit `governing` declaration — the nucleotide file for all kinetics,
  and for each initial condition the file `registry.jl` records as that species'
  source. Verify every parameter reports its chosen file and its rejected
  alternative; that `governing_choices` returns all 35; and that dropping any one
  `governing` declaration makes the load throw naming the identifier and both
  values.

- [ ] 2.4 Verify the governing choice for initial conditions matches the registry's
  own `source_file` field for all eight species — `M_atp_c`, `M_adp_c`, `M_pi_c`
  and `M_ppi_c` to `"central_balanced"`, `M_amp_c`, `M_gtp_c`, `M_gdp_c` and
  `M_gmp_c` to `"nucleotide_balanced"` — and that a deliberately mismatched
  declaration fails rather than importing from a file the registry rejects.

- [ ] 2.5 Verify informedness is derived, not asserted: the central extract's rows
  for these five reactions all classify as `:prior_default` (their geometric
  standard deviations run from 1e33 to 1e63, and all 17 of their K_M sit at prior
  width); the nucleotide file's forward constants classify as `:balanced` at gstd
  1.0513; `kcatR_R_PYK3`, `kcatR_R_GK1` and `kcatR_R_PPA` classify as
  `:prior_default` per design.md D4; `M_ppi_c` stays uninformed; and the composed
  model's asserted-prior enumeration contains none of this module's parameters.

- [ ] 2.6 Add the five enzyme concentrations of design.md D5 and expose the
  locus-to-concentration map. Verify each equals `copies / 20180` to six decimals
  (ADK1 0.010555, PPA 0.009415, GK1 0.009217, PGK3 0.020367, PYK3 0.027304), and
  that `PGK3` and `PYK3` report loci `JCVISYN3A_0606` and `JCVISYN3A_0221` as
  shared while the three recycling reactions report none.

## 3. Boundary declarations

- [ ] 3.1 Implement `coupling` returning the four `MassEdge`s and five
  `CurrencyEdge`s of design.md D7, every one with `peer` unnamed. Verify a test
  asserting exactly nine edges with the tabulated species, kinds and directions;
  that no edge names `M_atp_c` or `M_adp_c`; and that no deferred-counter,
  catalytic, rate-constant, volume or clamped edge appears.

- [ ] 3.2 Implement `inputs` returning exactly `M_13dpg_c` and `M_pep_c`. Verify
  `resolve_coupling` on the module alone raises no inputs/edges disagreement; that
  the five inbound currency edges on owned states require no input; and that
  removing `M_13dpg_c` from `inputs` makes it throw naming that species.

- [ ] 3.3 Verify cross-module edge consistency against the two sibling changes:
  compose this module's declarations with those of `add-central-glycolysis` and
  `add-pts-transport` and assert no species and direction is described as both
  mass and currency, per the table in design.md D7. Skip with a recorded reason if
  neither sibling has landed yet.

- [ ] 3.4 Implement `reduction_notes` carrying the held inbound metabolites in any
  composition that does not execute the coupling, and — the negative case —
  verify `reduction_declarations` contains **no** Michaelis-constant column
  deviation, because all 17 agree between the nucleotide file's two tables.

- [ ] 3.5 Verify standalone resolution reports rather than fails: composing the
  module alone asserts `resolve_coupling` succeeds, lists `M_13dpg_c` and
  `M_pep_c` among the unowned states, and distinguishes this deliberate partial
  composition from an incomplete one.

## 4. Dynamics, and the reason each reaction exists

- [ ] 4.1 Implement `dynamics`, assembling `du` for the eight owned states from the
  five reaction rates, reading `M_13dpg_c` and `M_pep_c` from `u_inputs`. Verify a
  hand-checked unit test at one state vector: each species' derivative equals the
  signed sum of the rates that touch it with its stoichiometric coefficient — two
  for ADK1's ADP and PPA's phosphate — and no derivative is returned for
  `M_13dpg_c`, `M_3pg_c`, `M_pep_c` or `M_pyr_c`.

- [ ] 4.2 Implement `reduction_notes` entries recording why each recycling
  reaction is required, with its magnitude: ADK1 against 553 AMP/s and a ~79,800
  particle pool exhausted in ~144 s of 6,300; PPA against 173 mM of unhydrolysed
  pyrophosphate; GK1 against ~24,000 stranded GMP on a ~39,800 particle pool.
  Verify the composed model returns all three, and that PGK3 and PYK3 are absent
  from that list.

- [ ] 4.3 Add `ChargingDrain <: AbstractSubModel` and `HeldGlycolytic` doubles to
  `test/test_corea_nucleotide_recycling.jl` — the first draining ATP to AMP and
  PPi at 553.1 per second, the second owning `M_13dpg_c` and `M_pep_c` with zero
  derivatives. Verify the three-module composition builds and every input resolves.

- [ ] 4.4 Integrate the composition over a full 6,300 s cycle with a stiff solver
  and explicit `abstol`/`reltol`. Verify a trajectory is produced, all eight owned
  states stay non-negative, `M_atp_c` stays positive throughout, and `M_ppi_c`
  settles to a bounded plateau rather than growing. Record the wall-clock, which
  design.md's second open question needs.

## 5. Conservation, and the record

- [ ] 5.1 **Adenylate and guanylate conservation — this module's own check, over a
  full cycle.** Assert `M_atp_c + M_adp_c + M_amp_c` and
  `M_gtp_c + M_gdp_c + M_gmp_c` each equal their initial values at every saved
  point of the 4.4 trajectory, as two independent bounds derived from the solver
  tolerances. Verify both pass, that the interval is at least one full cycle, and
  that a failure names which moiety drifted.

- [ ] 5.2 **Phosphate closure — this module's second conservation check.** Assert
  `Pi + 3·ATP + 2·ADP + AMP + 3·GTP + 2·GDP + GMP + 2·PPi` is exactly invariant
  with PGK3 and PYK3 inactive, and invariant after subtracting the phosphate
  delivered across the two inbound mass edges with all five active. Verify both
  forms pass and that pyrophosphate is counted as two throughout.

- [ ] 5.3 **Reproduce the dead end of record.** Verify that removing ADK1 under the
  same drain exhausts `M_atp_c` and `M_adp_c` within a few minutes — of the order
  of the recorded 141–144 s — and that removing PPA makes `M_ppi_c` rise without
  bound. Both must fail the 5.1 and 5.2 checks respectively, so a passing module
  is evidence rather than a tautology.

- [ ] 5.4 Verify the governing machinery is load-bearing: a test that loads only
  the nucleotide extract asserts the `governing` declarations no longer fail on
  anything, demonstrating that the rival extract is what makes them meaningful,
  and that all 35 identifiers are present in both files.

- [ ] 5.5 Correct `dev/notes/reduced-syn3a-scoping.md`, which records the
  cross-file inconsistency for PGK3 and PYK3 only, in both its open-questions list
  and its GTP-regeneration section. Add ADK1, GK1 and PPA with their ratios
  (1.4×, 54×, 902×) and the observation that the central file's geometric standard
  deviations — 1e33 to 1e63, against 1.0513 — identify the out-of-scope values
  mechanically. Verify the note and `design.md` D3 state the same five ratios.

- [ ] 5.6 Run the full suite and record the count. Verify `sbatch
  test/run_tests.slurm` passes with no failures, then update `docs/handoff.md`
  with what changed, the job id, the measured full-cycle wall-clock, and what wave
  2 inherits: that ATP and ADP deliberately carry no edge from this side, and that
  the shared PGK/PYK enzyme nominals disappear once live protein counts arrive.
