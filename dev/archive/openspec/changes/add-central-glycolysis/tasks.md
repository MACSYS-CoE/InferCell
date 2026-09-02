## 1. The vendored parameter extract

- [ ] 1.1 Write `dev/scripts/extract_central_glycolysis.jl`, taking a path to a
  `Minimal_Cell` checkout, and run it to produce
  `src/organisms/coreA/data/central_glycolysis.tsv`. Verify the output has exactly
  65 data rows — 10 `kcatF_R_*`, 10 `kcatR_R_*`, 32 `km_R_*_M_*`, 13
  `conc_M_*` — and that four spot values match the upstream file byte for byte:
  `kcatF_R_PGI = 804.3384` (gstd 1.0513), `kcatF_R_FBA = 59.7` (gstd 1.7466, the
  loosest forward constant in the core), `km_R_PGI_M_g6p_c = 22.9419`, and
  `conc_M_atp_c = 3.6529` (gstd 1.2825).

- [ ] 1.2 Write `src/organisms/coreA/data/README.md` recording the upstream
  filename `Central_AA_Zane_Balanced_direction_fixed_nounqATP.tsv`, the commit
  `db048ac`, which SBtab tables each row class came from, the fact that the
  reshape renames identifiers and changes no value, and the exact command that
  regenerates the extract. Verify by re-running that command from the README and
  getting an unchanged file.

## 2. The sub-model, its rate law and its parameters

- [ ] 2.1 Create `src/organisms/coreA/central_glycolysis.jl` with the
  `CentralGlycolysis <: AbstractSubModel` type, `formalism = :ode`, and
  `states` returning the 11 `:glycolytic` plus 2 `:redox` registry species. Wire
  one `include` into `src/InferCell.jl` after the registry and one into
  `test/runtests.jl` for `test/test_corea_central_glycolysis.jl`. Verify
  `using InferCell` loads clean and a test asserting `length(states(m)) == 13`,
  that the set equals `species_in_group(:glycolytic) ∪ species_in_group(:redox)`,
  and that `species_index.(states(m))` is strictly increasing, passes under
  `sbatch test/run_tests.slurm`.

- [ ] 2.2 Implement the modular rate law of design.md D1, generic over substrate
  and product counts, and the ten reactions' stoichiometry. Verify unit tests: a
  reaction whose forward and reverse catalytic terms are equal has zero net rate;
  raising products far above substrates makes the net rate negative; and each
  reaction reports the substrate/product term counts of D1 (PGI 1,1 … GAPD 3,2 …
  LDH_L 2,2), summing to 32 Michaelis constants.

- [ ] 2.3 Import all 65 values through `load_parameter` against the vendored
  extract, with `LogNormal(log(mode), log(gstd))` priors and `fixed = true` by
  default plus a `free` constructor keyword. Verify every returned
  `InferParameter` reports `"central_balanced"` as its source and its upstream
  identifier; `governing_choices` is empty (one file, no ambiguity yet);
  `M_g3p_c` and `M_lac__L_c` are marked uninformed and keep gstd 10; no parameter
  of this module appears in the composed model's asserted-prior enumeration; and
  `CentralGlycolysis(free = [:kcatF_R_FBA])` frees exactly that one.

- [ ] 2.4 Add the ten nominal enzyme concentrations of design.md D6, derived from
  copy number at 20,180 particles per mM, marked nominal rather than measured, and
  overridable. Verify each equals `copies / 20180` to six decimals (PGI 0.013181 …
  GAPD 0.067146 … LDH_L 0.054509); that overriding one scales exactly the rates
  that enzyme catalyses and no others; and that none of them is the published
  no-GPR default of 0.001 mM.

- [ ] 2.5 Confirm the registry-agreement check actually fires on the 13
  `conc_M_*` rows. Verify by a test that mutates one extract row in a temporary
  copy and asserts `load_parameter` throws naming the species, the imported value
  and the registry's — a passing load must be evidence, not a skipped check.

## 3. Boundary declarations

- [ ] 3.1 Implement `coupling` returning the five `CurrencyEdge`s and four
  `MassEdge`s of design.md D8, every one with `peer` unnamed. Verify a test
  asserting exactly nine edges with the tabulated species, kinds and directions;
  that no edge names `M_nad_c` or `M_nadh_c`; that no deferred-counter,
  rate-constant, volume or clamped edge is present; and that every edge species
  resolves to a registry position.

- [ ] 3.2 Implement `inputs` returning the three energy currencies and the ten
  enzyme protein counts. Verify `resolve_coupling` on the module alone raises no
  inputs/edges disagreement in either direction, and that removing `M_pi_c` from
  `inputs` makes it throw naming that species — the drift check must be live.

- [ ] 3.3 Implement `reduction_notes` carrying four entries: `NOX` dropped with
  its three reasons; the Michaelis-constant column choice naming all eight
  differing constants with both values each; the energy currencies held rather
  than integrated in any composition that does not execute the coupling; and the
  enzyme concentrations being nominal stand-ins for live counts. Verify
  `reduction_declarations` on a composition containing this module returns all
  four, and that the `NOX` and column-choice entries read as sentences naming the
  affected reactions and constants.

- [ ] 3.4 Verify standalone resolution reports rather than fails: a test composing
  the module alone asserts `resolve_coupling` succeeds, lists `M_atp_c`,
  `M_adp_c` and `M_pi_c` among the unowned states, and distinguishes this
  deliberate partial composition from an incomplete one.

## 4. Dynamics and standalone integration

- [ ] 4.1 Implement `dynamics`, assembling `du` for the 13 owned states from the
  ten reaction fluxes, reading `M_atp_c`, `M_adp_c`, `M_pi_c` and the ten enzyme
  counts from `u_inputs`. Verify a hand-checked unit test at one state vector:
  each species' derivative equals the signed sum of the fluxes of the reactions
  the D8 table says touch it, and no derivative is returned for a species outside
  the thirteen.

- [ ] 4.2 Add `HeldEnergyPool <: AbstractSubModel` to
  `test/test_corea_central_glycolysis.jl`, owning `M_atp_c`, `M_adp_c` and
  `M_pi_c` with zero derivatives at their registry values. Verify composing it
  with `CentralGlycolysis` builds a problem and that every input resolves.

- [ ] 4.3 Integrate the two-module composition from the registry's initial
  conditions at published parameters, with a stiff solver and explicit
  `abstol`/`reltol`. Verify a trajectory is produced over the requested interval
  and that all thirteen owned states stay non-negative at every saved point.

## 5. Conservation, and the record

- [ ] 5.1 **Redox balance — this module's own conservation check.** Assert
  `M_nad_c + M_nadh_c` equals its initial value at every saved point of the 4.3
  trajectory, within a bound derived from the solver tolerances rather than a
  hand-picked constant. Verify the assertion passes, and that tightening the
  tolerances tightens the bound the test enforces.

- [ ] 5.2 Show the redox check can fail. Verify a mutation test that perturbs
  `GAPD`'s NAD⁺ stoichiometry so it no longer mirrors `LDH_L`, and asserts the 5.1
  check fails naming the conserved pool and the drift — a conservation test that
  cannot fail is not evidence.

- [ ] 5.3 Add the correction-of-record row for the Michaelis-constant column to
  the inline table in `dev/notes/reduced-syn3a-scoping.md`, replacing the note's
  "~35 K_M from the balanced file (with priors)" with what the source actually
  does: the balanced column carries the priors, the `Quantity` column is what the
  published simulator reads, and the two disagree on 8 of these 32. Verify the
  note and `proposal.md` state the same eight identifiers and values.

- [ ] 5.4 Run the full suite and record the count. Verify `sbatch
  test/run_tests.slurm` passes with no failures and the total exceeds the
  pre-change 829 by the number of tests this change adds, then update
  `docs/handoff.md` with what changed, the job id, and what wave 2 inherits —
  in particular that the outbound currency edges are declared but unexecuted
  because `_build_rhs` lets a module write only its own states.
