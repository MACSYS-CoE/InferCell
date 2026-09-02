## 1. The per-gene sequence extract

- [ ] 1.1 Write `dev/scripts/extract_transcription_genes.jl`, taking a path to a
  `Minimal_Cell` checkout, and run it to produce
  `src/organisms/coreA/data/transcription_genes.tsv`: one row per locus with
  transcript length, the four base counts, protein copy number and measured mean
  transcript count. Verify 17 rows; that the four base counts sum to the length for
  every gene; and four spot rows — `JCVISYN3A_0445` (1284; A 521, C 125, G 197,
  U 441; 266 copies; 0.4403 mRNA), `JCVISYN3A_0607` (1017; 1355 copies; 2.1781),
  `JCVISYN3A_0779` (2238; 831 copies), `JCVISYN3A_0694` (270; 290 copies).

- [ ] 1.2 Verify the generator fails loudly rather than defaulting: a locus absent
  from the genome record, the proteomics table or the mRNA counts must abort with
  its name, since a silently missing gene would surface only as a model with
  sixteen transcripts.

- [ ] 1.3 Append the three upstream sources (`syn3A.gb`, `proteomics.xlsx`,
  `mRNA_counts.csv`), commit `db048ac`, the complement-strand handling, and the
  regeneration command to `src/organisms/coreA/data/README.md`, creating it if no
  sibling change has. Verify re-running that command leaves the extract unchanged.

## 2. The sub-model and its reactions

- [ ] 2.1 Create `src/organisms/coreA/transcription.jl` with
  `CoreATranscription <: AbstractSubModel`, `formalism = :jump`,
  `inference_mode = :simulation`, and `states` returning the 17 transcripts and 5
  counters. Wire one `include` into `src/InferCell.jl` after the registry and one
  into `test/runtests.jl`. Verify `using InferCell` loads clean, that
  `length(states(m)) == 22`, that no state is a registry species
  (`!any(is_registered, states(m))`), and that the suite passes under `sbatch
  test/run_tests.slurm`.

- [ ] 2.2 Implement the 17 reactions per design.md D1 as constant-rate jumps, each
  catalytic in its gene and producing one transcript plus the five counter
  increments. Verify a unit test that firing one gene's reaction leaves its gene
  quantity unchanged, raises its transcript by one, raises `ATP_trsc` by that
  gene's length and each monomer counter by that gene's base count, and touches no
  other gene's transcript.

- [ ] 2.3 Carry gene copy number as a fixed per-gene quantity rather than a state,
  per design.md D1, and record that the published model carries genes as species.
  Verify no gene appears in `states(m)`, that every gene quantity is 1, and that
  the composed model returns the recorded note.

- [ ] 2.4 Implement seeded initialisation drawing each transcript's initial count
  from its measured mean using the published distribution. Verify two runs with one
  seed are identical, two runs with different seeds may differ, and the per-gene
  mean over many draws is within sampling error of the recorded measured value.

- [ ] 2.5 Verify the gene set matches the metabolic modules: a test asserting this
  module's 17 loci equal the union of the loci `add-central-glycolysis`,
  `add-pts-transport` and `add-nucleotide-recycling` declare, and that a mismatch
  reports the loci on each side. Skip with a recorded reason if none of the three
  has landed.

## 3. Rate constants

- [ ] 3.1 Implement the rate-constant formula of design.md D2 with
  `base_mapping = :corrected` as the default and `:published` selectable. Verify
  each base's count enters only its own NTP term under `:corrected`; that
  `:published` reproduces the C→UTP, G→CTP, U→GTP permutation; and that the
  selection is recorded on the model so a trajectory is identifiable.

- [ ] 3.2 Take all four NTP concentrations from balanced values per design.md D3 —
  ATP 3.6529, GTP 1.6627, CTP 0.6874, UTP 2.7681 — each with its gstd and source
  file. Verify none equals the `MinCell_CMEODE.py` setup constants (1.04, 0.68,
  0.34, 0.68), and that ATP and GTP agree with the registry's values for the pools
  `add-nucleotide-recycling` owns.

- [ ] 3.3 Implement `recompute_rate_constants!` and `rate_constants` per design.md
  D4. Verify recomputation with different NTP concentrations updates all 17 and
  changes nothing else about the sub-model; that simulating past the declared 60 s
  cadence leaves the constants unchanged, because this module performs no rebuild;
  and that `rate_constants` returns all 17 keyed by locus.

- [ ] 3.4 Add the `/180` promoter proxy per design.md D5. Verify the 17 values
  match the table there to three decimals (PGI 1.478, GAPD 7.528, ptsG 4.617, GK1
  1.033); that two genes' promoter strengths are in the ratio of their protein
  counts; and that the copy numbers agree with those the metabolic modules derive
  enzyme concentrations from.

- [ ] 3.5 Verify the turnover cap is non-binding: assert every gene's
  promoter-scaled turnover is below the published ceiling of 180 nt/s, the largest
  being GAPD's ~8.85, and that the sub-model reports this rather than assuming it.

## 4. Costs, coupling and the caveats

- [ ] 4.1 Implement `coupling` returning the eleven edges of design.md D7 — four
  rate-constant edges at 60 s piecewise-constant, five deferred counters with the
  published clamped policy, two clamped edges at `origin = :ours` holding CTP at
  0.6874 and UTP at 2.7681 mM with provenance. Verify exactly eleven edges with the
  tabulated species, kinds and directions; that `inputs(m)` is empty; and that three
  kinds coexisting on `M_ctp_c` inbound resolves rather than conflicting.

- [ ] 4.2 Verify the drain semantics are recorded per design.md D6: each counter
  reports the registry species it debits and what that drain produces — `ATP_trsc`
  yielding ADP and Pi, the four monomer counters yielding PPi. Verify the composed
  model can be queried for this and returns all five.

- [ ] 4.3 Implement `reduction_notes` with five entries: the corrected base mapping,
  naming both mappings and the ~1.9x GTP-sensitivity effect; the CTP and UTP
  chemostats; the promoter-strength proxy; the proteomics circularity, naming the 17
  affected parameters; and the ATP double charge, stated as the published model's
  accounting rather than an error. Verify `reduction_declarations` returns all five,
  and — the negative case — that the five deferred counters do **not** appear,
  because they take the published clamped policy.

- [ ] 4.4 Implement the constitutive-motif record per the reaction spec. Verify
  exactly one rate constant governs each gene, that no promoter switching rate
  exists, and that the composed model returns the recorded reason the two-state
  identifiability result does not transfer.

- [ ] 4.5 Verify the chemostat exemption is recorded and the live pools are not
  exempt: resolving this module alone reports `M_atp_c` and `M_gtp_c` as dead ends
  naming this module as consumer and identifying their moieties, while `M_ctp_c`
  and `M_utp_c` are exempt with the exemption recorded.

- [ ] 4.6 Verify the gradient obstruction fires at the join and not here: composing
  this module alone reports none, because its inference mode is simulation-based;
  composing it with a differentiable sub-model reports the five clamped counters
  and the species they clip, as a diagnosable warning.

## 5. Simulation, calibration and the record

- [ ] 5.1 Simulate the module alone over a full cell cycle. Verify a trajectory is
  produced, every transcript count is a non-negative integer throughout, and the
  cost counters accumulate monotonically because no hook debits them.

- [ ] 5.2 **Calibration check against measured data.** Verify each gene's
  time-averaged transcript count over the trajectory lands within a factor of two of
  the measured mean in the extract, and that the computed rate constants fall in
  1.26e-3 to 8.29e-3 s⁻¹ per design.md D9. This is the check that the formula, the
  extract and the promoter proxy are wired together correctly.

- [ ] 5.3 **Measure channel 4's gain.** Implement an elasticity diagnostic
  reporting `∂ln k_g / ∂ln[NTP]` per gene, and verify it reproduces design.md D9:
  0.044–0.051 for all four pools together, 0.0079–0.0117 for GTP under
  `:corrected`, and 0.0156–0.0196 under `:published`. Record the result where wave 3
  will find it.

- [ ] 5.4 Verify undebited counters are labelled: a trajectory from a composition
  with no hook records that the five cost counters accrued without being debited,
  and names them.

- [ ] 5.5 Add to `dev/notes/reduced-syn3a-scoping.md` the base-mapping bug with its
  code locations and effect, the channel-4 elasticity measurement — which bears on
  the note's argument that channel 4 satisfies design criterion 2 — and the fact
  that transcription is a second source of pyrophosphate, which the note's PPA
  argument attributes to charging alone. Verify the note and `design.md` D9 state
  the same elasticity ranges.

- [ ] 5.6 Run the full suite and record the count. Verify `sbatch
  test/run_tests.slurm` passes with no failures, then update `docs/handoff.md` with
  what changed, the job id, and what the later waves inherit: that channel 4's gain
  is a few percent, that transcription contributes pyrophosphate to
  `add-nucleotide-recycling`'s phosphate closure, and that the proteomics
  circularity constrains wave 3's observable design.
