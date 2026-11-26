# Table 1's evidence, cell by cell

Read 2026-11-06, by a Claude Code agent from the primary sources below. Each cell's verdict is
that reading. A person in the field should spot-check them before submission.

"0 hits" means a sweep of the full text for `bayes|posterior|prior distribution|credible|likelihood`
found nothing.

| Framework | Priors | Provenance in model | Inference on dynamics | Declared reductions |
|---|---|---|---|---|
| **Karr 2012** | **No.** `Parameter.m` has `defaultValue`, `units`. 0 hits in the paper and WholeCellKB. | **Partial, a separate KB.** `Parameter` inherits `references` (`KnowledgeBaseObject.m`), and the simulation reads the KB through `initializeConstants`. NAR: "we tracked the species, experimental conditions and citation of each piece of evidence". | **Heuristic fit, no posterior.** "refined the values of the sub-model parameters to make the sub-models mutually consistent"; `FitConstants.m`: "use the MATLAB fmincon routine" (the code was not re-checked to say whether fmincon is the abandoned formal attempt or a step inside the heuristic). Data S1 §1.3 (p. 3): the constrained optimisation found no feasible solution, so "we ... developed a heuristic procedure that uses the constraints to calculate a consistent set of parameter values from a subset of the parameters". Fit targets are aggregate observables (doubling time, composition, RNA decay and expression; §A.1, p. 110). No uncertainty on any fitted parameter. | **Partial, in prose.** There is no whole-model reference, but some submodels state how they simplify the published models they adapt: FtsZ "simplified to not include annealing and cyclization" (§3.7, p. 50); cytokinesis filament bending is "modified from a model proposed by Li et al." (§3.3; a stated change to a published model, not stated as a simplification); the rRNA cleavage scheme "adapted and simplified" (§3.23, p. 89). Changed from "Not applicable" in the PR #202 review, to match the syn3A row, where file notes count as Partial. Nearest objects: `fittedConstantNames` (`Process.m`), `experimentallyConstrained` (`Parameter.m`). |
| **Macklin 2020** | **No.** `parameters.tsv` holds single values; the five-fold perturbation is a sensitivity analysis. 0 hits. | **Partial, annotated flat files.** `parameters.tsv` has `source` and `comments`: `cellDensity` cites "Baldwin et al."; `KcatExoRNase` is blank. | **Fit, no posterior.** "iterative parameter estimation approach based on ODEs"; "weighting factor, which we set to optimize a trade-off". 0 hits. | **Partial.** At 899df45: hard-coded dicts, `RNA_EXPRESSION_ADJUSTMENTS` in `fit_sim_data_1.py`. `*_removed.tsv` and `REMOVED_DATA` came later (47e1864, 2022) and are not the paper's version. |
| **Vivarium 2022** | **No.** "At present, there are no Vivarium processes to represent model uncertainty" (this concerns mechanism and initial state, not parameters specifically). | **No.** "Processes are configured by parameters"; plain dicts in `process.py`. The downstream `vivarium-ecoli` has `Parameter(value, source, ..., note)`, which is not credited to Vivarium here. | **No.** No fitting described. Core files: 0 hits. | **No.** `schema_override` and parameter overrides are configuration. |
| **Thornburg 2022** | **Partial, metabolic block.** The SBtab has `!UnconstrainedStd` and a `Parameter prior` table. The simulation uses "the mean values used as parameters in kinetic simulations". | **Partial, flat files.** The pre-balancing TSV has e.g. `R_PGI ... 441.96 ... Cocks-Finch (1985)`. The balanced `Quantity` table has only `!Value`, `!Unit`, `!Type`. | **Balancing, one block, literature data.** "Initial estimates... determined using Parameter Balancing... using Bayesian estimation and prior distributions derived from user-inputted surveys of kinetic data". Later: "kinetic parameters were adjusted such that the differences in the free energies... were minimized". STAR Methods: balancing "calculates standard deviations for each estimated parameter", but Table S2 has values and units only, with no SD column. Refinement also compared time-averaged fluxes with FBA (model against model, Table S3). One hand-tune outside metabolism: "Rather than the 33% calculated, 16% of RNAP was selected to be active" to match the observed doubling time. No fit to time-series data. | **Partial.** `PPP_patch.py`: "add patch changes"; `iMB155_NoH2O.xml` "modified from Breuer et al". Code and file notes, not enumerated objects. |
| **Thornburg 2026** | **No.** `kinetic_params.xlsx` has one `Value` column on every parameter sheet. | **Largely no.** Only the `Membrane Surface Area` sheet has `Source`, `Reference`. | **Manual, no posterior.** "adjustment of the ratio of the binding rates"; metabolism "mostly unchanged from our previous model". | **Partial at most.** `parameter_comparison.py` compares old and new parameter files. |

Table 1's one row for syn3A combines the last two.

## Supplements (read 2026-11-26, KIM-168)

Read by Claude Code agents; the author should spot-check them.

- **Karr 2012: read.** Data S1 (122 pp.) and the three supplementary tables, from the Elsevier CDN
  (`1-s2.0-S0092867412007763-mmc1…6`). All four verdicts hold. The inference cell now says the fit is
  heuristic, as Data S1 §1.3 says, and the reductions cell notes submodel-level simplifications.
- **Thornburg 2022: read in part.** STAR Methods (PMC9985924) and Tables S1–S4 plus Data S1 from the
  Elsevier CDN (`1-s2.0-S0092867421014884-mmc1…5`). Table S2 (kinetic parameters) has values and units only;
  sources appear only in Table S1's proteomics headers and the Membrane Surface Area sheet. All four
  verdicts hold. Document S1 (Figures S1–S7) could not be downloaded.
- **Macklin 2020: not obtained.** PMC served a bot check, Science returned 403, and Europe PMC lists the
  article as not open access. The main text (Covert lab PDF) was re-read; nothing in it changes a verdict.
- **Thornburg 2026: not obtained.** The published paper is Cell 189(9):2582–2597.e27, doi
  10.1016/j.cell.2026.02.009, with an erratum (10.1016/j.cell.2026.04.003). Cell, ScienceDirect and the
  Elsevier API refused access, and there is no PMC copy. Its cells still rest on the preprint and the
  repository; the published abstract differs from the preprint's only in copy-editing.

## Limits of this reading

- **Thornburg 2026:** this reads bioRxiv 10.1101/2025.06.10.658899, which may differ from the published
  text, and the repository Minimal_Cell_4DWCM. Its supplement (Table S2) and erratum are unread.
- **Supplements still unread:** Macklin's and Vivarium's, Thornburg 2022's Document S1, and Thornburg
  2026's. Macklin's holds its fitting method, so it limits that row's inference cell.
- **WholeCellKB's live database** was not queried; only `public/models.py` was read.
- **Macklin's commit:** 899df45 is the last commit before publication, not a confirmed paper tag.
- **odecell's** `paropt` offers Differential Evolution; Thornburg 2022's main text does not say it was used.

## Sources

- Papers, as PMC copies: PMC3413483, PMC3531061, PMC7990026, PMC8963310, PMC9985924.
- github.com/CovertLab/WholeCell at 6cdee6b; github.com/CovertLab/WholeCellKB at master;
  github.com/CovertLab/WholeCellEcoliRelease at 899df45 and 47e1864;
  github.com/vivarium-collective/vivarium-core at 9f534d6;
  github.com/Luthey-Schulten-Lab/Minimal_Cell at db048ac;
  github.com/Luthey-Schulten-Lab/Minimal_Cell_4DWCM at main, pushed 2026-09-30.
