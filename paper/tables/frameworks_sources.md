# Table 1's evidence, cell by cell

Read 2026-11-06, by a Claude Code agent from the primary sources below. Each cell's verdict is
that reading. A person in the field should spot-check them before submission.

"0 hits" means a sweep of the full text for `bayes|posterior|prior distribution|credible|likelihood`
found nothing.

| Framework | Priors | Provenance in model | Inference on dynamics | Declared reductions |
|---|---|---|---|---|
| **Karr 2012** | **No.** `Parameter.m` has `defaultValue`, `units`. 0 hits in the paper and WholeCellKB. | **Partial, a separate KB.** `Parameter` inherits `references` (`KnowledgeBaseObject.m`), and the simulation reads the KB through `initializeConstants`. NAR: "we tracked the species, experimental conditions and citation of each piece of evidence". | **Fit, no posterior.** "refined the values of the sub-model parameters to make the sub-models mutually consistent"; `FitConstants.m`: "use the MATLAB fmincon routine". | **None stated.** There is no reference model to depart from. Nearest: `fittedConstantNames` (`Process.m`), `experimentallyConstrained` (`Parameter.m`). |
| **Macklin 2020** | **No.** `parameters.tsv` holds single values; the five-fold perturbation is a sensitivity analysis. 0 hits. | **Partial, annotated flat files.** `parameters.tsv` has `source` and `comments`: `cellDensity` cites "Baldwin et al."; `KcatExoRNase` is blank. | **Fit, no posterior.** "iterative parameter estimation approach based on ODEs"; "weighting factor, which we set to optimize a trade-off". 0 hits. | **Partial.** At 899df45: hard-coded dicts, `RNA_EXPRESSION_ADJUSTMENTS` in `fit_sim_data_1.py`. `*_removed.tsv` and `REMOVED_DATA` came later (47e1864, 2022) and are not the paper's version. |
| **Vivarium 2022** | **No.** "At present, there are no Vivarium processes to represent model uncertainty" (this concerns mechanism and initial state, not parameters specifically). | **No.** "Processes are configured by parameters"; plain dicts in `process.py`. The downstream `vivarium-ecoli` has `Parameter(value, source, ..., note)`, which is not credited to Vivarium here. | **No.** No fitting described. Core files: 0 hits. | **No.** `schema_override` and parameter overrides are configuration. |
| **Thornburg 2022** | **Partial, metabolic block.** The SBtab has `!UnconstrainedStd` and a `Parameter prior` table. The simulation uses "the mean values used as parameters in kinetic simulations". | **Partial, flat files.** The pre-balancing TSV has e.g. `R_PGI ... 441.96 ... Cocks-Finch (1985)`. The balanced `Quantity` table has only `!Value`, `!Unit`, `!Type`. | **Balancing, one block, literature data.** "Initial estimates... determined using Parameter Balancing... using Bayesian estimation and prior distributions derived from user-inputted surveys of kinetic data". Later: "kinetic parameters were adjusted such that the differences in the free energies... were minimized". | **Partial.** `PPP_patch.py`: "add patch changes"; `iMB155_NoH2O.xml` "modified from Breuer et al". Code and file notes, not enumerated objects. |
| **Thornburg 2026** | **No.** `kinetic_params.xlsx` has one `Value` column on every parameter sheet. | **Largely no.** Only the `Membrane Surface Area` sheet has `Source`, `Reference`. | **Manual, no posterior.** "adjustment of the ratio of the binding rates"; metabolism "mostly unchanged from our previous model". | **Partial at most.** `parameter_comparison.py` compares old and new parameter files. |

Table 1's one row for syn3A combines the last two.

## Limits of this reading

- **Thornburg 2026:** the Cell text was out of reach (403). This reads bioRxiv 10.1101/2025.06.10.658899,
  which may differ, and the repository Minimal_Cell_4DWCM.
- **Supplements, unread:** Karr Data S1 and tables, Macklin's, Thornburg 2022 Tables S1–S3, Vivarium's,
  and Thornburg 2026's Table S2. Karr's and Macklin's hold their fitting methods, so these limit their
  provenance and inference cells.
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
