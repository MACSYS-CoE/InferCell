# 4DWCM: Simplified Pipeline

```
┌─────────────────────────────────────────────────────────────────────┐
│                        INITIALISATION                               │
├─────────────────────────────────────────────────────────────────────┤
│                                                                     │
│  1. PARSE GENOME                                                    │
│     syn3A.gb (GenBank) ──► 493 genes                                │
│     For each gene, extract:                                         │
│       - RNA sequence (all genes)                                    │
│       - AA sequence (protein-coding genes)                          │
│       - promoter strength (from kinetic_params.xlsx)                │
│     Register RDME species per gene:                                 │
│       mRNA (R_XXXX), protein (P_XXXX), ribosome-bound states,      │
│       degradation intermediates, membrane precursors                │
│                                                                     │
│  2. BUILD 3D LATTICE                                                │
│     64×64×128 grid, 10 nm spacing ──► ~400 nm diameter sphere       │
│     Assign site types from geometry:                                │
│       - extracellular (outside sphere)                              │
│       - membrane (shell)                                            │
│       - outer_cytoplasm (sub-membrane layer)                        │
│       - cytoplasm (interior)                                        │
│       - DNA (sites occupied by chromosome beads)                    │
│     Source: RegionsAndComplexes.py + oneParamMulder-local_min.json   │
│                                                                     │
│  3. POPULATE LATTICE                                                │
│     Place initial particle counts from initial_concentrations.xlsx   │
│     Generate initial chromosome config via sc_chain_generation       │
│     Map chromosome beads → DNA-type lattice sites                   │
│     Place ribosomes, RNAPs, metabolites in appropriate regions      │
│                                                                     │
│  4. REGISTER REACTIONS                                              │
│     RDME reactions: transcription, translation, degradation,        │
│       ribosome biogenesis, RNAP assembly, SecY translocation        │
│     Each reaction is region-specific (e.g. transcription only       │
│       at DNA sites, translation in cytoplasm)                       │
│     Diffusion rules per species (size-dependent rates)              │
│                                                                     │
└──────────────────────────────┬──────────────────────────────────────┘
                               │
                               ▼
┌─────────────────────────────────────────────────────────────────────┐
│                        MAIN SIMULATION LOOP                         │
│                     (≈ 105 min of biology)                          │
├─────────────────────────────────────────────────────────────────────┤
│                                                                     │
│  5. RDME (Lattice Microbes, GPU 1)              Δt = 50 μs         │
│     Stochastic reaction-diffusion on the lattice.                   │
│     Runs continuously; all other steps interrupt it via Hook.       │
│                                                                     │
│       every 12.5 ms (Hook fires) ──┐                                │
│                                    ▼                                │
│  6. RIBOSOME POSITIONS            every hook                        │
│     Move ribosome centre-of-mass particles.                         │
│     Update excluded volume (ribos span multiple sites).             │
│     Update polysome structures every 8th call.                      │
│                                                                     │
│       every 1.0 s ──┐                                               │
│                     ▼                                               │
│  7. CME (Lattice Microbes, subprocess)          Δt_bio = 1 s       │
│     Well-stirred stochastic sim for:                                │
│       - tRNA charging (20 tRNA species)                             │
│       - transcription NTP consumption                               │
│     Runs as separate process, writes .lm file to disk.              │
│     Results read back and applied to RDME lattice.                  │
│                                                                     │
│       every 1.0 s (after CME) ──┐                                   │
│                                 ▼                                   │
│  8. ODE METABOLISM (odecell + LSODA)            Δt_bio = 1 s       │
│     Deterministic integration of full metabolic network:            │
│       - glycolysis, PPP, nucleotide/lipid/AA synthesis              │
│       - enzyme concentrations from current RDME counts              │
│       - monomer demands from translation/transcription as sinks     │
│     Integrates 1s, writes updated metabolite counts back.           │
│                                                                     │
│       every 4.0 s ──┐                                               │
│                     ▼                                               │
│  9. CHROMOSOME BD (btree_chromo, LAMMPS, GPU 2) Δt_bio = 4 s      │
│     Brownian dynamics of coarse-grained bead-spring DNA polymer.    │
│     Includes SMC loop extrusion, topoisomerase activity.            │
│     Tracks replication fork progression (if initiated).             │
│     Runs async on 2nd GPU; result read at next cycle.               │
│     New bead positions mapped back onto RDME lattice.               │
│                                                                     │
│       every 4.0 s (if SA threshold met) ──┐                         │
│                                           ▼                         │
│  10. MORPHOLOGY UPDATE                                              │
│      Pre-division: inflate sphere (add lattice sites to membrane)   │
│      Division: form septum, partition to daughter cells              │
│      Triggered by accumulated lipid + membrane protein synthesis    │
│                                                                     │
└──────────────────────────────┬──────────────────────────────────────┘
                               │
                               ▼
┌─────────────────────────────────────────────────────────────────────┐
│                           OUTPUTS                                   │
│                     (saved every ~1 s)                              │
├─────────────────────────────────────────────────────────────────────┤
│                                                                     │
│  11. PER-SPECIES COUNTS                                             │
│      Copy number of every species at each save point.               │
│      Covers all 493 gene products + metabolites + intermediates.    │
│                                                                     │
│  12. METABOLIC FLUXES                                               │
│      Reaction rates from the ODE integrator at each save point.     │
│      Stored in working_directory/fluxes/                            │
│                                                                     │
│  13. LATTICE SNAPSHOTS                                              │
│      Full 3D particle lattice + site lattice (region assignments).  │
│      Enables spatial analysis: where are ribosomes, mRNAs, etc.     │
│                                                                     │
│  14. CHROMOSOME CONFIGURATIONS                                      │
│      3D bead positions at each DNA update step.                     │
│      Stored in working_directory/DNA/                               │
│                                                                     │
│  15. RESTART FILES                                                  │
│      Full sim_properties dict + lattice state.                      │
│      Allows resuming via Restart_Whole_Cell_Minimal_Cell.py         │
│                                                                     │
│  16. CME TRAJECTORIES                                               │
│      Per-second .lm files from each global CME run.                 │
│      Stored in working_directory/CME/                               │
│                                                                     │
└─────────────────────────────────────────────────────────────────────┘
```

## One-line summary

```
GenBank ──► species + reactions on 3D lattice ──► RDME (continuous)
  interrupted by Hook every 12.5ms to run { ribosomes | CME | ODE | BD | growth }
  ──► time-series of counts, fluxes, 3D snapshots over full cell cycle
```