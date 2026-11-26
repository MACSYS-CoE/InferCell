# KIM-184 `check`

- job 21666251 on dave5, AMD EPYC 7543 32-Core Processor; Julia 1.10.5
- threads 1, BLAS threads 1, GC threads 1 (sweep 0)
- window 0–3150.0 s, 52 saves; Rodas5P + KLUFactorization, abstol 1.0e-10, reltol 1.0e-8; chunk 10

## The system

| quantity | value |
|---|---|
| CME species in the table | 2661 |
| states (species a first-order reaction reads or changes) | 2112 |
| first-order reactions (LM type 1) | 2002 |
| second-order reactions (type 2), left out | 363 |
| other types | 0 |
| Jacobian non-zeros from the table | 20044 |
| … with the structural diagonal | 21095 |
| reactions with flux over the window | 1511 |

## Checks

- RHS at x0 vs the table's propensities summed in BigFloat: 677 non-zero entries, max relative error 2.642e-15 (bound 1e-10); 1435 entries exactly zero in both (0 mismatches). **PASS**
- J·x0 vs the RHS: max relative difference 1.527e-16; zeros kept: true
- Dual linear solve (the solver's linear cache under a dual state): LinearCache / LinearSolveFunction / W SparseMatrixCSC{Dual}
- DualSparseKLU vs LinearSolve's DualLinearCache, one pass of 10 partials solved to 30 s: 0.353 s (82 steps) against 297.110 s (82 steps); partials differ by 0.000e+00 of their largest, values by 0.000e+00 relative
- AV3, B = 20, abstol 1e-14, reltol 1e-12, h = 1e-04:

| θ | reaction | dual | central FD | relative |
|---|---|---|---|---|
| θ_1 | 1164 | -2.607077870429e+04 | -2.607078230273e+04 | 1.380e-07 |
| θ_10 | 699 | -9.774521756389e+03 | -9.774522260676e+03 | 5.159e-08 |
| θ_20 | 25 | -1.601730910017e+01 | -1.601731171263e+01 | 1.631e-07 |

  **PASS** (bound 1.0e-6)
- Pass-wise gradient vs `ForwardDiff.gradient` with Chunk{10}, timing tolerances: max |Δ| = 0.000e+00

## The block (ranked by flux integrated over the window at θ = 0)

| rank | reaction | reactant | k⁰ (1/s) | ∫ flux dt |
|---|---|---|---|---|
| 1 | 1164 | M_RNA_JCVISYN3A_0542_c | 0.0701226 | 505.28 |
| 2 | 1206 | M_RNA_JCVISYN3A_0607_c | 0.0700633 | 459.715 |
| 3 | 1065 | M_RNA_JCVISYN3A_0475_c | 0.0744698 | 401.165 |
| 4 | 717 | M_RNA_JCVISYN3A_0213_c | 0.0721994 | 322.073 |
| 5 | 501 | M_RNA_JCVISYN3A_0910_c | 0.121098 | 305.609 |
| 6 | 696 | M_RNA_JCVISYN3A_0150_c | 0.0687414 | 271.734 |
| 7 | 615 | M_RNA_JCVISYN3A_0076_c | 0.0717223 | 264.388 |
| 8 | 450 | M_RNA_JCVISYN3A_0482_c | 0.0986722 | 255.324 |
| 9 | 666 | M_RNA_JCVISYN3A_0131_c | 0.0697684 | 253.015 |
| 10 | 699 | M_RNA_JCVISYN3A_0151_c | 0.0749411 | 234.616 |
| 11 | 1149 | M_RNA_JCVISYN3A_0536_c | 0.0650587 | 230.783 |
| 12 | 1137 | M_RNA_JCVISYN3A_0528_c | 0.0596509 | 228.491 |
| 13 | 531 | M_RNA_JCVISYN3A_0526_c | 0.0903102 | 227.055 |
| 14 | 537 | M_RNA_JCVISYN3A_0644_c | 0.0746262 | 223.151 |
| 15 | 438 | M_RNA_JCVISYN3A_0658_c | 0.0873492 | 219.802 |
| 16 | 1323 | M_RNA_JCVISYN3A_0747_c | 0.0682068 | 219.395 |
| 17 | 1029 | M_RNA_JCVISYN3A_0437_c | 0.0726423 | 216.153 |
| 18 | 24 | M_RNA_JCVISYN3A_0779_c | 0.0635742 | 214.71 |
| 19 | 534 | M_RNA_JCVISYN3A_0199_c | 0.0845762 | 214.464 |
| 20 | 25 | ptsg_cyto | 0.0670241 | 213.446 |
| 21 | 453 | M_RNA_JCVISYN3A_0660_c | 0.082227 | 212.458 |
| 22 | 498 | M_RNA_JCVISYN3A_0653_c | 0.0816599 | 210.189 |
| 23 | 396 | M_RNA_JCVISYN3A_0238_c | 0.071158 | 209.788 |
| 24 | 1173 | M_RNA_JCVISYN3A_0545_c | 0.0664275 | 207.244 |
| 25 | 528 | M_RNA_JCVISYN3A_0422_c | 0.0819739 | 206.614 |
| 26 | 522 | M_RNA_JCVISYN3A_0833_c | 0.0805489 | 198.99 |
| 27 | 495 | M_RNA_JCVISYN3A_0666_c | 0.0800046 | 198.069 |
| 28 | 411 | M_RNA_JCVISYN3A_0149_c | 0.0763915 | 197.434 |
| 29 | 441 | M_RNA_JCVISYN3A_0665_c | 0.0762276 | 197.409 |
| 30 | 732 | M_RNA_JCVISYN3A_0220_c | 0.0726423 | 194.556 |
| 31 | 417 | M_RNA_JCVISYN3A_0025_c | 0.0710441 | 194.021 |
| 32 | 480 | M_RNA_JCVISYN3A_0638_c | 0.0784152 | 192.441 |
| 33 | 456 | M_RNA_JCVISYN3A_0655_c | 0.0765562 | 189.687 |
| 34 | 510 | M_RNA_JCVISYN3A_0661_c | 0.0727911 | 186.552 |
| 35 | 1074 | M_RNA_JCVISYN3A_0493_c | 0.072521 | 186.535 |
| 36 | 462 | M_RNA_JCVISYN3A_0198_c | 0.0746262 | 184.978 |
| 37 | 519 | M_RNA_JCVISYN3A_0807_c | 0.0717618 | 184.901 |
| 38 | 942 | M_RNA_JCVISYN3A_0394_c | 0.060258 | 183.844 |
| 39 | 423 | M_RNA_JCVISYN3A_0647_c | 0.0733927 | 181.372 |
| 40 | 435 | M_RNA_JCVISYN3A_0654_c | 0.0699254 | 181.187 |

```
Julia Version 1.10.5
Commit 6f3fdf7b362 (2024-08-27 14:19 UTC)
Build Info:
  Official https://julialang.org/ release
Platform Info:
  OS: Linux (x86_64-linux-gnu)
  CPU: 64 × AMD EPYC 7543 32-Core Processor
  WORD_SIZE: 64
  LIBM: libopenlibm
  LLVM: libLLVM-15.0.7 (ORCJIT, znver3)
Threads: 1 default, 0 interactive, 1 GC (on 64 virtual cores)
Environment:
  LD_LIBRARY_PATH = /apps/modules/software/Julia/1.10.5/lib:/apps/modules/software/Julia/1.10.5/lib/julia:/apps/slurm/latest/lib:/apps/slurm/latest/lib/slurm:/opt/nvidia/latest/usr/lib64
  JULIA_PKG_OFFLINE = 1
```
