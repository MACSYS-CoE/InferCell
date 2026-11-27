# KIM-195 `threaded`

- job 21785047 on dave19, AMD EPYC 7543 32-Core Processor; Julia 1.11.9; LinearSolve 3.87.0
- threads 4, BLAS threads 1, GC threads 5
- window 0–3150.0 s in minute segments, saves at the whole minutes; Rodas5P, KLUFactorization (primal) and DualSparseKLU (dual, unless a row says the library path); abstol 1.0e-10, reltol 1.0e-8; chunk 10

primal solve (GE B = 10 at θ = 0): median 9.904 s; steps 6415 (rejected 416), f 51426, J 5999, W 6415, solves 51320

| placement | B | passes | gradient median (s) | repeats (s) | per pass median (s) | pass 1 stats |
|---|---|---|---|---|---|---|
| GE | 20 | 2 | 77.002 | 76.723, 77.281 | 76.630 | steps 5325 (rejected 465), f 42706, J 4860, W 5325, solves 42600 |
| GE | 40 | 4 | 92.412 | 93.288, 91.535 | 90.519 | steps 5326 (rejected 466), f 42714, J 4860, W 5326, solves 42608 |
| MET | 20 | 2 | 89.569 | 88.518, 90.619 | 87.861 | steps 6348 (rejected 446), f 50890, J 5902, W 6348, solves 50784 |
| MET | 40 | 4 | 108.474 | 107.977, 108.972 | 99.243 | steps 6460 (rejected 442), f 51786, J 6018, W 6460, solves 51680 |
| BOTH | 20 | 2 | 85.904 | 85.711, 86.097 | 79.657 | steps 5155 (rejected 460), f 41346, J 4695, W 5155, solves 41240 |
| BOTH | 40 | 4 | 104.102 | 103.187, 105.018 | 94.148 | steps 5228 (rejected 454), f 41930, J 4774, W 5228, solves 41824 |
