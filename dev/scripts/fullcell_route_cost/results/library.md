# KIM-195 `library`

- job 21785046 on dave9, AMD EPYC 7543 32-Core Processor; Julia 1.11.9; LinearSolve 3.87.0
- threads 1, BLAS threads 1, GC threads 1
- window 0–3150.0 s in minute segments, saves at the whole minutes; Rodas5P, KLUFactorization (primal) and DualSparseKLU (dual, unless a row says the library path); abstol 1.0e-10, reltol 1.0e-8; chunk 10

| placement | LinearSolve dual path, one pass (s) | its stats | DualSparseKLU, median (s) | its stats | partials rel |
|---|---|---|---|---|---|
| GE | 83.3 | steps 5325 (rejected 465), f 42706, J 4860, W 5325, solves 42600 | 69.135 | steps 5325 (rejected 465), f 42706, J 4860, W 5325, solves 42600 | 1.41e-08 |
| MET | 103.6 | steps 6348 (rejected 446), f 50890, J 5902, W 6348, solves 50784 | 83.384 | steps 6348 (rejected 446), f 50890, J 5902, W 6348, solves 50784 | 1.23e-12 |
