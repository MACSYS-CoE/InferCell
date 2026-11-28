# KIM-195 `library`

- job 21842558 on dave11, AMD EPYC 7543 32-Core Processor; Julia 1.11.9; LinearSolve 3.87.0
- InferCell 491cb900ed25f6b851af8aa2d12aa7112fab8261; syn3a-ensemble 2075809cbb825140cb31a7bec8bc3529770c8e13 (data/model unmodified)
- threads 1, BLAS threads 1, GC threads 1
- window 0–3150.0 s in minute segments, saves at the whole minutes; Rodas5P, KLUFactorization (primal) and DualSparseKLU (dual, unless a row says the library path); abstol 1.0e-10, reltol 1.0e-8; chunk 10

| placement | LinearSolve dual path, one pass (s) | its stats | DualSparseKLU, median (s) | its stats | partials rel |
|---|---|---|---|---|---|
| GE | 246.4 | steps 14296 (rejected 1244), f 114474, J 13052, W 14296, solves 114368 | 199.643 | steps 14296 (rejected 1244), f 114474, J 13052, W 14296, solves 114368 | 7.85e-09 |
| MET | 304.0 | steps 17129 (rejected 1293), f 137138, J 15836, W 17129, solves 137032 | 257.218 | steps 17129 (rejected 1293), f 137138, J 15836, W 17129, solves 137032 | 5.07e-11 |
