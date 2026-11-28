# KIM-195 `threaded`

- job 21842559 on dave5, AMD EPYC 7543 32-Core Processor; Julia 1.11.9; LinearSolve 3.87.0
- InferCell 491cb900ed25f6b851af8aa2d12aa7112fab8261; syn3a-ensemble 2075809cbb825140cb31a7bec8bc3529770c8e13 (data/model unmodified)
- threads 4, BLAS threads 1, GC threads 5
- window 0–3150.0 s in minute segments, saves at the whole minutes; Rodas5P, KLUFactorization (primal) and DualSparseKLU (dual, unless a row says the library path); abstol 1.0e-10, reltol 1.0e-8; chunk 10

primal solve (GE B = 10 at θ = 0): median 26.497 s; steps 16487 (rejected 1161), f 132002, J 15326, W 16487, solves 131896

| placement | B | passes | gradient median (s) | repeats (s) | per pass median (s) | pass 1 stats |
|---|---|---|---|---|---|---|
| GE | 20 | 2 | 202.772 | 199.075, 206.468 | 200.106 | steps 14296 (rejected 1244), f 114474, J 13052, W 14296, solves 114368 |
| GE | 40 | 4 | 246.185 | 249.263, 243.106 | 240.580 | steps 14288 (rejected 1249), f 114410, J 13039, W 14288, solves 114304 |
| MET | 20 | 2 | 242.808 | 241.743, 243.873 | 242.474 | steps 17129 (rejected 1293), f 137138, J 15836, W 17129, solves 137032 |
| MET | 40 | 4 | 275.041 | 285.078, 265.004 | 263.571 | steps 16824 (rejected 1300), f 134698, J 15524, W 16824, solves 134592 |
| BOTH | 20 | 2 | 242.437 | 245.723, 239.151 | 218.856 | steps 13939 (rejected 1190), f 111618, J 12749, W 13939, solves 111512 |
| BOTH | 40 | 4 | 267.827 | 267.494, 268.160 | 235.513 | steps 13118 (rejected 1137), f 105050, J 11981, W 13118, solves 104944 |
