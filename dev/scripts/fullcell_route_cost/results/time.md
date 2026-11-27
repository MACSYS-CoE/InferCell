# KIM-195 `time`

- job 21785045 on dave9, AMD EPYC 7543 32-Core Processor; Julia 1.11.9; LinearSolve 3.87.0
- threads 1, BLAS threads 1, GC threads 1
- window 0–3150.0 s in minute segments, saves at the whole minutes; Rodas5P, KLUFactorization (primal) and DualSparseKLU (dual, unless a row says the library path); abstol 1.0e-10, reltol 1.0e-8; chunk 10

primal solve (GE B = 10 at θ = 0): median 10.456 s; steps 6415 (rejected 416), f 51426, J 5999, W 6415, solves 51320

| placement | B | passes | gradient median (s) | repeats (s) | per pass median (s) | pass 1 stats |
|---|---|---|---|---|---|---|
| GE | 10 | 1 | 69.310 | 67.277, 69.310, 69.925 | 69.310 | steps 5327 (rejected 465), f 42722, J 4862, W 5327, solves 42616 |
| GE | 20 | 2 | 139.944 | 137.554, 139.944, 140.046 | 69.684 | steps 5325 (rejected 465), f 42706, J 4860, W 5325, solves 42600 |
| GE | 40 | 4 | 274.571 | 270.815, 274.571, 275.143 | 68.435 | steps 5326 (rejected 466), f 42714, J 4860, W 5326, solves 42608 |
| MET | 10 | 1 | 79.370 | 78.544, 79.776, 79.370 | 79.370 | steps 6159 (rejected 449), f 49378, J 5710, W 6159, solves 49272 |
| MET | 20 | 2 | 158.905 | 156.999, 159.096, 158.905 | 79.104 | steps 6348 (rejected 446), f 50890, J 5902, W 6348, solves 50784 |
| MET | 40 | 4 | 297.739 | 297.739, 298.908, 292.802 | 74.808 | steps 6460 (rejected 442), f 51786, J 6018, W 6460, solves 51680 |
| BOTH | 10 | 1 | 72.869 | 73.188, 72.869, 70.747 | 72.869 | steps 6084 (rejected 466), f 48778, J 5618, W 6084, solves 48672 |
| BOTH | 20 | 2 | 131.021 | 130.159, 131.021, 131.491 | 65.445 | steps 5155 (rejected 460), f 41346, J 4695, W 5155, solves 41240 |
| BOTH | 40 | 4 | 264.503 | 262.443, 264.503, 265.876 | 64.815 | steps 5228 (rejected 454), f 41930, J 4774, W 5228, solves 41824 |

Looser tolerances (abstol 1e-06, reltol 1e-06), GE B = 20, one repeat: primal 3.437 s (steps 2301 (rejected 56), f 18514, J 2245, W 2301, solves 18408); gradient 49.177 s; pass 1 steps 1925 (rejected 34), f 15506, J 1891, W 1925, solves 15400
