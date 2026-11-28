# KIM-195 `time`

- job 21842557 on dave7, AMD EPYC 7543 32-Core Processor; Julia 1.11.9; LinearSolve 3.87.0
- InferCell 491cb900ed25f6b851af8aa2d12aa7112fab8261; syn3a-ensemble 2075809cbb825140cb31a7bec8bc3529770c8e13 (data/model unmodified)
- threads 1, BLAS threads 1, GC threads 1
- window 0–3150.0 s in minute segments, saves at the whole minutes; Rodas5P, KLUFactorization (primal) and DualSparseKLU (dual, unless a row says the library path); abstol 1.0e-10, reltol 1.0e-8; chunk 10

primal solve (GE B = 10 at θ = 0): median 25.305 s; steps 16487 (rejected 1161), f 132002, J 15326, W 16487, solves 131896

| placement | B | passes | gradient median (s) | repeats (s) | per pass median (s) | pass 1 stats |
|---|---|---|---|---|---|---|
| GE | 10 | 1 | 180.830 | 180.728, 180.830, 180.858 | 180.830 | steps 14316 (rejected 1247), f 114634, J 13069, W 14316, solves 114528 |
| GE | 20 | 2 | 358.796 | 360.309, 357.227, 358.796 | 179.278 | steps 14296 (rejected 1244), f 114474, J 13052, W 14296, solves 114368 |
| GE | 40 | 4 | 703.564 | 703.564, 707.479, 700.170 | 177.419 | steps 14288 (rejected 1249), f 114410, J 13039, W 14288, solves 114304 |
| MET | 10 | 1 | 212.348 | 210.892, 212.348, 212.407 | 212.348 | steps 16729 (rejected 1255), f 133938, J 15474, W 16729, solves 133832 |
| MET | 20 | 2 | 426.150 | 430.596, 426.150, 425.088 | 214.312 | steps 17129 (rejected 1293), f 137138, J 15836, W 17129, solves 137032 |
| MET | 40 | 4 | 811.064 | 813.575, 811.064, 810.832 | 204.398 | steps 16824 (rejected 1300), f 134698, J 15524, W 16824, solves 134592 |
| BOTH | 10 | 1 | 229.794 | 231.135, 227.359, 229.794 | 229.794 | steps 18141 (rejected 1508), f 145234, J 16633, W 18141, solves 145128 |
| BOTH | 20 | 2 | 397.281 | 390.518, 399.284, 397.281 | 198.892 | steps 13939 (rejected 1190), f 111618, J 12749, W 13939, solves 111512 |
| BOTH | 40 | 4 | 755.578 | 762.676, 754.163, 755.578 | 189.454 | steps 13118 (rejected 1137), f 105050, J 11981, W 13118, solves 104944 |

Looser tolerances (abstol 1e-06, reltol 1e-06), GE B = 20, one repeat: primal 7.725 s (steps 4948 (rejected 506), f 39690, J 4442, W 4948, solves 39584); gradient 118.413 s; pass 1 steps 4557 (rejected 487), f 36562, J 4070, W 4557, solves 36456
