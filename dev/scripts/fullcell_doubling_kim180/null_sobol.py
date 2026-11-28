# KIM-180: a synthetic check of prereg §4 at KIM-58's size (N = 10^3 draws, 838 constants, 20 bins).
# It shows why §4.3 needs a family-wise threshold: under a pure null the top 25 binned S_i
# sum to about 0.4. Run from the repo root, with numpy (e.g. /fred/oz022/tkimpson/miniconda3,
# numpy 2.2.6): python3 -I dev/scripts/fullcell_doubling_kim180/null_sobol.py (under 1 min).
# Output of 2026-11-28 (seed 20261128): pure null, top-25 sum 0.41, permutation threshold
# 0.0326, 0 eligible; diffuse (5 at 0.1, 50 at 0.004), threshold 0.0327, the 5 strong eligible
# and no truly zero one, eligible sum 0.501. The bootstrap SE (0.011) overstates the noise
# SD of a null S_i (0.006), so a rule on k SE is conservative and its error rate unknown
# (2 SE admits 1 false constant here); the permutation threshold is calibrated.
import numpy as np
rng = np.random.default_rng(20261128)
N, P, B = 1000, 838, 20
def S_bin(u, T):
    V = T.var()
    order = np.argsort(u, axis=0)
    out = np.empty(u.shape[1])
    for i in range(u.shape[1]):
        Ts = T[order[:, i]].reshape(B, -1)
        mb = Ts.mean(1); W = Ts.var(1, ddof=1).mean()
        out[i] = (((mb - T.mean())**2).mean() - (B - 1) / N * W) / V
    return out
def boot_se(u, T, reps=100):
    S = np.array([S_bin(u[idx], T[idx]) for idx in (rng.integers(0, N, N) for _ in range(reps))])
    return S.std(0)
for case in ["pure null", "diffuse: 5 strong (0.1 each) + 50 weak (0.004 each), rest noise"]:
    u = rng.standard_normal((N, P))
    if case == "pure null":
        T = rng.standard_normal(N); truth = np.zeros(P)
    else:
        c = np.zeros(P); c[:5] = np.sqrt(0.1); c[5:55] = np.sqrt(0.004)
        T = u @ c + rng.standard_normal(N) * np.sqrt(1 - c @ c); truth = c**2
    S = S_bin(u, T); se = boot_se(u, T)
    elig = S > 2 * se
    print(case, "| sum S over all", S.sum().round(3), "| true sum", truth.sum().round(3))
    print("  eligible at 2 SE:", elig.sum(), "of which truly zero:", (elig & (truth == 0)).sum())
    print("  noise S: sd", S[truth == 0].std().round(4), "max", S[truth == 0].max().round(4), "median SE", np.median(se).round(4))
    top = np.argsort(-S)[:25]; print("  top-25 truly nonzero:", (truth[top] > 0).sum(), " top-25 sum", S[top].sum().round(3))
    for z in [3, 4]:
        e = S > z * se; print(f"  eligible at {z} SE:", e.sum(), "false:", (e & (truth == 0)).sum())

# The permutation max-null threshold of prereg §4.3 (200 permutations here, 1,000 there).
for name in ["null", "diffuse"]:
    u = rng.standard_normal((N, P))
    c = np.zeros(P)
    if name == "diffuse":
        c[:5] = np.sqrt(0.1); c[5:55] = np.sqrt(0.004)
    T = u @ c + rng.standard_normal(N) * np.sqrt(1 - c @ c)
    S = S_bin(u, T)
    mx = [S_bin(u, rng.permutation(T)).max() for _ in range(200)]
    thr = np.quantile(mx, 0.95); e = S > thr
    print(name, "| permutation threshold", round(thr, 4), "| eligible", e.sum(), "| truly zero among them",
          (e & (c == 0)).sum(), "| eligible sum", S[e].sum().round(3))
