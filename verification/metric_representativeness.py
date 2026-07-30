#!/usr/bin/env python3
"""
Coherence: at what N does register.py's Procrustes-theta metric become trustworthy?

Take the KNOWN-GOOD codec<->v9-deployed pair (n=23043 -> 16.8 deg) and subsample it.
If 16.8 only emerges at large N and small subsamples inflate toward random (~90),
then the metric is a LARGE-N aggregate estimator and every small-N cross-version
number (my pair2=65 @250, pair3=39 @51; the old build_crossver 63/39.5) is noise,
not signal. Also test first-N (biased slice) vs random-N to separate N from
representativeness. Pure numpy, no GPU, no embedding.
"""
import numpy as np
from scipy.spatial import procrustes
KAPPA = 1.25; BALLR = 1/KAPPA**0.5
def logmap0(x):
    n = np.linalg.norm(x, axis=1, keepdims=True).clip(1e-9, BALLR*(1-1e-6))
    return (2/KAPPA**0.5)*np.arctanh(KAPPA**0.5*n)*(x/n)
def pca2(X):
    Xc = X - X.mean(0); _,_,Vt = np.linalg.svd(Xc, full_matrices=False); return Xc @ Vt[:2].T
def reg_angle(Za, Zb):
    Pa, Pb = pca2(logmap0(Za)), pca2(logmap0(Zb))
    m1, m2, _ = procrustes(Pa, Pb)
    th1 = np.arctan2(m1[:,1],m1[:,0]); th2 = np.arctan2(m2[:,1],m2[:,0])
    return np.degrees(np.median(np.abs((th1-th2+np.pi)%(2*np.pi)-np.pi)))

c = np.load("/home/rohit/codec_cache/coords_final.npz", allow_pickle=True)
Zc = c["z"].astype(np.float64); acc_c = np.array([str(a) for a in c["acc"]])
v = np.load("/fast/sentrybio/v9_karcher/karcher_v7.npz", allow_pickle=True)
Zv = v["coords"].astype(np.float64); acc_v = np.array([str(a) for a in v["accessions"]])
mv = {}
for i,a in enumerate(acc_v): mv.setdefault(a,i); mv.setdefault(a.rsplit(".",1)[0],i)
ci, vi = [], []
for i,a in enumerate(acc_c):
    j = mv.get(a, mv.get(a.rsplit(".",1)[0]))
    if j is not None: ci.append(i); vi.append(j)
ci, vi = np.array(ci), np.array(vi)
A, B = Zc[ci], Zv[vi]
print(f"[load] {len(ci):,} shared codec<->v9-deployed pairs; full-set angular = {reg_angle(A,B):.1f} deg\n")

rng = np.random.default_rng(0)
print(f"{'N':>7} | {'random-N (mean+-std over 8 draws)':>34} | {'first-N (biased slice)':>22}")
for N in [51, 100, 250, 500, 1000, 2000, 5000, 10000, 20000]:
    if N > len(ci): continue
    vals = []
    for _ in range(8):
        idx = rng.choice(len(ci), N, replace=False)
        vals.append(reg_angle(A[idx], B[idx]))
    vals = np.array(vals)
    firstN = reg_angle(A[:N], B[:N])
    print(f"{N:>7} | {vals.mean():>14.1f} +- {vals.std():>5.1f}  ({vals.min():.0f}-{vals.max():.0f})  | {firstN:>20.1f}")

print("\nREAD: if random-N stays ~16-20 down to small N, the metric is N-robust and my")
print("      pair2/3 failures are REPRESENTATIVENESS (biased first-250). If random-N ALSO")
print("      inflates at small N, the metric itself needs large N -> any cross-ver test must be large+random.")
