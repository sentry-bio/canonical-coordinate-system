#!/usr/bin/env python3
"""
δ ⊥ κ orthogonality, v2 — confound-controlled.

Fixes from v1: (a) CONSTANT leaf count (subsample all trees to N0=500) so the
τ-vs-b effect isn't a leaf-count artifact; (b) reticulation strength scaled to each
tree's median distance (comparable HGT across b); (c) sparse-focused f-grid plus a
heavy tail, to separate the clean sparse regime from the heavy-reticulation breakdown.

Report orthogonality on the SPARSE subset (where the decomposition should hold) and
FULL (to show graceful breakdown). Tree floor (f=0) must be exactly 0 for all b.
"""
import sys, numpy as np
sys.path.insert(0, "/Users/rohitfenn/Golden_500_genomes_broken/deploy/static/active-geometry/experiments/universality")
from volume_entropy import volume_b, delta_hyperbolicity
from scipy.sparse.csgraph import shortest_path
rng = np.random.default_rng(0)
N0 = 500

def bary_leafD(b, L):
    n = b**L
    paths = np.array([[(i//(b**l)) % b for l in range(L)] for i in range(n)], dtype=np.int16)
    cpl = np.zeros((n, n), dtype=np.int16); match = np.ones((n, n), bool)
    for l in range(L):
        match &= (paths[:, None, l] == paths[None, :, l]); cpl += match
    D = 2.0 * (L - cpl); np.fill_diagonal(D, 0.0)
    return D

def prep(b, L):
    D = bary_leafD(b, L)
    if D.shape[0] > N0:
        idx = rng.choice(D.shape[0], N0, replace=False); D = D[np.ix_(idx, idx)]
    return D

def reticulate(D, f):
    D = D.copy(); n = D.shape[0]
    w = 0.25 * np.median(D[D > 0])           # HGT strength ∝ tree scale (comparable across b)
    n_ret = int(round(f * n))
    for _ in range(n_ret):
        i, j = rng.choice(n, 2, replace=False)
        D[i, j] = D[j, i] = min(D[i, j], w)
    if n_ret:                                 # re-close the metric (shortcuts propagate)
        D = shortest_path(D, method='D', directed=False)
    return D

PAIRS = [(0,1),(2,3),(0,2),(1,3),(0,3),(1,2)]
def defect(D, M=30000):
    n = D.shape[0]; Q = rng.integers(0, n, size=(M, 4))
    good = np.ones(M, bool)
    for a in range(4):
        for bb in range(a+1, 4): good &= Q[:,a] != Q[:,bb]
    Q = Q[good]; g = lambda a,bb: D[Q[:,a], Q[:,bb]]
    S = np.sort(np.stack([g(0,1)+g(2,3), g(0,2)+g(1,3), g(0,3)+g(1,2)], 1), 1)
    L3,L2,L1 = S[:,0],S[:,1],S[:,2]; br = L1-L3; ok = br > 1e-9
    rho = (L1[ok]-L2[ok])/br[ok]
    return float((rho > 0.05).mean())

bs = [2,3,4,5]; fs = [0.0, 0.01, 0.02, 0.05, 0.10, 0.25]; Ls = {2:9,3:6,4:5,5:4}
KV = np.full((len(bs), len(fs)), np.nan); TAU = np.zeros_like(KV)
base = {b: prep(b, Ls[b]) for b in bs}
for ib, b in enumerate(bs):
    for jf, f in enumerate(fs):
        D = reticulate(base[b], f)
        vb = volume_b(D, unit=1.0); KV[ib,jf] = np.log(vb) if (vb==vb and vb>1) else np.nan
        TAU[ib,jf] = defect(D)
    print(f"  b={b} done ({N0} leaves, constant)")

def show(name, M):
    print(f"\n{name}\n        " + "".join(f"f={f:<6.2f}" for f in fs))
    for ib,b in enumerate(bs): print(f"  b={b}   " + "".join(f"{M[ib,jf]:<8.3f}" for jf in range(len(fs))))
show("κ_vol = ln(volume_b)  [curvature axis: predict tracks b, flat in f]", KV)
show("τ = frac(ρ>0.05)      [tree-defect axis: predict tracks f, flat in b]", TAU)

bb = np.array(bs,float)[:,None]*np.ones((1,len(fs))); ff = np.ones((len(bs),1))*np.array(fs,float)[None,:]
def corr(x,y):
    x,y=x.ravel(),y.ravel(); m=np.isfinite(x)&np.isfinite(y); x,y=x[m]-x[m].mean(),y[m]-y[m].mean()
    d=np.sqrt((x*x).sum()*(y*y).sum()); return float((x*y).sum()/d) if d>0 else np.nan
sp = [i for i,f in enumerate(fs) if f <= 0.05]     # sparse subset
print("\n" + "="*68)
print("ORTHOGONALITY  (own-axis high, cross-axis ~0)")
print(f"  FULL   κ_vol: corr b {corr(KV,bb):+.3f}  corr f {corr(KV,ff):+.3f}  |  τ: corr f {corr(TAU,ff):+.3f}  corr b {corr(TAU,bb):+.3f}")
print(f"  SPARSE κ_vol: corr b {corr(KV[:,sp],bb[:,sp]):+.3f}  corr f {corr(KV[:,sp],ff[:,sp]):+.3f}  |  "
      f"τ: corr f {corr(TAU[:,sp],ff[:,sp]):+.3f}  corr b {corr(TAU[:,sp],bb[:,sp]):+.3f}   (f≤0.05)")
print(f"\n[tree floor] τ(f=0) all b: {np.round(TAU[:,0],4)}   (must be 0 — Buneman exact on trees)")
print(f"[κ_vol flat in sparse f?] per-b spread over f≤0.05: {np.round(np.nanstd(KV[:,sp],axis=1),4)}  (small = flat)")
print("\nREAD: clean orthogonality in the SPARSE regime = the four-point object splits into")
print("      two independent invariants; heavy-f cross-talk = expected breakdown when HGT")
print("      is dense enough to reshape the space. Tree floor exact = Buneman ⟂ curvature.")
np.savez("/private/tmp/claude-501/-Users-rohitfenn-Golden-500-genomes-broken-deploy-static-active-geometry/692f28a8-314b-4a7e-a40c-89da6feae999/scratchpad/kernel_ortho_v2.npz", KV=KV, TAU=TAU, bs=bs, fs=fs)
