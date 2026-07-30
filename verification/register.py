"""
register.py — the datum REGISTRATION measurement: how tightly do the codec (dim-16, low-bias witness)
and v9 atlas (dim-129, operational realization) share ONE 2D coordinate backbone?
Not topology agreement (had that: 0.79) and not grounding — this is the actual COORDINATE alignment
that makes wiring possible. Extract each instrument's 2D backbone (log-map -> PCA-2), Procrustes-align
on shared genomes, report the residual, and DECOMPOSE it: reducible (determinate tree = the precision
target) vs overlay (v9 sees non-tree structure the 2D datum correctly omits). Anchors = gauge sanity.
"""
import numpy as np, csv
from scipy.spatial import procrustes
from scipy.stats import spearmanr
KAPPA = 1.25; BALLR = 1 / KAPPA ** 0.5
ANCHORS = {"GCF_000005845": "E.coli(meridian)", "GCF_000091665": "M.jann(chirality)"}


def logmap0(x):
    n = np.linalg.norm(x, axis=1, keepdims=True).clip(1e-9, BALLR * (1 - 1e-6))
    return (2 / KAPPA ** 0.5) * np.arctanh(KAPPA ** 0.5 * n) * (x / n)


def pca2(X):
    Xc = X - X.mean(0); _, _, Vt = np.linalg.svd(Xc, full_matrices=False); return Xc @ Vt[:2].T


c = np.load("/home/rohit/codec_cache/coords_final.npz", allow_pickle=True)
zc_all = c["z"]; acc_c = c["acc"].astype(str); tax = c["tax"]
v = np.load("/fast/sentrybio/v9_karcher/karcher_v7.npz", allow_pickle=True)
zv_all = v["coords"].astype(np.float64); acc_v = v["accessions"].astype(str)
vmap = {}
for i, a in enumerate(acc_v): vmap.setdefault(a, i); vmap.setdefault(a.rsplit(".", 1)[0], i)
pairs = [(i, vmap[acc_c[i]] if acc_c[i] in vmap else vmap[acc_c[i].rsplit(".", 1)[0]])
         for i in range(len(acc_c)) if acc_c[i] in vmap or acc_c[i].rsplit(".", 1)[0] in vmap]
ci = np.array([p[0] for p in pairs]); vi = np.array([p[1] for p in pairs])
zc = zc_all[ci].astype(np.float64); zv = zv_all[vi]; taxs = tax[ci]; accs = acc_c[ci]
print(f"[reg] {len(ci):,} shared genomes (codec {len(acc_c):,} / v9 {len(acc_v):,})")

Pc = pca2(logmap0(zc)); Pv = pca2(logmap0(zv))
rc = np.linalg.norm(zc, axis=1); rv = np.linalg.norm(zv, axis=1)
m1, m2, disp = procrustes(Pv, Pc)                    # m1=v9(std), m2=codec aligned to v9; disp in [0,~1]
resid = np.linalg.norm(m1 - m2, axis=1)              # per-genome registration residual
th1 = np.arctan2(m1[:, 1], m1[:, 0]); th2 = np.arctan2(m2[:, 1], m2[:, 0])
dth = np.abs((th1 - th2 + np.pi) % (2 * np.pi) - np.pi)
print(f"\n[reg] === REGISTRATION PRECISION (codec 2D-backbone <-> v9 2D-backbone) ===")
print(f"  Procrustes disparity        = {disp:.4f}   (0=identical frame, ~1=unrelated)")
print(f"  radial agreement (Spearman) = {spearmanr(rc, rv).correlation:+.3f}   (depth axis)")
print(f"  angular residual: median={np.degrees(np.median(dth)):.1f}deg  within-30deg={100*(dth<np.pi/6).mean():.0f}%  within-45deg={100*(dth<np.pi/4).mean():.0f}%")

# anchor gauge sanity
print(f"\n[reg] anchor gauge check (should have small residual if the anchored frame is consistent):")
base = np.array([a.rsplit(".", 1)[0] for a in accs])
for a, name in ANCHORS.items():
    hit = np.where(base == a)[0]
    if len(hit): print(f"  {name:18s} residual={resid[hit[0]]:.3f}  (median all={np.median(resid):.3f})")
    else: print(f"  {name:18s} NOT in shared set")

# decompose: reducible (determinate tree, dense clade) vs overlay (sparse / where v9 sees non-tree)
from collections import Counter
gen = taxs[:, 5]; gc = Counter(gen.tolist())
dense = np.array([gc[int(g)] >= 5 for g in gen])
print(f"\n[reg] === DECOMPOSITION (is the gap reducible noise or the overlay?) ===")
print(f"  determinate (genus>=5 members, {dense.sum():,}): median residual = {np.median(resid[dense]):.3f}")
print(f"  sparse/indeterminate ({(~dense).sum():,}):        median residual = {np.median(resid[~dense]):.3f}")
print(f"  ratio sparse/dense = {np.median(resid[~dense])/max(np.median(resid[dense]),1e-9):.2f}x")
print(f"\n[reg] READ: low disparity + high radial corr + small angular residual = wireable shared frame.")
print(f"      residual much larger in sparse than dense => gap is REDUCIBLE (tighten w/ data), not overlay.")
