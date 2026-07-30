#!/usr/bin/env python3
"""
Q1 : is within-clade fine structure SHARED across instruments (gauge-only theta
scramble) or IDIOSYNCRATIC (objective-driven all the way down)?

Instrument-independent test: the hyperbolic distance matrix is invariant under the
full isometry group (Mobius translations + rotations) = exactly the gauge that
scrambles theta inside a clade. So we compare codec vs v9 WITHIN-CLADE distance
matrices directly (no anchoring, no per-set PCA -> immune to the fragilities that
broke the 16.8 metric at small/biased N).

Statistics (both gauge-free):
  - Mantel (Spearman) on the within-clade distance matrices (global metric shape)
  - within-clade kNN overlap (Jaccard, local neighbor agreement)
Null: per-clade label permutation (does the codec<->v9 correspondence beat random).
Decay: within-phylum -> within-family -> within-genus. Partial-Mantel at family
(control same-genus) isolates sub-genus agreement. Aggregate over all clades of a
rank with bootstrap CI. PRE-REGISTERED verdict at the bottom.
"""
import numpy as np, csv, collections
KC = 1.25                       # codec kappa (datum 5/4, as register.py)
CAP = 120                       # cap members/clade (subsample larger) -> bounds O(n^2)+perm cost
NPERM = 80
KNN = 5
rng = np.random.default_rng(0)

c = np.load("/home/rohit/codec_cache/coords_final.npz", allow_pickle=True)
Zc_all = c["z"].astype(np.float64); acc_c = np.array([str(a) for a in c["acc"]])
vv = np.load("/fast/sentrybio/v9_karcher/karcher_v7.npz", allow_pickle=True)
Zv_all = vv["coords"].astype(np.float64); acc_v = np.array([str(a) for a in vv["accessions"]])
KV = float(vv["kappa"])
vmap = {}
for i, a in enumerate(acc_v): vmap.setdefault(a, i); vmap.setdefault(a.rsplit(".", 1)[0], i)

# clamp coords strictly inside each ball (numerical safety for arccosh formula)
def clamp(X, kappa):
    r = np.linalg.norm(X, axis=1, keepdims=True); mx = (1/np.sqrt(kappa))*(1-1e-4)
    s = np.minimum(1.0, mx/np.clip(r, 1e-12, None)); return X*s
Zc_all = clamp(Zc_all, KC); Zv_all = clamp(Zv_all, KV)

# shared genomes + taxonomy
man = {}
with open("/fast/sentrybio/v9_karcher/reference_manifest.csv") as f:
    for r in csv.DictReader(f): man[int(r["idx"])] = r
BAD = {"", "NA", "nan", "None", "unknown", "Unknown", "s__", "g__", "f__", "p__"}
recs = []   # (codec_idx, v9_idx, phylum, family, genus, species)
for i, a in enumerate(acc_c):
    j = vmap.get(a, vmap.get(a.rsplit(".", 1)[0]))
    if j is None: continue
    r = man.get(j)
    if not r: continue
    recs.append((i, j, r["phylum"], r["family"], r["genus"], r["species"]))
print(f"[data] {len(recs):,} shared+taxon'd genomes | codec kappa {KC} v9 kappa {KV:.4f}")

def hyp_dmat(X, kappa):
    sq = (X*X).sum(1); den = np.clip(1 - kappa*sq, 1e-9, None)
    d2 = np.clip(sq[:, None] + sq[None, :] - 2*(X@X.T), 0, None)
    arg = np.clip(1 + 2*kappa*d2/(den[:, None]*den[None, :]), 1.0, None)
    return np.arccosh(arg)/np.sqrt(kappa)
IU = {}
def upper(M):
    n = M.shape[0]
    if n not in IU: IU[n] = np.triu_indices(n, 1)
    return M[IU[n]]
def rankv(x): return np.argsort(np.argsort(x)).astype(np.float64)
def pear(a, b):
    a = a - a.mean(); b = b - b.mean()
    d = np.sqrt((a*a).sum()*(b*b).sum())
    return float((a*b).sum()/d) if d > 0 else np.nan

def analyze(idx_c, idx_v, genus_of=None):
    n = len(idx_c)
    if n > CAP:
        sel = rng.choice(n, CAP, replace=False); idx_c = idx_c[sel]; idx_v = idx_v[sel]
        if genus_of is not None: genus_of = genus_of[sel]
        n = CAP
    Dc = hyp_dmat(Zc_all[idx_c], KC); Dv = hyp_dmat(Zv_all[idx_v], KV)
    uc = upper(Dc); uv = upper(Dv)
    if uc.std() < 1e-9 or uv.std() < 1e-9: return None
    ruc = rankv(uc); r_obs = pear(ruc, rankv(uv))
    # kNN overlap
    kk = min(KNN, n-1)
    scn = np.argsort(Dc, axis=1)[:, 1:kk+1]; svn = np.argsort(Dv, axis=1)[:, 1:kk+1]
    sc = [set(row) for row in scn]
    jac = np.mean([len(sc[i] & set(svn[i]))/len(sc[i] | set(svn[i])) for i in range(n)])
    # partial Mantel controlling same-genus (only if genus labels vary)
    partial = np.nan
    if genus_of is not None:
        S = (genus_of[:, None] != genus_of[None, :]).astype(np.float64)
        us = upper(S)
        if us.std() > 1e-9:
            rcs = pear(ruc, rankv(us)); rvs = pear(rankv(uv), rankv(us))
            den = np.sqrt(max(1e-12, (1-rcs**2)*(1-rvs**2)))
            partial = (r_obs - rcs*rvs)/den
    # permutation null
    rp = np.empty(NPERM); jp = np.empty(NPERM)
    for t in range(NPERM):
        p = rng.permutation(n)
        Dvp = Dv[np.ix_(p, p)]
        rp[t] = pear(ruc, rankv(upper(Dvp)))
        svnp = np.argsort(Dvp, axis=1)[:, 1:kk+1]
        jp[t] = np.mean([len(sc[i] & set(svnp[i]))/len(sc[i] | set(svnp[i])) for i in range(n)])
    p_r = (1 + (rp >= r_obs).sum())/(1 + NPERM)
    return dict(n=n, r=r_obs, jac=jac, r_null=float(rp.mean()), jac_null=float(jp.mean()),
                p=p_r, partial=partial)

def group(level):
    g = collections.defaultdict(lambda: ([], []))
    genusidx = {"phylum": 2, "family": 3, "genus": 4}[level]
    for (ci, vi, ph, fa, ge, sp) in recs:
        key = (ph, fa, ge)[{"phylum": 0, "family": 1, "genus": 2}[level]]
        if key in BAD: continue
        g[key][0].append(ci); g[key][1].append(vi)
    return g

def run_level(level, min_members, extra_filter=None, with_partial=False):
    g = group(level)
    genmap = {}  # for partial: genus label per member
    if with_partial:
        gm = collections.defaultdict(list)
        for (ci, vi, ph, fa, ge, sp) in recs:
            if fa not in BAD: gm[fa].append(ge)
    out = []
    for key, (cis, vis) in g.items():
        if len(cis) < min_members: continue
        idx_c = np.array(cis); idx_v = np.array(vis)
        genus_of = None
        if with_partial:
            genus_of = np.array(gm[key])
        res = analyze(idx_c, idx_v, genus_of=genus_of)
        if res: res["clade"] = key; out.append(res)
    return out

def summ(name, out):
    if not out: print(f"  {name:26s} (no clades)"); return
    r = np.array([o["r"] for o in out]); j = np.array([o["jac"] for o in out])
    rn = np.array([o["r_null"] for o in out]); jn = np.array([o["jac_null"] for o in out])
    p = np.array([o["p"] for o in out]); frac_sig = float((p < 0.05).mean())
    # bootstrap CI on median r
    bs = [np.median(rng.choice(r, len(r), replace=True)) for _ in range(2000)]
    lo, hi = np.percentile(bs, [2.5, 97.5])
    part = np.array([o["partial"] for o in out if not np.isnan(o["partial"])])
    pstr = f" | partial(vs genus) med {np.median(part):+.3f}" if len(part) else ""
    print(f"  {name:26s} nclades={len(out):3d}  Mantel r med {np.median(r):+.3f} [95%CI {lo:+.2f},{hi:+.2f}] "
          f"(null {np.median(rn):+.3f})  kNN {np.median(j):.3f} (null {np.median(jn):.3f})  "
          f"sig@5% {100*frac_sig:3.0f}%{pstr}")
    return dict(r_med=float(np.median(r)), r_ci=[float(lo), float(hi)], knn=float(np.median(j)),
                r_null=float(np.median(rn)), knn_null=float(np.median(jn)), frac_sig=frac_sig,
                nclades=len(out))

print("\n" + "="*100)
print("Q1  within-clade DISTANCE-STRUCTURE agreement (codec dim-16  vs  v9 dim-129), gauge-free")
print("="*100)
res = {}
print("[decay curve: coarse -> fine]")
res["phylum"] = summ("within-PHYLUM", run_level("phylum", 50))
res["family"] = summ("within-FAMILY", run_level("family", 20, with_partial=True))
res["genus"]  = summ("within-GENUS (finest)", run_level("genus", 12))

print("\n[pre-registered verdict]")
gm = res["genus"]["r_med"]; gnull = res["genus"]["r_null"]; gci = res["genus"]["r_ci"]
gknn = res["genus"]["knn"]; gknnnull = res["genus"]["knn_null"]
if gci[0] <= gnull + 0.03 and gknn <= gknnnull*1.5:
    print(f"  IDIOSYNCRATIC at the tips: within-genus Mantel {gm:+.3f} (CI lo {gci[0]:+.3f}) ~ null {gnull:+.3f};")
    print(f"  kNN {gknn:.3f} ~ null {gknnnull:.3f}. Fine scale is objective-driven, NOT shared -> datum is coarse-only.")
elif gm > gnull + 0.1 and gci[0] > gnull:
    print(f"  SHARED fine structure (theta was gauge): within-genus Mantel {gm:+.3f} >> null {gnull:+.3f}, CI excludes null.")
    print(f"  Fine canonicity exists IN DISTANCE -> leaf-contrastive must be tested for whether it DEGRADES it (Q2).")
else:
    print(f"  PARTIAL/graded: within-genus Mantel {gm:+.3f} vs null {gnull:+.3f} (CI {gci}) — real but weak; report the gradient.")
print(f"\n  resolution gradient  phylum {res['phylum']['r_med']:+.3f} -> family {res['family']['r_med']:+.3f} "
      f"-> genus {res['genus']['r_med']:+.3f}   (null ~{gnull:+.3f})")

np.savez("/fast/sentrybio/v9_karcher/q1_gauge_results.npz",
         **{f"{lvl}_{k}": v for lvl, d in res.items() if d for k, v in d.items() if not isinstance(v, list)})
print("[saved] q1_gauge_results.npz")
