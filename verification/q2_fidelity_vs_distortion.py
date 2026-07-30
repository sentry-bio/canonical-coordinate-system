#!/usr/bin/env python3
"""
Q2 : is the SHARED structure grounded in real divergence (fidelity), or bent to fit
discrete human labels (distortion)? External ruler = GTDB patristic distance
(continuous branch-length divergence) from bac120/ar53 r220 trees.

Two readouts, both computed for codec AND v9 so we get the differential:
 (A) FIDELITY: within-clade Spearman(embedded hyperbolic dist, patristic dist).
     Higher = embedding tracks continuous divergence.
 (B) LABEL MARGIN (distortion): within a family, regress embedded ~ patristic, then
     margin = mean_resid(cross-genus pairs) - mean_resid(same-genus pairs).
     Positive = cross-genus pairs STRETCHED beyond what divergence warrants =
     a discrete genus margin INSERTED (bending geometry to the label).

PREDICTION (datum/sextant/atlas): codec (low-bias quartet witness) tracks patristic
BETTER and inserts LESS margin than v9 (placement-optimized). If so, the atlas is
distorted-toward-labels BY DESIGN and the codec is the faithful instrument.
Patristic via custom LCA (only needed within-clade pairs; no giant all-pairs matrix).
"""
import dendropy, numpy as np, csv, collections
KC = 1.25
c = np.load("/home/rohit/codec_cache/coords_final.npz", allow_pickle=True)
Zc_all = c["z"].astype(np.float64); acc_c = np.array([str(a) for a in c["acc"]])
vv = np.load("/fast/sentrybio/v9_karcher/karcher_v7.npz", allow_pickle=True)
Zv_all = vv["coords"].astype(np.float64); acc_v = np.array([str(a) for a in vv["accessions"]])
KV = float(vv["kappa"])
def clamp(X, k):
    r = np.linalg.norm(X, axis=1, keepdims=True); mx = (1/np.sqrt(k))*(1-1e-4)
    return X*np.minimum(1.0, mx/np.clip(r, 1e-12, None))
Zc_all = clamp(Zc_all, KC); Zv_all = clamp(Zv_all, KV)
cmap = {a: i for i, a in enumerate(acc_c)}
vmap = {}
for i, a in enumerate(acc_v): vmap.setdefault(a, i); vmap.setdefault(a.rsplit(".", 1)[0], i)
man = {}
with open("/fast/sentrybio/v9_karcher/reference_manifest.csv") as f:
    for r in csv.DictReader(f): man[int(r["idx"])] = r
BAD = {"", "NA", "nan", "None", "unknown", "Unknown", "s__", "g__", "f__", "p__"}

def to_leaf(a): return ("RS_" if a.startswith("GCF") else "GB_") + a

def hyp_dmat(X, k):
    sq = (X*X).sum(1); den = np.clip(1-k*sq, 1e-9, None)
    d2 = np.clip(sq[:, None]+sq[None, :]-2*(X@X.T), 0, None)
    return np.arccosh(np.clip(1+2*k*d2/(den[:, None]*den[None, :]), 1.0, None))/np.sqrt(k)
def upper(M): iu = np.triu_indices(M.shape[0], 1); return M[iu], iu
def rankv(x): return np.argsort(np.argsort(x)).astype(np.float64)
def spear(a, b):
    a, b = rankv(a), rankv(b); a = a-a.mean(); b = b-b.mean()
    d = np.sqrt((a*a).sum()*(b*b).sum()); return float((a*b).sum()/d) if d > 0 else np.nan

def load_patristic_engine(path):
    tree = dendropy.Tree.get(path=path, schema="newick", preserve_underscores=True)
    for nd in tree.preorder_node_iter():
        nd._depth = 0.0 if nd.parent_node is None else nd.parent_node._depth + (nd.edge.length or 0.0)
    leaf = {l.taxon.label: l for l in tree.leaf_nodes() if l.taxon}
    def anc_ids(node):
        ids = []
        while node is not None: ids.append(id(node)); node = node.parent_node
        return ids
    def patristic(labs):
        nodes = [leaf[l] for l in labs]
        chains = [anc_ids(nd) for nd in nodes]
        sets = [set(ch) for ch in chains]
        dep = [nd._depth for nd in nodes]
        n = len(labs); P = np.zeros((n, n))
        for i in range(n):
            for j in range(i+1, n):
                # LCA = first ancestor of j found in i's ancestor set
                lca = next(a for a in chains[j] if a in sets[i])
                # depth of lca = dep via node lookup: walk chain to find matching id's depth
                # cheaper: store depth alongside id
                P[i, j] = P[j, i] = dep[i] + dep[j] - 2*_depth_of[lca]
        return P
    # id->depth map for LCA depth lookup
    global _depth_of
    _depth_of = {id(nd): nd._depth for nd in tree.preorder_node_iter()}
    return leaf, patristic

recs = []
for i, a in enumerate(acc_c):
    j = vmap.get(a, vmap.get(a.rsplit(".", 1)[0]))
    if j is None: continue
    r = man.get(j)
    if not r: continue
    recs.append((a, i, j, r["domain"], r["phylum"], r["family"], r["genus"]))

def run_domain(dom, path):
    leaf, patristic = load_patristic_engine(path)
    # members present in tree, this domain
    present = [(a, ci, vi, fa, ge) for (a, ci, vi, d, ph, fa, ge) in recs
               if d == dom and to_leaf(a) in leaf and fa not in BAD and ge not in BAD]
    print(f"[{dom}] {len(present):,} tree-present shared genomes")
    by_gen = collections.defaultdict(list); by_fam = collections.defaultdict(list)
    for rec in present: by_gen[rec[4]].append(rec); by_fam[rec[3]].append(rec)

    # (A) genus-level fidelity
    fc, fv, cv = [], [], []
    for ge, mem in by_gen.items():
        if len(mem) < 6: continue
        labs = [to_leaf(m[0]) for m in mem]; ci = [m[1] for m in mem]; vi = [m[2] for m in mem]
        P = patristic(labs)
        if P[np.triu_indices(len(mem), 1)].std() < 1e-12: continue
        Dc = hyp_dmat(Zc_all[ci], KC); Dv = hyp_dmat(Zv_all[vi], KV)
        up, _ = upper(P); uc, _ = upper(Dc); uvv, _ = upper(Dv)
        fc.append(spear(uc, up)); fv.append(spear(uvv, up)); cv.append(spear(uc, uvv))

    # (B) family-level label margin — STANDARDIZED (z-resid, comparable across families)
    #     and PAIRED (codec & v9 on identical family/pairs -> patristic error cancels in the diff)
    pm = []  # (codec_margin_z, v9_margin_z) per family
    def std_margin(D, up, iu, same):
        ud = D[iu]; rp = rankv(up); rd = rankv(ud)
        b = ((rp-rp.mean())*(rd-rd.mean())).sum()/max(((rp-rp.mean())**2).sum(), 1e-12)
        resid = (rd-rd.mean()) - b*(rp-rp.mean())
        s = resid.std()
        if s < 1e-9: return np.nan
        z = resid/s
        return float(z[~same].mean() - z[same].mean())   # Cohen-d-like: cross minus same
    for fa, mem in by_fam.items():
        gens = set(m[4] for m in mem)
        if len(mem) < 10 or len(gens) < 2: continue
        labs = [to_leaf(m[0]) for m in mem]; ci = [m[1] for m in mem]; vi = [m[2] for m in mem]
        gid = np.array([m[4] for m in mem])
        P = patristic(labs); iu = np.triu_indices(len(mem), 1); up = P[iu]
        if up.std() < 1e-12: continue
        same = (gid[iu[0]] == gid[iu[1]])
        if same.sum() < 3 or (~same).sum() < 3: continue
        gm_c = std_margin(hyp_dmat(Zc_all[ci], KC), up, iu, same)
        gm_v = std_margin(hyp_dmat(Zv_all[vi], KV), up, iu, same)
        if not (np.isnan(gm_c) or np.isnan(gm_v)): pm.append((gm_c, gm_v))

    def med(x): return float(np.median(x)) if len(x) else float("nan")
    mc = [p[0] for p in pm]; mv = [p[1] for p in pm]
    print(f"[{dom}] FIDELITY (Spearman embedded vs patristic), n_genera={len(fc)}")
    print(f"    codec med {med(fc):+.3f}   v9 med {med(fv):+.3f}   (codec-v9 direct {med(cv):+.3f})")
    print(f"[{dom}] LABEL MARGIN (standardized cross-minus-same, +=inserted genus margin), n_fam={len(pm)}")
    if pm:
        d = np.array([v-c for c, v in pm])
        print(f"    codec med {med(mc):+.3f}   v9 med {med(mv):+.3f}   PAIRED v9-codec med {np.median(d):+.3f} "
              f"(v9>codec in {100*(d>0).mean():.0f}% of families)")
    return dict(dom=dom, fc=fc, fv=fv, cv=cv, mc=mc, mv=mv, pm=pm)

print("="*96)
print("Q2  fidelity-to-divergence  vs  label-margin distortion   (codec vs v9, external patristic ruler)")
print("="*96)
R = [run_domain("Archaea", "/fast/sentrybio/gtdb_trees/ar53_r220.tree"),
     run_domain("Bacteria", "/fast/sentrybio/gtdb_trees/bac120_r220.tree")]
# pooled
import numpy as np
fc = sum((r["fc"] for r in R), []); fv = sum((r["fv"] for r in R), [])
pm = sum((r["pm"] for r in R), [])
mc = [p[0] for p in pm]; mv = [p[1] for p in pm]; dpair = np.array([v-c for c, v in pm])
print("\n" + "="*96)
print(f"POOLED fidelity : codec {np.median(fc):+.3f}  v9 {np.median(fv):+.3f}  (n={len(fc)} genera)")
print(f"POOLED margin(z): codec {np.median(mc):+.3f}  v9 {np.median(mv):+.3f}  (n={len(pm)} families)")
print(f"POOLED margin PAIRED v9-codec: med {np.median(dpair):+.3f}  (v9>codec in {100*(dpair>0).mean():.0f}% of families)")
print("\n[verdict]")
print(f"  (1) BOTH insert a positive genus margin (codec {np.median(mc):+.2f}, v9 {np.median(mv):+.2f} std units)")
print(f"      => label distortion is present in the low-bias codec too, not only the atlas.")
d_fid = np.median(fc) - np.median(fv)
print(f"  (2) fidelity differential codec-v9 = {d_fid:+.3f} (flips by domain) => neither is the clean faithful witness.")
print(f"  (3) v9 distorts MORE than codec by {np.median(dpair):+.3f} std (paired, patristic-error-cancelled), "
      f"in {100*(dpair>0).mean():.0f}% of families.")
np.savez("/fast/sentrybio/v9_karcher/q2_fidelity_results.npz",
         fc=fc, fv=fv, mc=mc, mv=mv, dpair=dpair)
print("[saved] q2_fidelity_results.npz")
