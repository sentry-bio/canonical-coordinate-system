#!/usr/bin/env python3
"""
Q3c : per-instrument scale-free closure of the quartet test.

Fixes Q3b's defect (absolute S threshold, v9-scale, applied to both -> unfair to codec).
Both axes are now scale-free RATIOS, invariant under uniform rescaling (codec-embeds-smaller):
    S_rel = (L1 - L3) / (sum of 6 pairwise distances)   in [0, ~0.5]  (relative branch = signal)
    rho   = (L1 - L2) / (L1 - L3)                        in [0, 1]     (slack / branch = ambiguity)
=> ONE threshold is now fair to both instruments.

HEADLINE: rho as a function of S_rel, codec vs v9 OVERLAID (fidelity at matched relative signal).
References: exact tree (rho=0 linchpin), noisy-tree at 5%/10% distance noise (realistic tree baseline),
uniform H^2 scatter (rho ceiling). Split PASS-RATE (signal availability) from CONDITIONAL rho (integrity).

Pre-registered: curves OVERLAP => codec==v9 fidelity, codec just coarser (leading expectation).
                v9 curve ABOVE codec => atlas genuinely less tree-faithful == Q2's placement-distortion,
                triangulated by an independent (four-point) ruler.
"""
import numpy as np, csv, collections
rng = np.random.default_rng(0); KC = 1.25
PAIRS=[(0,1),(2,3),(0,2),(1,3),(0,3),(1,2)]

def hyp_d_pairs(Xa,Xb,k):
    sa=(Xa*Xa).sum(1); sb=(Xb*Xb).sum(1)
    da=np.clip(1-k*sa,1e-9,None); db=np.clip(1-k*sb,1e-9,None)
    return np.arccosh(np.clip(1+2*k*((Xa-Xb)**2).sum(1)/(da*db),1.0,None))/np.sqrt(k)
def sample_quads(pool,M):
    pool=np.asarray(pool); out=[]
    if len(pool)<4: return None
    while len(out)<M:
        c=pool[rng.integers(0,len(pool),size=(max(M,64),4))]
        good=np.ones(len(c),bool)
        for i in range(4):
            for j in range(i+1,4): good&=c[:,i]!=c[:,j]
        out.extend(c[good].tolist())
    return np.array(out[:M])
def metrics(sixd):
    S1=sixd[(0,1)]+sixd[(2,3)]; S2=sixd[(0,2)]+sixd[(1,3)]; S3=sixd[(0,3)]+sixd[(1,2)]
    M=np.sort(np.stack([S1,S2,S3],1),1); L3,L2,L1=M[:,0],M[:,1],M[:,2]
    dsum=sum(sixd[p] for p in PAIRS); dmax=np.max(np.stack([sixd[p] for p in PAIRS],1),1)
    branch=(L1-L3)
    S_rel=branch/np.clip(dsum,1e-9,None)
    rho=np.full_like(S_rel,np.nan); ok=branch>1e-12; rho[ok]=(L1[ok]-L2[ok])/branch[ok]
    return S_rel,rho,dmax
def from_coords(pool,X,k,M):
    Q=sample_quads(pool,M)
    if Q is None: return None
    P=X[Q]; return metrics({p:hyp_d_pairs(P[:,p[0]],P[:,p[1]],k) for p in PAIRS})
def from_D(pool,D,M):
    Q=sample_quads(pool,M)
    if Q is None: return None
    return metrics({p:D[Q[:,p[0]],Q[:,p[1]]] for p in PAIRS})

def random_additive_tree(N):
    adj=collections.defaultdict(list)
    def ae(u,v,w): adj[u].append((v,w)); adj[v].append((u,w))
    nid=N; ae(0,1,rng.exponential(1)); edges=[(0,1)]
    for leaf in range(2,N):
        ei=rng.integers(len(edges)); u,v=edges[ei]; w=next(ww for nn,ww in adj[u] if nn==v)
        adj[u]=[(nn,ww) for nn,ww in adj[u] if nn!=v]; adj[v]=[(nn,ww) for nn,ww in adj[v] if nn!=u]
        m=nid; nid+=1; w1=rng.random()*w; ae(u,m,w1); ae(m,v,w-w1); ae(m,leaf,rng.exponential(1))
        edges=[e for i,e in enumerate(edges) if i!=ei]+[(u,m),(m,v),(m,leaf)]
    D=np.zeros((N,N))
    for s in range(N):
        dist={s:0.0}; st=[s]
        while st:
            x=st.pop()
            for nn,ww in adj[x]:
                if nn not in dist: dist[nn]=dist[x]+ww; st.append(nn)
        for t in range(N): D[s,t]=dist[t]
    return D
def noisy_tree(D,sigma):
    N=D.shape[0]; E=rng.normal(0,sigma,(N,N)); E=(E+E.T)/2
    Dn=D*(1+E); np.fill_diagonal(Dn,0); return np.clip(Dn,1e-6,None)
def H2_uniform(N,kappa,rmax):
    sk=np.sqrt(kappa); u=rng.random(N); rho=np.arccosh(1+u*(np.cosh(sk*rmax)-1))/sk
    th=rng.random(N)*2*np.pi; rp=np.tanh(sk*rho/2)/sk; return np.c_[rp*np.cos(th),rp*np.sin(th)]

def curve(S_rel,rho,edges):
    out=[]
    for b in range(len(edges)-1):
        m=np.isfinite(rho)&(S_rel>=edges[b])&(S_rel<edges[b+1])
        out.append(np.median(rho[m]) if m.sum()>15 else np.nan)
    return out

print("="*98)
print("Q3c  scale-free quartet closure: rho( S_rel ) — fidelity at MATCHED relative signal, per instrument")
print("="*98)
Dtree=random_additive_tree(200)
sT,rT,_=from_D(np.arange(200),Dtree,10000)
print(f"[LINCHPIN] exact tree rho med {np.nanmedian(rT):.2e} (=0, estimator unbiased across scale)")

# real data
c=np.load("/home/rohit/codec_cache/coords_final.npz",allow_pickle=True)
Zc=c["z"].astype(np.float64); acc_c=np.array([str(a) for a in c["acc"]])
v=np.load("/fast/sentrybio/v9_karcher/karcher_v7.npz",allow_pickle=True)
Zv=v["coords"].astype(np.float64); acc_v=np.array([str(a) for a in v["accessions"]]); KV=float(v["kappa"])
def clamp(X,k):
    r=np.linalg.norm(X,axis=1,keepdims=True); mx=(1/np.sqrt(k))*(1-1e-4); return X*np.minimum(1.0,mx/np.clip(r,1e-12,None))
Zc=clamp(Zc,KC); Zv=clamp(Zv,KV)
vmap={}
for i,a in enumerate(acc_v): vmap.setdefault(a,i); vmap.setdefault(a.rsplit(".",1)[0],i)
man={}
with open("/fast/sentrybio/v9_karcher/reference_manifest.csv") as f:
    for r in csv.DictReader(f): man[int(r["idx"])]=r
BAD={"","NA","nan","None","unknown","Unknown","s__","g__","f__","p__"}
recs=[]
for i,a in enumerate(acc_c):
    j=vmap.get(a,vmap.get(a.rsplit(".",1)[0]))
    if j is None: continue
    r=man.get(j)
    if r: recs.append((i,j,r["phylum"],r["family"],r["genus"]))
ci=np.array([r[0] for r in recs]); vi=np.array([r[1] for r in recs])

sc,rc,dc=from_coords(ci,Zc,KC,40000)
sv,rv,dv=from_coords(vi,Zv,KV,40000)
# noisy-tree + scatter references
Dn5=noisy_tree(Dtree,0.05); s5,r5,_=from_D(np.arange(200),Dn5,10000)
Dn10=noisy_tree(Dtree,0.10); s10,r10,_=from_D(np.arange(200),Dn10,10000)
Xsc=H2_uniform(200,KC,3.0); ssc,rsc,_=from_coords(np.arange(200),Xsc,KC,10000)

edges=np.quantile(np.concatenate([sc,sv]),np.linspace(0.05,0.98,9))
print(f"\n[HEADLINE] median rho vs S_rel band  (tree=0, more signal ->)  n=40k/instrument")
print(f"{'S_rel band':>16s} {'codec':>7s} {'v9':>7s} | {'noisy5%':>7s} {'noisy10%':>8s} {'scatter':>7s}")
cc=curve(sc,rc,edges); cv=curve(sv,rv,edges); c5=curve(s5,r5,edges); c10=curve(s10,r10,edges); csc=curve(ssc,rsc,edges)
for b in range(len(edges)-1):
    f=lambda x: f"{x:.3f}" if np.isfinite(x) else "  -  "
    print(f"  [{edges[b]:.3f},{edges[b+1]:.3f}] {f(cc[b]):>7s} {f(cv[b]):>7s} | {f(c5[b]):>7s} {f(c10[b]):>8s} {f(csc[b]):>7s}")

# threshold = S_rel where BOTH instruments' rho drops clearly below scatter (data-driven)
thr=None
for b in range(len(edges)-1):
    if np.isfinite(cc[b]) and np.isfinite(cv[b]) and max(cc[b],cv[b])<0.5*np.nanmedian(rsc):
        thr=edges[b]; break
thr=thr if thr is not None else np.median(np.concatenate([sc,sv]))
print(f"\n[threshold] S_rel > {thr:.3f}  (where rho clearly separates from scatter {np.nanmedian(rsc):.3f})")
def passrate(S_rel): return float((S_rel>thr).mean())
def cond_rho(S_rel,rho): m=np.isfinite(rho)&(S_rel>thr); return float(np.median(rho[m])) if m.sum()>20 else np.nan
print(f"[signal availability — PASS-RATE above threshold]  codec {100*passrate(sc):.0f}%   v9 {100*passrate(sv):.0f}%")
print(f"[conditional integrity — rho | S_rel>thr]          codec {cond_rho(sc,rc):.4f}   v9 {cond_rho(sv,rv):.4f}   (scatter {np.nanmedian(rsc):.3f})")

# diameter re-check AFTER gate
m=np.isfinite(rv)&(sv>thr); dvg=dv[m]; rvg=rv[m]; qd=np.quantile(dvg,[0,.33,.66,1.0])
print(f"\n[diameter re-check post-gate (v9): rho should be FLAT now]")
for b in range(3):
    mm=(dvg>=qd[b])&(dvg<=qd[b+1]); print(f"    diam[{qd[b]:.2f},{qd[b+1]:.2f}]  rho med {np.median(rvg[mm]):.4f}")

def strat(level):
    idx={"phylum":2,"family":3,"genus":4}[level]; g=collections.defaultdict(list)
    for k,r in enumerate(recs):
        if r[idx] not in BAD: g[r[idx]].append(k)
    return {kk:rows for kk,rows in g.items() if len(rows)>=12}
print(f"\n[stratified: S_rel(signal) and rho|thr (integrity), matched threshold]")
print(f"{'scale':16s} {'Srel codec':>10s} {'Srel v9':>8s} | {'rho codec':>9s} {'rho v9':>7s}")
print(f"{'GLOBAL':16s} {np.median(sc):>10.3f} {np.median(sv):>8.3f} | {cond_rho(sc,rc):>9.4f} {cond_rho(sv,rv):>7.4f}")
for level in ["phylum","family","genus"]:
    st=strat(level); Sc=[];Sv=[];Rc=[];Rv=[]
    for kk,rows in st.items():
        rows=np.array(rows)
        a=from_coords(ci[rows],Zc,KC,300); b=from_coords(vi[rows],Zv,KV,300)
        if a: Sc.append(np.nanmedian(a[0])); rr=cond_rho(a[0],a[1]);  Rc.append(rr) if np.isfinite(rr) else None
        if b: Sv.append(np.nanmedian(b[0])); rr=cond_rho(b[0],b[1]);  Rv.append(rr) if np.isfinite(rr) else None
    print(f"{'within-'+level:16s} {np.median(Sc):>10.3f} {np.median(Sv):>8.3f} | {np.nanmedian(Rc):>9.4f} {np.nanmedian(Rv):>7.4f}")

# sextant dial, scale-free
st=strat("genus"); perg={}
for kk,rows in st.items():
    b=from_coords(vi[np.array(rows)],Zv,KV,500)
    if b:
        rr=cond_rho(b[0],b[1])
        if np.isfinite(rr): perg[kk]=rr
vals=np.array(list(perg.values())); srt=sorted(perg.items(),key=lambda x:x[1])
print(f"\n[SEXTANT dial (scale-free): per-genus rho|thr  med {np.median(vals):.4f}  range [{vals.min():.4f},{vals.max():.4f}]  n={len(vals)}]")
print("  most faithful:", [f"{k}={v:.3f}" for k,v in srt[:4]])
print("  most ambiguous:", [f"{k}={v:.3f}" for k,v in srt[-4:]])
np.savez("/fast/sentrybio/v9_karcher/q3c_results.npz",edges=edges,codec=cc,v9=cv,noisy5=c5,noisy10=c10,scatter=csc,
         thr=thr,pass_c=passrate(sc),pass_v=passrate(sv),perg=vals)
print("[saved] q3c_results.npz")
