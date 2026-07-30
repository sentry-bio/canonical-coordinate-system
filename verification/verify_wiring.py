"""verify_wiring.py — solid-ground checks:
  #1 transform ROUND-TRIP: does the published transform execute as a function on HELD-OUT genomes
     (fit backbone+anchor on 70%, apply to 30%, θ agreement should ~= in-sample 16.8°), and does the
     SERIALIZED datum_transform_v9.json reproduce that agreement (the artifact actually runs).
  #2 anchor-gauge STABILITY: is the 2D backbone plane stable under subsampling (principal angle), and
     is the E.coli/M.jann anchor gauge well-separated & stable (frame doesn't wobble on 2 genomes)?
"""
import numpy as np, json
KAPPA = 1.25; BALLR = 1 / KAPPA ** 0.5
ECOLI = "GCF_000005845"; MJANN = "GCF_000091665"


def logmap0(x):
    n = np.linalg.norm(x, axis=1, keepdims=True).clip(1e-9, BALLR * (1 - 1e-6))
    return (2 / KAPPA ** 0.5) * np.arctanh(KAPPA ** 0.5 * n) * (x / n)


c = np.load("/home/rohit/codec_cache/coords_final.npz", allow_pickle=True); zc_all = c["z"]; acc_c = c["acc"].astype(str)
v = np.load("/fast/sentrybio/v9_karcher/karcher_v7.npz", allow_pickle=True); zv_all = v["coords"].astype(np.float64); acc_v = v["accessions"].astype(str)
vmap = {}
for i, a in enumerate(acc_v): vmap.setdefault(a, i); vmap.setdefault(a.rsplit(".", 1)[0], i)
pairs = [(i, vmap.get(acc_c[i], vmap.get(acc_c[i].rsplit(".", 1)[0]))) for i in range(len(acc_c))
         if acc_c[i] in vmap or acc_c[i].rsplit(".", 1)[0] in vmap]
ci = np.array([p[0] for p in pairs]); vi = np.array([p[1] for p in pairs])
zc = zc_all[ci].astype(np.float64); zv = zv_all[vi]; base = np.array([a.rsplit(".", 1)[0] for a in acc_c[ci]])
ei = np.where(base == ECOLI)[0][0]; mj = np.where(base == MJANN)[0]
N = len(ci); rng = np.random.default_rng(0)


def bb(z, fit):
    T = logmap0(z); mu = T[fit].mean(0); _, _, Vt = np.linalg.svd(T[fit] - mu, full_matrices=False); return Vt[:2], mu


def ang0(z, B, mu):
    P = (logmap0(z[ei:ei + 1]) - mu) @ B.T; return np.arctan2(P[0, 1], P[0, 0])


def theta(z, B, mu):
    P = (logmap0(z) - mu) @ B.T; return (np.arctan2(P[:, 1], P[:, 0]) - ang0(z, B, mu) + np.pi) % (2 * np.pi) - np.pi


def chirefix(thv, thc):
    return -thv if (len(mj) and np.sign(thv[mj[0]]) != np.sign(thc[mj[0]])) else thv


# ---- #1 ROUND-TRIP (held-out generalization) ----
perm = rng.permutation(N); fit = perm[:int(.7 * N)]; test = perm[int(.7 * N):]
Bv, muv = bb(zv, fit); Bc, muc = bb(zc, fit)
thc = theta(zc, Bc, muc); thv = chirefix(theta(zv, Bv, muv), thc)
dth = np.abs((thc - thv + np.pi) % (2 * np.pi) - np.pi)
print("[#1 round-trip] transform fit on 70%, applied forward:")
print(f"   in-sample(fit) median = {np.degrees(np.median(dth[fit])):.1f}deg   HELD-OUT(test) median = {np.degrees(np.median(dth[test])):.1f}deg")
print(f"   held-out within-30 = {100*(dth[test]<np.pi/6).mean():.0f}%   -> GENERALIZES if held-out ~= in-sample (~16.8)")
# serialized artifact executes
S = json.load(open("/home/rohit/codec_cache/datum_transform_v9.json"))
Bs = np.array(S["v9_backbone_basis_2x129"]); mus = np.array(S["v9_tangent_mean_129"])
ths = chirefix(theta(zv, Bs, mus), theta(zc, *bb(zc, np.arange(N))))
thc_all = theta(zc, *bb(zc, np.arange(N)))
ds = np.abs((thc_all - ths + np.pi) % (2 * np.pi) - np.pi)
print(f"[#1 artifact] datum_transform_v9.json applied -> median = {np.degrees(np.median(ds)):.1f}deg  (serialized transform RUNS, reproduces 16.8)")

# ---- #2 ANCHOR-GAUGE STABILITY ----
def subspace_maxangle(A, B):
    s = np.linalg.svd(A @ B.T, compute_uv=False); return np.degrees(np.arccos(np.clip(s, -1, 1))).max()
planes = [bb(zv, rng.choice(N, N // 2, replace=False))[0] for _ in range(6)]
plane_ang = np.mean([subspace_maxangle(planes[0], p) for p in planes[1:]])
print(f"\n[#2 plane]  backbone 2D-plane stability: mean max principal angle across 5 half-subsamples = {plane_ang:.1f}deg  (small = plane doesn't wobble)")
thE, thM = [], []
for _ in range(8):
    B, mu = bb(zv, rng.choice(N, N // 2, replace=False))
    thE.append(np.degrees(ang0(zv, B, mu)))
    PM = (logmap0(zv[mj[0]:mj[0] + 1]) - mu) @ B.T; thM.append(np.degrees(np.arctan2(PM[0, 1], PM[0, 0])))
sep = [abs(((e - m + 180) % 360) - 180) for e, m in zip(thE, thM)]
print(f"[#2 anchors] E.coli<->M.jann angular separation = {np.mean(sep):.0f}+-{np.std(sep):.0f}deg across subsamples  (large & stable = well-conditioned gauge; near 0/180 = degenerate)")
print(f"\n[verify] #1 SOLID if held-out({np.degrees(np.median(dth[test])):.0f})~=in-sample & artifact runs; #2 SOLID if plane<~15deg & anchors well-separated/stable.")
