"""
register_transform.py — FINISH the wiring: produce the actual published transform T: v9 -> datum-canonical.
The codec is the coarse canonical WITNESS; v9 is the operational realization. The wiring registers v9's
ANGULAR coordinate onto the anchored datum frame (θ shared, certifiable-candidate) and carries RADIAL as
advisory. v9's FINE structure (genus-resolution the dim-16 codec cannot hold) is v9's own job, not reconciled.
Emits datum_transform_v9.json: the 2D backbone basis + anchor-fixed θ0 + measured precision. Candidate, not
certified (16.8° median > freeze tolerance). This is the "published transform" the datum architecture specified.
"""
import numpy as np, csv, json
KAPPA = 1.25; BALLR = 1 / KAPPA ** 0.5
ECOLI = "GCF_000005845"; MJANN = "GCF_000091665"


def logmap0(x):
    n = np.linalg.norm(x, axis=1, keepdims=True).clip(1e-9, BALLR * (1 - 1e-6))
    return (2 / KAPPA ** 0.5) * np.arctanh(KAPPA ** 0.5 * n) * (x / n)


c = np.load("/home/rohit/codec_cache/coords_final.npz", allow_pickle=True)
zc_all = c["z"]; acc_c = c["acc"].astype(str)
v = np.load("/fast/sentrybio/v9_karcher/karcher_v7.npz", allow_pickle=True)
zv_all = v["coords"].astype(np.float64); acc_v = v["accessions"].astype(str)
vmap = {}
for i, a in enumerate(acc_v): vmap.setdefault(a, i); vmap.setdefault(a.rsplit(".", 1)[0], i)
pairs = [(i, vmap.get(acc_c[i], vmap.get(acc_c[i].rsplit(".", 1)[0]))) for i in range(len(acc_c))
         if acc_c[i] in vmap or acc_c[i].rsplit(".", 1)[0] in vmap]
ci = np.array([p[0] for p in pairs]); vi = np.array([p[1] for p in pairs])
zc = zc_all[ci].astype(np.float64); zv = zv_all[vi]; accs = acc_c[ci]
base = np.array([a.rsplit(".", 1)[0] for a in accs])


def backbone(z):                                          # -> 2D coords, PCA-2 basis, mean
    T = logmap0(z); mu = T.mean(0); _, _, Vt = np.linalg.svd(T - mu, full_matrices=False)
    return (T - mu) @ Vt[:2].T, Vt[:2], mu


Pc, _, _ = backbone(zc); Pv, Bv, muv = backbone(zv)


def anchored_theta(P, ecoli_i):                           # rotate so E.coli is at θ=0
    th = np.arctan2(P[:, 1], P[:, 0]); return (th - th[ecoli_i] + np.pi) % (2 * np.pi) - np.pi


ei = np.where(base == ECOLI)[0]
if not len(ei): raise SystemExit("E. coli anchor not in shared set")
ei = ei[0]
thc = anchored_theta(Pc, ei); thv = anchored_theta(Pv, ei)
# chirality: align handedness using M. jannaschii (flip v9 if it lands on the opposite side)
mj = np.where(base == MJANN)[0]
if len(mj) and np.sign(thv[mj[0]]) != np.sign(thc[mj[0]]):
    Pv[:, 1] *= -1; Bv[1] *= -1; thv = anchored_theta(Pv, ei)
dth = np.abs((thc - thv + np.pi) % (2 * np.pi) - np.pi)
rc = np.linalg.norm(zc, axis=1); rv = np.linalg.norm(zv, axis=1)
from scipy.stats import spearmanr
rho_r = spearmanr(rc, rv).correlation

print(f"[wire] {len(ci):,} shared genomes; E.coli=θ0 gauge, M.jann chirality-fixed")
print(f"[wire] datum-canonical θ agreement codec<->v9: median={np.degrees(np.median(dth)):.1f}deg  "
      f"within-30={100*(dth<np.pi/6).mean():.0f}%  within-45={100*(dth<np.pi/4).mean():.0f}%")
print(f"[wire] radial (advisory) Spearman = {rho_r:+.3f}")

T = {"transform": "v9 -> datum-canonical (angular)", "certified": False, "precision_note": "VALIDATED-CANDIDATE",
     "gauge": {"prime_meridian_anchor": ECOLI, "chirality_anchor": MJANN, "theta0_convention": "E.coli=0"},
     "v9_backbone_basis_2x129": Bv.tolist(), "v9_tangent_mean_129": muv.tolist(),
     "kappa": KAPPA, "ball_radius": BALLR,
     "measured_precision": {"angular_median_deg": float(np.degrees(np.median(dth))),
                            "within_30deg_frac": float((dth < np.pi / 6).mean()),
                            "radial_advisory_spearman": float(rho_r), "n_shared": int(len(ci))},
     "apply": "z_v9(ball,129) -> logmap0 -> minus v9_tangent_mean -> @ basis.T -> 2D -> atan2 - E.coli_angle = datum θ; r = ||z_v9||(advisory)",
     "scope": "CANONICAL BACKBONE ONLY (θ). Fine/genus resolution is v9's operational job, NOT registered.",
     "freeze_gate": "certify when angular median < ~10deg and determinate-quartet agreement >= 0.90"}
json.dump(T, open("/home/rohit/codec_cache/datum_transform_v9.json", "w"), indent=2)
print(f"[wire] wrote datum_transform_v9.json — the published v9->datum transform (candidate precision).")
print(f"[wire] FINISHED: codec=coarse witness / v9=operational; wired on θ (advisory r); fine=v9's own job.")
