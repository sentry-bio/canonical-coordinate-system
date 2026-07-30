"""
eff_dim.py — EFFECTIVE / INTRINSIC dimensionality of a codec coordinate cloud.
The direct practical test of the n=2 thesis: how many dimensions does the learned
point cloud actually occupy, regardless of the ambient latent size?

Reports three complementary reads:
  1. PCA participation ratio + variance thresholds  (linear effective dim, tangent space)
  2. TwoNN intrinsic dim (Facco et al. 2017)         (nonlinear, distance-ratio MLE)
     computed with BOTH Euclidean and the proper hyperbolic (Poincaré, κ=5/4) metric,
     since the coords live in the Poincaré ball and Euclidean NN would misread curvature.
If PR and TwoNN both land ~2-3, that IS the n=2 thesis in the geometry the codec learned.
numpy-only; runs anywhere (CPU, minutes).

  python eff_dim.py coords.npz
"""
import sys, numpy as np
KAPPA = 1.25
BALLR = 1.0 / KAPPA ** 0.5


def load_z(path):
    d = np.load(path, allow_pickle=True)
    for k in ("z", "coords", "zQ", "reps", "embeddings"):
        if k in d:
            return np.asarray(d[k], float)
    raise SystemExit(f"no coordinate key in {path}; keys={list(d.keys())}")


def log_map_zero(x):
    """Poincaré log map at origin -> tangent (Euclidean) space, so PCA is geometrically honest."""
    n = np.linalg.norm(x, axis=1, keepdims=True).clip(1e-9, BALLR * (1 - 1e-6))
    return (2.0 / KAPPA ** 0.5) * np.arctanh(KAPPA ** 0.5 * n) * (x / n)


def pca_report(Z):
    X = log_map_zero(Z); X = X - X.mean(0)
    ev = np.sort(np.linalg.eigvalsh(np.cov(X.T)))[::-1]
    ev = ev[ev > ev.max() * 1e-12]
    cum = np.cumsum(ev) / ev.sum()
    pr = float((ev.sum() ** 2) / (ev ** 2).sum())          # participation ratio = effective #dims
    thr = {p: int(np.searchsorted(cum, p) + 1) for p in (0.90, 0.95, 0.99)}
    return ev, cum, pr, thr


def _poincare(A, B):
    diff = ((A - B) ** 2).sum(-1)
    na = 1 - KAPPA * (A ** 2).sum(-1); nb = 1 - KAPPA * (B ** 2).sum(-1)
    return (1 / KAPPA ** 0.5) * np.arccosh(np.clip(1 + 2 * KAPPA * diff / (na * nb).clip(1e-12, None), 1, None))


def twonn(Z, metric="euclid", n=4000, seed=0, chunk=64):
    rng = np.random.default_rng(seed)
    P = Z[rng.permutation(len(Z))[:min(n, len(Z))]]
    mus = []
    for i in range(0, len(P), chunk):
        c = P[i:i + chunk]
        D = np.sqrt(((c[:, None, :] - Z[None, :, :]) ** 2).sum(-1)) if metric == "euclid" \
            else _poincare(c[:, None, :], Z[None, :, :])
        D.sort(axis=1)
        m = D[:, 2] / np.maximum(D[:, 1], 1e-12)            # r2/r1 (skip self at col 0)
        mus.append(m[m > 1 + 1e-9])
    mu = np.concatenate(mus)
    return float(len(mu) / np.log(mu).sum())               # TwoNN MLE: d = N / Σ ln μ


if __name__ == "__main__":
    path = sys.argv[1] if len(sys.argv) > 1 else "coords_final.npz"
    Z = load_z(path)
    print(f"[effdim] {path}: {len(Z)} points in ambient dim {Z.shape[1]} (max||z||={np.linalg.norm(Z,axis=1).max():.3f}, ballR={BALLR:.3f})")
    ev, cum, pr, thr = pca_report(Z)
    print(f"[effdim] PCA (tangent): participation ratio = {pr:.2f}  |  dims for var: 90%={thr[0.90]} 95%={thr[0.95]} 99%={thr[0.99]}")
    print(f"[effdim]   top-8 variance share: {np.round(ev[:8]/ev.sum(),3)}")
    de = twonn(Z, "euclid"); dh = twonn(Z, "poincare")
    print(f"[effdim] TwoNN intrinsic dim:  Euclidean={de:.2f}   Poincaré(κ=5/4)={dh:.2f}")
    print(f"[effdim] READ: effective dim ≈ {np.mean([pr, dh]):.1f}  (PR={pr:.1f}, TwoNN-hyp={dh:.1f}); n=2 thesis wants ~2-3")
