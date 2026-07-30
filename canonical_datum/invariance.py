"""
Invariance — the operational definition of canonicity.

Canonical = INVARIANT ACROSS INSTRUMENTS, not the output of any one model. The canonical backbone is the
angular coordinate that two maximally-different instruments both recover. This module measures it: the
median datum-θ residual between two realizations of the same genomes, in the anchored frame.

Measured for the dim-16 low-bias codec vs the dim-129 v9 atlas: 16.8° median, 67% within 30°
(random ~90°). That agreement is what lets you *trust* the frame is real, not a v9 artifact.
"""
import numpy as np
from .datum import BALL_RADIUS, KAPPA


def _logmap0(x):
    n = np.linalg.norm(x, axis=1, keepdims=True).clip(1e-9, BALL_RADIUS * (1 - 1e-6))
    return (2.0 / np.sqrt(KAPPA)) * np.arctanh(np.sqrt(KAPPA) * n) * (x / n)


def _datum_theta(coords, ecoli_i, mjann_i):
    """coords: [N, D] aligned so row ecoli_i is E. coli, mjann_i is M. jannaschii. -> anchored θ."""
    T = _logmap0(np.asarray(coords, float)); mu = T.mean(0)
    _, _, Vt = np.linalg.svd(T - mu, full_matrices=False)        # 2D backbone = top-2 tangent components
    P = (T - mu) @ Vt[:2].T
    th = np.arctan2(P[:, 1], P[:, 0]); th = (th - th[ecoli_i] + np.pi) % (2 * np.pi) - np.pi
    if th[mjann_i] < 0: th = -th
    return th


def cross_instrument_agreement(coords_a, coords_b, ecoli_i, mjann_i):
    """
    coords_a, coords_b : [N, D_a], [N, D_b] — the SAME N genomes, in row order, from two instruments
                         (dimensions may differ). ecoli_i / mjann_i index the two anchors.
    Returns dict: angular median residual (deg), fraction within 30°, and radial-axis Spearman (advisory).
    """
    th_a = _datum_theta(coords_a, ecoli_i, mjann_i)
    th_b = _datum_theta(coords_b, ecoli_i, mjann_i)
    d = np.abs((th_a - th_b + np.pi) % (2 * np.pi) - np.pi)
    ra = np.linalg.norm(np.asarray(coords_a, float), axis=1)
    rb = np.linalg.norm(np.asarray(coords_b, float), axis=1)
    order = lambda x: np.argsort(np.argsort(x))
    rho = float(np.corrcoef(order(ra), order(rb))[0, 1])
    return {"angular_median_deg": float(np.degrees(np.median(d))),
            "within_30deg_frac": float((d < np.pi / 6).mean()),
            "radial_advisory_spearman": rho, "n": int(len(th_a)),
            "note": "canonical iff angular residual << random (90deg); radial is advisory only"}
