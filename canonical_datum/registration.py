"""
Registration — express any atlas realization's coordinate in the datum-canonical angular frame.

This is the published transform (the WGS84 "the meter is *this*"): a realization (v9, v10.x, or a future
encoder) registers to the SAME frozen form via the two anchors, so coordinates become comparable across
instruments and versions. The ANGULAR coordinate (θ = lineage) is the certifiable-candidate output;
the RADIAL coordinate (r = depth) is carried as ADVISORY.

Reference implementation — demonstrates the mechanism. Verified as a generalizing function
(held-out θ agreement 16.9° ≈ in-sample 16.8°; see docs/EVIDENCE.md #1).
"""
import json
import numpy as np
from .datum import KAPPA, BALL_RADIUS


def logmap0(x):
    """Poincaré ball -> tangent space at the origin (Euclidean), so the 2D backbone is geometrically honest."""
    x = np.atleast_2d(np.asarray(x, float))
    n = np.linalg.norm(x, axis=1, keepdims=True).clip(1e-9, BALL_RADIUS * (1 - 1e-6))
    return (2.0 / np.sqrt(KAPPA)) * np.arctanh(np.sqrt(KAPPA) * n) * (x / n)


def load_transform(path):
    return json.load(open(path))


def to_datum_theta(atlas_coords, transform, ecoli_coord, mjann_coord=None):
    """
    atlas_coords : [N, D] Poincaré-ball coords from an atlas realization.
    transform    : dict from datum_transform_v9.json (backbone basis + tangent mean).
    ecoli_coord  : [D] the prime-meridian anchor's coord in the SAME realization (fixes θ = 0).
    mjann_coord  : [D] the chirality anchor's coord (fixes handedness); optional.
    Returns datum-canonical θ in radians (E. coli = 0), advisory r = ||coord||.
    """
    B = np.asarray(transform["v9_backbone_basis_2x129"], float)   # [2, D]
    mu = np.asarray(transform["v9_tangent_mean_129"], float)      # [D]
    P = (logmap0(atlas_coords) - mu) @ B.T                        # [N, 2]
    th = np.arctan2(P[:, 1], P[:, 0])
    th0 = np.arctan2(*((logmap0(ecoli_coord) - mu) @ B.T)[0][::-1])
    th = (th - th0 + np.pi) % (2 * np.pi) - np.pi                 # anchor E. coli at 0
    if mjann_coord is not None:
        thm = np.arctan2(*((logmap0(mjann_coord) - mu) @ B.T)[0][::-1]) - th0
        if ((thm + np.pi) % (2 * np.pi) - np.pi) < 0:             # fix chirality via M. jannaschii
            th = -th
    r = np.linalg.norm(np.atleast_2d(np.asarray(atlas_coords, float)), axis=1)   # ADVISORY
    return th, r
