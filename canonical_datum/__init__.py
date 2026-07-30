"""
canonical_datum — the frozen coordinate FRAME for the vertical history of life.

The core of the canonical coordinate system: a theorem-grounded datum (curvature, dimension, anchors),
the cross-instrument invariance that defines canonicity, and the registration transform that expresses
any atlas realization in the shared frame. The retrieval engine in `canonical_hybrid` is the operational
layer that *consumes* this frame.

Status: VALIDATED CANDIDATE. See docs/EVIDENCE.md for the honest ledger (what is verified vs. open).
"""
from . import datum
from .datum import (KAPPA, INTRINSIC_DIM, BALL_RADIUS, PRIME_MERIDIAN, CHIRALITY_ANCHOR,
                    REFERENCE_TREE, DATUM_VERSION, CONFORMANCE_TOLERANCE, summary)
from .registration import load_transform, to_datum_theta, logmap0
from .invariance import cross_instrument_agreement

__all__ = ["datum", "KAPPA", "INTRINSIC_DIM", "BALL_RADIUS", "PRIME_MERIDIAN", "CHIRALITY_ANCHOR",
           "REFERENCE_TREE", "DATUM_VERSION", "CONFORMANCE_TOLERANCE", "summary",
           "load_transform", "to_datum_theta", "logmap0", "cross_instrument_agreement"]
