"""
canonical_hybrid — Multi-Resolution Canonical Coordinate System
================================================================

Provides three operating modes:

  1. KESTREL standalone: CPU-only canonical address via tangent-space
     Voronoi tessellation. Returns (cell_id, confidence, taxonomy).
     No Atlas required.

  2. Atlas standalone: Full-precision tangent-space embedding with
     genus-level resolution.

  3. Hybrid tandem: KESTREL assigns coarse cell → Atlas searches
     within cell for precise match. 76% R@1 at 32× speedup.

The canonical coordinate is defined by a tangent-space Voronoi
tessellation with K centers obtained by K-means clustering of
Atlas teacher embeddings. Each center carries taxonomic metadata
(dominant domain, family, genus) derived from the training set.

Usage:
    from canonical_hybrid import CanonicalHybrid

    hybrid = CanonicalHybrid.load("tessellation_K50.npz")

    # KESTREL standalone
    addr = hybrid.kestrel_address(tangent_vector)
    # → CanonicalAddress(cell=23, confidence=0.87, domain="Bacteria",
    #                    family="Bacillaceae", top3=[23, 17, 41])

    # Hybrid: KESTREL coarse → Atlas precise
    match = hybrid.hybrid_search(kestrel_tangent, atlas_tangent)
    # → HybridMatch(cell=23, nn_idx=14592, nn_dist=0.34,
    #               candidates_searched=1471, speedup=32x)
"""

from canonical_hybrid.core import (
    CanonicalHybrid,
    CanonicalAddress,
    HybridMatch,
    TessellationMeta,
)
from canonical_hybrid.build import build_tessellation
from canonical_hybrid.v2 import HybridEngine, CanonicalResult

__version__ = "3.0.0"
__all__ = [
    # v2 (production)
    "HybridEngine",
    "CanonicalResult",
    # v1 (experimental)
    "CanonicalHybrid",
    "CanonicalAddress",
    "HybridMatch",
    "TessellationMeta",
    "build_tessellation",
]
