"""
Core module: CanonicalHybrid — tangent-space Voronoi tessellation engine.
"""

from __future__ import annotations

import json
from dataclasses import dataclass, field
from pathlib import Path
from typing import Dict, List, Optional, Tuple, Union

import numpy as np


# ---------------------------------------------------------------------------
# Data classes
# ---------------------------------------------------------------------------
@dataclass
class TessellationMeta:
    """Metadata for a single Voronoi cell."""
    cell_id: int
    n_members: int                     # genomes in training set assigned here
    dominant_domain: str               # most common domain
    dominant_family: str               # most common family
    domain_purity: float               # fraction of members in dominant domain
    family_purity: float               # fraction of members in dominant family
    mean_radius: float                 # mean tangent-space norm of members
    center: np.ndarray                 # (D,) tangent-space Voronoi center

    def to_dict(self) -> dict:
        return {
            "cell_id": self.cell_id,
            "n_members": self.n_members,
            "dominant_domain": self.dominant_domain,
            "dominant_family": self.dominant_family,
            "domain_purity": round(self.domain_purity, 4),
            "family_purity": round(self.family_purity, 4),
            "mean_radius": round(self.mean_radius, 4),
        }


@dataclass
class CanonicalAddress:
    """Canonical coordinate for a genome embedding."""
    cell_id: int                       # primary Voronoi cell
    confidence: float                  # margin: (d_second - d_nearest) / d_nearest
    distance_to_center: float          # L2 to assigned center
    domain: str                        # taxonomy from cell metadata
    family: str                        # taxonomy from cell metadata
    top_k_cells: List[int] = field(default_factory=list)  # nearest K cells
    top_k_dists: List[float] = field(default_factory=list)
    within_cell_position: Optional[np.ndarray] = None  # offset from center

    def to_dict(self) -> dict:
        return {
            "cell_id": self.cell_id,
            "confidence": round(self.confidence, 4),
            "distance_to_center": round(self.distance_to_center, 4),
            "domain": self.domain,
            "family": self.family,
            "top_k_cells": self.top_k_cells[:5],
            "top_k_dists": [round(d, 4) for d in self.top_k_dists[:5]],
        }


@dataclass
class HybridMatch:
    """Result of hybrid coarse→precise search."""
    cell_id: int                       # KESTREL-assigned coarse cell
    nn_idx: int                        # index in reference set of nearest match
    nn_dist: float                     # L2 distance to nearest match
    nn_accession: Optional[str]        # accession of nearest match
    candidates_searched: int           # how many references Atlas checked
    total_references: int              # total reference set size
    speedup: float                     # total / candidates_searched
    r_at_k: Dict[int, bool] = field(default_factory=dict)  # R@1, R@5, R@10

    def to_dict(self) -> dict:
        return {
            "cell_id": self.cell_id,
            "nn_idx": self.nn_idx,
            "nn_dist": round(self.nn_dist, 4),
            "nn_accession": self.nn_accession,
            "candidates_searched": self.candidates_searched,
            "speedup": round(self.speedup, 1),
        }


# ---------------------------------------------------------------------------
# Main class
# ---------------------------------------------------------------------------
class CanonicalHybrid:
    """
    Multi-resolution canonical coordinate system.

    Built from a tangent-space Voronoi tessellation of Atlas teacher embeddings.
    Supports three modes:
      - kestrel_address(): KESTREL standalone canonical coordinate
      - atlas_address(): Atlas standalone coordinate (same API, higher precision)
      - hybrid_search(): KESTREL coarse → Atlas precise within cell
    """

    def __init__(
        self,
        centers: np.ndarray,               # (K, D) Voronoi centers in tangent space
        cell_meta: List[TessellationMeta],  # metadata per cell
        reference_tangents: np.ndarray,     # (N, D) full teacher tangent vectors
        reference_cells: np.ndarray,        # (N,) cell assignment per reference
        reference_domains: Optional[np.ndarray] = None,
        reference_families: Optional[np.ndarray] = None,
        reference_accessions: Optional[np.ndarray] = None,
        kappa: float = 1.25,
    ):
        self.centers = centers.astype(np.float32)
        self.cell_meta = cell_meta
        self.reference_tangents = reference_tangents.astype(np.float32)
        self.reference_cells = reference_cells.astype(np.int32)
        self.reference_domains = reference_domains
        self.reference_families = reference_families
        self.reference_accessions = reference_accessions
        self.kappa = kappa
        self.K = centers.shape[0]
        self.N = reference_tangents.shape[0]
        self.D = centers.shape[1]

        # Precompute per-cell member indices for fast lookup
        self._cell_members: Dict[int, np.ndarray] = {}
        for k in range(self.K):
            self._cell_members[k] = np.where(reference_cells == k)[0]

        # Precompute center norms for fast distance
        self._center_norms_sq = (centers ** 2).sum(axis=1)  # (K,)

    # -----------------------------------------------------------------
    # KESTREL standalone: canonical address
    # -----------------------------------------------------------------
    def kestrel_address(self, tangent: np.ndarray, top_k: int = 5) -> CanonicalAddress:
        """
        Compute canonical address for a single tangent vector (from KESTREL or Atlas).

        Args:
            tangent: (D,) tangent-space vector
            top_k: number of nearest cells to return
        """
        tangent = tangent.astype(np.float32).ravel()
        dists = self._distances_to_centers(tangent)
        sorted_idx = np.argsort(dists)
        nearest = sorted_idx[0]
        d_nearest = dists[nearest]
        d_second = dists[sorted_idx[1]]

        confidence = float((d_second - d_nearest) / (d_nearest + 1e-8))
        meta = self.cell_meta[nearest]

        return CanonicalAddress(
            cell_id=int(nearest),
            confidence=confidence,
            distance_to_center=float(d_nearest),
            domain=meta.dominant_domain,
            family=meta.dominant_family,
            top_k_cells=[int(i) for i in sorted_idx[:top_k]],
            top_k_dists=[float(dists[i]) for i in sorted_idx[:top_k]],
            within_cell_position=tangent - self.centers[nearest],
        )

    def kestrel_address_batch(self, tangents: np.ndarray, top_k: int = 5) -> List[CanonicalAddress]:
        """Batch canonical address computation."""
        tangents = tangents.astype(np.float32)
        N = tangents.shape[0]

        # Vectorized distances to all centers
        dists = self._distances_to_centers_batch(tangents)  # (N, K)
        sorted_idx = np.argsort(dists, axis=1)  # (N, K)

        results = []
        for i in range(N):
            nearest = sorted_idx[i, 0]
            d_nearest = dists[i, nearest]
            d_second = dists[i, sorted_idx[i, 1]]
            confidence = float((d_second - d_nearest) / (d_nearest + 1e-8))
            meta = self.cell_meta[nearest]

            results.append(CanonicalAddress(
                cell_id=int(nearest),
                confidence=confidence,
                distance_to_center=float(d_nearest),
                domain=meta.dominant_domain,
                family=meta.dominant_family,
                top_k_cells=[int(j) for j in sorted_idx[i, :top_k]],
                top_k_dists=[float(dists[i, j]) for j in sorted_idx[i, :top_k]],
                within_cell_position=tangents[i] - self.centers[nearest],
            ))
        return results

    # -----------------------------------------------------------------
    # Atlas standalone: same API, higher precision
    # -----------------------------------------------------------------
    def atlas_address(self, tangent: np.ndarray, top_k: int = 5) -> CanonicalAddress:
        """Identical to kestrel_address — same canonical frame, different input model."""
        return self.kestrel_address(tangent, top_k)

    # -----------------------------------------------------------------
    # Hybrid: KESTREL coarse → Atlas precise
    # -----------------------------------------------------------------
    def hybrid_search(
        self,
        kestrel_tangent: np.ndarray,
        atlas_tangent: np.ndarray,
        n_cells: int = 1,
        correct_idx: Optional[int] = None,
    ) -> HybridMatch:
        """
        KESTREL assigns coarse cell, Atlas searches within it.

        Args:
            kestrel_tangent: (D,) KESTREL tangent vector (for cell assignment)
            atlas_tangent:   (D,) Atlas tangent vector (for precise search)
            n_cells:         number of KESTREL cells to search (1=primary, 3=top-3)
            correct_idx:     if known, check whether correct match is found
        """
        kestrel_tangent = kestrel_tangent.astype(np.float32).ravel()
        atlas_tangent = atlas_tangent.astype(np.float32).ravel()

        # KESTREL assigns coarse cell(s)
        dists_to_centers = self._distances_to_centers(kestrel_tangent)
        coarse_cells = np.argsort(dists_to_centers)[:n_cells]

        # Gather all reference members in the selected cells
        candidate_indices = np.concatenate(
            [self._cell_members[int(c)] for c in coarse_cells]
        )
        n_candidates = len(candidate_indices)

        if n_candidates == 0:
            return HybridMatch(
                cell_id=int(coarse_cells[0]),
                nn_idx=-1, nn_dist=float("inf"), nn_accession=None,
                candidates_searched=0, total_references=self.N,
                speedup=float("inf"),
            )

        # Atlas L2-NN within candidates
        candidate_vecs = self.reference_tangents[candidate_indices]  # (M, D)
        dists = np.linalg.norm(candidate_vecs - atlas_tangent, axis=1)
        best_local = dists.argmin()
        nn_idx = int(candidate_indices[best_local])
        nn_dist = float(dists[best_local])

        nn_acc = None
        if self.reference_accessions is not None:
            nn_acc = str(self.reference_accessions[nn_idx])

        r_at_k = {}
        if correct_idx is not None:
            sorted_local = np.argsort(dists)
            global_sorted = candidate_indices[sorted_local]
            for k in [1, 5, 10]:
                r_at_k[k] = bool(correct_idx in global_sorted[:k])

        return HybridMatch(
            cell_id=int(coarse_cells[0]),
            nn_idx=nn_idx,
            nn_dist=nn_dist,
            nn_accession=nn_acc,
            candidates_searched=n_candidates,
            total_references=self.N,
            speedup=float(self.N / max(n_candidates, 1)),
            r_at_k=r_at_k,
        )

    def hybrid_search_batch(
        self,
        kestrel_tangents: np.ndarray,
        atlas_tangents: np.ndarray,
        n_cells: int = 1,
        correct_indices: Optional[np.ndarray] = None,
    ) -> List[HybridMatch]:
        """Batch hybrid search."""
        N = kestrel_tangents.shape[0]
        results = []
        for i in range(N):
            ci = int(correct_indices[i]) if correct_indices is not None else None
            results.append(self.hybrid_search(
                kestrel_tangents[i], atlas_tangents[i],
                n_cells=n_cells, correct_idx=ci,
            ))
        return results

    # -----------------------------------------------------------------
    # Distance utilities
    # -----------------------------------------------------------------
    def _distances_to_centers(self, v: np.ndarray) -> np.ndarray:
        """L2 distances from vector v to all K centers. Returns (K,)."""
        # ||v - c||^2 = ||v||^2 + ||c||^2 - 2*v.c
        v_sq = (v ** 2).sum()
        dots = self.centers @ v  # (K,)
        dists_sq = v_sq + self._center_norms_sq - 2 * dots
        return np.sqrt(np.maximum(dists_sq, 0))

    def _distances_to_centers_batch(self, V: np.ndarray) -> np.ndarray:
        """L2 distances from (N, D) to all K centers. Returns (N, K)."""
        V_sq = (V ** 2).sum(axis=1, keepdims=True)  # (N, 1)
        dots = V @ self.centers.T  # (N, K)
        dists_sq = V_sq + self._center_norms_sq[np.newaxis, :] - 2 * dots
        return np.sqrt(np.maximum(dists_sq, 0))

    # -----------------------------------------------------------------
    # Cell inspection
    # -----------------------------------------------------------------
    def cell_info(self, cell_id: int) -> Dict:
        """Inspect a Voronoi cell."""
        meta = self.cell_meta[cell_id]
        members = self._cell_members[cell_id]
        return {
            **meta.to_dict(),
            "member_indices": members.tolist()[:20],  # first 20
            "total_members_in_reference": len(members),
        }

    def cell_taxonomy_summary(self) -> List[Dict]:
        """Summary of all cells with taxonomy info."""
        return [m.to_dict() for m in self.cell_meta]

    # -----------------------------------------------------------------
    # Serialization
    # -----------------------------------------------------------------
    def save(self, path: Union[str, Path]) -> None:
        """Save tessellation to a .npz + .json pair."""
        path = Path(path)
        path.parent.mkdir(parents=True, exist_ok=True)

        # Numpy arrays
        np.savez_compressed(
            path.with_suffix(".npz"),
            centers=self.centers,
            reference_tangents=self.reference_tangents,
            reference_cells=self.reference_cells,
            reference_domains=self.reference_domains if self.reference_domains is not None else np.array([]),
            reference_families=self.reference_families if self.reference_families is not None else np.array([]),
        )

        # Metadata JSON
        meta_json = {
            "version": "0.1.0",
            "K": self.K,
            "N": self.N,
            "D": self.D,
            "kappa": self.kappa,
            "cells": [m.to_dict() for m in self.cell_meta],
        }
        if self.reference_accessions is not None:
            meta_json["accessions"] = self.reference_accessions.tolist()

        with open(path.with_suffix(".json"), "w") as f:
            json.dump(meta_json, f, indent=2)

    @classmethod
    def load(cls, path: Union[str, Path]) -> "CanonicalHybrid":
        """Load tessellation from .npz + .json pair."""
        path = Path(path)
        npz_path = path.with_suffix(".npz")
        json_path = path.with_suffix(".json")

        data = np.load(npz_path, allow_pickle=False)
        with open(json_path) as f:
            meta_json = json.load(f)

        centers = data["centers"]
        ref_tangents = data["reference_tangents"]
        ref_cells = data["reference_cells"]
        ref_domains = data["reference_domains"] if data["reference_domains"].size > 0 else None
        ref_families = data["reference_families"] if data["reference_families"].size > 0 else None

        accessions = None
        if "accessions" in meta_json:
            accessions = np.array(meta_json["accessions"])

        cell_meta = []
        for i, c in enumerate(meta_json["cells"]):
            cell_meta.append(TessellationMeta(
                cell_id=c["cell_id"],
                n_members=c["n_members"],
                dominant_domain=c["dominant_domain"],
                dominant_family=c["dominant_family"],
                domain_purity=c["domain_purity"],
                family_purity=c["family_purity"],
                mean_radius=c["mean_radius"],
                center=centers[i],
            ))

        return cls(
            centers=centers,
            cell_meta=cell_meta,
            reference_tangents=ref_tangents,
            reference_cells=ref_cells,
            reference_domains=ref_domains,
            reference_families=ref_families,
            reference_accessions=accessions,
            kappa=meta_json.get("kappa", 1.25),
        )

    def __repr__(self) -> str:
        return (f"CanonicalHybrid(K={self.K}, N={self.N}, D={self.D}, "
                f"kappa={self.kappa})")
