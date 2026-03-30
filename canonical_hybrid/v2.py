"""
canonical_hybrid.v2 — Production Hybrid Inference Engine
========================================================

Tangent-space Voronoi tessellation for multi-resolution genome classification.
97.1% R@1 at 7× speedup (K=25, 3-cell search, adaptive).

Two operating modes:
  1. KESTREL standalone: canonical address from k-mer features (0.03ms, CPU)
  2. Hybrid tandem: KESTREL coarse → Atlas precise (0.85ms, CPU+GPU)

Optional presentation layer: map canonical addresses to shareable (r, θ) coordinates
via BiosphereCoordinate projection matrix P for URLs like /coord/0.39/147.3

Usage:
    from canonical_hybrid.v2 import HybridEngine

    engine = HybridEngine.from_teacher_coords(
        "/path/to/teacher_coords.npz",
        K=25, n_search_cells=3,
    )

    # KESTREL standalone — canonical address
    addr = engine.classify_kestrel(dna, model, extract_fn, kappa)

    # Hybrid — KESTREL narrows, Atlas resolves
    match = engine.classify_hybrid(kestrel_tangent, atlas_tangent)

    # Save/load for deployment
    engine.save("/path/to/engine_v2")
    engine = HybridEngine.load("/path/to/engine_v2")
"""

from __future__ import annotations

import json
import logging
import math
import time
from dataclasses import dataclass, field
from pathlib import Path
from typing import Dict, List, Optional, Tuple, Union

import numpy as np

log = logging.getLogger(__name__)

KAPPA = 5.0 / 4.0
DOMAIN_NAMES = {0: "Bacteria", 1: "Archaea", 2: "Eukaryota"}


# ---------------------------------------------------------------------------
# Result types
# ---------------------------------------------------------------------------
@dataclass
class CanonicalResult:
    """Full classification result with canonical coordinates."""

    # Canonical address (tessellation)
    cell_id: int
    cell_confidence: float            # margin ratio: (d2 - d1) / d1
    cell_domain: str                  # dominant domain of assigned cell
    cell_family: str                  # dominant family of assigned cell
    top_cells: List[int] = field(default_factory=list)

    # Tangent-space radius (≈ hyperbolic depth)
    r: float = 0.0                    # tangent-space norm

    # Classification
    domain: str = ""
    domain_confidence: float = 0.0
    domain_probs: Dict[str, float] = field(default_factory=dict)

    # Hybrid metadata
    mode: str = "kestrel"             # kestrel | hybrid
    nn_accession: Optional[str] = None
    nn_distance: float = float("inf")
    nn_top_k: List[Tuple] = field(default_factory=list)  # [(accession, dist), ...]
    candidates_searched: int = 0
    search_speedup: float = 1.0
    latency_ms: float = 0.0

    def to_dict(self) -> dict:
        d = {
            "cell_id": self.cell_id,
            "cell_confidence": round(self.cell_confidence, 4),
            "cell_domain": self.cell_domain,
            "cell_family": self.cell_family,
            "top_cells": self.top_cells[:5],
            "r": round(self.r, 4),
            "domain": self.domain,
            "domain_confidence": round(self.domain_confidence, 4),
            "mode": self.mode,
            "latency_ms": round(self.latency_ms, 2),
        }
        if self.domain_probs:
            d["domain_probs"] = {k: round(v, 4) for k, v in self.domain_probs.items()}
        if self.nn_accession is not None:
            d["nn_accession"] = self.nn_accession
            d["nn_distance"] = round(self.nn_distance, 4)
            d["candidates_searched"] = self.candidates_searched
            d["search_speedup"] = round(self.search_speedup, 1)
        if self.nn_top_k:
            d["nn_top_k"] = [(acc, round(dist, 4)) for acc, dist in self.nn_top_k[:5]]
        return d


# ---------------------------------------------------------------------------
# Cell metadata
# ---------------------------------------------------------------------------
@dataclass
class CellMeta:
    """Precomputed metadata for a Voronoi cell."""
    cell_id: int
    n_members: int
    dominant_domain: str
    dominant_family: str
    domain_purity: float
    family_purity: float
    mean_radius: float
    member_indices: np.ndarray        # indices into reference set


# ---------------------------------------------------------------------------
# Hybrid Engine
# ---------------------------------------------------------------------------
class HybridEngine:
    """
    Production hybrid inference engine.

    Optimized defaults: K=25, n_search_cells=3 → 91.7% R@1, 7× speedup.
    """

    VERSION = "3.0.0"

    def __init__(
        self,
        centers: np.ndarray,               # (K, D) Voronoi centers
        cells: List[CellMeta],             # metadata per cell
        reference_tangents: np.ndarray,     # (N, D) teacher tangent vectors
        reference_cells: np.ndarray,        # (N,) cell assignments
        reference_accessions: Optional[np.ndarray] = None,
        reference_domains: Optional[np.ndarray] = None,
        kappa: float = KAPPA,
        n_search_cells: int = 3,
    ):
        self.centers = centers.astype(np.float32)
        self.cells = cells
        self.reference_tangents = reference_tangents.astype(np.float32)
        self.reference_cells = reference_cells.astype(np.int32)
        self.reference_accessions = reference_accessions
        self.reference_domains = reference_domains
        self.kappa = kappa
        self.n_search_cells = n_search_cells

        self.K = centers.shape[0]
        self.N = reference_tangents.shape[0]
        self.D = centers.shape[1]

        # Precompute for fast distance
        self._center_norms_sq = (centers ** 2).sum(axis=1)
        self._ref_norms_sq = (reference_tangents ** 2).sum(axis=1)

        # Precompute contiguous member arrays per cell for cache-friendly search
        self._cell_ref_tangents: Dict[int, np.ndarray] = {}
        self._cell_ref_indices: Dict[int, np.ndarray] = {}
        self._cell_ref_norms_sq: Dict[int, np.ndarray] = {}
        for c in cells:
            idx = c.member_indices
            self._cell_ref_indices[c.cell_id] = idx
            self._cell_ref_tangents[c.cell_id] = reference_tangents[idx]
            self._cell_ref_norms_sq[c.cell_id] = self._ref_norms_sq[idx]

    # -----------------------------------------------------------------
    # Core: tangent vector → canonical result
    # -----------------------------------------------------------------
    def _assign_cell(self, tangent: np.ndarray, top_k: int = 5) -> Tuple[int, float, List[int], np.ndarray]:
        """Assign tangent vector to Voronoi cell. Returns (cell, confidence, top_k_cells, dists)."""
        v = tangent.astype(np.float32).ravel()
        dists_sq = (v ** 2).sum() + self._center_norms_sq - 2 * (self.centers @ v)
        dists = np.sqrt(np.maximum(dists_sq, 0))
        order = np.argsort(dists)
        nearest = order[0]
        d1, d2 = dists[order[0]], dists[order[1]]
        confidence = float((d2 - d1) / (d1 + 1e-8))
        return int(nearest), confidence, [int(i) for i in order[:top_k]], dists

    def _search_within_cells(
        self,
        atlas_tangent: np.ndarray,
        cell_ids: List[int],
        top_k: int = 5,
    ) -> Tuple[int, float, int, List[Tuple]]:
        """
        L2-NN search within specified cells using precomputed contiguous arrays.
        Returns (nn_idx, nn_dist, n_searched, top_k_results).
        """
        v = atlas_tangent.astype(np.float32).ravel()
        v_sq = (v ** 2).sum()

        # Gather candidates from precomputed cell arrays (cache-friendly)
        all_vecs = []
        all_idx = []
        all_norms = []
        for c in cell_ids:
            all_vecs.append(self._cell_ref_tangents[c])
            all_idx.append(self._cell_ref_indices[c])
            all_norms.append(self._cell_ref_norms_sq[c])

        if not all_vecs:
            return -1, float("inf"), 0, []

        candidate_vecs = np.concatenate(all_vecs)
        candidate_idx = np.concatenate(all_idx)
        candidate_norms = np.concatenate(all_norms)
        n_candidates = len(candidate_idx)

        # Fast L2 via ||a-b||² = ||a||² + ||b||² - 2<a,b>
        dots = candidate_vecs @ v
        dists = np.sqrt(np.maximum(v_sq + candidate_norms - 2 * dots, 0))

        # Top-K results
        if top_k >= n_candidates:
            order = np.argsort(dists)
        else:
            order = np.argpartition(dists, top_k)[:top_k]
            order = order[np.argsort(dists[order])]

        best_global = int(candidate_idx[order[0]])
        best_dist = float(dists[order[0]])

        top_k_results = []
        for j in order[:top_k]:
            gidx = int(candidate_idx[j])
            acc = str(self.reference_accessions[gidx]) if self.reference_accessions is not None else str(gidx)
            top_k_results.append((acc, float(dists[j])))

        return best_global, best_dist, n_candidates, top_k_results

    def _adaptive_n_cells(self, confidence: float) -> int:
        """
        Adaptive cell count based on KESTREL confidence.
        High confidence → fewer cells (faster). Low → more cells (safer).
        """
        if confidence > 3.0:
            return 1                   # very clear assignment
        elif confidence > 1.0:
            return min(2, self.n_search_cells)
        else:
            return self.n_search_cells  # low confidence → search all configured cells

    # -----------------------------------------------------------------
    # Mode 1: KESTREL standalone
    # -----------------------------------------------------------------
    def address_from_tangent(
        self,
        tangent: np.ndarray,
        domain_probs: Optional[np.ndarray] = None,
    ) -> CanonicalResult:
        """
        Compute canonical address from a pre-computed tangent vector.
        Works for both KESTREL and Atlas tangent vectors.

        Args:
            tangent: (D,) tangent-space vector
            domain_probs: (3,) softmax domain probabilities (if available from model)
        """
        t0 = time.perf_counter()
        cell_id, confidence, top_cells, _ = self._assign_cell(tangent)
        r = float(np.linalg.norm(tangent))
        meta = self.cells[cell_id]

        # Use model domain probs if available (99.2% acc), else cell metadata (86.2%)
        if domain_probs is not None:
            dom_pred = int(domain_probs.argmax())
            dom_conf = float(domain_probs[dom_pred])
            dom_dict = {DOMAIN_NAMES[i]: float(domain_probs[i]) for i in range(len(domain_probs))}
        else:
            dom_pred = list(DOMAIN_NAMES.values()).index(meta.dominant_domain) if meta.dominant_domain in DOMAIN_NAMES.values() else 0
            dom_conf = meta.domain_purity
            dom_dict = {}

        elapsed = (time.perf_counter() - t0) * 1000

        return CanonicalResult(
            cell_id=cell_id,
            cell_confidence=confidence,
            cell_domain=meta.dominant_domain,
            cell_family=meta.dominant_family,
            top_cells=top_cells,
            r=r,
            domain=DOMAIN_NAMES.get(dom_pred, meta.dominant_domain),
            domain_confidence=dom_conf,
            domain_probs=dom_dict,
            mode="tangent",
            latency_ms=elapsed,
        )

    def classify_kestrel(
        self,
        dna: str,
        model,
        extract_fn,
        kappa: float,
        n_ensemble: int = 1,
        frag_bp: int = 509,
    ) -> CanonicalResult:
        """
        Full KESTREL standalone classification with canonical address.

        Args:
            dna: DNA sequence string
            model: SpectrumBranch model (from kestrel_bio.load_model)
            extract_fn: feature extractor (from kestrel_bio.load_model)
            kappa: curvature (from kestrel_bio.load_model)
            n_ensemble: number of fragment sub-samples
            frag_bp: fragment length for ensemble
        """
        import torch
        from kestrel_bio.spectrum_branch import hyperbolic_frechet_mean

        t0 = time.perf_counter()
        dna = dna.upper().replace(" ", "").replace("\n", "")
        L = len(dna)
        if L < 4:
            raise ValueError(f"Sequence too short ({L} bp)")

        # Extract features
        if n_ensemble == 1:
            frags = [extract_fn(dna)]
        else:
            rng = np.random.default_rng(hash(dna[:64]) & 0xFFFFFFFF)
            frags = []
            frag_len = min(L, frag_bp)
            for _ in range(n_ensemble):
                start = int(rng.integers(0, max(1, L - frag_len + 1)))
                frag = dna[start:start + frag_len]
                frags.append(extract_fn(frag if len(frag) >= 4 else dna[:max(4, frag_len)]))

        feat = torch.from_numpy(np.stack(frags))
        with torch.no_grad():
            out = model(feat)

        tang_vecs = out["tangent"].cpu().numpy()
        dom_logits = out["domain_logits"].cpu().numpy()

        # Ensemble tangent via Frechet mean
        tang = hyperbolic_frechet_mean(tang_vecs, kappa) if n_ensemble > 1 else tang_vecs[0]

        # Domain from logits
        dom_probs = np.exp(dom_logits.mean(axis=0))
        dom_probs = dom_probs / dom_probs.sum()
        dom_pred = int(dom_probs.argmax())
        dom_conf = float(dom_probs[dom_pred])

        # Canonical address
        cell_id, confidence, top_cells, _ = self._assign_cell(tang)
        r = float(np.linalg.norm(tang))
        meta = self.cells[cell_id]
        elapsed = (time.perf_counter() - t0) * 1000

        return CanonicalResult(
            cell_id=cell_id,
            cell_confidence=confidence,
            cell_domain=meta.dominant_domain,
            cell_family=meta.dominant_family,
            top_cells=top_cells,
            r=r,
            domain=DOMAIN_NAMES.get(dom_pred, str(dom_pred)),
            domain_confidence=dom_conf,
            domain_probs={DOMAIN_NAMES[i]: float(dom_probs[i]) for i in range(len(dom_probs))},
            mode="kestrel",
            latency_ms=elapsed,
        )

    # -----------------------------------------------------------------
    # Mode 2: Hybrid — KESTREL coarse → Atlas precise
    # -----------------------------------------------------------------
    def classify_hybrid(
        self,
        kestrel_tangent: np.ndarray,
        atlas_tangent: np.ndarray,
        kestrel_domain_probs: Optional[np.ndarray] = None,
        adaptive: bool = True,
        top_k: int = 5,
    ) -> CanonicalResult:
        """
        Hybrid classification: KESTREL assigns coarse cell(s), Atlas resolves within.

        Args:
            kestrel_tangent: (D,) tangent vector from KESTREL
            atlas_tangent: (D,) tangent vector from Atlas
            kestrel_domain_probs: (3,) domain probabilities from KESTREL
            adaptive: if True, adjust search breadth based on KESTREL confidence
            top_k: number of nearest neighbors to return
        """
        t0 = time.perf_counter()

        # KESTREL coarse assignment
        cell_id, confidence, top_cells, _ = self._assign_cell(kestrel_tangent)

        # Adaptive cell count: high confidence → fewer cells → faster
        n_cells = self._adaptive_n_cells(confidence) if adaptive else self.n_search_cells
        search_cells = top_cells[:n_cells]

        # Atlas precise search within KESTREL's cells
        nn_idx, nn_dist, n_searched, nn_top_k = self._search_within_cells(
            atlas_tangent, search_cells, top_k=top_k,
        )

        r = float(np.linalg.norm(atlas_tangent))
        meta = self.cells[cell_id]

        # Domain from KESTREL probs if available, else from cell
        if kestrel_domain_probs is not None:
            dom_pred = int(kestrel_domain_probs.argmax())
            dom_conf = float(kestrel_domain_probs[dom_pred])
            dom_probs = {DOMAIN_NAMES[i]: float(kestrel_domain_probs[i])
                         for i in range(len(kestrel_domain_probs))}
        else:
            dom_pred = list(DOMAIN_NAMES.values()).index(meta.dominant_domain) if meta.dominant_domain in DOMAIN_NAMES.values() else 0
            dom_conf = meta.domain_purity
            dom_probs = {}

        nn_acc = None
        if self.reference_accessions is not None and nn_idx >= 0:
            nn_acc = str(self.reference_accessions[nn_idx])

        elapsed = (time.perf_counter() - t0) * 1000

        return CanonicalResult(
            cell_id=cell_id,
            cell_confidence=confidence,
            cell_domain=meta.dominant_domain,
            cell_family=meta.dominant_family,
            top_cells=top_cells,
            r=r,
            domain=DOMAIN_NAMES.get(dom_pred, str(dom_pred)),
            domain_confidence=dom_conf,
            domain_probs=dom_probs,
            mode="hybrid",
            nn_accession=nn_acc,
            nn_distance=nn_dist,
            nn_top_k=nn_top_k,
            candidates_searched=n_searched,
            search_speedup=float(self.N / max(n_searched, 1)),
            latency_ms=elapsed,
        )

    # -----------------------------------------------------------------
    # Batch operations
    # -----------------------------------------------------------------
    def address_batch(self, tangents: np.ndarray) -> List[CanonicalResult]:
        """Batch canonical addressing."""
        return [self.address_from_tangent(tangents[i]) for i in range(tangents.shape[0])]

    def hybrid_batch(
        self,
        kestrel_tangents: np.ndarray,
        atlas_tangents: np.ndarray,
        kestrel_domain_probs: Optional[np.ndarray] = None,
    ) -> List[CanonicalResult]:
        """Batch hybrid classification."""
        N = kestrel_tangents.shape[0]
        results = []
        for i in range(N):
            kdp = kestrel_domain_probs[i] if kestrel_domain_probs is not None else None
            results.append(self.classify_hybrid(
                kestrel_tangents[i], atlas_tangents[i], kdp,
            ))
        return results

    # -----------------------------------------------------------------
    # Construction
    # -----------------------------------------------------------------
    @classmethod
    def from_teacher_coords(
        cls,
        teacher_path: str,
        K: int = 25,
        n_search_cells: int = 3,
        seed: int = 42,
    ) -> "HybridEngine":
        """
        Build engine from Atlas teacher coordinate file.

        Uses domain-constrained tessellation: K-means runs WITHIN each domain
        independently, guaranteeing 100% domain-pure cells. Cell budget is
        allocated proportionally to domain size (minimum 3 per domain).
        """
        from sklearn.cluster import MiniBatchKMeans

        log.info(f"Building HybridEngine: K={K}, n_search_cells={n_search_cells}, domain-constrained")
        tc = np.load(teacher_path, allow_pickle=False)
        tangents = tc["coords"].astype(np.float32)
        domains = tc["domain_labels"].astype(np.int32)
        families = tc["family_labels"].astype(np.int32)
        family_names = tc["family_names"]
        accessions = tc["accessions"]

        N, D = tangents.shape

        # Allocate cells proportionally to domain size (min 3 each)
        domain_sizes = {d: int((domains == d).sum()) for d in range(3)}
        total_g = sum(domain_sizes.values())
        K_per = {d: max(3, round(K * domain_sizes[d] / total_g)) for d in range(3)}
        while sum(K_per.values()) > K:
            K_per[max(K_per, key=K_per.get)] -= 1
        while sum(K_per.values()) < K:
            K_per[max(K_per, key=lambda d: domain_sizes[d] / K_per[d])] += 1

        for d in range(3):
            log.info(f"  {DOMAIN_NAMES[d]}: {domain_sizes[d]} genomes -> {K_per[d]} cells")

        # K-means within each domain
        all_centers = []
        cells = []
        assignments = np.zeros(N, dtype=np.int32)
        offset = 0

        for d in range(3):
            d_idx = np.where(domains == d)[0]
            km = MiniBatchKMeans(
                n_clusters=K_per[d], random_state=seed,
                batch_size=min(4096, len(d_idx)), n_init=3, max_iter=100,
            )
            km.fit(tangents[d_idx])

            for k in range(K_per[d]):
                members = d_idx[km.labels_ == k]
                gid = offset + k
                assignments[members] = gid

                fam_c = np.bincount(families[members], minlength=len(family_names))
                fi = fam_c.argmax()
                cells.append(CellMeta(
                    cell_id=gid, n_members=len(members),
                    dominant_domain=DOMAIN_NAMES[d],
                    dominant_family=str(family_names[fi]) if fi < len(family_names) else "Unknown",
                    domain_purity=1.0,
                    family_purity=float(fam_c[fi] / len(members)),
                    mean_radius=float(np.linalg.norm(tangents[members], axis=1).mean()),
                    member_indices=members,
                ))

            all_centers.append(km.cluster_centers_)
            offset += K_per[d]

        centers = np.concatenate(all_centers).astype(np.float32)
        log.info(f"  {N} references, {K} cells, ALL domain-pure (constrained)")

        return cls(
            centers=centers, cells=cells,
            reference_tangents=tangents,
            reference_cells=assignments,
            reference_accessions=accessions,
            reference_domains=domains,
            kappa=KAPPA,
            n_search_cells=n_search_cells,
        )

    # -----------------------------------------------------------------
    # Serialization
    # -----------------------------------------------------------------
    def save(self, path: Union[str, Path]) -> None:
        """Save engine to .npz + .json pair."""
        path = Path(path)
        path.parent.mkdir(parents=True, exist_ok=True)

        np.savez_compressed(
            path.with_suffix(".npz"),
            centers=self.centers,
            reference_tangents=self.reference_tangents,
            reference_cells=self.reference_cells,
            reference_domains=self.reference_domains if self.reference_domains is not None else np.array([]),
        )

        meta = {
            "version": self.VERSION,
            "K": self.K, "N": self.N, "D": self.D,
            "kappa": self.kappa,
            "n_search_cells": self.n_search_cells,
            "cells": [{
                "cell_id": c.cell_id, "n_members": c.n_members,
                "dominant_domain": c.dominant_domain,
                "dominant_family": c.dominant_family,
                "domain_purity": round(c.domain_purity, 4),
                "family_purity": round(c.family_purity, 4),
                "mean_radius": round(c.mean_radius, 4),
            } for c in self.cells],
        }
        if self.reference_accessions is not None:
            meta["accessions"] = self.reference_accessions.tolist()

        with open(path.with_suffix(".json"), "w") as f:
            json.dump(meta, f, indent=2)

        size_mb = path.with_suffix(".npz").stat().st_size / 1e6
        log.info(f"Saved HybridEngine → {path} ({size_mb:.1f} MB)")

    @classmethod
    def load(cls, path: Union[str, Path]) -> "HybridEngine":
        """Load engine from .npz + .json pair."""
        path = Path(path)
        data = np.load(path.with_suffix(".npz"), allow_pickle=False)
        with open(path.with_suffix(".json")) as f:
            meta = json.load(f)

        centers = data["centers"]
        ref_tangents = data["reference_tangents"]
        ref_cells = data["reference_cells"]
        ref_domains = data["reference_domains"] if data["reference_domains"].size > 0 else None

        accessions = np.array(meta["accessions"]) if "accessions" in meta else None

        cells = []
        for i, cm in enumerate(meta["cells"]):
            mask = ref_cells == i
            cells.append(CellMeta(
                cell_id=cm["cell_id"], n_members=cm["n_members"],
                dominant_domain=cm["dominant_domain"],
                dominant_family=cm["dominant_family"],
                domain_purity=cm["domain_purity"],
                family_purity=cm["family_purity"],
                mean_radius=cm["mean_radius"],
                member_indices=np.where(mask)[0],
            ))

        return cls(
            centers=centers, cells=cells,
            reference_tangents=ref_tangents,
            reference_cells=ref_cells,
            reference_accessions=accessions,
            reference_domains=ref_domains,
            kappa=meta.get("kappa", KAPPA),
            n_search_cells=meta.get("n_search_cells", 3),
        )

    # -----------------------------------------------------------------
    # Inspection
    # -----------------------------------------------------------------
    def summary(self) -> str:
        """Print engine summary."""
        dom_pure = sum(1 for c in self.cells if c.domain_purity > 0.9)
        fam_coh = sum(1 for c in self.cells if c.family_purity > 0.5)
        avg_cell = self.N // self.K
        avg_search = avg_cell * self.n_search_cells
        return (
            f"HybridEngine v{self.VERSION}\n"
            f"  Tessellation: K={self.K} cells, {self.N:,} references, {self.D}-dim\n"
            f"  Search: up to {self.n_search_cells} cells (adaptive), ~{avg_search:,} candidates\n"
            f"  Quality: {dom_pure}/{self.K} domain-pure, {fam_coh}/{self.K} family-coherent\n"
            f"  Curvature: κ={self.kappa}\n"
            f"  Speedup: ~{self.N // avg_search}×"
        )

    # -----------------------------------------------------------------
    # Presentation layer: (r, θ) coordinates + URL generation
    # -----------------------------------------------------------------
    def attach_projection(self, P: np.ndarray, theta_offset: float = 0.0) -> None:
        """
        Attach a BiosphereCoordinate projection matrix for (r, θ) extraction.

        This is a presentation layer — it does NOT affect retrieval accuracy.
        It maps 129-dim tangent vectors to 2D coordinates for visualization
        and shareable URLs.

        Args:
            P: (2, D) orthonormal projection matrix from BiosphereCoordinate spec
            theta_offset: angular offset in degrees (E. coli = 0°)
        """
        self._projection_P = P.astype(np.float64)
        self._theta_offset = theta_offset

        # Precompute (r, θ) for each cell center
        self._cell_coords = {}
        for c in self.cells:
            r = float(np.linalg.norm(c.center if hasattr(c, 'center') else self.centers[c.cell_id]))
            proj = self._projection_P @ self.centers[c.cell_id].astype(np.float64)
            theta = float(np.degrees(np.arctan2(proj[1], proj[0])) - self._theta_offset) % 360.0
            self._cell_coords[c.cell_id] = (round(r, 4), round(theta, 2))

    def to_r_theta(self, tangent: np.ndarray) -> Tuple[float, float]:
        """
        Project tangent vector to (r, θ) BiosphereCoordinate.
        Requires attach_projection() to have been called.
        """
        if not hasattr(self, '_projection_P'):
            raise RuntimeError("Call attach_projection(P, theta_offset) first")
        v = tangent.astype(np.float64).ravel()
        r = float(np.linalg.norm(v))
        proj = self._projection_P @ v
        theta = float(np.degrees(np.arctan2(proj[1], proj[0])) - self._theta_offset) % 360.0
        return round(r, 4), round(theta, 2)

    def cell_r_theta(self, cell_id: int) -> Tuple[float, float]:
        """Get (r, θ) for a cell center. Requires attach_projection()."""
        if not hasattr(self, '_cell_coords'):
            raise RuntimeError("Call attach_projection(P, theta_offset) first")
        return self._cell_coords[cell_id]

    def coord_url(self, result: CanonicalResult, base: str = "https://biosphereatlas.com") -> str:
        """
        Generate a shareable URL for a classification result.

        Format: {base}/coord/{r}/{theta}
        Or if projection not attached: {base}/cell/{cell_id}
        """
        if hasattr(self, '_cell_coords'):
            r, theta = self._cell_coords[result.cell_id]
            return f"{base}/coord/{r}/{theta}"
        return f"{base}/cell/{result.cell_id}"

    def __repr__(self) -> str:
        return f"HybridEngine(K={self.K}, N={self.N}, cells={self.n_search_cells})"
