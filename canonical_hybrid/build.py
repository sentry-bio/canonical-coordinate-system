"""
Build a CanonicalHybrid tessellation from Atlas teacher embeddings.

Usage:
    from canonical_hybrid import build_tessellation

    hybrid = build_tessellation(
        teacher_path="/home/rohit/e1_results/teacher_coords.npz",
        K=50,
        output_path="/home/rohit/tessellation_K50",
    )

Or from command line:
    python -m canonical_hybrid.build \
        --teacher /home/rohit/e1_results/teacher_coords.npz \
        --K 50 \
        --output /home/rohit/tessellation_K50
"""

from __future__ import annotations

import argparse
import logging
import time
from pathlib import Path
from typing import Optional

import numpy as np

log = logging.getLogger(__name__)

DOMAIN_NAMES = {0: "Bacteria", 1: "Archaea", 2: "Eukaryota"}


def build_tessellation(
    teacher_path: str,
    K: int = 50,
    output_path: Optional[str] = None,
    seed: int = 42,
    kappa: float = 1.25,
) -> "CanonicalHybrid":
    """
    Build a CanonicalHybrid from Atlas teacher embeddings.

    Args:
        teacher_path: path to teacher_coords.npz
        K: number of Voronoi cells
        output_path: if provided, save the tessellation here
        seed: random seed for K-means
        kappa: curvature constant

    Returns:
        CanonicalHybrid instance
    """
    from sklearn.cluster import MiniBatchKMeans
    from canonical_hybrid.core import CanonicalHybrid, TessellationMeta

    log.info(f"Building K={K} tessellation from {teacher_path}")

    # Load teacher data
    tc = np.load(teacher_path, allow_pickle=False)
    tangents = tc["coords"].astype(np.float32)          # (N, D)
    domains = tc["domain_labels"].astype(np.int32)       # (N,)
    families = tc["family_labels"].astype(np.int32)      # (N,)
    accessions = tc["accessions"]                         # (N,)
    family_names = tc["family_names"]                     # (F,)

    N, D = tangents.shape
    log.info(f"  {N} embeddings, {D}-dim, {len(np.unique(families))} families")

    # K-means
    t0 = time.time()
    km = MiniBatchKMeans(
        n_clusters=K, random_state=seed,
        batch_size=min(4096, N), n_init=3, max_iter=100,
    )
    km.fit(tangents)
    centers = km.cluster_centers_.astype(np.float32)
    assignments = km.labels_.astype(np.int32)
    log.info(f"  K-means: {time.time() - t0:.1f}s")

    # Build cell metadata
    cell_meta = []
    for k in range(K):
        mask = assignments == k
        n_members = int(mask.sum())

        if n_members == 0:
            cell_meta.append(TessellationMeta(
                cell_id=k, n_members=0,
                dominant_domain="Unknown", dominant_family="Unknown",
                domain_purity=0.0, family_purity=0.0,
                mean_radius=0.0, center=centers[k],
            ))
            continue

        # Domain stats
        cell_domains = domains[mask]
        dom_counts = np.bincount(cell_domains, minlength=3)
        dom_idx = dom_counts.argmax()
        dom_purity = float(dom_counts[dom_idx] / n_members)

        # Family stats
        cell_families = families[mask]
        fam_counts = np.bincount(cell_families, minlength=len(family_names))
        fam_idx = fam_counts.argmax()
        fam_purity = float(fam_counts[fam_idx] / n_members)
        fam_name = str(family_names[fam_idx]) if fam_idx < len(family_names) else "Unknown"

        # Radial stats
        cell_radii = np.linalg.norm(tangents[mask], axis=1)
        mean_r = float(cell_radii.mean())

        cell_meta.append(TessellationMeta(
            cell_id=k, n_members=n_members,
            dominant_domain=DOMAIN_NAMES.get(int(dom_idx), "Unknown"),
            dominant_family=fam_name,
            domain_purity=dom_purity,
            family_purity=fam_purity,
            mean_radius=mean_r,
            center=centers[k],
        ))

    # Summary
    dom_pure = sum(1 for m in cell_meta if m.domain_purity > 0.9)
    fam_pure = sum(1 for m in cell_meta if m.family_purity > 0.5)
    log.info(f"  Domain-pure cells (>90%): {dom_pure}/{K}")
    log.info(f"  Family-coherent cells (>50%): {fam_pure}/{K}")
    log.info(f"  Mean cell size: {N/K:.0f}")

    hybrid = CanonicalHybrid(
        centers=centers,
        cell_meta=cell_meta,
        reference_tangents=tangents,
        reference_cells=assignments,
        reference_domains=domains,
        reference_families=families,
        reference_accessions=accessions,
        kappa=kappa,
    )

    if output_path:
        hybrid.save(output_path)
        log.info(f"  Saved → {output_path}.{{npz,json}}")

    return hybrid


def main():
    logging.basicConfig(level=logging.INFO, format="[%(asctime)s] %(message)s")

    parser = argparse.ArgumentParser(description="Build canonical tessellation")
    parser.add_argument("--teacher", required=True, help="teacher_coords.npz")
    parser.add_argument("--K", type=int, default=50, help="Number of Voronoi cells")
    parser.add_argument("--output", required=True, help="Output path (without extension)")
    parser.add_argument("--seed", type=int, default=42)
    args = parser.parse_args()

    hybrid = build_tessellation(
        teacher_path=args.teacher,
        K=args.K,
        output_path=args.output,
        seed=args.seed,
    )
    print(f"\n{hybrid}")
    print(f"\nCell taxonomy summary (first 10):")
    for c in hybrid.cell_taxonomy_summary()[:10]:
        print(f"  Cell {c['cell_id']:3d}: {c['dominant_domain']:12s} "
              f"{c['dominant_family']:30s} "
              f"n={c['n_members']:5d}  dom_purity={c['domain_purity']:.0%}  "
              f"fam_purity={c['family_purity']:.0%}")


if __name__ == "__main__":
    main()
