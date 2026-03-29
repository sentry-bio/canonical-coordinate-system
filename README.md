# Canonical Coordinate System for the Tree of Life

A tangent-space Voronoi tessellation that provides a multi-resolution canonical addressing scheme for genomic sequences.

**97.1% R@1 retrieval accuracy at 7× speedup** over brute-force search through 47,000 reference genomes.

## What This Does

Every genome gets a **canonical address** — a cell in a Voronoi tessellation of the hyperbolic tangent space at curvature κ = 5/4. Each cell carries taxonomic metadata (domain, family, purity). Two independent encoders (KESTREL from k-mer spectra, Atlas from full sequences) converge to the same cell assignment 85.6% of the time at K=25 resolution.

```python
from canonical_hybrid import HybridEngine

engine = HybridEngine.load("tessellation_K25")

# From a tangent-space embedding (129-dim)
address = engine.address_from_tangent(tangent_vector)
print(address.cell_id)        # 17
print(address.cell_domain)    # "Bacteria"
print(address.cell_family)    # "Lactobacillaceae"
print(address.cell_confidence)# 2.3

# Hybrid search: KESTREL narrows, Atlas resolves
match = engine.classify_hybrid(kestrel_tangent, atlas_tangent)
print(match.nn_accession)     # "GCA_000005845.2"
print(match.nn_distance)      # 0.34
print(match.search_speedup)   # 7.0
```

## Architecture

```
                    ┌─────────────────────────────┐
                    │  Canonical Coordinate System │
                    │  K=25 Voronoi tessellation   │
                    │  47,000 reference genomes    │
                    └──────────┬──────────────────┘
                               │
              ┌────────────────┼────────────────┐
              │                │                │
     KESTREL (CPU)      Atlas (GPU)      Hybrid (CPU+GPU)
     0.03ms/query      10ms/query       0.85ms/query
     85.6% cell agree  100% self-R@1    97.1% R@1
     standalone addr   precise coords   7× speedup
```

## Benchmarks (47K genomes, K=25)

| Mode | R@1 | R@10 | Speedup | Latency |
|------|-----|------|---------|---------|
| Hybrid (3 cells, fixed) | 97.1% | 97.1% | 7× | 0.85ms |
| Hybrid (adaptive) | 91.2% | — | 14× | 0.39ms |
| KESTREL standalone | 85.6% cell agree | 97.1% top-3 | — | 0.03ms |

## Scientific Foundation

The coordinate system is grounded in the [BiosphereAtlas](https://biosphereatlas.com) research program:

- **Curvature κ = 5/4**: Derived from the entropy of the genetic code (h ≈ 1.61 bits/nucleotide), verified by machine-checked Lean 4 proof. Measured empirically at κ = 1.2505 across 120+ training runs.
- **Frame 2 validation**: 85%+ cross-domain transfer (E2), cos = 0.974 agreement between independent encoders, Manning entropy n ≈ 2.
- **Tangent-space stability**: L2 distances in tangent space have Spearman ρ = 0.918 correlation with true Poincaré geodesic distances, without the float32 precision collapse at the ball boundary.

## Installation

```bash
pip install canonical-coordinate-system
```

Building a tessellation requires `scikit-learn`:
```bash
pip install canonical-coordinate-system[build]
```

## Building a Tessellation

```python
from canonical_hybrid import HybridEngine

engine = HybridEngine.from_teacher_coords(
    "teacher_coords.npz",  # Atlas tangent vectors
    K=25,                   # number of Voronoi cells
    n_search_cells=3,       # cells to search in hybrid mode
)
engine.save("my_tessellation")
```

## API

The engine exposes three modes:

- `address_from_tangent(v)` — canonical address from any 129-dim tangent vector
- `classify_kestrel(dna, model, ...)` — full DNA → address pipeline (requires `kestrel-bio`)
- `classify_hybrid(kestrel_v, atlas_v)` — KESTREL coarse → Atlas precise search

## License

MIT. Created by [Sentry Bio](https://sentry.bio).
