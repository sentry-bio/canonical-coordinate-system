# Canonical Coordinate System for the Tree of Life

> **The curvature of life's coordinate manifold is the entropy of its alphabet — κ = (h·ln2)². This is the
> coordinate standard built on that frame, with every claim tagged by how much it has earned: useful as
> infrastructure whether or not the law proves true, and honest enough that finding out is worth doing.**

A **reference coordinate datum** for the vertical history of life — a minimal, theorem-grounded frame that
genomic tools register against, plus the operational retrieval engine that composes on it. *Canonical* here
means **information-tight**: the chart where coordinate = information at one bit of place per bit of biology,
which two maximally-different encoders both converge on — so cross-instrument agreement is the *symptom* that
the frame is real, not the definition. The package defines that frame, measures the agreement gauge-free, and
publishes the transform that registers a working atlas onto it.

**Status: validated *candidate* datum.** κ (two-instrument-validated, audit clean), n≈2 (measured 2.01±0.06),
and the coarse cross-instrument backbone are held; external grounding beyond GTDB and the extension past biology
are open **bets**, named as such. See [`docs/EVIDENCE.md`](docs/EVIDENCE.md) for the tagged ledger
(DERIVED / MEASURED / CIRCULAR / BET) and [`docs/SUBSTRATE_TABLE.md`](docs/SUBSTRATE_TABLE.md) for the
organizing principle. This README does not overclaim past them.

---

## How to read this

The package is a **descent from altitude** — claims get more certain and more concrete as you go down.
Enter where your temperament puts you:

- **Summit (the vision):** the law above → [`docs/SUBSTRATE_TABLE.md`](docs/SUBSTRATE_TABLE.md) (the organizing
  principle, with its predict-then-measure rows) → [`docs/CANONICITY.md`](docs/CANONICITY.md) (why it coheres).
- **Bedrock (the skeptic):** [`verification/`](verification/) (ten receipts that run) →
  [`docs/EVIDENCE.md`](docs/EVIDENCE.md) (every claim tagged DERIVED / MEASURED / CIRCULAR / BET) → then climb.
- **The rope** between the two is the ledger: at every altitude, `EVIDENCE.md` tells you how solid the ground is.
  Nothing in this package asks you to trust it above the confidence it declares.

The **view from altitude** — the full theoretical frame and what it would mean if the law holds — lives in
[`STATE_EQUATION_TIMECAPSULE.md`](STATE_EQUATION_TIMECAPSULE.md). This package is the *instrument*; the time
capsule is the *view*. The instrument stands alone; the view links down here for its grounding.

## The idea in one paragraph

A phylogenetic tree embeds in a 2-dimensional hyperbolic surface whose curvature is fixed by the entropy of
the molecular alphabet: **κ = (h·ln2)² = 5/4** for DNA (h ≈ 1.61 bits/nt) — machine-checked from a
description-length optimum, and confirmed across three alphabets (the [substrate table](docs/SUBSTRATE_TABLE.md):
DNA, protein, RNA-virus, each κ *predicted before measured*). **n = 2 is the hinge** — the one dimension where
the law is parameter-free — and it is *measured* (2.01±0.06), not assumed. Two coordinates exhaust the tree:
**angle** (which lineage) and **radius** (depth). A genome's canonical address is its angular coordinate θ,
anchored so *E. coli* = 0°; independent encoders recover the same coarse backbone (≈16.8° median at large
representative N) — the invariance that signals the frame is real, not an artifact of one network.

## The three layers — sextant · map · projection (different jobs — they do not compete)

| Layer | Geodesy role | What | Role | Module |
|---|---|---|---|---|
| **Datum** | **projection** | the frozen form — κ=5/4, n=2, quartet, anchors, tokenizer, tree | the shared *bounded* frame everything registers against | [`canonical_datum/datum.py`](canonical_datum/datum.py) |
| **Codec** (dim-16) | **sextant** | the minimal low-bias witness | reports (position, **confidence**); makes the map's distortion measurable | (reference realization) |
| **Atlas** (dim-129, live) | **map** | the operational placer | places genomes; distorts deliberately for use (characterized) | `canonical_hybrid/` |

A sextant reports its own error bar; a map distorts on purpose for usability; a projection is useful because it
is *shared and bounded*, not because it is *true*. The **codec↔atlas residual is the map's label-distortion**
(measured, +0.087 std) — the codec earns its keep as the ruler, not as a rival placer. See
[`docs/CANONICITY.md`](docs/CANONICITY.md) (synthesis) and [`docs/POSITIONING.md`](docs/POSITIONING.md) (ledger).

## Scientific foundation (single source of truth)

- **Curvature κ = 5/4**, derived from genetic-code entropy κ = (h·ln2)², verified across alphabets
  (DNA 1.246 pred / 1.247 meas; protein 3.85 / 3.80; RNA virus 0.955 / 0.95). Measured live ≈ 1.237–1.245.
- **Intrinsic dimension n ≈ 2** globally (PCA participation ratio 2.5; top-2 hold 72% of variance;
  tree-consistency saturates at dim-2). Local structure is ~11–12-dimensional — the biology the 2D tree
  cannot hold, resolved by the full-dim atlas.
- **Cross-instrument invariance** — the operational definition of canonical. A dim-16 quartet model and a
  dim-129 multi-loss model recover the same angular coordinate to 16.8° median (random ≈ 90°).

## Quickstart

```python
from canonical_datum import summary, load_transform, to_datum_theta, cross_instrument_agreement

summary()          # the frozen form: κ, n, anchors, version, certified/advisory axes

# register an atlas coordinate onto the datum-canonical angular frame:
T = load_transform("canonical_datum/transforms/datum_transform_v9.json")
theta, radius = to_datum_theta(atlas_coords, T, ecoli_coord, mjann_coord)   # θ certified, r advisory

# measure canonicity = angular agreement between two instruments on the same genomes:
cross_instrument_agreement(codec_coords, atlas_coords, ecoli_i, mjann_i)
# -> {'angular_median_deg': 16.8, 'within_30deg_frac': 0.67, 'radial_advisory_spearman': 0.46, ...}
```

The **operational retrieval engine** (coarse→fine hybrid search, Voronoi cell addressing) lives in
`canonical_hybrid/` and *consumes* the datum — it is the fast-nearest-neighbor layer, not the frame itself.

## What is verified, and what is open

**Held** (see [`docs/EVIDENCE.md`](docs/EVIDENCE.md)): κ from theorem, *two-instrument-validated* (geometry κ\*
⟂ entropy (h·ln2)², circularity audit clean); n≈2 measured; within-clade structure decays coarse→fine (Q1);
the embedded tree is *faithful and signal-limited* (Q3c); registration generalizes (16.9°≈16.8°); anchor gauge
stable (2.4°, 111°±1°). The atlas carries two *orthogonal* distortions — signal-loss (the fade) and label-margin.

**Open — the load-bearing BETs (named, not hidden):**
- **Radius = accumulated information** — the keystone; supported abductively (MDL predicts E10), not proven.
- **External grounding beyond GTDB** — every sequence ruler shares data with the reference tree; needs a
  non-sequence phenotype or an independent tree. Truly open.
- **The law past biology** — neural/linguistic are blank [substrate-table](docs/SUBSTRATE_TABLE.md) rows
  (prior attempts circular); filling them cleanly is a pre-registered, adversarial effort.
- **Certification** — precision candidate-grade; and *cross-version* harmonization is demoted (v10 warm-starts
  from v9, so agreement is largely inheritance — the codec, an independent lineage, is the real invariance test).

## Repository layout

```
canonical_datum/        # the FRAME (core): datum form, registration transform, invariance metric
  ├── datum.py          #   frozen constants (κ, n, anchors, tokenizer, resolution boundary)
  ├── registration.py   #   atlas coord -> datum-canonical θ (advisory r)
  ├── invariance.py     #   cross-instrument angular-agreement metric
  └── transforms/datum_transform_v9.json   # the published v9 -> datum transform
canonical_hybrid/       # the ENGINE (operational layer): Voronoi tessellation + hybrid coarse→fine retrieval
docs/                   # SUBSTRATE_TABLE.md (the organizing principle) · CANONICITY.md (synthesis)
                        #   POSITIONING.md (sextant/map/projection ledger) · EVIDENCE.md (tagged receipts)
verification/           # receipts: eff_dim · q1_within_clade_gauge · q2_fidelity_vs_distortion
                        #   q3c_scalefree · metric_representativeness · phase0_reproduce_pipeline · register
```

The full theoretical frame — the state equation, the MDL/Lyapunov spine, the Buneman≡Gromov quartet bridge,
the honest ledger of the three live bets — is in [`STATE_EQUATION_TIMECAPSULE.md`](STATE_EQUATION_TIMECAPSULE.md).

---
*A validated candidate reference frame, positioned as infrastructure — not a discovery. The discovery is the
κ/E-series program that grounds it; this is the coordinate standard that makes downstream work comparable.*
