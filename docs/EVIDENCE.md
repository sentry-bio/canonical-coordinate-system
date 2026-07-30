# Evidence Ledger

> **The curvature of life's coordinate manifold is the entropy of its alphabet — κ = (h·ln2)². This is the
> coordinate standard built on that frame, with every claim tagged by how much it has earned: useful as
> infrastructure whether or not the law proves true, and honest enough that finding out is worth doing.**

*The honest state of the datum. Every claim is tagged **DERIVED** (deductive, machine-checked) / **MEASURED**
(empirical, holds) / **CIRCULAR** (semi-circular, flagged) / **BET** (load-bearing, open). Written to survive
scrutiny — the load-bearing open bets are named, not buried. This ledger is what separates the substrate table
from a Bode's-law coincidence: the gaps are real predictions and the residuals are honestly seated.*

---

## The claims, by status

### DERIVED — follows deductively; the math is machine-checked
| Claim | Number | Source |
|---|---|---|
| State equation κ=(h·ln2)² from a description-length optimum | 9 theorems `#check`-verified | `KappaCurvature.lean` (Fenn & Fenn, Zenodo) |
| n=2 is optimal (embeddability floor + L-degeneracy above it) | `n_floor_optimal`, `n_degenerate` | Lean Part III/IV |
| κ\* is a Lyapunov-stable attractor (why a loss finds it) | `h_stable`, `lyapunov_zero_iff` | Lean Part III |

*What "machine-checked" buys: the theorems follow from the definitions. The empirical content lives in the
definitions (radius=information at rate h·ln2); Lean certifies the deduction, not that biology obeys the model.*

### MEASURED — empirical, holds under scrutiny
| Claim | Number | Source |
|---|---|---|
| κ from entropy, zero free parameters, across alphabets | DNA 1.246/1.247 · protein 3.85/3.80 · virus 0.955/0.95 | E-series; audit below |
| **Two-instrument, not tautology** (geometry κ\* ⟂ entropy (h·ln2)²) | both recover known truth synthetically | hinge audit (2026-07-29) |
| Intrinsic dimension n≈2 (global) | measured 2.01±0.06; PCA participation ratio 2.5 | `verification/eff_dim.py`, results.yaml R6 |
| Local structure higher-dim | TwoNN ≈ 11.6 (Euclidean & hyperbolic agree) | `verification/eff_dim.py` |
| Shared within-clade structure decays coarse→fine | Mantel phylum .55 → family .51 → genus .23 (null 0) | `verification/q1_within_clade_gauge.py` |
| Embedded tree is faithful & signal-limited (not noise-limited) | ρ 0.014 at good signal < 5%-noise tree; scatter 0.174 | `verification/q3c_scalefree.py` |
| Self-similar shape (relative branch flat across ranks) | S_rel ≈ 0.07 phylum→genus | `verification/q3c_scalefree.py` |
| Fine structure partially grounded in real divergence | Spearman embedded↔patristic ≈ 0.29 (both instruments) | `verification/q2_fidelity_vs_distortion.py` |
| **Four-point object splits: curvature ⊥ tree-defect** (the Buneman≡Gromov kernel) | κ_vol reticulation-blind (corr b +0.96, f +0.05); tree floor τ=0 exact for all κ | `verification/kernel_orthogonality.py` |
| Cross-instrument angular backbone (coarse) | 16.8° median codec↔v9, 67% within 30° (random ~90°); topology agreement 0.79 | `verification/register.py` |

### CIRCULAR — real but semi-circular; flagged, not leaned on
| Claim | Why flagged |
|---|---|
| Patristic "grounding" (~0.29) | GTDB-derived; separates continuous-divergence from discrete-rank, but not nature-vs-GTDB |
| The audit's residual | h and geometry are orthogonal statistics of the *same* sequences (not independent data) |
| Cross-instrument θ agreement (16.8°) | coarse aggregate, stable only at large *representative* N; anchored-θ is gauge-fragile |

### BET — load-bearing and open; the entire remaining risk
| Bet | Status |
|---|---|
| **Radius = accumulated information** | supported *abductively* (MDL predicts E10: radius=complexity-not-time) — not proven |
| **κ measured independently of h in the wild** | audit CLEAN on the pipeline; the same-sequences residual keeps it a bet |
| **Biology is a clean tree at sequence level** | Q3 is about the *embedding*; needs tier-2 raw-sequence four-point |
| **The law extends beyond biology** | neural/linguistic are blank substrate-table rows (prior attempts circular) |

---

## The two distortions the atlas carries (measured, orthogonal)

The operational atlas is a **map** — it distorts deliberately for usability, and we measured both distortions,
which turn out to be *independent axes*:

1. **Signal-loss (the fade).** Branch signal thins toward the tips (Q1 decay, Q3c S_rel). Not corruption —
   the tips are functionally flat, so there is less to resolve. The datum's resolution *fades*; it has no cliff.
2. **Label-margin (Q2).** Both instruments insert a genus boundary beyond what divergence warrants (codec +0.18,
   v9 +0.37 std; v9 > codec in 70% of families). Forcing cluster-separability is the map's "Mercator distortion."

**These two are orthogonal** (a retraction from this session): inflating a genus branch keeps the embedding a
*valid tree*, so four-point fidelity (Q3) is *blind* to the label-margin (Q2). The atlas is a **cleaner tree
that over-separates genera** (v9 ρ 0.013 < codec 0.020 — dimension buys fidelity). Q2 and Q3 measure different
wrongness; do not triangulate them.

## The four solid-ground checks — final status

1. **Round-trip — PASS.** Transform fit on 70%, applied to held-out 30%, θ 16.9° ≈ in-sample 16.8°.
2. **Anchor-gauge stability — PASS.** 2D plane stable to 2.4°; E.coli/M.jann gauge 111°±1°.
3. **Number/doc audit — PASS.** Docs match runs; the independent measurements tell one coherent story.
4. **Cross-version harmonization — REFRAMED (not a fixable pass-as-stated).** The confound was diagnosed:
   (a) the deployed reference was pinned to *stale tokens* (regenerated a month post-build → 12.7° drift), and
   (b) anchored-θ is *representativeness-fragile* (biased/single-clade samples → ~90°). Deeper reframe:
   cross-*version* agreement is largely *inheritance* (v10 warm-starts from v9), so the **codec — an independent
   lineage — is the real invariance test**, and it already fired at the coarse scale. Cross-version is demoted
   from "the central promise" to "a weak, ancestry-contaminated test."

## Retractions logged this session (the ledger self-corrects)

- "Codec is the *faithful* witness" → **least-distorted**, not faithful (it inserts +0.18 margin too).
- "The dial fingers HGT genera" → threshold artifact; the confidence dial survives, the HGT reading does not.
- "Q2×Q3 triangulate one distortion" → orthogonal axes; refuted.
- "θ injects ~40° windowing noise" → pipeline is bit-deterministic; the fragility is *representativeness*, not noise.

## What this licenses — and what it does not

**Earned:** a theorem-grounded, resolution-bounded, gauge-free coordinate frame; a two-instrument (non-tautological)
curvature law confirmed across three alphabets; a measured anatomy (faithful self-similar tree + two orthogonal
distortions); a self-correcting method with its negatives on the record.

**Not earned:** that the frame is grounded beyond GTDB; that the law extends past biology; that radius is proven
to be the information coordinate. These are BETs, and closing them is the frontier — the h-circularity audit
(done, clean), then tier-2 raw-sequence, then one non-biological substrate-table row under pre-registration.

## Reproduction

`verification/` holds the receipts: `eff_dim.py` (dimension), `q1_within_clade_gauge.py` (decay),
`q2_fidelity_vs_distortion.py` (grounding+distortion), `q3c_scalefree.py` (tree-fidelity closure),
`kernel_orthogonality.py` (the δ⊥κ decomposition — curvature ⊥ tree-defect, synthetic ground truth),
`metric_representativeness.py` (the θ-fragility floor), `phase0_reproduce_pipeline.py` (the token-drift
diagnosis), `register.py` (coarse agreement). The κ-measurement pipeline itself (Mash distance, empirical
entropy, geometry κ\*) was independently validated against known synthetic truth in the hinge audit. They run
against the reference realizations (dim-16 codec `coords_final.npz`, dim-129 v9 `karcher_v7.npz`).
