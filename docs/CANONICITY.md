# Canonicity — Codec, Atlas, Datum

*A north-star synthesis. What makes a coordinate "canonical," what the three artifacts are, why they
must not converge, and the vector of their mutual evolution. Written to guide, and to be falsified.*

> **The curvature of life's coordinate manifold is the entropy of its alphabet — κ = (h·ln2)². This is the
> coordinate standard built on that frame, with every claim tagged by how much it has earned: useful as
> infrastructure whether or not the law proves true, and honest enough that finding out is worth doing.**

---

## Thesis

> **Canonical = information-tight; cross-instrument invariance is the *symptom*, not the definition.**
> On the description-length optimum the map from genome to coordinate is an *isometry of information* — one
> bit of position per bit of biology (`n_degenerate`, machine-checked). Two maximally-different instruments
> agree because they both converge on *that* chart. The three artifacts are not competitors and must not
> converge; their disagreements are the signal, and their harmonization is a stable, decomposable transform.

---

## I. The reframe: canonical = information-tight (invariance is how we detect it)

We spent a long time treating "canonical" as something the datum *form* confers by declaration, then as mere
**invariance** across instruments. The deeper reading, from the description-length structure: there is a unique
chart in which coordinate = information at rate 1:1 (the addressing term is constant on the critical surface).
*That* is the canonical frame. Invariance is downstream of it — two honest instruments agree precisely because
they are both pulled toward the one information-tight chart.

The evidence for the invariance symptom: a 16-dimensional quartet-only codec and a 129-dimensional classifier,
sharing *nothing* but the data and the entropy theorem, recover the **same angular backbone** — 16.8° median at
large representative N (a *coarse aggregate*; on biased or single-clade samples anchored-θ collapses toward 90°,
which is gauge-fragility, not disagreement — measure it gauge-free). WGS84 is not "true"; it is invariant *and*
bounded. So is this — and the deeper "why" is that both instruments are approximating the information-tight chart.

**Consequence:** the codec's role is not "a weaker atlas." It is the **null instrument** — minimal,
low-bias, task-free — whose only job is to hold a line that no task pressure can move, so that the atlas's
movement becomes *measurable*. You cannot measure distortion without an undistorted reference. The codec
is that reference *because* it is too simple to have an opinion.

### The kernel — one four-point object carries both the tree and the curvature

Why is *information-tightness* the same thing as *tree geometry*? Because the atom of both is the **quartet**.
Buneman's four-point condition (a metric is a tree iff, for every four points, the two largest of the three
pairwise-sums are equal) is *identical* to Gromov's definition of δ-hyperbolicity with δ=0. The same four
points, read two ways:

- **δ (the defect) = deviation from tree = information flowing off the genealogy** (HGT, homoplasy). This is the
  data-processing inequality for Markov descent, written as a metric condition.
- **The Gromov product = shared descent information**, whose exponential decay into angular separation *is* the
  curvature: sin(θ/2) = e^{−√κ·(x·y)}. Capacity-matching — lineages produced (e^{h·ln2·r}) filling the angular
  slots the geometry provides (e^{√κ·r}) — forces √κ = h·ln2.

Proven as a **decomposition**, not asserted: on a synthetic grid (branching × reticulation, ground truth on
both axes), the four-point object splits cleanly — curvature κ is reticulation-blind (corr +0.96 with branching,
+0.05 with reticulation) while the tree-defect is *exactly zero* on pure trees of every curvature
([`verification/kernel_orthogonality.py`](../verification/kernel_orthogonality.py)). So **(κ, δ) are two
independent invariants of any distance matrix**: κ = the entropy/branching, δ = the non-tree flow. That the
*same* object is the atom of phylogenetics (Buneman) and of curvature (Gromov) is not an analogy — it is one
theorem, and it is why this frame is simultaneously a tree and an information geometry. *(One honest residual:
the defect-meter's gain scales with curvature — directionally the kernel's own δ_max ∝ 1/√κ, flagged not
claimed. And the angle-form of the test, fit on a built embedding, is circular by construction — retracted; the
non-circular content is δ.)*

---

## II. The three artifacts

| Artifact | Is | Role | Judged by |
|---|---|---|---|
| **Datum** | the frozen *form*: κ=5/4, n=2, quartet objective, anchors, reference tree | the ideal both realizations approximate | theorem-consistency |
| **Codec** (dim-16) | the minimal realization = **null/witness instrument** | hold the canonical line fixed; make drift measurable; anchor versioning | cross-instrument invariance |
| **Atlas / v9** (dim-129) | the operational realization, **wired to classify** | place genomes accurately (stays live) | placement fidelity |

The datum is a form. The codec is incorruptible but coarse. The atlas is powerful but opinionated. None is
"the truth." The truth is what is invariant across them.

---

## III. The discrepancy decomposes — into two MEASURED, orthogonal axes

We can now replace the old conjecture ("task-distortion vs biological-overlay") with what was actually measured
this iteration. The codec↔atlas discrepancy is two *independent* things:

1. **Signal-loss (the fade).** Shared branch structure thins coarse→fine (within-clade Mantel: phylum 0.55 →
   family 0.51 → genus 0.23; relative branch S_rel flat, so the *shape* is self-similar — the fade is
   magnitude). This is **not corruption**: the tips are functionally flat, so there is less to resolve. The
   embedding stays a *faithful tree wherever there is signal* (ρ≈0.014 at good signal vs 0.174 scatter). It is
   **signal-limited, not noise-limited.**
2. **Label-margin (the map's deliberate distortion).** Both instruments insert a genus boundary beyond what
   divergence warrants (codec +0.18, v9 +0.37 std; v9 > codec in 70% of families). This is the atlas's
   "Mercator distortion" — bought for discriminability.

**The retraction that matters:** these two are *orthogonal*, and four-point fidelity (Q3) is **blind** to the
label-margin (Q2) — because inflating a genus branch keeps the embedding a *valid tree*. So the atlas is a
**cleaner tree that over-separates genera** (v9 ρ 0.013 < codec 0.020: dimension buys tree-fidelity *and* costs
label-honesty). Do not fuse the two axes; they are different kinds of wrongness, measured by different rulers.

---

## IV. Registration is de-distortion

Registering the atlas to the datum — transforming its angle onto the anchored canonical frame — is, in
part, **undoing the atlas's discriminative warp** to re-expose the tree beneath its classifier training.
The datum is a corrective: it strips the task-specific opinion off the atlas to recover the shared geometry.
The certified axis is **θ (lineage)** — cross-instrument-shared. The **radial axis (depth) is advisory**
(cross-instrument Spearman only +0.46) — present, carried, labeled, but not certified.

---

## V. What the harmonized picture gets us

1. **Canonicity by triangulation, grounded in information-tightness.** No artifact is authoritative; the
   canonical is the information-tight chart both instruments approximate. A defensible epistemics, not a declared one.
2. **A distortion governor — the label-margin.** The codec↔atlas residual is a *measured* distortion budget
   (Q2: +0.087 std, v9>codec in 70% of families) — it flags when the map has over-separated genera beyond
   divergence. This is the map's characterized "Mercator distortion," not a mystery residual.
3. **Two orthogonal residual axes, both first-class** — signal-loss (the fade: how far the datum's resolution
   reaches) and label-margin (how much the map bends to labels). Measured separately (Q1/Q3c vs Q2); do not fuse.
4. **Comparability through the datum** — a v9 placement carries a datum-registered θ, so coordinates are
   comparable across labs and time. *Cross-version* comparability specifically is a weak case (v10 warm-starts
   from v9 → inheritance, not rediscovery); the **codec, an independent lineage, is the load-bearing witness**.

---

## VI. The vector of mutual evolution — they must NOT converge

The instinct is that harmonization means the three become one. **It does not.** Convergence would mean
either the codec grew opinions (losing its witness role) or the atlas got dumber (losing its fidelity).
Harmonization is a **stable, decomposable transform between them**, and the vector of evolution is toward
**legibility of the residual**, not identity:

- **Codec** → tightens as an *incorruptible witness* (more anchors, more precision). Evolves toward a
  sharper canonical line, **never toward capability**. Its power is its poverty, done precisely.
- **Atlas** → free to specialize and improve as a placer, **now governed**: the datum flags when
  discriminative sharpening is *drift* rather than *fidelity*. It can be anchored to the datum on θ (kept
  honest to the tree) while left free in the residual dimensions (where task-utility and biology live).
  Anchored where it must be, free where it should be.
- **Datum** → frozen form; its only motion is the reference tree improving across releases, with published
  transforms so all prior coordinates migrate.

The end-state is not one thing. It is a **governed ecosystem**: a reference too simple to be corrupted, a
realization powerful enough to be useful, and a contract between them that keeps the powerful one honest
about the difference between getting better and drifting.

---

## VII. Evidence status — and the prediction that was tested (and refuted)

**Held (measured):** κ from theorem, two-instrument-validated (audit clean); n≈2 global (measured 2.01±0.06);
coarse→fine decay (Q1); faithful self-similar signal-limited tree (Q3c); two orthogonal distortions (§III);
partial patristic grounding (~0.29). Anchored-θ is a fragile *chart* — measure gauge-free.

**A prior version of this section predicted** that v10.6 (trained harder on leaf-contrastive separation) would
show a *larger* fine-scale four-point residual than v9 — the "task-distortion signature." **This was tested and
refuted.** Q3c found the opposite ordering: the higher-dimensional, more-trained atlas is the *cleaner* tree
(v9 ρ 0.013 < codec 0.020). And the framing was wrong at the root — four-point fidelity is *orthogonal* to the
label-margin (inflating a genus branch stays a valid tree), so a distortion-via-four-point test could never have
worked. **The label-margin is real and measured (Q2); it simply is not visible to a tree-fidelity ruler.** This
correction — a refuted prediction, logged rather than hidden — is exactly the self-correction the ledger is for.

**The real open frontier** is no longer cross-version (that is ancestry-contaminated inheritance) but the three
BETs in `EVIDENCE.md`: radius=information, κ-independence-of-h in the wild (audit clean; residual remains), and
whether the law extends past biology (the blank substrate-table rows). Those are where evidence moves the needle.

---

## VIII. Guardrails (so the north star stays honest)

- The 16.8° backbone is **fuzzy** — "the canonical tree" is aspirational shorthand for "the
  cross-instrument-invariant backbone, our current best estimate." Precision refinement earns the crisper word.
- The two distortions (signal-loss, label-margin) are now **measured and orthogonal** (§III). Do not revive the
  old "task-distortion vs biological-overlay" conjecture or the refuted v10.6 prediction (§VII) as if live.
- Do not let the codec's minimality (fails function, loses at placement) read as failure — for a witness,
  those are the purity tests passing.
- Canonicity is **invariance**, not authority. If a future instrument disagrees with the backbone, the
  backbone is what must be re-examined — not the instrument dismissed.

---
*Companion: `CODEC_DATUM_POSITIONING.md` (the tagged ledger). This document is the why; that one is the what.
Both are candidates, not certainties — the honest state of a reference frame still being validated.*
