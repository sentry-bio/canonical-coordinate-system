# The Substrate Table — an organizing principle for information-generating hierarchies

*The Mendeleevian object. One law orders the rows; the filled cells are confirmations; the blank cells are
predictions. Read the layout, not the prose — the way a periodic table is read.*

---

## The law (one line, with its cause)

> **The curvature of an information-generating hierarchy's coordinate manifold is the entropy of its alphabet.**
>
> $$(n-1)\sqrt{\kappa} = h\ln 2 \quad\xrightarrow{\ n=2\ }\quad \kappa = (h\ln 2)^2 \qquad \text{(zero free parameters)}$$

Left side is the manifold's **volume entropy** (the rate at which room appears with radius — pure geometry).
Right side is the alphabet's **information rate** (typical-sequence growth — pure Shannon). The law sets them
equal: geometry holds exactly as many distinguishable states as the code generates. **n = 2 is the hinge** —
it is the only dimension where $(n-1)=1$ and the chemistry reads off *directly* as a curvature; it is also the
minimal faithful embedding of a tree (Sarkar). Machine-checked from a description-length optimum
(`theory/lean/…/KappaCurvature.lean`); the empirical bet is that biology sits at that optimum.

---

## The table

Rows are substrates ordered by alphabet entropy *h*. The law fills the κ column. **Filled rows are
measured confirmations; blank rows are live predictions — compute κ from *h*, then go measure it.**

```
 SUBSTRATE       alphabet        h (bits/sym)   n         κ = (h·ln2)²        measured κ      status
 ───────────────────────────────────────────────────────────────────────────────────────────────────
 RNA virus       4 (constrained)  ~1.4          2         0.955  (predict)    ~0.95           ✓ confirmed
 DNA genome      4  nt            ~1.61         2.01±.06   1.246  (predict)    ~1.247          ✓ confirmed
 Protein         20 aa            ~2.8*         2.03       3.85   (predict)    ~3.80           ✓ confirmed
 ───────────────────────────────────────────────────────────────────────────────────────────────────
 tRNA / rRNA     4 (structured)    ?            2          ← PREDICT then test   —             ○ open
 Codon           64                ?            2          ← PREDICT             —             ○ open
 Neural code     ?                 ?            2          ← PREDICT             —             ○ pending
 Linguistic      ?                 ?            2          ← PREDICT             —             △ flagged*
 ───────────────────────────────────────────────────────────────────────────────────────────────────
        n = 2 column is CONSTANT → the hinge.   (*protein h back-computed from κ; linguistic prior attempt
        was tree-circular — a blank row is a prediction, not a claim.)
```

Everything Mendeleevian is readable off the grid: the **ordering** (by *h*), the **law** (each filled row
obeys it), the **substrate-independence** (form constant, κ varies), the **confirmations** (predict≈measure,
noisily — DNA 0.1%, protein ~1.3%, virus ~0.5%; the scatter is the signature of *measurement*, not identity),
and — the decisive feature — the **blank rows are literal falsifiable predictions.**

## Why this is a law and not a curve-fit (the eka-silicon move)

Clade-breadth is not alphabet-breadth. Measuring across all clades of the tree pins *one* point on the κ curve
unshakably (same alphabet, same *h*) but adds no independent points. The **law** is tested by moving *h* and
watching κ follow — and that is the DNA/protein/virus rows: a *different* alphabet, a *different* entropy, κ
**predicted before measured**, and it lands. That is Mendeleev predicting eka-silicon, not Bode fitting a curve.

## The audit that licenses the table (2026-07-29)

The load-bearing risk was circularity: is *h* smuggled into the κ measurement? Audited and **cleared** — two
genuinely independent instruments that meet only at the comparison:

- **κ\* (geometry)** — the curvature that best embeds a **Mash/MinHash** distance matrix (alignment-free k-mer
  Jaccard; *no* substitution model, *no* base-frequency assumption, *no h* input). Validated on synthetic
  trees: δ calls trees hyperbolic and Euclidean nulls Euclidean; branching recovered monotonically.
- **κ_Manning (information)** — (h·ln2)² from *empirical* per-symbol entropy (conditional/block rate; the naive
  log₂α is explicitly held out as a broken control). Validated: recovers a known Markov rate to ~0.5%.

The prediction side (MDL/Lambert-W) is pure theory, a separate pipeline. **Documented residual:** *h* and the
geometry are orthogonal statistics of the *same* sequences — not independent data sources — so the law's
content is precisely that these two projections coincide.

## What filling the blank rows would mean

If neural/linguistic fill cleanly (independent *h*, parameter-free κ prediction, *h*-independent κ measurement,
a within-domain *h*→κ curve, and a discriminating null), the table extends past biology and becomes a
**periodic table of information-generating processes** — and biology's rows reclassify *upward*, from "a law of
biology" to "an instance of a law of information hierarchies." If they don't fill, they stay blank and biology's
rows lose nothing. **The table is honest either way — which is the whole point.** See `docs/EVIDENCE.md` for the
tagged ledger and `STATE_EQUATION_TIMECAPSULE.md` for the full frame.
