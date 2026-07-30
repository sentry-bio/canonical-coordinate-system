# The State Equation — A Time Capsule

> **The curvature of life's coordinate manifold is the entropy of its alphabet — κ = (h·ln2)². This is the
> coordinate standard built on that frame, with every claim tagged by how much it has earned: useful as
> infrastructure whether or not the law proves true, and honest enough that finding out is worth doing.**

*A record of the view, at its most information-dense — the elegant collapse, honestly ledgered.
Not a product claim; a snapshot of what the whole program looks like when it clicks into one equation.
Dated 2026-07-26. Companion: `KappaCurvature.lean` (machine-checked), Fenn & Fenn, "A Description-Length
Principle for Information-Generating Hierarchies" (Zenodo 10.5281/zenodo.19381558).*

---

## I. The one equation

$$(n-1)\sqrt{\kappa} \;=\; h\ln 2 \qquad\xrightarrow{\;n=2\;}\qquad \boxed{\kappa = (h\ln 2)^2}$$

Both sides are exponential growth rates, and the equation just sets them equal:

- **Left — the geometry's capacity.** A ball of radius $r$ in $\mathbb{H}^n$ of curvature $-\kappa$ has volume $\sim e^{(n-1)\sqrt{\kappa}\,r}$. The **volume entropy** $(n-1)\sqrt{\kappa}$ is the rate at which room appears as you move outward. (Pure geometry.)
- **Right — the data's demand.** The number of distinguishable typical sequences at divergence $r$ grows as $e^{(h\ln 2)\,r}$, $h$ = Shannon entropy of the molecular alphabet (DNA $h\approx 1.61$ bits/nt). (Pure information theory.)
- **The match.** Demand equal volumes hold equal numbers of distinguishable organisms → the exponents equate.

**n = 2 is the hinge, not a fitted parameter.** At $n=2$, $(n-1)=1$ and the law is *parameter-free*: the alphabet entropy **is** the curvature, undressed. At $n=3$ a factor 4 appears and the magic dies. And $n=2$ is forced from the other side too: a tree embeds in $\mathbb{H}^2$ with vanishing distortion (Sarkar 2011) and needs no more. **DNA and protein share $n=2$; only $\kappa$ moves with the alphabet** ($\kappa_{\text{DNA}}\!\approx\!1.25$, $\kappa_{\text{protein}}\!\approx\!3.8$, both $(h\ln2)^2$; protein $n=2.03$ across 14 families). $n$ is geometric, $\kappa$ is informational.

## II. Why it is forced, not fit — the description-length spine (machine-checked)

The state equation is the stationary point of a **minimum-description-length** optimization, and the nine load-bearing steps are verified in Lean 4:

- Description length per transmitted bit splits into **transmission** ("what it is") + **addressing** ("where it is"): $L = (h+c)/I + (n-1)\sqrt{\kappa}/(I\ln 2)$.
- $L$ strictly increases in $\kappa$ (`L_monotone_kappa`) → the optimizer wants the *flattest* space.
- Losslessness requires capacity $\ge$ information: $(n-1)\sqrt{\kappa}\ge h\ln2$. The constraint **binds** at the optimum (`constraint_active`: the rate-matching residual is zero *iff* $\kappa=\kappa_{\text{critical}}$). → the state equation is the tightest lossless curvature. **$\kappa$ is a close-packing value.**
- On the critical surface the **addressing term is a constant — one bit of position per bit of biology** (`n_degenerate`). This is the deepest line in the whole theory:

  > **Canonical = information-tight.** The canonical frame is the unique chart in which coordinate = information at rate 1:1. Cross-instrument invariance (two encoders agreeing) is the *symptom*; information-tightness is the *cause*. Two honest instruments agree because they both converge on the one 1-bit-per-bit chart.

- $n=2$ is optimal twice over: **embeddability floor** (trees need $\mathbb{H}^2$) + **degeneracy above it** ($L$ is $n$-independent on the critical surface, `n_degenerate`), so parsimony takes the floor.
- $\kappa^\*$ is a **Lyapunov-stable attractor** (`h_stable`): $V=(\sqrt\kappa-\sqrt{\kappa^\*})^2$ is a Lyapunov function. This is *why* a contrastive loss locks $\kappa\to$ the value in ~5 epochs regardless of embedding dimension (E6) — training descends $U$ to its unique zero.

*What "machine-checked" buys, exactly:* it certifies the theorems **follow from the definitions**. The empirical content lives in the definitions ($I=h\ln2$ per radius; the description-length objective). Lean closes the deduction; it cannot and does not certify that biology obeys the model.

## III. The quartet — atom of both the tree and the curvature

Four is the minimal unit of branching (3 taxa = a star, no decision; 4 = the first topology choice), and a tree is *determined* by its quartets. **Buneman's four-point condition** (a metric is a tree iff, for every 4 points, the two largest of the three pairwise-sums are equal) is *literally* **Gromov's definition of $\delta$-hyperbolicity** with $\delta=0$. So the quartet is simultaneously the atom of phylogenetics and the definitional probe of curvature — not by analogy, by theorem.

**But the two faces are orthogonal on leaf data, and conflating them is a trap:**

| face | quantity | what it measures on data | our instrument |
|---|---|---|---|
| **volume** | $\kappa=(n-1)^{-2}(h\ln2)^2$ | curvature ↔ proliferation rate | volume entropy (count growth) |
| **thinness** | four-point $\delta$ | *deviation from tree* (0 for any tree, regardless of $\kappa$) | quartet test |

A perfect tree has $\delta=0$ **whatever the ambient curvature** — so $\delta$ on genome points measures reticulation/scatter, *not* $\kappa$. $\kappa$ is read from how the *count* of organisms grows with radius. Two faces of one hyperbolic geometry; two different measurements.

## IV. The instruments, and what they actually found

- **$\kappa$-meter (volume entropy).** Certified-collapse: $h_{\text{vol}}=h_{\text{eff}}\ln2$, genomic **on the line**, slope 1.008, $r$ 0.9997. The curvature side is the most closed part of the program.
- **Cross-instrument invariance (16.8° codec↔v9).** Real but a **coarse aggregate**: stable only at large *representative* $N$ (16.8±0.2 at 20k; a single-clade or biased slice → ~90°). It is the shadow of information-tightness seen through two lenses, not the definition.
- **Q1 — within-clade distance structure (gauge-free Mantel).** Shared structure **decays coarse→fine**: phylum 0.55 → family 0.51 → genus 0.23. The angular coordinate $\theta$ had *overstated* the idiosyncrasy (single-clade $\theta\!\approx\!90°$) because it conflates internal shape with rotational **gauge**; the distance matrix recovers ~half of it. The datum's certified resolution *fades*; there is no cliff.
- **Q2 — grounded vs distorted (patristic ruler).** Fine structure is **partially grounded** (Spearman embedded-vs-patristic ~0.3, both instruments) **and distorted toward labels** (positive genus-margin in *both* — codec +0.18, v9 +0.37 std; v9 > codec in 70% of families). Answer to *"native representation or human-label distortion?"* = **measurably both**; the placement-optimized atlas distorts more. In MDL currency, **label distortion = excess description length above the optimum.**
- **Q3 — quartet tree-fidelity (three-iteration closure to scale-free ρ(S_rel)).** Two-axis decomposition per quadruple: **Signal S** = internal branch (is there structure to resolve?), **Ambiguity ρ** = slack/branch ∈ [0,1] (given structure, is it tree-faithful?). Linchpin: exact tree → ρ=0 at all scales (Buneman), estimator unbiased. Closure finding: at matched relative signal the embedding is a **faithful hyperbolic tree — cleaner than a 5%-distance-noise tree at good signal (ρ 0.014), far below scatter (0.174).** It is **signal-limited, not noise-limited**: the coarse→fine fade is *magnitude* (bits thinning at the tips), while the *shape is self-similar* (relative branch S_rel ≈ flat across ranks). **v9 (dim-129) is more tree-faithful than codec (dim-16) at matched signal** (ρ 0.013 vs 0.020) — dimension buys fidelity, compression injects ambiguity. **Two retractions the cleanup forced:** (i) Q2 and Q3 are *orthogonal*, not one distortion — inflating genus branches keeps it a valid tree, so ρ is *blind* to Q2's label-margin (v9 is a *cleaner tree that over-separates genera*); (ii) the earlier "dial fingers HGT genera" was a threshold artifact — the confidence dial survives (per-genus ρ 0.00–0.16) but its HGT-specific reading does not. **Still the embedding, not the biology** — the clean-tree claim for raw sequence needs tier-2 (k-mer four-point). *This bullet is the record of a self-correcting instrument: two decorations dissolved under the fair test; the load-bearing structure came out harder.*

- **Kernel test — δ ⊥ κ (the Buneman≡Gromov unification, made concrete).** On a synthetic grid (branching *b* × reticulation *f*, ground truth on both axes), the four-point object splits into two invariants: **curvature κ_vol is reticulation-blind** (corr with *b* +0.96, with *f* +0.05; flat to 4 decimals in sparse *f*) and the **tree-defect τ is exactly 0 on pure trees of every curvature** (Buneman holds independent of κ) and tracks *f*. So you can read a *curvature* and a *tree-defect* as separate coordinates of any distance matrix — the operational form of the unification. **Honest residual:** the defect-meter's *gain* depends on curvature (τ↔*b* = −0.35) — directionally exactly the kernel's *second* relation δ_max ∝ 1/√κ (higher curvature permits less slack), but flagged-not-claimed (ρ is scale-free, so τ∝1/√κ is not cleanly isolated). And the prior "angle–depth" version of this test was **retracted as circular** — sin(θ/2)=e^{−√κ·GP} is a theorem true-by-construction in any ℍ²_κ, so on a built embedding it just returns the planted κ. *Confirms the decomposition as mechanism; not a claim about biology (needs native-sequence δ vs an independent HGT signal).*

## V. The honest ledger

**Derived (deductive, machine-checked):** the state equation from MDL; $n=2$ floor+degeneracy; Lyapunov stability. *Certified — but of the model, not of nature.*

**Measured (empirical, holds):** $\kappa$ from volume entropy (genomic on the line), now **two-instrument-validated** (geometry κ\* ⟂ entropy (h·ln2)²; circularity audit clean — Mash distance is entropy-agnostic, $h$ empirical, prediction a separate pipeline; both instruments recover known synthetic truth); the coarse→fine decay of shared structure (Q1 distance-Mantel, Q2 patristic, Q3 ρ); label distortion in both instruments; the δ⊥κ decomposition (curvature blind to reticulation, tree floor exact).

**Circular / semi-circular (flagged, not leaned on):** patristic is GTDB-derived; embedded coords carry the net *and* $\kappa$ baked in (so the angle–depth kernel test is circular on them); $h$ and the geometry are orthogonal statistics of the *same* sequences (the audit's residual — real agreement, but one dataset viewed twice).

**The three live empirical bets — the entire remaining risk:**
1. **Radius is the information/addressing coordinate.** *Supported abductively:* MDL *predicts* E10 (radius = functional complexity, not time), because bits accrue irregularly in time — a former caveat that flipped into a prediction. Now sharpened (via MDL+quartet): radius = addressing cost from the root is near-*definitional* given a divergence metric + a root; the live residual is *which* variable owns the radius when time/divergence/complexity diverge (E10), plus the instrument gap (radial cross-instrument Spearman only 0.46). Not proven.
2. **$\kappa$ is measured independently of $h$.** **AUDIT DONE (clean).** Mash distance is entropy-agnostic; $h$ empirical (log₂α held out as broken control); geometry κ\* and entropy (h·ln2)² are separate pipelines; both recover known synthetic truth. Residual: $h$ and geometry are orthogonal statistics of the *same* sequences — real agreement, not two independent data sources.
3. **The biology is a clean tree at the sequence level.** Needs tier-2: raw k-mer four-point $\delta$, net-free — the one substrate free of both the network and the reference tree. (δ⊥κ decomposition proven as *mechanism* on synthetic data; biological δ=HGT still needs native distances vs an independent transfer signal.)

## VI. The elegance, stated once

> Life is descent-with-modification over a finite alphabet → a branching process → a tree → and a tree embeds in the hyperbolic plane and nowhere smaller. Demand the shortest lossless description, and the flattest hyperbolic curvature that still holds the proliferation is forced: **$\kappa=(h\ln2)^2$**, the alphabet's entropy written as a curvature, at the one dimension ($n=2$) where that law is free of parameters. The coordinate is **canonical where it is information-tight** — one bit of place per bit of biology — it **fades** where function saturates and evolution reticulates, and it **distorts** exactly to the degree we force it to honor discrete labels the data does not natively contain. One equation; the E-series and volume-entropy collapse as its confirmations; Q1/Q2/Q3 as its honest edges.

## VII. Frontier

- **Tier-2 quartet:** raw k-mer four-point $\delta$ (net-free, tree-free) — does the *biology* satisfy the four-point condition, separating it from embedding distortion.
- **$\kappa$-independence audit:** confirm the volume-entropy distances don't smuggle in $h$.
- **The sextant:** Gromov-product placement + local $\delta$ as a production confidence/novelty dial (the operational form of "know when you don't know").

*Filed as a record of the view on 2026-07-26. If it reads as too clean later, §V is where to look first.*
