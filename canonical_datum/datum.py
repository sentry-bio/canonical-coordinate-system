"""
The DATUM — the frozen form. Not fitted; theorem-grounded constants.

A geodetic datum (WGS84) is minimal, stable, and invariant — it defines the frame everything else
composes on, and it does NOT know terrain, weather, or biology. This is the genomic analogue: the
low-bias coordinate frame for the vertical history of life. Its power is parsimony done precisely.

Status: VALIDATED CANDIDATE (not certified). See docs/EVIDENCE.md — κ, n=2, and cross-instrument
angular invariance are verified; cross-version harmonization and external grounding are OPEN.
"""
import numpy as np

KAPPA = 5.0 / 4.0                      # curvature = (h·ln2)^2 from molecular-alphabet entropy h≈1.61 bits/nt.
                                       #   Two-instrument-validated (geometry κ* ≈ entropy (h·ln2)^2); live ≈1.237.
INTRINSIC_DIM = 2                      # n=2: MEASURED (2.01±0.06), and the HINGE — at n=2, (n-1)=1, so the law is
                                       #   parameter-free. Also the minimal faithful tree embedding (Sarkar).
BALL_RADIUS = 1.0 / np.sqrt(KAPPA)     # Poincaré-ball radius at this curvature

OBJECTIVE = "quartet"                  # the four-point condition — the ATOM of both the tree (Buneman) and the
                                       #   curvature (Gromov). The datum is specified INTRINSICALLY on quartets.
REFERENCE_TREE = "GTDB_r214"           # the versioned CONTENT the form is realized on
TOKENIZER = "BPE-4096"                 # PART OF THE FROZEN FRAME: the tokenizer moves θ ~12.7° (Phase 0), so it is
                                       #   a datum parameter, not an implementation detail. Freeze it with κ.

PRIME_MERIDIAN = "GCF_000005845.2"     # E. coli K-12  -> canonical θ = 0
CHIRALITY_ANCHOR = "GCF_000091665.1"   # M. jannaschii -> fixes handedness (the SO(2) gauge)

CONFORMANCE_TOLERANCE = 0.90           # determinate-quartet agreement required to CERTIFY (freeze) the frame
DATUM_VERSION = "candidate-1.0"        # not frozen until the freeze-gate is met (see docs/EVIDENCE.md)

# what the datum certifies vs. carries vs. flags -- the honest tiering:
CANONICAL_PRIMITIVE = "distance"       # the canonical object is the INTRINSIC (gauge-free) distance/quartet
                                       #   structure. θ below is a DERIVED chart with known gauge fragility.
CERTIFIED_AXIS = "theta"               # angular (lineage): cross-instrument-shared — but a FRAGILE anchored chart
                                       #   (per-set-SVD θ collapses on biased samples; use intrinsic distances).
ADVISORY_AXIS = "radius"               # radial (depth): clade-coupled, not an independent canonical law
# RESOLUTION BOUNDARY (measured this iteration): the datum is a coarse-shared frame whose certified resolution
# FADES with scale — it does not have a clean cliff. Within-clade distance structure is shared strongly at
# phylum/family (Mantel ~0.55/0.51) and thins through genus (~0.23); the embedded tree is faithful (ρ≈0 vs
# scatter 0.18) wherever there is branch SIGNAL, and goes SIGNAL-limited (not corrupt) toward the tips.
CERTIFIED_RESOLUTION = "phylum→family (shared, tree-faithful); fades through genus; advisory at species/strain"
FREEZE_GATE = "angular median < ~10deg AND determinate-quartet agreement >= 0.90"


def summary():
    return {"kappa": KAPPA, "intrinsic_dim": INTRINSIC_DIM, "ball_radius": round(BALL_RADIUS, 4),
            "objective": OBJECTIVE, "reference_tree": REFERENCE_TREE, "tokenizer": TOKENIZER,
            "prime_meridian": PRIME_MERIDIAN, "chirality_anchor": CHIRALITY_ANCHOR,
            "version": DATUM_VERSION, "canonical_primitive": CANONICAL_PRIMITIVE,
            "certified_axis": CERTIFIED_AXIS, "advisory_axis": ADVISORY_AXIS,
            "certified_resolution": CERTIFIED_RESOLUTION,
            "status": "VALIDATED CANDIDATE — see docs/EVIDENCE.md"}
