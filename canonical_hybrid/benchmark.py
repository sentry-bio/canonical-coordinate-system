"""
Benchmark the v1 hybrid system end-to-end.

Tests all three modes:
  1. KESTREL standalone canonical addressing
  2. Atlas standalone canonical addressing
  3. Hybrid tandem (KESTREL coarse → Atlas precise)

Reports R@1, R@10, cell agreement, latency, and speedup.

Run on Inference:
    source /fast/sentrybio/venv/bin/activate
    python -m canonical_hybrid.benchmark \
        --teacher    /home/rohit/e1_results/teacher_coords.npz \
        --phase0_dir /home/rohit/canonical_experiment/ \
        --output_dir /home/rohit/hybrid_v1_benchmark/
"""

from __future__ import annotations

import argparse
import json
import logging
import os
import time
from pathlib import Path

import numpy as np
import torch.nn.functional as F
import torch

from canonical_hybrid.build import build_tessellation
from canonical_hybrid.core import CanonicalHybrid

logging.basicConfig(level=logging.INFO, format="[%(asctime)s] %(message)s",
                    datefmt="%Y-%m-%d %H:%M:%S")
log = logging.getLogger(__name__)

DOMAIN_NAMES = {0: "Bacteria", 1: "Archaea", 2: "Eukaryota"}


def run_benchmark(
    teacher_path: str,
    phase0_dir: str,
    output_dir: str,
    K_values: list = None,
    seed: int = 42,
):
    if K_values is None:
        K_values = [25, 50, 100]

    os.makedirs(output_dir, exist_ok=True)
    output_path = Path(output_dir)

    # ------------------------------------------------------------------
    # Load test data
    # ------------------------------------------------------------------
    log.info("Loading test data...")
    dp = np.load(os.path.join(phase0_dir, "distance_profiles.npz"))
    A_tang = dp["atlas_tangents"].astype(np.float32)
    K_tang = dp["kestrel_tangents"].astype(np.float32)
    N_test = A_tang.shape[0]

    # Match test to teacher for ground-truth indices
    tc = np.load(teacher_path, allow_pickle=False)
    all_tangents = tc["coords"].astype(np.float32)
    all_domains = tc["domain_labels"].astype(np.int32)

    A_t = torch.from_numpy(A_tang)
    all_t = torch.from_numpy(all_tangents)
    test_to_teacher = torch.zeros(N_test, dtype=torch.long)
    for start in range(0, N_test, 200):
        end = min(start + 200, N_test)
        sim = F.cosine_similarity(A_t[start:end].unsqueeze(1), all_t.unsqueeze(0), dim=-1)
        test_to_teacher[start:end] = sim.argmax(dim=1)
    test_to_teacher = test_to_teacher.numpy()
    test_domains = all_domains[test_to_teacher]

    log.info(f"  {N_test} test pairs, {all_tangents.shape[0]} references")

    results = {}

    # ------------------------------------------------------------------
    # Baseline: direct tangent-space NN
    # ------------------------------------------------------------------
    log.info("\n" + "=" * 60)
    log.info("Baseline: Direct tangent-space L2-NN")
    log.info("=" * 60)

    t0 = time.time()
    K_sq = torch.from_numpy(K_tang).pow(2).sum(dim=1)
    A_sq = torch.from_numpy(all_tangents).pow(2).sum(dim=1)
    dot = torch.from_numpy(K_tang) @ torch.from_numpy(all_tangents).T
    l2_mat = (K_sq.unsqueeze(1) + A_sq.unsqueeze(0) - 2 * dot).clamp(min=0).sqrt()
    baseline_nn = l2_mat.argmin(dim=1).numpy()
    baseline_time = time.time() - t0

    baseline_r1 = float((baseline_nn == test_to_teacher).mean())
    # R@10
    top10 = l2_mat.topk(10, dim=1, largest=False).indices.numpy()
    baseline_r10 = float(np.mean([test_to_teacher[i] in top10[i] for i in range(N_test)]))

    log.info(f"  R@1:  {baseline_r1:.1%}")
    log.info(f"  R@10: {baseline_r10:.1%}")
    log.info(f"  Time: {baseline_time:.2f}s ({N_test} queries × {all_tangents.shape[0]} refs)")
    results["baseline"] = {
        "r_at_1": baseline_r1, "r_at_10": baseline_r10,
        "time_s": baseline_time, "mode": "kestrel-direct-NN",
    }

    # Also: Atlas direct NN (upper bound)
    t0 = time.time()
    A_sq_test = torch.from_numpy(A_tang).pow(2).sum(dim=1)
    dot_atlas = torch.from_numpy(A_tang) @ torch.from_numpy(all_tangents).T
    l2_atlas = (A_sq_test.unsqueeze(1) + A_sq.unsqueeze(0) - 2 * dot_atlas).clamp(min=0).sqrt()
    # Exclude self-match
    for i in range(N_test):
        l2_atlas[i, test_to_teacher[i]] = float("inf")
    atlas_nn = l2_atlas.argmin(dim=1).numpy()
    atlas_time = time.time() - t0

    # For Atlas, R@1 = 0% by design (self-excluded), so check R@1 with self included
    l2_atlas_with_self = (A_sq_test.unsqueeze(1) + A_sq.unsqueeze(0) - 2 * dot_atlas).clamp(min=0).sqrt()
    atlas_nn_self = l2_atlas_with_self.argmin(dim=1).numpy()
    atlas_r1_self = float((atlas_nn_self == test_to_teacher).mean())
    log.info(f"\n  Atlas self-retrieval R@1: {atlas_r1_self:.1%} (should be ~100%)")
    results["atlas_self"] = {"r_at_1": atlas_r1_self, "time_s": atlas_time}

    # ------------------------------------------------------------------
    # For each K: build tessellation → benchmark all modes
    # ------------------------------------------------------------------
    for K in K_values:
        log.info(f"\n{'=' * 60}")
        log.info(f"K = {K}")
        log.info(f"{'=' * 60}")

        # Build tessellation
        save_path = str(output_path / f"tessellation_K{K}")
        hybrid = build_tessellation(
            teacher_path=teacher_path, K=K,
            output_path=save_path, seed=seed,
        )

        # --- Mode 1: KESTREL standalone ---
        log.info(f"\n  Mode 1: KESTREL standalone canonical addressing")
        t0 = time.time()
        k_addresses = hybrid.kestrel_address_batch(K_tang)
        k_time = time.time() - t0

        # Cell agreement with Atlas
        a_addresses = hybrid.kestrel_address_batch(A_tang)
        cell_agree = float(np.mean([
            k_addresses[i].cell_id == a_addresses[i].cell_id
            for i in range(N_test)
        ]))
        top3_agree = float(np.mean([
            a_addresses[i].cell_id in k_addresses[i].top_k_cells[:3]
            for i in range(N_test)
        ]))

        # Confidence statistics
        k_conf = [a.confidence for a in k_addresses]
        log.info(f"    Cell agreement:      {cell_agree:.1%}")
        log.info(f"    Top-3 agreement:     {top3_agree:.1%}")
        log.info(f"    Confidence (mean):   {np.mean(k_conf):.3f}")
        log.info(f"    Confidence (median): {np.median(k_conf):.3f}")
        log.info(f"    Time: {k_time:.3f}s ({k_time/N_test*1000:.2f}ms/query)")

        # Domain accuracy from canonical address
        k_dom_correct = float(np.mean([
            k_addresses[i].domain == DOMAIN_NAMES.get(test_domains[i], "?")
            for i in range(N_test)
        ]))
        log.info(f"    Domain accuracy:     {k_dom_correct:.1%}")

        results[f"kestrel_K{K}"] = {
            "mode": "kestrel-standalone",
            "K": K,
            "cell_agreement": cell_agree,
            "top3_agreement": top3_agree,
            "confidence_mean": float(np.mean(k_conf)),
            "confidence_median": float(np.median(k_conf)),
            "domain_accuracy": k_dom_correct,
            "time_s": k_time,
            "ms_per_query": k_time / N_test * 1000,
        }

        # --- Mode 2: Hybrid (1 cell) ---
        for n_cells in [1, 3]:
            log.info(f"\n  Mode 2: Hybrid (n_cells={n_cells})")
            t0 = time.time()
            matches = hybrid.hybrid_search_batch(
                K_tang, A_tang, n_cells=n_cells,
                correct_indices=test_to_teacher,
            )
            h_time = time.time() - t0

            h_r1 = float(np.mean([m.r_at_k.get(1, False) for m in matches]))
            h_r5 = float(np.mean([m.r_at_k.get(5, False) for m in matches]))
            h_r10 = float(np.mean([m.r_at_k.get(10, False) for m in matches]))
            mean_searched = float(np.mean([m.candidates_searched for m in matches]))
            mean_speedup = float(np.mean([m.speedup for m in matches]))

            log.info(f"    R@1:                 {h_r1:.1%}")
            log.info(f"    R@5:                 {h_r5:.1%}")
            log.info(f"    R@10:                {h_r10:.1%}")
            log.info(f"    Candidates searched:  {mean_searched:.0f} / {hybrid.N}")
            log.info(f"    Speedup:             {mean_speedup:.0f}x")
            log.info(f"    Time: {h_time:.2f}s ({h_time/N_test*1000:.2f}ms/query)")

            results[f"hybrid_K{K}_cells{n_cells}"] = {
                "mode": f"hybrid-{n_cells}cell",
                "K": K,
                "n_cells": n_cells,
                "r_at_1": h_r1,
                "r_at_5": h_r5,
                "r_at_10": h_r10,
                "mean_candidates": mean_searched,
                "mean_speedup": mean_speedup,
                "time_s": h_time,
                "ms_per_query": h_time / N_test * 1000,
            }

        # --- Per-domain hybrid breakdown (1 cell) ---
        log.info(f"\n  Per-domain hybrid (K={K}, 1 cell):")
        matches_1cell = hybrid.hybrid_search_batch(
            K_tang, A_tang, n_cells=1, correct_indices=test_to_teacher,
        )
        for d in range(3):
            mask = test_domains == d
            if mask.sum() < 10:
                continue
            d_matches = [matches_1cell[i] for i in range(N_test) if mask[i]]
            d_r1 = float(np.mean([m.r_at_k.get(1, False) for m in d_matches]))
            d_searched = float(np.mean([m.candidates_searched for m in d_matches]))
            log.info(f"    {DOMAIN_NAMES[d]:12s}: R@1={d_r1:.1%}  "
                     f"searched={d_searched:.0f}")

    # ------------------------------------------------------------------
    # Summary table
    # ------------------------------------------------------------------
    log.info(f"\n{'=' * 70}")
    log.info(f"V1 HYBRID BENCHMARK SUMMARY")
    log.info(f"{'=' * 70}")

    log.info(f"\n  {'Mode':<30s} {'R@1':>8s} {'R@10':>8s} {'Searched':>10s} "
             f"{'Speedup':>8s} {'ms/q':>8s}")
    log.info(f"  {'-' * 72}")

    log.info(f"  {'KESTREL direct NN':<30s} {baseline_r1:>8.1%} {baseline_r10:>8.1%} "
             f"{all_tangents.shape[0]:>10d} {'1x':>8s} "
             f"{baseline_time/N_test*1000:>8.1f}")

    for K in K_values:
        for nc in [1, 3]:
            key = f"hybrid_K{K}_cells{nc}"
            if key not in results:
                continue
            r = results[key]
            label = f"Hybrid K={K} ({nc} cell{'s' if nc > 1 else ''})"
            log.info(f"  {label:<30s} {r['r_at_1']:>8.1%} {r['r_at_10']:>8.1%} "
                     f"{r['mean_candidates']:>10.0f} "
                     f"{r['mean_speedup']:>7.0f}x "
                     f"{r['ms_per_query']:>8.1f}")

    log.info(f"\n  Atlas self-retrieval: {atlas_r1_self:.1%}")

    # ------------------------------------------------------------------
    # Save
    # ------------------------------------------------------------------
    def to_serializable(obj):
        if isinstance(obj, (np.floating, np.integer)):
            return float(obj)
        if isinstance(obj, np.ndarray):
            return obj.tolist()
        if isinstance(obj, dict):
            return {k: to_serializable(v) for k, v in obj.items()}
        return obj

    with open(output_path / "benchmark_results.json", "w") as f:
        json.dump(to_serializable(results), f, indent=2)
    log.info(f"\nResults → {output_path / 'benchmark_results.json'}")


def main():
    parser = argparse.ArgumentParser(description="Benchmark v1 hybrid system")
    parser.add_argument("--teacher", required=True)
    parser.add_argument("--phase0_dir", required=True)
    parser.add_argument("--output_dir", default="hybrid_v1_benchmark")
    parser.add_argument("--K", type=int, nargs="+", default=[25, 50, 100])
    parser.add_argument("--seed", type=int, default=42)
    args = parser.parse_args()

    run_benchmark(
        teacher_path=args.teacher,
        phase0_dir=args.phase0_dir,
        output_dir=args.output_dir,
        K_values=args.K,
        seed=args.seed,
    )


if __name__ == "__main__":
    main()
