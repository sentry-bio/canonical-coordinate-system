#!/usr/bin/env python3
"""
Phase 0 GATE for solid-ground check #4 (cross-version harmonization).

Reproduce the DEPLOYED v9 embedding pipeline faithfully, then verify that
re-embedded v9 lands reference genomes back onto the deployed karcher_v7 coords
at near-zero angular residual.

Recipe (copied EXACTLY from extract_karcher.py, whose geometry we IMPORT):
    windows -> encode_angular_only  (ball points) -> Karcher/Frechet mean on ball.
The DirectRadialHead never touches coords (advisory radius only), so it is omitted.

Also computes the OLD Euclidean mean-pool aggregation as a diagnostic: if that
lands ~30-40 deg off karcher_v7 while Karcher lands ~0, we have PROVEN the 39 deg
cross-version confound was the aggregation swap, not a real embedding instability.

Runs on the inference box. GPU is hard-capped to protect live production.
"""
import sys, os, numpy as np, torch
sys.path.insert(0, "/fast/sentrybio/v9_training")
import extract_karcher as EK                      # exact geometry + karcher_mean
from model_v15_5 import load_v15_5_model

CKPT = "/fast/sentrybio/v9_training/checkpoints/best.pt"
TOK  = "/fast/sentrybio/v9_training/tokenized"
REF  = "/fast/sentrybio/v9_karcher/karcher_v7.npz"
N        = 250
ENC_BS   = 32                                      # sub-batch for encode (bounds GPU mem)
MAXWIN   = 100                                     # match recipe
ANCHORS  = ["GCF_000005845.2", "GCF_000091665.1"]  # E.coli (prime meridian), M.jannaschii (chirality)

device = "cuda" if torch.cuda.is_available() else "cpu"
if device == "cuda":
    torch.cuda.set_per_process_memory_fraction(0.30, 0)   # ~1.8GB hard cap -> protect production

ref       = np.load(REF, allow_pickle=True)
ref_acc   = [str(a) for a in ref["accessions"]]
ref_coords= ref["coords"]
kappa     = float(ref["kappa"])
acc2idx   = {a: i for i, a in enumerate(ref_acc)}
print(f"[load] karcher_v7 coords {ref_coords.shape}  kappa {kappa:.6f}  ball_radius {EK.ball_radius(kappa):.4f}")

def tok_path(acc):
    p = os.path.join(TOK, acc + ".npy")
    if os.path.exists(p): return p
    b = acc.rsplit(".", 1)[0]
    p2 = os.path.join(TOK, b + ".npy")
    return p2 if os.path.exists(p2) else None

# select: anchors first (needed for datum-theta), then first N-2 with tokens present
sel = []
for a in ANCHORS:
    if a in acc2idx and tok_path(a): sel.append(a)
anchor_ok = len(sel) == len(ANCHORS)
for a in ref_acc:
    if len(sel) >= N: break
    if a in sel: continue
    if tok_path(a): sel.append(a)
print(f"[select] {len(sel)} genomes ({'anchors PRESENT' if anchor_ok else 'ANCHORS MISSING'})")

model = load_v15_5_model(CKPT, device=device); model.eval()
print(f"[model] loaded live_kappa={getattr(model,'live_kappa',None)}")

def embed_genome(acc, rng):
    arr = np.load(tok_path(acc)).astype(np.int32)
    if arr.ndim == 1: arr = arr[None, :]
    n_win = min(arr.shape[0], MAXWIN)
    wt = []
    for wi in range(n_win):
        tokens = arr[wi]
        if len(tokens) > 512:
            start = int(rng.integers(0, len(tokens) - 512 + 1))
            tokens = tokens[start:start + 512]
        else:
            tokens = tokens[:512]
        if len(tokens) < 512:
            tokens = np.pad(tokens, (0, 512 - len(tokens)))
        wt.append(tokens)
    stacked = np.stack(wt)
    zs = []
    for j in range(0, len(stacked), ENC_BS):
        batch = torch.from_numpy(stacked[j:j+ENC_BS]).long().to(device)
        with torch.no_grad(), torch.amp.autocast("cuda"):
            zs.append(model.encode_angular_only(batch).float().cpu())
    z = torch.cat(zs, 0)                                  # (n_win, 129) ball points, on CPU
    if z.shape[0] > 1:
        mu, _ = EK.karcher_mean(z, kappa=kappa)           # EXACT deployed aggregation
        eucl  = z.mean(0)                                 # OLD confounded aggregation
    else:
        mu = z[0]; eucl = z[0]
    win_len = arr.shape[1] if arr.ndim == 2 else len(arr)
    return mu.numpy(), eucl.numpy(), int(win_len)

rngA = np.random.default_rng(0)   # draw A (random 512-slices)
rngB = np.random.default_rng(1)   # draw B (independent slices) -> noise floor
CA, CB, CE, CR, WL = [], [], [], [], []
for k, a in enumerate(sel):
    muA, euA, wl = embed_genome(a, rngA)
    muB, _,  _   = embed_genome(a, rngB)
    CA.append(muA); CB.append(muB); CE.append(euA); CR.append(ref_coords[acc2idx[a]]); WL.append(wl)
    if (k + 1) % 50 == 0: print(f"  embedded {k+1}/{len(sel)}")
CA, CB, CE, CR = map(np.array, (CA, CB, CE, CR))
WL = np.array(WL)

def dir_angle(X, Y):
    xn = X / np.linalg.norm(X, axis=1, keepdims=True)
    yn = Y / np.linalg.norm(Y, axis=1, keepdims=True)
    c = (xn * yn).sum(1).clip(-1, 1)
    return np.degrees(np.arccos(c))

ang_repro = dir_angle(CA, CR)   # Karcher reproduction vs deployed
ang_noise = dir_angle(CA, CB)   # slice-noise floor (Karcher A vs Karcher B)
ang_eucl  = dir_angle(CE, CR)   # OLD Euclidean-mean vs deployed (the confound)

print("\n" + "=" * 74)
print("PHASE 0 - 129-d direction agreement (degrees)")
print("=" * 74)
def rpt(name, v):
    print(f"{name:46s} median {np.median(v):6.2f}  p90 {np.percentile(v,90):6.2f}  "
          f"max {v.max():6.2f}  <2deg {100*(v<2).mean():4.0f}%  <5deg {100*(v<5).mean():4.0f}%")
rpt("Karcher-repro  <->  deployed karcher_v7", ang_repro)
rpt("slice-noise floor  (Karcher A <-> B)", ang_noise)
rpt("Euclidean-mean <-> karcher_v7 [OLD/confound]", ang_eucl)
print(f"\nwindow-length: {int((WL>512).sum())}/{len(WL)} genomes have >512-token windows "
      f"(random-slice active); unique lens {sorted(set(WL.tolist()))[:6]}")

# datum-theta residual (the claim's own metric), in the anchored 2D backbone frame
if anchor_ok:
    ei, mi = sel.index(ANCHORS[0]), sel.index(ANCHORS[1])
    def logmap0(x, K):
        n = np.linalg.norm(x, axis=1, keepdims=True).clip(1e-9, EK.ball_radius(K) * (1 - 1e-6))
        return (2 / np.sqrt(K)) * np.arctanh(np.sqrt(K) * n) * (x / n)
    def datum_theta(C, K):
        T = logmap0(C, K); mu = T.mean(0)
        _, _, Vt = np.linalg.svd(T - mu, full_matrices=False)
        P = (T - mu) @ Vt[:2].T
        th = np.arctan2(P[:, 1], P[:, 0]); th = (th - th[ei] + np.pi) % (2*np.pi) - np.pi
        if th[mi] < 0: th = -th
        return th
    tR, tA, tE = datum_theta(CR, kappa), datum_theta(CA, kappa), datum_theta(CE, kappa)
    d_ra = np.degrees(np.abs((tA - tR + np.pi) % (2*np.pi) - np.pi))
    d_re = np.degrees(np.abs((tE - tR + np.pi) % (2*np.pi) - np.pi))
    print(f"\ndatum-theta residual vs karcher_v7:  Karcher-repro median {np.median(d_ra):5.2f}deg   "
          f"Euclidean[OLD] median {np.median(d_re):5.2f}deg")

# verdict
med = np.median(ang_repro)
print("\n" + "=" * 74)
if med < 5:
    print(f"GATE PASS: Karcher reproduction lands on deployed karcher_v7 (median {med:.2f} deg < 5).")
    print("Pipeline is faithful -> proceed to Phase 2 (apply SAME harness to v10.6).")
elif med < 12:
    print(f"GATE MARGINAL: median {med:.2f} deg. Reproduced but with residual - inspect noise floor.")
else:
    print(f"GATE FAIL: median {med:.2f} deg. Reproduction does not match deployed - deeper diff remains.")
print("=" * 74)

np.savez("/fast/sentrybio/v9_karcher/phase0_repro.npz",
         sel=np.array(sel), CA=CA, CB=CB, CE=CE, CR=CR, WL=WL,
         ang_repro=ang_repro, ang_noise=ang_noise, ang_eucl=ang_eucl, kappa=kappa)
print("[saved] /fast/sentrybio/v9_karcher/phase0_repro.npz")
