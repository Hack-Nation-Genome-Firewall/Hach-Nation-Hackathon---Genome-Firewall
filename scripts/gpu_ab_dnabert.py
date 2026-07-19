#!/usr/bin/env python
"""
Stretch experiment: DNABERT-2 genome-embedding baseline vs. AMRFinderPlus
features on ONE drug (ciprofloxacin), on the SAME homology-grouped split.

Purpose is a single defensible A/B number for the Tech video:
"we tested a genome language model; our interpretable calibrated feature model
matched/beat it." It is NOT meant to replace the pipeline.

Method (deliberately simple, bounded to run overnight on one GPU):
  * for each genome, download its FASTA, tile into WINDOW-bp non-overlapping
    windows, sample up to N_WIN windows, embed each with DNABERT-2 (mean of
    token embeddings), then mean-pool windows -> one genome vector;
  * train LogisticRegression on the TRAIN split, evaluate balanced accuracy on
    the TEST split;
  * print DNABERT-2 balanced acc next to the feature-model balanced acc for the
    same drug/split (pass --feature-balacc, or it reads results/pitch_metrics.csv).

Caveats (state these in the video, don't hide them): windows are random tiles,
not gene-aware; a subset of genomes is used for speed; mean-pooling is a crude
aggregator. This is a baseline probe, not a tuned model.
"""
from __future__ import annotations
import argparse, io, sys, time, urllib.request
from pathlib import Path
import numpy as np
import pandas as pd


def log(*a): print(*a, flush=True)


def fetch_fasta(gid, timeout=90):
    url = f"https://www.bv-brc.org/api/genome_sequence/?eq(genome_id,{gid})&http_accept=application/dna+fasta"
    for attempt in range(3):
        try:
            return urllib.request.urlopen(urllib.request.Request(url), timeout=timeout).read().decode()
        except Exception:
            if attempt == 2: return None
            time.sleep(2)


def concat_seq(fasta_text):
    return "".join(l.strip() for l in fasta_text.splitlines() if l and not l.startswith(">")).upper()


def windows(seq, w, n, rng):
    idx = list(range(0, max(1, len(seq) - w), w))
    if not idx: return [seq[:w]]
    if len(idx) > n: idx = list(rng.choice(idx, size=n, replace=False))
    return [seq[i:i + w] for i in idx]


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--labels", default="data/manifests/labels.csv")
    ap.add_argument("--splits", default="data/manifests/split_manifest.csv")
    ap.add_argument("--drug", default="ciprofloxacin")
    ap.add_argument("--window", type=int, default=10000)
    ap.add_argument("--n-windows", type=int, default=24)
    ap.add_argument("--max-train", type=int, default=400)
    ap.add_argument("--max-test", type=int, default=150)
    ap.add_argument("--model", default="zhihan1996/DNABERT-2-117M")
    ap.add_argument("--batch", type=int, default=16)
    ap.add_argument("--feature-balacc", type=float, default=None,
                    help="feature-model balanced acc for this drug/split, for the printout")
    ap.add_argument("--out", default="results/gpu_ab_ciprofloxacin.txt")
    args = ap.parse_args()

    import torch
    from transformers import AutoTokenizer, AutoModel
    from sklearn.linear_model import LogisticRegression
    from sklearn.metrics import balanced_accuracy_score

    dev = "cuda" if torch.cuda.is_available() else "cpu"
    log(f"[ab] device={dev}  drug={args.drug}  model={args.model}")
    if dev == "cpu":
        log("[ab] WARNING: no GPU visible; this will be very slow.")

    lab = pd.read_csv(args.labels, dtype={"genome_id": str})
    lab = lab[(lab.antibiotic == args.drug) & (lab.phenotype.isin(["Resistant", "Susceptible"]))]
    lab = lab[["genome_id", "label"]].drop_duplicates("genome_id")
    spl = pd.read_csv(args.splits, dtype={"genome_id": str})[["genome_id", "split"]]
    df = lab.merge(spl, on="genome_id")
    rng = np.random.default_rng(42)
    tr = df[df.split == "train"].sample(min(args.max_train, (df.split == "train").sum()), random_state=1)
    te = df[df.split == "test"].sample(min(args.max_test, (df.split == "test").sum()), random_state=1)
    log(f"[ab] train={len(tr)}  test={len(te)}  (subsampled)")

    tok = AutoTokenizer.from_pretrained(args.model, trust_remote_code=True)
    mdl = AutoModel.from_pretrained(args.model, trust_remote_code=True).to(dev).eval()

    def embed_windows(wins):
        vecs = []
        for i in range(0, len(wins), args.batch):
            chunk = wins[i:i + args.batch]
            enc = tok(chunk, return_tensors="pt", padding=True, truncation=True, max_length=512).to(dev)
            with torch.no_grad():
                out = mdl(**enc)[0]                      # (B, T, H)
                mask = enc["attention_mask"].unsqueeze(-1)
                pooled = (out * mask).sum(1) / mask.sum(1).clamp(min=1)
            vecs.append(pooled.cpu().numpy())
        return np.vstack(vecs).mean(0)                   # mean over windows -> (H,)

    def genome_vecs(frame, tag):
        X, y, keep = [], [], 0
        for gid, lb in zip(frame.genome_id, frame.label):
            fa = fetch_fasta(gid)
            if not fa: continue
            seq = concat_seq(fa)
            if len(seq) < args.window: continue
            X.append(embed_windows(windows(seq, args.window, args.n_windows, rng)))
            y.append(int(lb)); keep += 1
            if keep % 25 == 0: log(f"[ab] {tag} embedded {keep}/{len(frame)}")
        return np.array(X), np.array(y)

    t0 = time.time()
    Xtr, ytr = genome_vecs(tr, "train")
    Xte, yte = genome_vecs(te, "test")
    log(f"[ab] embedding done in {(time.time()-t0)/60:.1f} min  Xtr={Xtr.shape} Xte={Xte.shape}")

    clf = LogisticRegression(max_iter=2000, class_weight="balanced").fit(Xtr, ytr)
    dnabert_ba = balanced_accuracy_score(yte, clf.predict(Xte))

    feat = args.feature_balacc
    if feat is None:
        mp = Path("results/pitch_metrics.csv")
        if mp.exists():
            m = pd.read_csv(mp)
            row = m[(m.get("drug") == args.drug) & (m.get("split") == "grouped")]
            if len(row) and "balanced_accuracy" in row:
                feat = float(row.balanced_accuracy.iloc[0])

    lines = [
        f"A/B on {args.drug} (grouped test split, n_test={len(yte)})",
        f"  DNABERT-2 genome-embedding + LogReg : balanced accuracy = {dnabert_ba:.3f}",
        f"  AMRFinderPlus features (our model)  : balanced accuracy = "
        + (f"{feat:.3f}" if feat is not None else "(run make_pitch_results.py first)"),
        "",
        "Caveat: random-tile windows, subsampled genomes, mean-pool aggregation;",
        "baseline probe, not a tuned model.",
    ]
    Path(args.out).parent.mkdir(parents=True, exist_ok=True)
    Path(args.out).write_text("\n".join(lines) + "\n")
    log("\n" + "\n".join(lines))
    log(f"[ab] wrote {args.out}")


if __name__ == "__main__":
    main()
