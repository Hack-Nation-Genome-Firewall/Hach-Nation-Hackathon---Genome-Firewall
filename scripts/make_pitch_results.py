#!/usr/bin/env python
"""
One-shot pitch-results runner. Fire this the moment real features.csv exists.

Trains + evaluates the Track B model on BOTH splits (homology-grouped and
temporal), collects the rubric metrics into one table, and writes a
grouped-vs-temporal comparison figure. Reuses Track B's own train.py /
evaluate.py via subprocess (no coupling to internals).

  python scripts/make_pitch_results.py \
      --features data/manifests/features.csv \
      --labels   data/manifests/labels.csv \
      --spec     data/manifests/feature_spec.json \
      --grouped  data/manifests/split_manifest.csv \
      --temporal data/manifests/split_manifest_temporal.csv \
      --out results/

Outputs (under --out):
  models/kp_grouped.joblib, models/kp_temporal.joblib
  eval_grouped/*, eval_temporal/*      (Track B evaluate.py output)
  pitch_metrics.csv                    (both splits, all drugs, side by side)
  fig_grouped_vs_temporal.png          (balanced-accuracy comparison)
"""
from __future__ import annotations
import argparse, json, subprocess, sys
from pathlib import Path

import pandas as pd


def run(cmd):
    print("::", " ".join(str(c) for c in cmd), flush=True)
    subprocess.run([str(c) for c in cmd], check=True)


def train_eval(tag, splits, args, outdir):
    bundle = outdir / "models" / f"kp_{tag}.joblib"
    bundle.parent.mkdir(parents=True, exist_ok=True)
    evdir = outdir / f"eval_{tag}"
    run([sys.executable, "-m", "module2_predictor.train",
         "--features", args.features, "--labels", args.labels,
         "--splits", splits, "--spec", args.spec,
         "--output", bundle, "--model-version", f"kp-real-{tag}-v1"])
    run([sys.executable, "-m", "module2_predictor.evaluate",
         "--features", args.features, "--labels", args.labels,
         "--splits", splits, "--spec", args.spec,
         "--bundle", bundle, "--output-dir", evdir])
    # evaluate.py writes a per-drug metrics file; find the csv/json it produced
    for cand in ("metrics.csv", "overall_metrics.csv", "metrics.json"):
        p = evdir / cand
        if p.exists():
            df = pd.read_csv(p) if p.suffix == ".csv" else pd.DataFrame(json.load(open(p)))
            df["split"] = tag
            return df
    print(f"[warn] no metrics file found under {evdir}; check its contents")
    return pd.DataFrame()


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--features", required=True)
    ap.add_argument("--labels", required=True)
    ap.add_argument("--spec", required=True)
    ap.add_argument("--grouped", required=True)
    ap.add_argument("--temporal", required=True)
    ap.add_argument("--out", default="results")
    args = ap.parse_args()
    outdir = Path(args.out); outdir.mkdir(parents=True, exist_ok=True)

    frames = []
    for tag, splits in (("grouped", args.grouped), ("temporal", args.temporal)):
        try:
            frames.append(train_eval(tag, splits, args, outdir))
        except subprocess.CalledProcessError as e:
            print(f"[error] {tag} split failed: {e}")

    if frames:
        allm = pd.concat(frames, ignore_index=True)
        allm.to_csv(outdir / "pitch_metrics.csv", index=False)
        print("\n=== PITCH METRICS (both splits) ===")
        print(allm.to_string(index=False))

        # comparison figure if a balanced-accuracy column is present
        bcol = next((c for c in allm.columns if "balanced" in c.lower()), None)
        dcol = next((c for c in allm.columns if c.lower() in ("drug", "antibiotic")), None)
        if bcol and dcol:
            import matplotlib
            matplotlib.use("Agg")
            import matplotlib.pyplot as plt
            import numpy as np
            piv = allm.pivot_table(index=dcol, columns="split", values=bcol)
            fig, ax = plt.subplots(figsize=(7, 4))
            x = np.arange(len(piv)); w = 0.38
            ax.bar(x - w/2, piv.get("grouped", 0), w, label="homology-grouped")
            ax.bar(x + w/2, piv.get("temporal", 0), w, label="temporal (future isolates)")
            ax.set_xticks(x); ax.set_xticklabels(piv.index, rotation=15, ha="right")
            ax.set_ylabel("balanced accuracy"); ax.set_ylim(0, 1)
            ax.set_title("Held-out performance: homology-grouped vs temporal split")
            ax.legend(frameon=False, fontsize=9)
            fig.tight_layout()
            fig.savefig(outdir / "fig_grouped_vs_temporal.png", dpi=200, bbox_inches="tight")
            print(f"wrote {outdir/'fig_grouped_vs_temporal.png'}")
    print("\nDONE. See", outdir / "pitch_metrics.csv")


if __name__ == "__main__":
    main()
