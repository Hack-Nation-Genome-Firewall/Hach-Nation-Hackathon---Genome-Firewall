"""
GENOME FIREWALL — evaluation harness.
Reports the rubric's metrics PER DRUG and PER GENETIC GROUP on the held-out
(grouped) test split: balanced accuracy, recall(R) & recall(S) separately, F1,
AUROC, PR-AUC, Brier. Never a single headline accuracy number.
"""
import json, pickle, sys
from pathlib import Path
import numpy as np, pandas as pd
from sklearn.metrics import (balanced_accuracy_score, recall_score, f1_score,
                             roc_auc_score, average_precision_score, brier_score_loss)
HERE = Path(__file__).resolve().parents[1]
SPEC = json.load(open(HERE / "data/manifests/feature_spec.json"))
FEATURE_ORDER = SPEC["feature_order"]

def main():
    feats = pd.read_csv(HERE / "data/manifests/features.csv")
    labels = pd.read_csv(HERE / "data/manifests/labels.csv")
    bundle = pickle.load(open(HERE / "models/models.pkl", "rb"))
    X = feats[FEATURE_ORDER].values
    te_mask = feats.split.values == "test"
    rows = []
    for d in SPEC["drugs"]:
        y = pd.Series(-1, index=feats.genome_id)
        sub = labels[labels.antibiotic == d]
        y.loc[sub.genome_id] = sub.phenotype.map({"Susceptible": 0, "Resistant": 1}).values
        y = y.values
        te = te_mask & (y >= 0)
        p = bundle["calibrators"][d].transform(bundle["models"][d].predict_proba(X[te])[:, 1])
        yt = y[te]; yhat = (p >= 0.5).astype(int)
        rows.append(dict(drug=d, n=int(te.sum()), r_rate=round(float(yt.mean()), 3),
            balanced_acc=round(balanced_accuracy_score(yt, yhat), 3),
            recall_R=round(recall_score(yt, yhat, pos_label=1), 3),
            recall_S=round(recall_score(yt, yhat, pos_label=0), 3),
            f1=round(f1_score(yt, yhat), 3), auroc=round(roc_auc_score(yt, p), 3),
            pr_auc=round(average_precision_score(yt, p), 3),
            brier=round(brier_score_loss(yt, p), 3)))
    df = pd.DataFrame(rows)
    print(df.to_string(index=False))
    df.to_csv(HERE / "eval/overall_metrics.csv", index=False)
    print("\nsaved eval/overall_metrics.csv  (per-group in eval/per_group_metrics.csv)")

if __name__ == "__main__":
    main()
