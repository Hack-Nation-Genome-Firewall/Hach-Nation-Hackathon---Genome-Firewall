"""
GENOME FIREWALL — Module 2 training (REFERENCE).

Per-drug regularized logistic regression on AMR-marker features, trained under a
HOMOLOGY-GROUPED split (whole clusters in one split), calibrated with isotonic
regression on a dedicated calibration split. Emits models/models.pkl consumed by
predict.py and the app.

Real-data swap: replace features.csv / labels.csv with AMRFinderPlus output +
BV-BRC lab-measured labels, and split_manifest.csv with a Mash/sourmash grouped
split. Nothing else changes.
"""
import json, pickle, warnings
from pathlib import Path
import numpy as np, pandas as pd
from sklearn.linear_model import LogisticRegression
from sklearn.isotonic import IsotonicRegression
warnings.filterwarnings("ignore")

HERE = Path(__file__).resolve().parents[1]
SPEC = json.load(open(HERE / "data/manifests/feature_spec.json"))
DRUGS = SPEC["drugs"]
FEATURE_ORDER = SPEC["feature_order"]   # markers + target-gene flags — same order predict.py uses
NOCALL_BAND = 0.15   # |p-0.5|*2 < BAND -> no-call; tune on calibration split

def main():
    feats = pd.read_csv(HERE / "data/manifests/features.csv")
    labels = pd.read_csv(HERE / "data/manifests/labels.csv")
    X = feats[FEATURE_ORDER].values
    split = feats.split.values
    models, calibrators = {}, {}
    for d in DRUGS:
        y = pd.Series(-1, index=feats.genome_id)
        sub = labels[labels.antibiotic == d]
        y.loc[sub.genome_id] = sub.phenotype.map({"Susceptible": 0, "Resistant": 1}).values
        y = y.values
        tested = y >= 0
        tr = tested & (split == "train")
        ca = tested & (split == "calibration")
        clf = LogisticRegression(C=1.0, class_weight="balanced", max_iter=3000).fit(X[tr], y[tr])
        iso = IsotonicRegression(out_of_bounds="clip").fit(clf.predict_proba(X[ca])[:, 1], y[ca])
        models[d], calibrators[d] = clf, iso
        print(f"trained {d}: n_train={tr.sum()} n_cal={ca.sum()}")
    (HERE / "models").mkdir(exist_ok=True)
    with open(HERE / "models/models.pkl", "wb") as f:
        pickle.dump({"models": models, "calibrators": calibrators,
                     "feature_order": SPEC["feature_order"],
                     "drug_targets": SPEC["drug_targets"], "nocall_band": NOCALL_BAND}, f)
    print("saved models/models.pkl")

if __name__ == "__main__":
    main()
