"""
GENOME FIREWALL — reproduce the evaluation figures.

Produces:
  eval/fig_leakage.png     random-split vs homology-grouped-split balanced accuracy
  eval/fig_reliability.png per-drug calibration reliability curves (held-out test)

Run after train.py:  python eval/make_figures.py
Reads the frozen data/manifests/* and models/models.pkl. No seaborn/global-state
plotting — uses explicit Figure objects so results are reproducible.
"""
import json, pickle
from pathlib import Path
import numpy as np, pandas as pd
import matplotlib.pyplot as plt
from sklearn.linear_model import LogisticRegression
from sklearn.metrics import balanced_accuracy_score
from sklearn.calibration import calibration_curve
from sklearn.metrics import brier_score_loss

HERE = Path(__file__).resolve().parents[1]
SPEC = json.load(open(HERE / "data/manifests/feature_spec.json"))
DRUGS = SPEC["drugs"]
FEATURE_ORDER = SPEC["feature_order"]

META_GREY = "#8a8a8a"; RED = "#c44e52"; BLUE = "#4c72b0"


def load():
    feats = pd.read_csv(HERE / "data/manifests/features.csv")
    labels = pd.read_csv(HERE / "data/manifests/labels.csv")
    bundle = pickle.load(open(HERE / "models/models.pkl", "rb"))
    return feats, labels, bundle


def y_for(feats, labels, drug):
    y = pd.Series(-1, index=feats.genome_id)
    sub = labels[labels.antibiotic == drug]
    y.loc[sub.genome_id] = sub.phenotype.map({"Susceptible": 0, "Resistant": 1}).values
    return y.values


def fig_leakage(feats, labels):
    """Compare a random row split (leaks) vs the frozen grouped split (honest)."""
    X = feats[FEATURE_ORDER].values
    cluster = feats.cluster_id.values
    split = feats.split.values
    train_cal = np.isin(split, ["train", "calibration"])
    test = split == "test"
    rand_res, grp_res = {}, {}
    for d in DRUGS:
        y = y_for(feats, labels, d); tested = y >= 0
        idx = np.where(tested)[0]
        # random split (WRONG): rows shuffled, clusters ignored
        rr = np.random.default_rng(0); perm = rr.permutation(idx); cut = int(0.75 * len(perm))
        def fe(tr, te):
            m = LogisticRegression(C=1.0, class_weight="balanced", max_iter=2000).fit(X[tr], y[tr])
            return balanced_accuracy_score(y[te], (m.predict_proba(X[te])[:, 1] >= 0.5).astype(int))
        rand_res[d] = fe(perm[:cut], perm[cut:])
        # grouped split (RIGHT): test = held-out clusters
        tr = idx[train_cal[idx]]; te = idx[test[idx]]
        grp_res[d] = fe(tr, te)

    deltas = {d: (rand_res[d] - grp_res[d]) * 100 for d in DRUGS}
    n_inflated = sum(v > 0 for v in deltas.values())
    mean_delta = np.mean(list(deltas.values()))

    fig, ax = plt.subplots(figsize=(5.4, 3.4))
    x = np.arange(len(DRUGS)); w = 0.36
    ax.bar(x - w/2, [rand_res[d] for d in DRUGS], w, label="Random split (leaks)", color=RED)
    ax.bar(x + w/2, [grp_res[d] for d in DRUGS], w, label="Grouped split (honest)", color=BLUE)
    for i, d in enumerate(DRUGS):
        dv = deltas[d]
        lbl = f"+{dv:.0f}pt" if dv >= 0 else f"\u2212{abs(dv):.0f}pt"
        ax.annotate(lbl, (i, max(rand_res[d], grp_res[d]) + 0.02), ha="center", fontsize=6, color=RED)
    ax.set_xticks(x); ax.set_xticklabels(DRUGS, fontsize=6.5); ax.set_ylim(0, 1)
    ax.set_ylabel("Balanced accuracy")
    ax.set_title(f"Random split vs. homology-grouped split\n"
                 f"random inflates on {n_inflated} of {len(DRUGS)} drugs; mean {mean_delta:+.1f} pt",
                 fontsize=9)
    ax.axhline(0.5, ls=":", c=META_GREY, lw=0.8)
    ax.legend(frameon=False, fontsize=6, loc="upper right")
    fig.tight_layout()
    fig.savefig(HERE / "eval/fig_leakage.png", dpi=200)
    print(f"fig_leakage.png  deltas={ {d: round(v,1) for d,v in deltas.items()} }")


def fig_reliability(feats, labels, bundle):
    X = feats[FEATURE_ORDER].values
    test = feats.split.values == "test"
    fig, axes = plt.subplots(1, len(DRUGS), figsize=(2.8 * len(DRUGS), 3.0))
    for ax, d in zip(np.atleast_1d(axes), DRUGS):
        y = y_for(feats, labels, d); te = test & (y >= 0)
        p = bundle["calibrators"][d].transform(bundle["models"][d].predict_proba(X[te])[:, 1])
        yt = y[te]
        frac, mean_pred = calibration_curve(yt, p, n_bins=8, strategy="quantile")
        brier = brier_score_loss(yt, p)
        ax.plot([0, 1], [0, 1], ":", c=META_GREY, lw=0.8)
        ax.plot(mean_pred, frac, "o-", color=BLUE, ms=4)
        ax.set_title(f"{d[:10]}  (Brier={brier:.2f})")
        ax.set_xlim(0, 1); ax.set_ylim(0, 1); ax.set_xlabel("Predicted P(fail)")
    np.atleast_1d(axes)[0].set_ylabel("Observed fraction resistant")
    fig.suptitle("Calibration on held-out grouped test set", y=1.02)
    fig.tight_layout()
    fig.savefig(HERE / "eval/fig_reliability.png", dpi=200)
    print("fig_reliability.png  saved")


def main():
    feats, labels, bundle = load()
    fig_leakage(feats, labels)
    fig_reliability(feats, labels, bundle)


if __name__ == "__main__":
    main()
