"""
GENOME FIREWALL — reproduce the evaluation figures (Track C).

Produces:
  eval/fig_leakage.png     random-split vs homology-grouped-split balanced accuracy
  eval/fig_reliability.png per-drug calibration reliability curves (held-out test)

Reads the frozen data/synthetic/* fixture and models/synthetic_bundle.joblib.
Calibrated P(fail) is recomputed straight from the bundle (classifier ->
calibrator), matching module2_predictor.predict, so the figures reflect the
deployed model without duplicating gate logic.

Run after training the bundle:  python TrackC/make_figures.py
"""
import json
import sys
from pathlib import Path

import joblib
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
from sklearn.linear_model import LogisticRegression
from sklearn.metrics import balanced_accuracy_score, brier_score_loss
from sklearn.calibration import calibration_curve

HERE = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(HERE))  # bundle unpickles a calibrator from module2_predictor
SYNTHETIC = HERE / "data/synthetic"
SPEC = json.load(open(SYNTHETIC / "feature_spec.json"))
DRUGS = SPEC["drugs"]
MODEL_FEATURES = SPEC["model_features"]

META_GREY = "#8a8a8a"; RED = "#c44e52"; BLUE = "#4c72b0"


def load():
    feats = pd.read_csv(SYNTHETIC / "features.csv", dtype={"genome_id": str})
    labels = pd.read_csv(SYNTHETIC / "labels.csv", dtype={"genome_id": str})
    splits = pd.read_csv(SYNTHETIC / "split_manifest.csv",
                         dtype={"genome_id": str, "cluster_id": str})
    bundle = joblib.load(HERE / "models/synthetic_bundle.joblib")
    feats = feats.merge(splits, on="genome_id", validate="one_to_one")
    return feats, labels, bundle


def y_for(feats, labels, drug):
    """Lab/fixture phenotype aligned to feats rows; NaN where untested."""
    sub = labels[labels.antibiotic == drug][["genome_id", "phenotype"]]
    m = feats[["genome_id"]].merge(sub, on="genome_id", how="left")
    return m.phenotype.map({"Susceptible": 0, "Resistant": 1}).values


def calibrated_p_fail(bundle, drug, X):
    """classifier -> calibrator, identical to module2_predictor.predict.p_fail."""
    model = bundle["drug_models"][drug]
    raw = model["classifier"].predict_proba(X)[:, 1]
    return model["calibrator"].predict(raw)


def fig_leakage(feats, labels):
    """Random row split (leaks) vs frozen homology-grouped split (honest)."""
    X = feats[MODEL_FEATURES].values.astype(np.uint8)
    split = feats.split.values
    train_cal = np.isin(split, ["train", "calibration"])
    test = split == "test"
    rand_res, grp_res = {}, {}
    for d in DRUGS:
        y = y_for(feats, labels, d)
        tested = ~np.isnan(y)
        idx = np.where(tested)[0]
        yv = y.astype(float)

        def fit_eval(tr, te):
            m = LogisticRegression(C=1.0, class_weight="balanced", max_iter=2000)
            m.fit(X[tr], yv[tr].astype(int))
            pred = (m.predict_proba(X[te])[:, 1] >= 0.5).astype(int)
            return balanced_accuracy_score(yv[te].astype(int), pred)

        rr = np.random.default_rng(0)
        perm = rr.permutation(idx)
        cut = int(0.75 * len(perm))
        rand_res[d] = fit_eval(perm[:cut], perm[cut:])
        tr = idx[train_cal[idx]]; te = idx[test[idx]]
        grp_res[d] = fit_eval(tr, te)

    deltas = {d: (rand_res[d] - grp_res[d]) * 100 for d in DRUGS}
    n_inflated = sum(v > 0 for v in deltas.values())
    mean_delta = float(np.mean(list(deltas.values())))

    fig, ax = plt.subplots(figsize=(6.0, 3.4))
    x = np.arange(len(DRUGS)); w = 0.36
    ax.bar(x - w/2, [rand_res[d] for d in DRUGS], w, label="Random split (leaks)", color=RED)
    ax.bar(x + w/2, [grp_res[d] for d in DRUGS], w, label="Grouped split (honest)", color=BLUE)
    for i, d in enumerate(DRUGS):
        dv = deltas[d]
        lbl = f"+{dv:.0f}pt" if dv >= 0 else f"−{abs(dv):.0f}pt"
        ax.annotate(lbl, (i, max(rand_res[d], grp_res[d]) + 0.02), ha="center", fontsize=6, color=RED)
    ax.set_xticks(x); ax.set_xticklabels(DRUGS, fontsize=6.5, rotation=15); ax.set_ylim(0, 1.05)
    ax.set_ylabel("Balanced accuracy")
    ax.set_title(f"Random split vs. homology-grouped split\n"
                 f"random inflates on {n_inflated} of {len(DRUGS)} drugs; mean {mean_delta:+.1f} pt",
                 fontsize=9)
    ax.axhline(0.5, ls=":", c=META_GREY, lw=0.8)
    ax.legend(frameon=False, fontsize=6, loc="lower right")
    fig.tight_layout()
    (HERE / "eval").mkdir(exist_ok=True)
    fig.savefig(HERE / "eval/fig_leakage.png", dpi=200, bbox_inches="tight", pad_inches=0.25)
    plt.close(fig)
    print(f"fig_leakage.png  deltas={ {d: round(v,1) for d,v in deltas.items()} }")


def fig_reliability(feats, labels, bundle):
    X_all = feats[MODEL_FEATURES].values.astype(np.uint8)
    test = feats.split.values == "test"
    fig, axes = plt.subplots(1, len(DRUGS), figsize=(2.8 * len(DRUGS), 3.0))
    axes = np.atleast_1d(axes)
    for ax, d in zip(axes, DRUGS):
        y = y_for(feats, labels, d)
        te = test & (~np.isnan(y))
        p = calibrated_p_fail(bundle, d, X_all[te])
        yt = y[te].astype(int)
        brier = brier_score_loss(yt, p)
        ax.plot([0, 1], [0, 1], ":", c=META_GREY, lw=0.8)
        if len(np.unique(yt)) == 2:
            n_bins = min(8, max(2, len(yt) // 4))
            frac, mean_pred = calibration_curve(yt, p, n_bins=n_bins, strategy="quantile")
            ax.plot(mean_pred, frac, "o-", color=BLUE, ms=4)
        else:
            ax.text(0.5, 0.5, "single-class\nheld-out set", ha="center", va="center",
                    fontsize=7, color=META_GREY, transform=ax.transAxes)
        ax.set_title(f"{d[:12]}  (Brier={brier:.2f})", fontsize=8)
        ax.set_xlim(0, 1); ax.set_ylim(0, 1); ax.set_xlabel("Predicted P(fail)", fontsize=7)
    axes[0].set_ylabel("Observed fraction resistant", fontsize=7)
    fig.suptitle("Calibration on held-out grouped test set", y=0.99, fontsize=10)
    fig.tight_layout(rect=(0, 0, 1, 0.93))
    fig.savefig(HERE / "eval/fig_reliability.png", dpi=200, bbox_inches="tight", pad_inches=0.25)
    plt.close(fig)
    print("fig_reliability.png  saved")


def main():
    feats, labels, bundle = load()
    fig_leakage(feats, labels)
    fig_reliability(feats, labels, bundle)


if __name__ == "__main__":
    main()
