"""Evaluate calibrated predictions and abstention on a grouped held-out split."""

from __future__ import annotations

import argparse
import json
from collections import Counter
from pathlib import Path
from typing import Any

import numpy as np
import pandas as pd
from sklearn.metrics import (
    accuracy_score,
    average_precision_score,
    balanced_accuracy_score,
    brier_score_loss,
    f1_score,
    recall_score,
    roc_auc_score,
)

from data.bvbrc import write_csv_atomic, write_json_atomic
from module2_predictor.contracts import load_feature_spec, load_training_frames
from module2_predictor.predict import load_bundle, predict_drug


ROOT = Path(__file__).resolve().parents[1]
SYNTHETIC = ROOT / "data/synthetic"


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--features", type=Path, default=SYNTHETIC / "features.csv")
    parser.add_argument("--labels", type=Path, default=SYNTHETIC / "labels.csv")
    parser.add_argument("--splits", type=Path, default=SYNTHETIC / "split_manifest.csv")
    parser.add_argument("--spec", type=Path, default=SYNTHETIC / "feature_spec.json")
    parser.add_argument("--bundle", type=Path, default=ROOT / "models/synthetic_bundle.joblib")
    parser.add_argument("--output-dir", type=Path, default=ROOT / "eval")
    return parser.parse_args()


def _classification_metrics(y: np.ndarray, probabilities: np.ndarray) -> dict[str, float | None]:
    predictions = (probabilities >= 0.5).astype(int)
    result: dict[str, float | None] = {
        "balanced_accuracy": None,
        "recall_resistant": None,
        "recall_susceptible": None,
        "f1": None,
        "auroc": None,
        "pr_auc": None,
        "brier": None,
    }
    if not len(y):
        return result
    result["brier"] = float(brier_score_loss(y, probabilities))
    result["f1"] = float(f1_score(y, predictions, zero_division=0))
    if len(np.unique(y)) == 2:
        result.update(
            {
                "balanced_accuracy": float(balanced_accuracy_score(y, predictions)),
                "recall_resistant": float(recall_score(y, predictions, pos_label=1)),
                "recall_susceptible": float(recall_score(y, predictions, pos_label=0)),
                "auroc": float(roc_auc_score(y, probabilities)),
                "pr_auc": float(average_precision_score(y, probabilities)),
            }
        )
    return result


def evaluate_bundle(
    features: pd.DataFrame,
    labels: pd.DataFrame,
    splits: pd.DataFrame,
    spec: dict[str, Any],
    bundle: dict[str, Any],
) -> tuple[pd.DataFrame, pd.DataFrame, pd.DataFrame, dict[str, Any]]:
    held_out_features = features.merge(splits, on="genome_id", validate="one_to_one")
    held_out_features = held_out_features.loc[held_out_features["split"] == "test"]
    overall_rows: list[dict[str, Any]] = []
    group_rows: list[dict[str, Any]] = []
    prediction_rows: list[dict[str, Any]] = []

    for drug in spec["drugs"]:
        drug_labels = labels.loc[labels["antibiotic"] == drug, ["genome_id", "phenotype"]]
        data = held_out_features.merge(drug_labels, on="genome_id", how="inner", validate="one_to_one")
        if data.empty:
            raise ValueError(f"No held-out labels are available for {drug}")
        records = []
        for row in data.to_dict(orient="records"):
            prediction = predict_drug(row, bundle, drug, spec)
            y_true = 1 if row["phenotype"] == "Resistant" else 0
            records.append(
                {
                    "genome_id": row["genome_id"],
                    "cluster_id": row["cluster_id"],
                    "drug": drug,
                    "y_true": y_true,
                    "phenotype": row["phenotype"],
                    "verdict": prediction["verdict"],
                    "p_fail": prediction["p_fail"],
                    "verdict_probability": prediction["verdict_probability"],
                    "evidence_tier": prediction["evidence_tier"],
                    "no_call_reason": prediction["no_call_reason"],
                    "no_call_reasons": ";".join(prediction["no_call_reasons"]),
                    "ood_distance": prediction["ood"]["distance"],
                }
            )
        predictions = pd.DataFrame(records)
        prediction_rows.extend(records)
        y = predictions["y_true"].to_numpy(dtype=int)
        probabilities = predictions["p_fail"].to_numpy(dtype=float)
        called = predictions["verdict"].ne("no_call").to_numpy()
        all_metrics = _classification_metrics(y, probabilities)
        called_metrics = _classification_metrics(y[called], probabilities[called])
        called_accuracy = (
            float(accuracy_score(y[called], (probabilities[called] >= 0.5).astype(int))) if called.any() else None
        )
        reason_counts = Counter(
            reason
            for value in predictions["no_call_reasons"]
            for reason in str(value).split(";")
            if reason
        )
        overall_rows.append(
            {
                "drug": drug,
                "n": int(len(predictions)),
                "resistant": int(y.sum()),
                "susceptible": int(len(y) - y.sum()),
                **all_metrics,
                "no_call_rate": float(1.0 - called.mean()),
                "called_coverage": float(called.mean()),
                "called_n": int(called.sum()),
                "called_accuracy": called_accuracy,
                "called_balanced_accuracy": called_metrics["balanced_accuracy"],
                "no_call_reasons": json.dumps(reason_counts, sort_keys=True),
            }
        )

        for cluster_id, group in predictions.groupby("cluster_id", sort=True):
            group_y = group["y_true"].to_numpy(dtype=int)
            group_probabilities = group["p_fail"].to_numpy(dtype=float)
            group_called = group["verdict"].ne("no_call").to_numpy()
            metrics = _classification_metrics(group_y, group_probabilities)
            group_rows.append(
                {
                    "drug": drug,
                    "cluster_id": cluster_id,
                    "n": int(len(group)),
                    "resistant": int(group_y.sum()),
                    "susceptible": int(len(group_y) - group_y.sum()),
                    **metrics,
                    "called_coverage": float(group_called.mean()),
                    "called_n": int(group_called.sum()),
                }
            )

    overall = pd.DataFrame(overall_rows)
    per_group = pd.DataFrame(group_rows)
    predictions = pd.DataFrame(prediction_rows)
    report = {
        "schema_version": 1,
        "model_version": bundle["model_version"],
        "feature_spec_sha256": bundle["feature_spec_sha256"],
        "split": "test",
        "grouping": "homology cluster",
        "notes": [
            "All probability metrics use every held-out example.",
            "Called metrics exclude explicit no-call predictions.",
            "Per-group class-dependent metrics are null when a group has only one class.",
            "Synthetic results are integration tests, not scientific performance claims."
            if spec.get("synthetic")
            else "Results remain research-only and require laboratory confirmation.",
        ],
        "overall": overall.to_dict(orient="records"),
    }
    return overall, per_group, predictions, report


def main() -> None:
    args = parse_args()
    spec = load_feature_spec(args.spec)
    features, labels, splits = load_training_frames(args.features, args.labels, args.splits, spec)
    bundle = load_bundle(args.bundle)
    overall, per_group, predictions, report = evaluate_bundle(features, labels, splits, spec, bundle)
    args.output_dir.mkdir(parents=True, exist_ok=True)
    write_csv_atomic(args.output_dir / "overall_metrics.csv", overall)
    write_csv_atomic(args.output_dir / "per_group_metrics.csv", per_group)
    write_csv_atomic(args.output_dir / "held_out_predictions.csv", predictions)
    write_json_atomic(args.output_dir / "evaluation_report.json", report)
    print(overall.to_string(index=False))
    print(f"wrote evaluation artifacts to {args.output_dir}")


if __name__ == "__main__":
    main()
