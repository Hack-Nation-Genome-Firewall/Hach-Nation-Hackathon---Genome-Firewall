"""Train one calibrated, abstaining logistic-regression model per antibiotic."""

from __future__ import annotations

import argparse
import json
import os
import platform
from datetime import datetime, timezone
from pathlib import Path
from typing import Any

import joblib
import numpy as np
import pandas as pd
import sklearn
from sklearn.linear_model import LogisticRegression
from sklearn.metrics import balanced_accuracy_score, brier_score_loss

from module2_predictor.calibration import choose_call_threshold, fit_calibrator
from module2_predictor.contracts import (
    canonical_json_sha256,
    file_sha256,
    load_feature_spec,
    load_training_frames,
)


ROOT = Path(__file__).resolve().parents[1]
SYNTHETIC = ROOT / "data/synthetic"


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--features", type=Path, default=SYNTHETIC / "features.csv")
    parser.add_argument("--labels", type=Path, default=SYNTHETIC / "labels.csv")
    parser.add_argument("--splits", type=Path, default=SYNTHETIC / "split_manifest.csv")
    parser.add_argument("--spec", type=Path, default=SYNTHETIC / "feature_spec.json")
    parser.add_argument("--output", type=Path, default=ROOT / "models/synthetic_bundle.joblib")
    parser.add_argument("--model-version", default="track-b-v1")
    parser.add_argument("--calibration-method", choices=("auto", "isotonic", "sigmoid"), default="auto")
    parser.add_argument("--target-called-balanced-accuracy", type=float, default=0.80)
    parser.add_argument("--minimum-called-coverage", type=float, default=0.50)
    parser.add_argument("--ood-quantile", type=float, default=0.99)
    parser.add_argument("--random-state", type=int, default=20260718)
    return parser.parse_args()


def _minimum_hamming_distances(reference: np.ndarray, query: np.ndarray) -> np.ndarray:
    reference = np.asarray(reference, dtype=np.uint8)
    query = np.asarray(query, dtype=np.uint8)
    if reference.ndim != 2 or query.ndim != 2 or reference.shape[1] != query.shape[1]:
        raise ValueError("Reference and query feature matrices must be aligned two-dimensional arrays")
    if not len(reference):
        raise ValueError("OOD reference matrix cannot be empty")
    distances = []
    for row in query:
        distances.append(float(np.mean(reference != row, axis=1).min()))
    return np.asarray(distances, dtype=float)


def train_bundle(
    features: pd.DataFrame,
    labels: pd.DataFrame,
    splits: pd.DataFrame,
    spec: dict[str, Any],
    *,
    model_version: str,
    calibration_method: str = "auto",
    target_called_balanced_accuracy: float = 0.80,
    minimum_called_coverage: float = 0.50,
    ood_quantile: float = 0.99,
    random_state: int = 20260718,
    source_hashes: dict[str, str] | None = None,
) -> tuple[dict[str, Any], dict[str, Any]]:
    if not 0.5 <= target_called_balanced_accuracy <= 1.0:
        raise ValueError("target_called_balanced_accuracy must be between 0.5 and 1")
    if not 0.0 < minimum_called_coverage <= 1.0:
        raise ValueError("minimum_called_coverage must be in (0, 1]")
    if not 0.5 <= ood_quantile <= 1.0:
        raise ValueError("ood_quantile must be between 0.5 and 1")

    merged_features = features.merge(splits, on="genome_id", how="inner", validate="one_to_one")
    model_features = list(spec["model_features"])
    drug_models: dict[str, Any] = {}
    summaries: dict[str, Any] = {}

    for drug in spec["drugs"]:
        drug_labels = labels.loc[labels["antibiotic"] == drug, ["genome_id", "phenotype"]].copy()
        if drug_labels.empty:
            raise ValueError(f"No labels are available for {drug}")
        data = merged_features.merge(drug_labels, on="genome_id", how="inner", validate="one_to_one")
        data["label"] = data["phenotype"].map({"Susceptible": 0, "Resistant": 1}).astype(int)
        masks = {split: data["split"].eq(split).to_numpy() for split in ("train", "calibration", "test")}
        for split, mask in masks.items():
            _require_binary_split(data.loc[mask, "label"].to_numpy(), drug, split)

        x = data[model_features].to_numpy(dtype=np.uint8)
        y = data["label"].to_numpy(dtype=int)
        train_mask = masks["train"]
        calibration_mask = masks["calibration"]
        classifier = LogisticRegression(
            C=1.0,
            class_weight="balanced",
            max_iter=3000,
            solver="liblinear",
            random_state=random_state,
        )
        classifier.fit(x[train_mask], y[train_mask])
        raw_calibration = classifier.predict_proba(x[calibration_mask])[:, 1]
        calibrator = fit_calibrator(raw_calibration, y[calibration_mask], calibration_method)
        calibrated = calibrator.predict(raw_calibration)
        threshold = choose_call_threshold(
            calibrated,
            y[calibration_mask],
            target_balanced_accuracy=target_called_balanced_accuracy,
            minimum_coverage=minimum_called_coverage,
        )

        reference = x[train_mask].astype(np.uint8)
        calibration_distances = _minimum_hamming_distances(reference, x[calibration_mask])
        distance_floor = 1.0 / max(1, len(model_features))
        ood_threshold = max(float(np.quantile(calibration_distances, ood_quantile)), distance_floor)
        ood_threshold = min(1.0, ood_threshold)

        split_counts = {
            split: {
                "n": int(mask.sum()),
                "resistant": int(y[mask].sum()),
                "susceptible": int(mask.sum() - y[mask].sum()),
            }
            for split, mask in masks.items()
        }
        drug_models[drug] = {
            "classifier": classifier,
            "calibrator": calibrator,
            "call_threshold": float(threshold["threshold"]),
            "ood_reference": reference,
            "ood_threshold": ood_threshold,
            "split_counts": split_counts,
        }
        summaries[drug] = {
            "calibration_method": calibrator.method,
            "calibration_brier": float(brier_score_loss(y[calibration_mask], calibrated)),
            "calibration_balanced_accuracy_at_0_5": float(
                balanced_accuracy_score(y[calibration_mask], calibrated >= 0.5)
            ),
            "call_threshold": threshold,
            "ood_quantile": ood_quantile,
            "ood_threshold": ood_threshold,
            "split_counts": split_counts,
        }

    bundle = {
        "bundle_schema_version": 1,
        "model_version": model_version,
        "created_at": datetime.now(timezone.utc).replace(microsecond=0).isoformat(),
        "species": spec["species"],
        "drugs": list(spec["drugs"]),
        "model_features": model_features,
        "feature_spec": spec,
        "feature_spec_sha256": canonical_json_sha256(spec),
        "source_hashes": source_hashes or {},
        "runtime": {
            "python": platform.python_version(),
            "numpy": np.__version__,
            "pandas": pd.__version__,
            "scikit_learn": sklearn.__version__,
        },
        "drug_models": drug_models,
    }
    summary = {
        "bundle_schema_version": 1,
        "model_version": model_version,
        "feature_spec_sha256": bundle["feature_spec_sha256"],
        "source_hashes": source_hashes or {},
        "training": summaries,
    }
    return bundle, summary


def _require_binary_split(labels: np.ndarray, drug: str, split: str) -> None:
    counts = {int(value): int((labels == value).sum()) for value in np.unique(labels)}
    if set(counts) != {0, 1}:
        raise ValueError(f"{drug} {split} split must contain both classes, got {counts}")
    if min(counts.values()) < 5:
        raise ValueError(f"{drug} {split} split needs at least five examples per class, got {counts}")


def save_bundle(path: Path, bundle: dict[str, Any], summary: dict[str, Any]) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    temporary = path.with_name(f".{path.name}.tmp")
    joblib.dump(bundle, temporary, compress=3)
    os.replace(temporary, path)
    summary_path = path.with_suffix(".training.json")
    summary_path.write_text(json.dumps(summary, indent=2, sort_keys=True) + "\n", encoding="utf-8")


def main() -> None:
    args = parse_args()
    spec = load_feature_spec(args.spec)
    features, labels, splits = load_training_frames(args.features, args.labels, args.splits, spec)
    source_hashes = {
        "features": file_sha256(args.features),
        "labels": file_sha256(args.labels),
        "splits": file_sha256(args.splits),
        "feature_spec_file": file_sha256(args.spec),
    }
    bundle, summary = train_bundle(
        features,
        labels,
        splits,
        spec,
        model_version=args.model_version,
        calibration_method=args.calibration_method,
        target_called_balanced_accuracy=args.target_called_balanced_accuracy,
        minimum_called_coverage=args.minimum_called_coverage,
        ood_quantile=args.ood_quantile,
        random_state=args.random_state,
        source_hashes=source_hashes,
    )
    save_bundle(args.output, bundle, summary)
    for drug, details in summary["training"].items():
        counts = details["split_counts"]
        print(
            f"trained {drug}: train={counts['train']['n']} calibration={counts['calibration']['n']} "
            f"test={counts['test']['n']} calibrator={details['calibration_method']} "
            f"call_threshold={details['call_threshold']['threshold']:.2f}"
        )
    print(f"saved {args.output}")


if __name__ == "__main__":
    main()
