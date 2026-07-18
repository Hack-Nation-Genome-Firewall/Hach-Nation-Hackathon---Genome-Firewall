"""Validated single-genome inference for calibrated Track B model bundles."""

from __future__ import annotations

import argparse
import json
import math
from pathlib import Path
from typing import Any, Mapping

import joblib
import numpy as np
import pandas as pd

from module2_predictor.contracts import canonical_json_sha256, load_feature_spec, validate_inference_row


ROOT = Path(__file__).resolve().parents[1]
SYNTHETIC = ROOT / "data/synthetic"
DEFAULT_BUNDLE = ROOT / "models/synthetic_bundle.joblib"
DEFAULT_SPEC = SYNTHETIC / "feature_spec.json"


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--bundle", type=Path, default=DEFAULT_BUNDLE)
    parser.add_argument("--spec", type=Path, default=DEFAULT_SPEC)
    parser.add_argument("--features", type=Path, default=SYNTHETIC / "features.csv")
    parser.add_argument("--splits", type=Path, default=SYNTHETIC / "split_manifest.csv")
    parser.add_argument("--genome-id")
    return parser.parse_args()


def load_bundle(path: Path | str = DEFAULT_BUNDLE) -> dict[str, Any]:
    bundle = joblib.load(path)
    required = {
        "bundle_schema_version",
        "model_version",
        "drugs",
        "model_features",
        "feature_spec",
        "feature_spec_sha256",
        "drug_models",
    }
    missing = sorted(required - bundle.keys())
    if missing:
        raise ValueError(f"Model bundle is missing keys: {missing}")
    actual_sha = canonical_json_sha256(bundle["feature_spec"])
    if actual_sha != bundle["feature_spec_sha256"]:
        raise ValueError("Model bundle feature specification checksum is invalid")
    return bundle


def _missing(value: Any) -> bool:
    return value is None or value is pd.NA or (isinstance(value, float) and math.isnan(value))


def _target_gate(row: Mapping[str, Any], spec: Mapping[str, Any], drug: str) -> dict[str, Any]:
    target_features = list(spec["drug_targets"][drug])
    missing = [feature for feature in target_features if feature not in row or _missing(row.get(feature))]
    invalid = [
        feature
        for feature in target_features
        if feature in row and not _missing(row.get(feature)) and row.get(feature) not in (0, 1, False, True)
    ]
    if missing or invalid:
        return {
            "target_features": target_features,
            "status": "unknown",
            "present": None,
            "action": "route_to_no_call",
            "details": {"missing": missing, "invalid": invalid},
        }
    absent = [feature for feature in target_features if int(row[feature]) == 0]
    if absent:
        return {
            "target_features": target_features,
            "status": "absent_or_disrupted",
            "present": False,
            "action": "route_to_no_call",
            "details": {"absent_or_disrupted": absent},
        }
    return {
        "target_features": target_features,
        "status": "present",
        "present": True,
        "action": "proceed",
        "details": {},
    }


def _quality_gate(row: Mapping[str, Any], spec: Mapping[str, Any]) -> dict[str, Any]:
    features = spec["quality_features"]
    policy = spec["quality_policy"]
    missing = [column for column in features.values() if column not in row or _missing(row.get(column))]
    if missing:
        return {"status": "unknown", "action": "route_to_no_call", "reasons": [f"missing:{','.join(missing)}"]}
    reasons: list[str] = []
    completeness = float(row[features["completeness"]])
    contamination = float(row[features["contamination"]])
    contigs = int(row[features["contigs"]])
    if completeness < float(policy["minimum_completeness"]):
        reasons.append("completeness_below_minimum")
    if contamination > float(policy["maximum_contamination"]):
        reasons.append("contamination_above_maximum")
    if contigs > int(policy["maximum_contigs"]):
        reasons.append("contigs_above_maximum")
    return {
        "status": "fail" if reasons else "pass",
        "action": "route_to_no_call" if reasons else "proceed",
        "reasons": reasons,
        "observed": {
            "completeness": completeness,
            "contamination": contamination,
            "contigs": contigs,
        },
    }


def _minimum_hamming_distance(reference: np.ndarray, row: np.ndarray) -> float:
    matrix = np.asarray(reference, dtype=np.uint8)
    vector = np.asarray(row, dtype=np.uint8).reshape(-1)
    if matrix.ndim != 2 or matrix.shape[1] != len(vector) or not len(matrix):
        raise ValueError("Invalid OOD reference matrix in model bundle")
    return float(np.mean(matrix != vector, axis=1).min())


def _evidence(
    row: Mapping[str, Any],
    spec: Mapping[str, Any],
    drug: str,
    coefficients: np.ndarray,
) -> tuple[str, list[dict[str, Any]], list[dict[str, Any]]]:
    known: list[dict[str, Any]] = []
    known_names: set[str] = set()
    for marker, metadata in spec["marker_evidence"].items():
        if drug in metadata.get("drugs", []) and int(row.get(marker, 0)) == 1:
            known_names.add(marker)
            known.append(
                {
                    "marker": marker,
                    "type": metadata.get("type", "unknown"),
                    "source": metadata.get("source"),
                    "amrfinder_method": metadata.get("amrfinder_method"),
                }
            )
    statistical: list[dict[str, Any]] = []
    for feature, coefficient in zip(spec["model_features"], coefficients, strict=True):
        if int(row[feature]) == 1 and feature not in known_names and abs(float(coefficient)) > 1e-9:
            statistical.append({"feature": feature, "coefficient": round(float(coefficient), 6)})
    statistical.sort(key=lambda item: abs(item["coefficient"]), reverse=True)
    statistical = statistical[:5]
    if known:
        tier = "known_marker"
    elif statistical:
        tier = "statistical_only"
    else:
        tier = "no_signal"
    return tier, known, statistical


def predict_drug(
    row: Mapping[str, Any],
    bundle: Mapping[str, Any],
    drug: str,
    spec: Mapping[str, Any] | None = None,
) -> dict[str, Any]:
    spec = spec or bundle["feature_spec"]
    validate_inference_row(row, spec)
    if canonical_json_sha256(spec) != bundle["feature_spec_sha256"]:
        raise ValueError("Inference feature specification does not match the trained bundle")
    if drug not in bundle["drug_models"]:
        raise ValueError(f"Unsupported antibiotic: {drug}")
    model = bundle["drug_models"][drug]
    vector = np.asarray([[int(row[feature]) for feature in bundle["model_features"]]], dtype=np.uint8)
    raw_probability = model["classifier"].predict_proba(vector)[:, 1]
    p_fail = float(model["calibrator"].predict(raw_probability)[0])
    verdict_probability = max(p_fail, 1.0 - p_fail)
    statistical_verdict = "likely_to_fail" if p_fail >= 0.5 else "likely_to_work"

    target_gate = _target_gate(row, spec, drug)
    quality_gate = _quality_gate(row, spec)
    ood_distance = _minimum_hamming_distance(model["ood_reference"], vector[0])
    ood = ood_distance > float(model["ood_threshold"]) + 1e-12
    tier, known_markers, statistical_features = _evidence(
        row,
        spec,
        drug,
        model["classifier"].coef_[0],
    )
    known_marker_conflict = bool(known_markers) and statistical_verdict == "likely_to_work"

    no_call_reasons: list[str] = []
    if target_gate["status"] == "unknown":
        no_call_reasons.append("target_status_unknown")
    elif target_gate["status"] != "present":
        no_call_reasons.append("drug_target_absent_or_disrupted")
    if quality_gate["status"] == "unknown":
        no_call_reasons.append("quality_status_unknown")
    elif quality_gate["status"] == "fail":
        no_call_reasons.append("low_assembly_quality")
    if ood:
        no_call_reasons.append("out_of_distribution")
    if known_marker_conflict:
        no_call_reasons.append("known_marker_conflicts_with_model")
    if verdict_probability < float(model["call_threshold"]):
        no_call_reasons.append("low_confidence")
    verdict = "no_call" if no_call_reasons else statistical_verdict

    return {
        "record_schema_version": 1,
        "drug": drug,
        "verdict": verdict,
        "p_fail": round(p_fail, 6),
        "verdict_probability": round(verdict_probability, 6),
        "confidence": round(verdict_probability, 6),
        "calibrated": True,
        "evidence_tier": tier,
        "supporting_markers": known_markers,
        "statistical_features": statistical_features,
        "target_gate": target_gate,
        "quality_gate": quality_gate,
        "ood": {
            "distance": round(ood_distance, 6),
            "threshold": round(float(model["ood_threshold"]), 6),
            "is_out_of_distribution": ood,
        },
        "call_threshold": round(float(model["call_threshold"]), 6),
        "no_call_reason": no_call_reasons[0] if no_call_reasons else None,
        "no_call_reasons": no_call_reasons,
        "model_version": bundle["model_version"],
        "feature_spec_sha256": bundle["feature_spec_sha256"],
    }


def predict_genome(
    row: Mapping[str, Any],
    bundle: Mapping[str, Any] | None = None,
    spec: Mapping[str, Any] | None = None,
) -> list[dict[str, Any]]:
    bundle = bundle or load_bundle()
    spec = spec or bundle["feature_spec"]
    return [predict_drug(row, bundle, drug, spec) for drug in bundle["drugs"]]


def main() -> None:
    args = parse_args()
    bundle = load_bundle(args.bundle)
    spec = load_feature_spec(args.spec)
    features = pd.read_csv(args.features, dtype={"genome_id": str})
    splits = pd.read_csv(args.splits, dtype={"genome_id": str})
    held_out = features.merge(splits, on="genome_id", validate="one_to_one")
    held_out = held_out.loc[held_out["split"] == "test"]
    if held_out.empty:
        raise ValueError("No test genomes are available for the demo")
    if args.genome_id:
        held_out = held_out.loc[held_out["genome_id"] == args.genome_id]
        if held_out.empty:
            raise ValueError(f"Genome {args.genome_id} is not in the test split")
    row = held_out.iloc[0].to_dict()
    print(json.dumps({"genome_id": row["genome_id"], "predictions": predict_genome(row, bundle, spec)}, indent=2))


if __name__ == "__main__":
    main()
