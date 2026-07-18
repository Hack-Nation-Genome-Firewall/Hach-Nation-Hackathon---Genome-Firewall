"""Validation and I/O contracts shared by Track B training and inference."""

from __future__ import annotations

import hashlib
import json
import math
from pathlib import Path
from typing import Any, Mapping

import pandas as pd


VALID_SPLITS = {"train", "calibration", "test"}
VALID_PHENOTYPES = {"Susceptible", "Resistant"}


def canonical_json_sha256(value: Any) -> str:
    encoded = json.dumps(value, sort_keys=True, separators=(",", ":")).encode("utf-8")
    return hashlib.sha256(encoded).hexdigest()


def file_sha256(path: Path) -> str:
    digest = hashlib.sha256()
    with path.open("rb") as handle:
        for block in iter(lambda: handle.read(1024 * 1024), b""):
            digest.update(block)
    return digest.hexdigest()


def load_feature_spec(path: Path) -> dict[str, Any]:
    with path.open("r", encoding="utf-8") as handle:
        spec = json.load(handle)
    validate_feature_spec(spec)
    return spec


def validate_feature_spec(spec: Mapping[str, Any]) -> None:
    required = {
        "schema_version",
        "species",
        "drugs",
        "model_features",
        "marker_evidence",
        "drug_targets",
        "quality_features",
        "quality_policy",
        "expected_label_evidence",
    }
    missing = sorted(required - spec.keys())
    if missing:
        raise ValueError(f"Feature specification is missing keys: {missing}")
    drugs = list(spec["drugs"])
    if not drugs or len(drugs) != len(set(drugs)):
        raise ValueError("Feature specification drugs must be non-empty and unique")
    features = list(spec["model_features"])
    if not features or len(features) != len(set(features)):
        raise ValueError("model_features must be non-empty and unique")
    targets = spec["drug_targets"]
    for drug in drugs:
        if drug not in targets or not targets[drug]:
            raise ValueError(f"No target features are configured for {drug}")
    evidence = spec["marker_evidence"]
    unknown_markers = sorted(set(evidence) - set(features))
    if unknown_markers:
        raise ValueError(f"marker_evidence references unknown model features: {unknown_markers}")
    for marker, record in evidence.items():
        marker_drugs = set(record.get("drugs", []))
        if not marker_drugs <= set(drugs):
            raise ValueError(f"Marker {marker} references unsupported drugs: {sorted(marker_drugs - set(drugs))}")


def load_training_frames(
    features_path: Path,
    labels_path: Path,
    splits_path: Path,
    spec: Mapping[str, Any],
) -> tuple[pd.DataFrame, pd.DataFrame, pd.DataFrame]:
    features = pd.read_csv(features_path, dtype={"genome_id": str})
    labels = pd.read_csv(labels_path, dtype={"genome_id": str})
    splits = pd.read_csv(splits_path, dtype={"genome_id": str, "cluster_id": str})
    validate_training_frames(features, labels, splits, spec)
    return features, labels, splits


def validate_training_frames(
    features: pd.DataFrame,
    labels: pd.DataFrame,
    splits: pd.DataFrame,
    spec: Mapping[str, Any],
) -> None:
    _require_columns(features, {"genome_id", *spec["model_features"]}, "features")
    _require_columns(labels, {"genome_id", "antibiotic", "phenotype", "evidence"}, "labels")
    _require_columns(splits, {"genome_id", "cluster_id", "split"}, "split manifest")
    if features["genome_id"].duplicated().any():
        raise ValueError("Feature matrix contains duplicate genome_id rows")
    if splits["genome_id"].duplicated().any():
        raise ValueError("Split manifest contains duplicate genome_id rows")
    if labels.duplicated(["genome_id", "antibiotic"]).any():
        raise ValueError("Labels contain duplicate genome_id/antibiotic pairs")
    split_values = set(splits["split"].dropna().astype(str))
    if split_values != VALID_SPLITS:
        raise ValueError(f"Split manifest must contain exactly {sorted(VALID_SPLITS)}, got {sorted(split_values)}")
    cluster_split_counts = splits.groupby("cluster_id")["split"].nunique()
    leaking = cluster_split_counts[cluster_split_counts > 1]
    if not leaking.empty:
        raise ValueError(f"Homology clusters cross splits: {leaking.index.tolist()[:10]}")
    if set(features["genome_id"]) != set(splits["genome_id"]):
        missing_splits = sorted(set(features["genome_id"]) - set(splits["genome_id"]))[:10]
        missing_features = sorted(set(splits["genome_id"]) - set(features["genome_id"]))[:10]
        raise ValueError(
            f"Feature/split genome IDs differ; missing splits={missing_splits}, missing features={missing_features}"
        )
    unsupported_drugs = sorted(set(labels["antibiotic"]) - set(spec["drugs"]))
    if unsupported_drugs:
        raise ValueError(f"Labels contain unsupported antibiotics: {unsupported_drugs}")
    invalid_phenotypes = sorted(set(labels["phenotype"]) - VALID_PHENOTYPES)
    if invalid_phenotypes:
        raise ValueError(f"Labels contain unsupported phenotypes: {invalid_phenotypes}")
    expected_evidence = spec["expected_label_evidence"]
    invalid_evidence = labels.loc[labels["evidence"] != expected_evidence, "evidence"].drop_duplicates().tolist()
    if invalid_evidence:
        raise ValueError(
            f"Labels violate expected evidence {expected_evidence!r}: {invalid_evidence[:10]}"
        )
    unknown_genomes = sorted(set(labels["genome_id"]) - set(features["genome_id"]))[:10]
    if unknown_genomes:
        raise ValueError(f"Labels reference genomes without features: {unknown_genomes}")
    for feature in spec["model_features"]:
        values = features[feature]
        if values.isna().any() or not set(values.unique()) <= {0, 1, False, True}:
            raise ValueError(f"Model feature {feature} must be complete and binary")


def validate_inference_row(row: Mapping[str, Any], spec: Mapping[str, Any]) -> None:
    missing_model_features = [feature for feature in spec["model_features"] if feature not in row]
    if missing_model_features:
        raise ValueError(f"Inference row is missing model features: {missing_model_features}")
    for feature in spec["model_features"]:
        value = row[feature]
        if _missing(value) or value not in (0, 1, False, True):
            raise ValueError(f"Inference feature {feature} must be binary, got {value!r}")


def _missing(value: Any) -> bool:
    return value is None or value is pd.NA or (isinstance(value, float) and math.isnan(value))


def _require_columns(frame: pd.DataFrame, required: set[str], name: str) -> None:
    missing = sorted(required - set(frame.columns))
    if missing:
        raise ValueError(f"{name} is missing columns: {missing}")
