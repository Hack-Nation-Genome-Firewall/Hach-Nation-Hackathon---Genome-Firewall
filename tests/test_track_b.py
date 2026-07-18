from __future__ import annotations

from copy import deepcopy

import pandas as pd
import pytest

from module2_predictor.contracts import validate_inference_row, validate_training_frames
from module2_predictor.evaluate import evaluate_bundle
from module2_predictor.predict import predict_drug
from module2_predictor.train import train_bundle
from scripts.make_synthetic_fixture import build


@pytest.fixture(scope="module")
def trained_fixture():
    features, labels, splits, spec = build()
    validate_training_frames(features, labels, splits, spec)
    bundle, summary = train_bundle(
        features,
        labels,
        splits,
        spec,
        model_version="test-v1",
        random_state=11,
    )
    return features, labels, splits, spec, bundle, summary


def test_group_leakage_is_rejected():
    features, labels, splits, spec = build()
    leaking = splits.copy()
    first_cluster = leaking.loc[0, "cluster_id"]
    indices = leaking.index[leaking["cluster_id"] == first_cluster]
    leaking.loc[indices[0], "split"] = "test"

    with pytest.raises(ValueError, match="Homology clusters cross splits"):
        validate_training_frames(features, labels, leaking, spec)


def test_nonconforming_label_evidence_is_rejected():
    features, labels, splits, spec = build()
    invalid = labels.copy()
    invalid.loc[0, "evidence"] = "Computational Method"

    with pytest.raises(ValueError, match="violate expected evidence"):
        validate_training_frames(features, invalid, splits, spec)


def test_missing_model_feature_is_an_error(trained_fixture):
    features, _labels, _splits, spec, _bundle, _summary = trained_fixture
    row = features.iloc[0].to_dict()
    row.pop(spec["model_features"][0])

    with pytest.raises(ValueError, match="missing model features"):
        validate_inference_row(row, spec)


def test_missing_target_is_a_no_call_not_assumed_present(trained_fixture):
    features, _labels, _splits, spec, bundle, _summary = trained_fixture
    row = features.iloc[0].to_dict()
    row.pop("target__meropenem")

    result = predict_drug(row, bundle, "meropenem", spec)

    assert result["verdict"] == "no_call"
    assert "target_status_unknown" in result["no_call_reasons"]
    assert result["target_gate"]["present"] is None


def test_target_absence_and_low_quality_produce_explicit_no_calls(trained_fixture):
    features, _labels, _splits, spec, bundle, _summary = trained_fixture
    target_absent = features.loc[features["genome_id"] == "synthetic_050_0"].iloc[0].to_dict()
    low_quality = features.loc[features["genome_id"] == "synthetic_051_0"].iloc[0].to_dict()

    target_result = predict_drug(target_absent, bundle, "meropenem", spec)
    quality_result = predict_drug(low_quality, bundle, "ciprofloxacin", spec)

    assert "drug_target_absent_or_disrupted" in target_result["no_call_reasons"]
    assert "low_assembly_quality" in quality_result["no_call_reasons"]


def test_feature_spec_mismatch_is_rejected(trained_fixture):
    features, _labels, _splits, spec, bundle, _summary = trained_fixture
    altered = deepcopy(spec)
    altered["status"] = "tampered"

    with pytest.raises(ValueError, match="does not match"):
        predict_drug(features.iloc[0].to_dict(), bundle, "meropenem", altered)


def test_evaluation_uses_prediction_path_and_reports_no_calls(trained_fixture):
    features, labels, splits, spec, bundle, _summary = trained_fixture

    overall, per_group, predictions, report = evaluate_bundle(features, labels, splits, spec, bundle)

    assert set(overall["drug"]) == set(spec["drugs"])
    assert (overall["no_call_rate"] > 0).all()
    assert not per_group.empty
    assert predictions["verdict"].eq("no_call").any()
    assert report["split"] == "test"
