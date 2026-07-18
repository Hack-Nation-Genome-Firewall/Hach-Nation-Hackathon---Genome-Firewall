"""Probability calibration and abstention-threshold selection."""

from __future__ import annotations

from dataclasses import dataclass

import numpy as np
from sklearn.isotonic import IsotonicRegression
from sklearn.linear_model import LogisticRegression
from sklearn.metrics import balanced_accuracy_score


@dataclass
class ProbabilityCalibrator:
    method: str
    model: object

    def predict(self, raw_probabilities: np.ndarray) -> np.ndarray:
        raw = np.asarray(raw_probabilities, dtype=float).reshape(-1)
        if self.method == "isotonic":
            values = self.model.predict(raw)
        elif self.method == "sigmoid":
            scores = _logit(raw).reshape(-1, 1)
            values = self.model.predict_proba(scores)[:, 1]
        else:
            raise ValueError(f"Unknown calibration method: {self.method}")
        return np.clip(np.asarray(values, dtype=float), 0.0, 1.0)


def fit_calibrator(
    raw_probabilities: np.ndarray,
    labels: np.ndarray,
    method: str = "auto",
) -> ProbabilityCalibrator:
    raw = np.asarray(raw_probabilities, dtype=float).reshape(-1)
    y = np.asarray(labels, dtype=int).reshape(-1)
    if len(raw) != len(y) or len(raw) < 10:
        raise ValueError("Calibration requires at least 10 aligned examples")
    classes, counts = np.unique(y, return_counts=True)
    if set(classes) != {0, 1}:
        raise ValueError("Calibration data must contain both classes")
    selected = method
    if method == "auto":
        selected = "isotonic" if len(y) >= 500 and counts.min() >= 50 else "sigmoid"
    if selected == "isotonic":
        model = IsotonicRegression(out_of_bounds="clip", y_min=0.0, y_max=1.0)
        model.fit(raw, y)
    elif selected == "sigmoid":
        model = LogisticRegression(C=1_000_000.0, solver="lbfgs", max_iter=2000)
        model.fit(_logit(raw).reshape(-1, 1), y)
    else:
        raise ValueError("Calibration method must be auto, isotonic, or sigmoid")
    return ProbabilityCalibrator(method=selected, model=model)


def choose_call_threshold(
    calibrated_probabilities: np.ndarray,
    labels: np.ndarray,
    *,
    target_balanced_accuracy: float = 0.80,
    minimum_coverage: float = 0.50,
) -> dict[str, float]:
    probabilities = np.asarray(calibrated_probabilities, dtype=float).reshape(-1)
    y = np.asarray(labels, dtype=int).reshape(-1)
    if len(probabilities) != len(y) or not len(y):
        raise ValueError("Threshold tuning requires aligned non-empty arrays")
    verdict_probability = np.maximum(probabilities, 1.0 - probabilities)
    predictions = (probabilities >= 0.5).astype(int)
    candidates: list[dict[str, float]] = []
    for threshold in np.linspace(0.50, 0.99, 50):
        called = verdict_probability >= threshold
        coverage = float(called.mean())
        if coverage < minimum_coverage or len(np.unique(y[called])) < 2:
            continue
        score = float(balanced_accuracy_score(y[called], predictions[called]))
        candidates.append(
            {
                "threshold": float(threshold),
                "coverage": coverage,
                "balanced_accuracy": score,
            }
        )
    if not candidates:
        return {"threshold": 0.50, "coverage": 1.0, "balanced_accuracy": 0.5}
    meeting_target = [item for item in candidates if item["balanced_accuracy"] >= target_balanced_accuracy]
    if meeting_target:
        return min(meeting_target, key=lambda item: item["threshold"])
    return max(candidates, key=lambda item: (item["balanced_accuracy"], item["coverage"], -item["threshold"]))


def _logit(probabilities: np.ndarray) -> np.ndarray:
    clipped = np.clip(np.asarray(probabilities, dtype=float), 1e-6, 1.0 - 1e-6)
    return np.log(clipped / (1.0 - clipped))
