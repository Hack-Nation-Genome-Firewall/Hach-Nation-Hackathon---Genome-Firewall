"""
GENOME FIREWALL — Module 2 inference.

Turns one genome's feature row into a per-drug decision record with:
  - deterministic target-presence gate
  - calibrated probability
  - no-call abstention
  - honest evidence tier + supporting markers

This is REFERENCE code trained on SYNTHETIC data so the whole team can build
against a working example. Swap in the real AMRFinderPlus features + BV-BRC
labels and retrain (module2_predictor/train.py) — the contract does not change.
"""
import json, pickle
from pathlib import Path

HERE = Path(__file__).resolve().parents[1]
SPEC = json.load(open(HERE / "data/manifests/feature_spec.json"))
MARKER_DRUG = SPEC["marker_drug"]     # marker -> drug it drives
MARKER_TYPE = SPEC["marker_types"]    # marker -> gene/point/efflux


def load_bundle(path=None):
    path = path or (HERE / "models/models.pkl")
    with open(path, "rb") as f:
        return pickle.load(f)


def _evidence(row, drug):
    """(i) known resistance marker present -> known_marker; else decided later."""
    present = [m for m, d in MARKER_DRUG.items() if d == drug and row.get(m, 0) == 1]
    return ("known_marker", present) if present else (None, [])


def predict_drug(row, bundle, drug):
    clf, iso = bundle["models"][drug], bundle["calibrators"][drug]
    order = bundle["feature_order"]
    import numpy as np
    x = np.array([[row.get(m, 0) for m in order]])
    p = float(iso.transform(clf.predict_proba(x)[:, 1])[0])   # calibrated P(fail)
    conf = abs(p - 0.5) * 2
    band = bundle["nocall_band"]

    # deterministic target-presence gate (independent of the ML model)
    targets = bundle["drug_targets"][drug]
    tgt_present = all(row.get(f"target__{g}", 1) == 1 for g in targets)
    gate = {"target_genes": targets, "present": bool(tgt_present)}

    tier, markers = _evidence(row, drug)

    reason = None
    if not tgt_present:
        verdict, reason = "no_call", "drug target absent/disrupted — cannot assert susceptibility"
        gate["action"] = "route_to_no_call"
    elif row.get("qc_complete", 1.0) < 0.90 or row.get("qc_contigs", 0) > 500:
        verdict, reason, gate["action"] = "no_call", "low assembly quality", "proceed"
    elif conf < band:
        verdict, reason, gate["action"] = "no_call", f"low confidence (conf={conf:.2f} < {band})", "proceed"
    else:
        verdict = "likely_to_fail" if p >= 0.5 else "likely_to_work"
        gate["action"] = "proceed"

    if tier is None:
        tier = "statistical_only" if verdict in ("likely_to_fail", "likely_to_work") else "no_signal"

    return {
        "drug": drug, "verdict": verdict, "p_fail": round(p, 3),
        "confidence": round(conf, 3), "calibrated": True, "evidence_tier": tier,
        "supporting_markers": [{"marker": m, "type": MARKER_TYPE[m]} for m in markers],
        "target_gate": gate, "no_call_reason": reason,
    }


def predict_genome(row, bundle=None):
    bundle = bundle or load_bundle()
    return [predict_drug(row, bundle, d) for d in SPEC["drugs"]]


if __name__ == "__main__":
    import pandas as pd
    feats = pd.read_csv(HERE / "data/manifests/features.csv")
    b = load_bundle()
    demo = feats[feats.split == "test"].iloc[0].to_dict()
    print(f"Demo genome {demo['genome_id']} (held-out cluster {demo['cluster_id']}):")
    for rec in predict_genome(demo, b):
        print(json.dumps(rec, indent=2))
