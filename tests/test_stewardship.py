"""Tests for the stewardship recommendation layer."""
from module2_predictor.stewardship import recommend


def _rec(drug, verdict, conf, tier="known_marker", markers=None):
    return {"drug": drug, "verdict": verdict, "confidence": conf,
            "evidence_tier": tier,
            "supporting_markers": [{"marker": m} for m in (markers or [])]}


def test_prefers_narrowest_effective_agent():
    ab = [
        _rec("meropenem", "likely_to_work", 0.95),
        _rec("gentamicin", "likely_to_work", 0.90),   # narrower -> should win
        _rec("ceftazidime", "likely_to_fail", 0.80),
    ]
    r = recommend(ab)
    assert r["escalate_to_lab"] is False
    assert r["primary_choice"] == "gentamicin"           # narrowest workable
    assert "meropenem" in r["alternatives_narrow_to_broad"]


def test_all_fail_escalates():
    ab = [
        _rec("meropenem", "likely_to_fail", 0.99, markers=["marker__blaKPC"]),
        _rec("ciprofloxacin", "likely_to_fail", 0.97),
        _rec("gentamicin", "likely_to_fail", 0.85),
        _rec("ceftazidime", "likely_to_fail", 0.9),
    ]
    r = recommend(ab)
    assert r["escalate_to_lab"] is True
    assert r["primary_choice"] is None
    assert set(r["predicted_resistant"]) == {"meropenem", "ciprofloxacin", "gentamicin", "ceftazidime"}


def test_no_call_never_recommended():
    ab = [
        _rec("meropenem", "no_call", 0.55),
        _rec("gentamicin", "no_call", 0.52),
    ]
    r = recommend(ab)
    assert r["escalate_to_lab"] is True
    assert set(r["abstained"]) == {"meropenem", "gentamicin"}


def test_min_confidence_tightens_only():
    ab = [_rec("gentamicin", "likely_to_work", 0.60)]
    assert recommend(ab)["primary_choice"] == "gentamicin"
    # raising the floor above the call confidence forces escalation
    assert recommend(ab, min_confidence=0.8)["escalate_to_lab"] is True


def test_safety_notice_always_present():
    for ab in ([_rec("gentamicin", "likely_to_work", 0.9)],
               [_rec("gentamicin", "no_call", 0.5)]):
        assert "laboratory" in recommend(ab)["safety_notice"].lower()
