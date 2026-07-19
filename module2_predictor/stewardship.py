"""
Stewardship recommendation layer (Track A/B add-on).

Turns a per-drug antibiogram (list of predict_drug records) into a single,
defensible clinical-support recommendation:

  * prefer the NARROWEST-spectrum agent predicted to work (likely_to_work),
    supporting antibiotic stewardship — narrow-spectrum use slows the very
    resistance this tool exists to detect;
  * never override a no_call: if every agent is no_call or likely_to_fail,
    the recommendation is to escalate to full laboratory AST;
  * carry the mechanism + calibrated confidence through, so the recommendation
    is explainable, not a black box.

This is DECISION SUPPORT ONLY. Every recommendation must be confirmed by
standard laboratory antimicrobial susceptibility testing before use.
"""
from __future__ import annotations
from typing import Any, Mapping, Sequence

# Narrow -> broad spectrum ordering for the four in-scope drugs.
# Lower rank = narrower spectrum = preferred when efficacy is equal.
# (gentamicin: narrow aminoglycoside; ciprofloxacin: fluoroquinolone;
#  ceftazidime: 3rd-gen cephalosporin; meropenem: last-resort carbapenem,
#  deliberately ranked broadest so it is spared unless nothing narrower works.)
SPECTRUM_RANK = {
    "gentamicin": 1,
    "ciprofloxacin": 2,
    "ceftazidime": 3,
    "meropenem": 4,
}
DEFAULT_RANK = 99  # unknown drug -> treated as broad, deprioritised


def _rank(drug: str) -> int:
    return SPECTRUM_RANK.get(drug, DEFAULT_RANK)


def recommend(
    antibiogram: Sequence[Mapping[str, Any]],
    *,
    min_confidence: float = 0.0,
) -> dict[str, Any]:
    """
    antibiogram: list of predict_drug() records (each has drug/verdict/confidence/
                 supporting_markers/evidence_tier/no_call_reasons).
    min_confidence: optional extra floor on top of the model's own calibrated
                    call_threshold (the model has already abstained below its
                    threshold; this only *tightens*, never loosens).

    Returns a recommendation dict; safe by construction (no_call -> escalate).
    """
    by_drug = {r["drug"]: r for r in antibiogram}

    workable = [
        r for r in antibiogram
        if r.get("verdict") == "likely_to_work"
        and float(r.get("confidence", 0.0)) >= min_confidence
    ]
    workable.sort(key=lambda r: (_rank(r["drug"]), -float(r.get("confidence", 0.0))))

    resistant = sorted(
        (r for r in antibiogram if r.get("verdict") == "likely_to_fail"),
        key=lambda r: _rank(r["drug"]),
    )
    no_call = [r for r in antibiogram if r.get("verdict") == "no_call"]

    if workable:
        top = workable[0]
        alternatives = [w["drug"] for w in workable[1:]]
        rec = {
            "recommendation": "use_narrowest_effective_agent",
            "primary_choice": top["drug"],
            "primary_confidence": round(float(top.get("confidence", 0.0)), 4),
            "primary_mechanism_tier": top.get("evidence_tier"),
            "primary_supporting_markers": [
                m.get("marker") for m in top.get("supporting_markers", [])
            ],
            "alternatives_narrow_to_broad": alternatives,
            "spectrum_rationale": (
                f"{top['drug']} is the narrowest-spectrum agent predicted to work "
                f"(spectrum rank {_rank(top['drug'])}); reserving broader agents "
                f"supports stewardship."
            ),
            "escalate_to_lab": False,
        }
    else:
        rec = {
            "recommendation": "escalate_to_laboratory_ast",
            "primary_choice": None,
            "reason": (
                "no in-scope agent is confidently predicted to work"
                + (" (all predicted to fail)" if resistant and not no_call else "")
                + (" (insufficient evidence / abstained)" if no_call else "")
            ),
            "predicted_resistant": [r["drug"] for r in resistant],
            "abstained": [r["drug"] for r in no_call],
            "escalate_to_lab": True,
        }

    rec["antibiogram_summary"] = {
        r["drug"]: {
            "verdict": r.get("verdict"),
            "confidence": round(float(r.get("confidence", 0.0)), 4),
            "evidence_tier": r.get("evidence_tier"),
        }
        for r in antibiogram
    }
    rec["safety_notice"] = (
        "Decision support only. Confirm with standard laboratory antimicrobial "
        "susceptibility testing before any treatment decision."
    )
    return rec


def recommend_from_genome(row, bundle=None, spec=None, **kw):
    """Convenience: run inference then recommend, in one call."""
    from module2_predictor.predict import predict_genome
    antibiogram = predict_genome(row, bundle, spec)
    return {"predictions": antibiogram, "stewardship": recommend(antibiogram, **kw)}
