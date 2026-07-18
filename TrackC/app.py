"""
GENOME FIREWALL — Module 3: The Decision Report (Streamlit).

Run:  streamlit run module3_app/app.py

Reference demo on SYNTHETIC data. In production, wire the FASTA uploader to
Module 1 (AMRFinderPlus) to build the feature row, then call predict_genome().
"""
import os
import sys
from pathlib import Path
import pandas as pd
import streamlit as st

HERE = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(HERE))
from module2_predictor.contracts import load_feature_spec  # noqa: E402
from module2_predictor.predict import load_bundle, predict_genome  # noqa: E402

SPEC_PATH = Path(os.environ.get("GENOME_FIREWALL_SPEC", HERE / "data/synthetic/feature_spec.json"))
FEATURES_PATH = Path(os.environ.get("GENOME_FIREWALL_FEATURES", HERE / "data/synthetic/features.csv"))
SPLITS_PATH = Path(os.environ.get("GENOME_FIREWALL_SPLITS", HERE / "data/synthetic/split_manifest.csv"))
BUNDLE_PATH = Path(os.environ.get("GENOME_FIREWALL_BUNDLE", HERE / "models/synthetic_bundle.joblib"))
SPEC = load_feature_spec(SPEC_PATH)

st.set_page_config(page_title="Genome Firewall", layout="wide")

# ---- mandatory safety banner (non-negotiable per brief) ----
st.error(
    "⚠️ RESEARCH PROTOTYPE — Every antibiotic-response report **must be confirmed "
    "with standard laboratory testing**. This tool is decision support only and "
    "must never make a treatment decision on its own. Not for clinical use."
)

st.title("🧬 Genome Firewall")
st.caption("Defensive decision support: predicts and explains antibiotic resistance "
           "that already exists. It never designs, modifies, or optimizes an organism.")

if SPEC.get("synthetic"):
    st.warning(
        "Synthetic integration fixture: these results test the software only and "
        "must not be presented as biological performance."
    )

# ---- coverage statement (honest scope) ----
with st.sidebar:
    st.header("Coverage")
    st.write(f"**Species:** *{SPEC['species']}*")
    st.write("**Antibiotics:** " + ", ".join(SPEC["drugs"]))
    st.info("Outside this species/antibiotic set the system returns **no-call**.")
    st.divider()
    st.caption(f"Contract status: {SPEC.get('status', 'unspecified')}")

VERDICT_STYLE = {
    "likely_to_fail": ("🔴 Likely to FAIL", "The antibiotic is predicted NOT to work."),
    "likely_to_work": ("🟢 Likely to WORK", "The antibiotic is predicted to work."),
    "no_call":        ("⚪ NO-CALL", "Evidence too weak / conflicting / out-of-distribution."),
}
TIER_LABEL = {
    "known_marker": "🧬 Known resistance gene / DNA change detected",
    "statistical_only": "📊 Statistical association only (NOT proof of biological cause)",
    "no_signal": "— No known resistance signal found",
}

feats = pd.read_csv(FEATURES_PATH, dtype={"genome_id": str})
splits = pd.read_csv(SPLITS_PATH, dtype={"genome_id": str, "cluster_id": str})
held = feats.merge(splits, on="genome_id", validate="one_to_one")
held = held[held.split == "test"]
gid = st.selectbox("Choose a held-out demo genome (or wire the FASTA uploader to Module 1):",
                   held.genome_id.tolist())
row = held[held.genome_id == gid].iloc[0].to_dict()

if not BUNDLE_PATH.exists():
    st.error("Model bundle is missing. Run `python -m module2_predictor.train` first.")
    st.stop()
bundle = load_bundle(BUNDLE_PATH)
recs = predict_genome(row, bundle, SPEC)

st.subheader(f"Antibiotic-response report — genome `{gid}`")
for rec in recs:
    head, expl = VERDICT_STYLE[rec["verdict"]]
    with st.container(border=True):
        c1, c2, c3 = st.columns([2, 1, 3])
        c1.markdown(f"### {rec['drug']}")
        c1.markdown(f"**{head}**")
        c2.metric("Confidence", f"{rec['confidence']*100:.0f}%")
        c2.caption(f"Calibrated P(fail)={rec['p_fail']}")
        c3.markdown(f"**Evidence:** {TIER_LABEL[rec['evidence_tier']]}")
        if rec["supporting_markers"]:
            c3.markdown("**Supporting markers:** " +
                        ", ".join(f"`{m['marker']}` ({m['type']})" for m in rec["supporting_markers"]))
        gate = rec["target_gate"]
        target_status = gate["status"].replace("_", " ")
        c3.caption(f"Target gate: {', '.join(gate['target_features'])} — "
                   f"{target_status} → {gate['action']}")
        if rec["no_call_reason"]:
            c3.warning(f"No-call reason: {rec['no_call_reason']}")

with st.expander("Why you can trust the *uncertainty* (calibration & generalization)"):
    st.write("- Models evaluated on a **homology-grouped** held-out split "
             "(near-identical genomes never span train/test).")
    st.write("- Confidence scores are **isotonic-calibrated**; see reliability plots in `eval/`.")
    st.write("- **No-call** is returned for weak, conflicting, or out-of-distribution evidence.")
    reliability = HERE / "eval/fig_reliability.png"
    if reliability.exists():
        st.image(str(reliability))
    else:
        st.caption("Reliability figure has not been generated for this bundle.")

st.divider()
st.caption("Human oversight required: a trained healthcare or laboratory professional "
           "must confirm every result before any treatment decision.")
