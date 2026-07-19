"""
GENOME FIREWALL — Module 3: The Decision Report (Streamlit).

Run:  streamlit run TrackC/app.py

Reference demo. Runs on a SYNTHETIC fixture today (disclosed prominently in the
UI); the FASTA-upload path is wired and waiting for Track A's genome reader
(AMRFinderPlus -> feature row) — see `build_feature_row_from_fasta` below.

Visual theme: clinical-blue (see .streamlit/config.toml), adapted from the
"healthcare" theme in github.com/jmedia65/awesome-streamlit-themes.
"""
import os
import sys
from collections import Counter
from html import escape
from pathlib import Path

import pandas as pd
import streamlit as st

HERE = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(HERE))
sys.path.insert(0, str(Path(__file__).resolve().parent))  # for local charts.py
from module2_predictor.contracts import load_feature_spec  # noqa: E402
from module2_predictor.predict import load_bundle, predict_genome  # noqa: E402
from charts import performance_figure, reliability_figure  # noqa: E402
from chat_assistant import render_floating_assistant  # noqa: E402
from report_generator import render_report_section  # noqa: E402


def _flatten(html_str: str) -> str:
    """Strip per-line indentation so Streamlit markdown does not treat the HTML
    as an indented code block (which would render raw <div> tags as text)."""
    return " ".join(line.strip() for line in html_str.splitlines() if line.strip())


SPEC_PATH = Path(os.environ.get("GENOME_FIREWALL_SPEC", HERE / "data/synthetic/feature_spec.json"))
FEATURES_PATH = Path(os.environ.get("GENOME_FIREWALL_FEATURES", HERE / "data/synthetic/features.csv"))
SPLITS_PATH = Path(os.environ.get("GENOME_FIREWALL_SPLITS", HERE / "data/synthetic/split_manifest.csv"))
BUNDLE_PATH = Path(os.environ.get("GENOME_FIREWALL_BUNDLE", HERE / "models/synthetic_bundle.joblib"))
EVAL_DIR = Path(os.environ.get("GENOME_FIREWALL_EVAL", HERE / "eval"))
SPEC = load_feature_spec(SPEC_PATH)
IS_SYNTHETIC = bool(SPEC.get("synthetic"))

_FAVICON = Path(__file__).resolve().parent / "assets" / "shield.svg"
st.set_page_config(
    page_title="Genome Firewall",
    page_icon=str(_FAVICON) if _FAVICON.exists() else None,
    layout="wide",
)

# Monochrome shield mark reused for the hero + assistant panel (no emoji anywhere).
GF_MARK = (
    '<svg viewBox="0 0 24 24" fill="none" xmlns="http://www.w3.org/2000/svg" '
    'style="width:1.6rem;height:1.6rem;flex:none;vertical-align:middle">'
    '<path d="M12 2 4 5v6c0 5 3.4 8.5 8 11 4.6-2.5 8-6 8-11V5l-8-3Z" fill="#e7f0fb" '
    'stroke="#0052a3" stroke-width="1.5" stroke-linejoin="round"/>'
    '<path d="M8 12h8M12 8.5v7" stroke="#0052a3" stroke-width="1.6" stroke-linecap="round"/></svg>'
)

# ---------------------------------------------------------------------------
# Presentation layer — IBM Plex fonts + clinical-blue card system.
# ---------------------------------------------------------------------------
st.markdown(
    """
    <style>
    @import url('https://fonts.googleapis.com/css2?family=IBM+Plex+Sans:wght@400;500;600;700&family=IBM+Plex+Mono:wght@400;600&display=swap');

    html, body, [class*="css"], .stMarkdown, p, label { font-family:'IBM Plex Sans',-apple-system,sans-serif; }
    h1,h2,h3,h4 { font-family:'IBM Plex Sans',sans-serif; font-weight:600; letter-spacing:-0.01em; }
    code,pre,kbd { font-family:'IBM Plex Mono',monospace; }
    .block-container { padding-top:2rem; padding-bottom:3rem; max-width:1140px; }
    #MainMenu, footer { visibility:hidden; }
    [data-testid="stToolbar"], [data-testid="stDecoration"] { display:none; }

    .gf-hero { display:flex; align-items:center; gap:.6rem; margin:.2rem 0 .1rem; }
    .gf-hero-title { font-size:1.9rem; font-weight:700; color:#1a2b3c; }
    .gf-hero-sub { color:#5b6b7b; font-size:.95rem; }
    .gf-tagchip { display:inline-block; background:#e7f0fb; color:#0052a3; font-size:.72rem;
        font-weight:600; padding:.18rem .55rem; border-radius:1rem; letter-spacing:.02em; margin-right:.25rem; }

    /* Synthetic-mode disclosure ribbon (honest by design) */
    .gf-ribbon { border:1px solid #cfe0f5; background:linear-gradient(90deg,#eef5fd,#f7fbff);
        border-left:5px solid #eda100; border-radius:.55rem; padding:.75rem 1rem; margin:.4rem 0 .2rem;
        color:#33404d; font-size:.9rem; }
    .gf-ribbon b { color:#8a5a00; }

    /* Decision cards */
    .gf-card { background:#fff; border:1px solid #e5e9ee; border-left:5px solid #b8c2cc;
        border-radius:.6rem; padding:1rem 1.15rem; margin-bottom:.85rem; box-shadow:0 1px 2px rgba(30,50,80,.05); }
    .gf-card.fail { border-left-color:#d64545; }
    .gf-card.work { border-left-color:#1e9e63; }
    .gf-card.nocall { border-left-color:#9aa6b2; background:#fbfcfd; }
    .gf-card-head { display:flex; align-items:center; justify-content:space-between; margin-bottom:.55rem; }
    .gf-drug { font-size:1.15rem; font-weight:600; color:#1a2b3c; text-transform:capitalize; }
    .gf-badge { font-size:.74rem; font-weight:700; padding:.24rem .62rem; border-radius:1rem; letter-spacing:.03em; white-space:nowrap; }
    .gf-badge.fail { background:#fdeaea; color:#b52b2b; }
    .gf-badge.work { background:#e7f6ed; color:#137a48; }
    .gf-badge.nocall { background:#eef1f4; color:#5b6b7b; }

    .gf-row { display:flex; gap:1.4rem; align-items:flex-start; flex-wrap:wrap; }
    .gf-conf { min-width:150px; }
    .gf-conf-val { font-size:1.5rem; font-weight:700; color:#1a2b3c; font-family:'IBM Plex Mono',monospace; }
    .gf-conf-lbl { font-size:.72rem; color:#7c8996; margin-top:-.15rem; }
    .gf-bar { height:6px; background:#eef1f4; border-radius:3px; margin-top:.4rem; overflow:hidden; }
    .gf-bar > span { display:block; height:100%; border-radius:3px; }
    .gf-bar > span.fail { background:#d64545; } .gf-bar > span.work { background:#1e9e63; } .gf-bar > span.nocall { background:#9aa6b2; }

    .gf-evi { flex:1; min-width:290px; }
    /* Evidence-tier pills — visually distinct grades of evidence strength */
    .gf-tierpill { display:inline-block; font-size:.72rem; font-weight:600; padding:.16rem .5rem;
        border-radius:.35rem; margin-bottom:.4rem; border:1px solid transparent; }
    .gf-tierpill.known { background:#e7f6ed; color:#137a48; border-color:#bfe6cf; }
    .gf-tierpill.stat  { background:#fdf4e3; color:#8a5a00; border-color:#f3e2bf; }
    .gf-tierpill.none  { background:#eef1f4; color:#5b6b7b; border-color:#dde3ea; }
    .gf-tierdesc { font-size:.82rem; color:#5b6b7b; margin-bottom:.3rem; }
    .gf-chip { display:inline-block; background:#f1f3f4; color:#33404d; font-family:'IBM Plex Mono',monospace;
        font-size:.74rem; padding:.15rem .45rem; border-radius:.35rem; margin:.12rem .28rem .12rem 0; border:1px solid #e5e9ee; }
    .gf-gate { font-size:.78rem; color:#7c8996; margin-top:.45rem; }
    .gf-nocall { font-size:.82rem; color:#8a5a00; background:#fdf4e3; border:1px solid #f3e2bf;
        padding:.4rem .6rem; border-radius:.4rem; margin-top:.5rem; }
    .gf-reason { display:inline-block; background:#fff; color:#8a5a00; font-size:.72rem; padding:.1rem .4rem;
        border-radius:.3rem; margin:.15rem .25rem 0 0; border:1px solid #f0dcae; }
    </style>
    """,
    unsafe_allow_html=True,
)

# ---- mandatory safety banner (non-negotiable per brief) ----
st.error(
    "RESEARCH PROTOTYPE — Every antibiotic-response report **must be confirmed "
    "with standard laboratory testing**. This tool is decision support only and "
    "must never make a treatment decision on its own. Not for clinical use."
)

# ---- synthetic-mode disclosure (prominent, honest) ----
if IS_SYNTHETIC:
    st.markdown(
        _flatten(
            """
            <div class="gf-ribbon">
              <b>SYNTHETIC INTEGRATION MODE — we are disclosing this openly.</b><br>
              This demo runs on a <b>synthetic fixture</b>, not real BV-BRC genomes. The
              <em>machinery</em> is real and reproducible — calibration, the deterministic
              target gate, no-call abstention, the homology-grouped split and every metric.
              Only the genomes are stand-ins until Track A's AMRFinderPlus features and the
              frozen BV-BRC lab labels are wired in. Showing this honestly beats pretending
              we already have real data.
            </div>
            """
        ),
        unsafe_allow_html=True,
    )

# ---- hero header ----
st.markdown(
    _flatten(
        """
        <div class="gf-hero">""" + GF_MARK + """
          <span class="gf-hero-title">Genome Firewall</span></div>
        <div class="gf-hero-sub">Defensive decision support: it predicts and <em>explains</em> antibiotic
          resistance that already exists — it never designs, modifies, or optimizes an organism.</div>
        """
    ),
    unsafe_allow_html=True,
)
st.markdown(
    '<div style="margin:.4rem 0 .2rem">'
    '<span class="gf-tagchip">DEFENSIVE BY CONSTRUCTION</span>'
    '<span class="gf-tagchip">CALIBRATED UNCERTAINTY</span>'
    '<span class="gf-tagchip">HUMAN OVERSIGHT REQUIRED</span></div>',
    unsafe_allow_html=True,
)
st.write("")

# ---- sidebar: coverage + honest status ----
with st.sidebar:
    st.header("Coverage")
    species = SPEC["species"]
    species_name = species["name"] if isinstance(species, dict) else species
    st.markdown(f"**Species**  \n*{species_name}*")
    st.markdown("**Antibiotics**  \n" + "".join(f"- {d}\n" for d in SPEC["drugs"]))
    st.info("Outside this species / antibiotic set the system returns **no-call**.")
    st.divider()
    if IS_SYNTHETIC:
        st.warning("**Data mode: SYNTHETIC**\n\nReplace with real BV-BRC + AMRFinderPlus features to go live.")
    else:
        st.success("**Data mode: REAL** (BV-BRC lab-measured)")
    st.caption(f"Contract status: `{SPEC.get('status', 'unspecified')}`")

VERDICT_META = {
    "likely_to_fail": ("fail", "● LIKELY TO FAIL"),
    "likely_to_work": ("work", "● LIKELY TO WORK"),
    "no_call":        ("nocall", "● NO-CALL"),
}
TIER_META = {
    "known_marker":     ("known", "KNOWN RESISTANCE MARKER",
                         "A known resistance gene / DNA change was detected — the strongest evidence tier."),
    "statistical_only": ("stat", "STATISTICAL ASSOCIATION ONLY",
                         "Model association only — <b>not</b> proof of a biological mechanism. Treat with care."),
    "no_signal":        ("none", "NO KNOWN SIGNAL",
                         "No known resistance signal was found for this drug."),
}
REASON_LABEL = {
    "drug_target_absent_or_disrupted": "drug target absent / disrupted",
    "target_status_unknown": "target status unknown",
    "low_assembly_quality": "low assembly quality",
    "quality_status_unknown": "quality status unknown",
    "out_of_distribution": "out-of-distribution genome",
    "known_marker_conflicts_with_model": "known marker conflicts with model",
    "low_confidence": "confidence below call threshold",
}


def render_card(rec: dict) -> str:
    cls, badge = VERDICT_META[rec["verdict"]]
    conf = int(round(rec["confidence"] * 100))
    tcls, tlabel, tdesc = TIER_META[rec["evidence_tier"]]
    markers = "".join(
        f'<span class="gf-chip">{escape(m["marker"])} · {escape(str(m.get("type", "")))}</span>'
        for m in rec["supporting_markers"]
    )
    markers_block = f"<div>{markers}</div>" if markers else ""
    gate = rec["target_gate"]
    gate_txt = (f'Target gate: {escape(", ".join(gate["target_features"]))} — '
                f'{escape(gate["status"].replace("_", " "))} → {escape(gate["action"])}')
    reasons = rec.get("no_call_reasons") or ([rec["no_call_reason"]] if rec.get("no_call_reason") else [])
    nocall = ""
    if reasons:
        chips = "".join(f'<span class="gf-reason">{escape(REASON_LABEL.get(r, r))}</span>' for r in reasons)
        nocall = (f'<div class="gf-nocall"><b>Why no-call:</b> not enough trustworthy evidence '
                  f'to make a call.<br>{chips}</div>')
    card = f"""
    <div class="gf-card {cls}">
      <div class="gf-card-head">
        <span class="gf-drug">{escape(rec["drug"])}</span>
        <span class="gf-badge {cls}">{badge}</span>
      </div>
      <div class="gf-row">
        <div class="gf-conf">
          <div class="gf-conf-val">{conf}%</div>
          <div class="gf-conf-lbl">confidence · calibrated P(fail)={rec["p_fail"]:.3f}</div>
          <div class="gf-bar"><span class="{cls}" style="width:{conf}%"></span></div>
        </div>
        <div class="gf-evi">
          <span class="gf-tierpill {tcls}">{tlabel}</span>
          <div class="gf-tierdesc">{tdesc}</div>
          {markers_block}
          <div class="gf-gate">{gate_txt}</div>
          {nocall}
        </div>
      </div>
    </div>
    """
    return _flatten(card)


MODULE1_DIR = HERE / "module1_reader"
SAMPLE_TSV = MODULE1_DIR / "fixtures" / "sample_amrfinder.tsv"


def _amrfinder_available() -> bool:
    """Is the AMRFinderPlus binary on PATH? Decides whether the upload tab defaults
    to real annotation or to the tool-free wiring demo."""
    import shutil
    return shutil.which("amrfinder") is not None


def build_feature_row_from_fasta(
    fasta_bytes: bytes, *, genome_id: str = "uploaded_genome",
    organism: str = "Klebsiella_pneumoniae", tsv_override=None,
) -> tuple[dict, list]:
    """SEAM for Track A — now wired to `module1_reader.run_genome_reader`.

    Writes the uploaded bytes to a temp FASTA and hands it to Track A's genome
    reader, pinning it to THIS app's feature spec so the row it returns already
    matches the contract the deployed bundle was trained on. Track A's files are
    used unmodified. Returns (validated feature row, unknown markers preserved by
    the reader). `tsv_override` runs the reader against a saved AMRFinderPlus TSV
    (Track A's bundled sample) so the wiring is demonstrable without the tool
    installed. Raises FileNotFoundError if AMRFinderPlus is not on PATH."""
    import tempfile
    if str(MODULE1_DIR) not in sys.path:
        sys.path.insert(0, str(MODULE1_DIR))  # Track A uses bare intra-module imports
    from build_features import run_genome_reader  # Track A entry point (untouched)

    unknown: list = []
    with tempfile.NamedTemporaryFile("wb", suffix=".fasta") as tmp:
        tmp.write(fasta_bytes)
        tmp.flush()
        row = run_genome_reader(
            tmp.name, genome_id=genome_id, backend="amrfinderplus", spec=SPEC,
            organism=organism, tsv_override=tsv_override, unknown_markers_out=unknown,
        )
    row.setdefault("genome_id", genome_id)
    return row, unknown


# ---------------------------------------------------------------------------
# Input: demo held-out genome  OR  upload a FASTA (ready for Track A).
# ---------------------------------------------------------------------------
feats = pd.read_csv(FEATURES_PATH, dtype={"genome_id": str})
splits = pd.read_csv(SPLITS_PATH, dtype={"genome_id": str, "cluster_id": str})
held = feats.merge(splits, on="genome_id", validate="one_to_one")
held = held[held.split == "test"]

if not BUNDLE_PATH.exists():
    st.error("Model bundle is missing. Run `python -m module2_predictor.train` first.")
    st.stop()
bundle = load_bundle(BUNDLE_PATH)


def render_genome_report(row: dict, gid: str, *, from_upload: bool = False) -> list:
    """Render the genome-specific report (verdict strip + cards + AI report) for one
    genome. Called INSIDE each input tab so every tab owns its own report — the Upload
    tab therefore stays blank until a FASTA is actually uploaded, instead of sharing
    (and overwriting) the demo genome's report. Returns the per-drug recommendations."""
    recs = predict_genome(row, bundle, SPEC)
    counts = Counter(r["verdict"] for r in recs)
    st.subheader(f"Antibiotic-response report — `{gid}`")

    if from_upload:
        _uploaded = st.session_state.get("_uploaded")
        if _uploaded and _uploaded[0] == gid:
            _name, _unknown, _via_sample = _uploaded
            _src = "Track A's bundled sample annotation" if _via_sample else "AMRFinderPlus"
            _found = ", ".join(f"`{m}`" for m in _unknown[:12]) + (" …" if len(_unknown) > 12 else "")
            st.info(
                f"**Live Track A wiring.** This report was built end to end from your "
                f"upload: **{_name}** → {_src} → feature row → prediction. The reader parsed "
                f"the genome and preserved **{len(_unknown)} marker(s)** it found"
                + (f": {_found}" if _unknown else "")
                + ". Because the *deployed model* is still the synthetic fixture (its "
                "vocabulary is `marker__known__*`, not real gene names), those real markers "
                "are carried as **unknown** rather than scored — so the verdicts below are a "
                "**pipeline demonstration, not a biological result**. They become meaningful "
                "once Phase 0 publishes the shared spec and Track B trains a bundle on it.")

    m1, m2, m3, m4 = st.columns(4)
    m1.metric("Antibiotics", len(recs))
    m2.metric("Likely to work", counts.get("likely_to_work", 0))
    m3.metric("Likely to fail", counts.get("likely_to_fail", 0))
    m4.metric("No-call", counts.get("no_call", 0))
    st.write("")
    st.markdown("".join(render_card(r) for r in recs), unsafe_allow_html=True)
    # Detailed AI report (grounded, one-shot) — download/print as PDF.
    render_report_section(recs, SPEC, gid, bundle)
    return recs


def _use_demo():
    st.session_state["_active_source"] = "demo"


def _use_upload():
    st.session_state["_active_source"] = "upload"


# Each tab renders its OWN report, so switching tabs is a pure client-side swap and
# the Upload tab shows nothing until a FASTA is uploaded. `_active_source` (set by the
# widgets' on_change) records the input the user last touched — it decides which genome
# the single floating assistant answers about (its widget keys are fixed, so it must be
# rendered exactly once, outside the tabs).
demo_recs, demo_gid = None, None
up_recs, up_gid = None, None

tab_demo, tab_upload = st.tabs(["Demo genome (held-out)", "Upload a genome (FASTA)"])
with tab_demo:
    gid = st.selectbox("Held-out demo genome:", held.genome_id.tolist(),
                       key="demo_gid", on_change=_use_demo)
    row = held[held.genome_id == gid].iloc[0].to_dict()
    demo_recs, demo_gid = render_genome_report(row, gid), gid
with tab_upload:
    up = st.file_uploader("Assembled genome — FASTA (.fasta / .fa / .fna)",
                          type=["fasta", "fa", "fna"], key="fasta_up",
                          on_change=_use_upload)
    st.caption("This runs the real pipeline: your genome → Track A's reader "
               "(`module1_reader.run_genome_reader`) → feature row → prediction. "
               "The report below then renders for **your** genome.")
    demo_wiring = st.checkbox(
        "AMRFinderPlus isn't installed here — prove the wiring with Track A's bundled "
        "sample annotation instead", value=not _amrfinder_available(),
        help="Runs the reader against module1_reader/fixtures/sample_amrfinder.tsv "
             "(read-only) so the FASTA → reader → prediction path is demonstrable "
             "without the tool.")
    if up is not None:
        tsv_override = str(SAMPLE_TSV) if (demo_wiring and SAMPLE_TSV.exists()) else None
        try:
            up_row, up_unknown = build_feature_row_from_fasta(
                up.getvalue(), genome_id=up.name, tsv_override=tsv_override)
            st.session_state["_uploaded"] = (up.name, up_unknown, bool(tsv_override))
            st.success(f"Track A's reader parsed **{up.name}** "
                       f"({up.size/1000:.0f} kB) into a contract-valid feature row.")
            up_recs, up_gid = render_genome_report(up_row, up.name, from_upload=True), up.name
        except FileNotFoundError:
            st.warning(
                "The reader is fully wired, but **AMRFinderPlus is not installed on "
                "this host**, so a raw FASTA can't be annotated here. Tick the box "
                "above to demonstrate the FASTA → reader → prediction path with Track "
                "A's bundled sample annotation, or install AMRFinderPlus to run for real.")
        except Exception as e:  # ContractError and friends — surface the real reason
            st.error(f"Track A's reader could not build a feature row: {e}")
    else:
        st.info("⬆️ Upload an assembled genome (FASTA) above to generate its "
                "antibiotic-response report here.")

# ---------------------------------------------------------------------------
# Held-out performance (interactive) + metrics table.  [model-level — always shown]
# ---------------------------------------------------------------------------
st.divider()
st.subheader("Held-out performance & calibration")
st.caption("Evaluated on a **homology-grouped** test split — near-identical genomes never "
           "span train/test, so these numbers are not inflated by leakage.")

overall_path = EVAL_DIR / "overall_metrics.csv"
pred_path = EVAL_DIR / "held_out_predictions.csv"
if overall_path.exists():
    odf = pd.read_csv(overall_path)
    st.plotly_chart(performance_figure(odf), use_container_width=True,
                    theme=None, config={"displayModeBar": False})
    show = ["drug", "n", "balanced_accuracy", "recall_resistant", "recall_susceptible",
            "f1", "auroc", "pr_auc", "brier", "no_call_rate"]
    show = [c for c in show if c in odf.columns]
    st.dataframe(
        odf[show].rename(columns={
            "balanced_accuracy": "bal_acc", "recall_resistant": "recall_R",
            "recall_susceptible": "recall_S", "no_call_rate": "no_call"}),
        use_container_width=True, hide_index=True,
    )
else:
    st.caption("Run `python -m module2_predictor.evaluate` to populate performance metrics.")

if pred_path.exists():
    with st.expander("Calibration reliability — predicted vs. observed (interactive)", expanded=True):
        st.caption("Perfect calibration follows the dotted diagonal. Hover any point for the "
                   "predicted probability vs. the observed resistant fraction in that bin.")
        pdf = pd.read_csv(pred_path)
        st.plotly_chart(reliability_figure(pdf, SPEC["drugs"]), use_container_width=True,
                        theme=None, config={"displayModeBar": False})
else:
    st.caption("Run `python -m module2_predictor.evaluate` to populate the reliability curves.")

st.divider()
st.caption("Human oversight required: a trained healthcare or laboratory professional "
           "must confirm every result before any treatment decision.")

# ---------------------------------------------------------------------------
# Floating "explain this report" assistant (grounded + safety-guarded).
# Rendered once, for whichever genome is the active input source (last touched).
# ---------------------------------------------------------------------------
if st.session_state.get("_active_source") == "upload" and up_recs is not None:
    fa_recs, fa_gid = up_recs, up_gid
else:
    fa_recs, fa_gid = demo_recs, demo_gid
render_floating_assistant(fa_recs, SPEC, fa_gid)
