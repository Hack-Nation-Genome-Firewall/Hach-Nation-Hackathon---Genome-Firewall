"""
Detailed AI report generator for the Genome Firewall app (Track C).

A one-shot, GROUNDED, one-click-PDF per-isolate report — the formal-document sibling
of the floating chat assistant. Design choices (deliberately leaner than the chat):

  * ONE-SHOT / stateless. No conversation history — the request is [system prompt,
    grounding block, "write the report"].
  * GROUNDED. Reuses `chat_assistant.build_context()` so the report and the chat are
    built from the SAME assembled facts. The model only PHRASES the pipeline output.
  * DETERMINISTIC WORDING. `temperature=0` (facts are fixed by grounding, not sampling).
  * HYBRID DOCUMENT. The header + color-coded results table are built DETERMINISTICALLY
    from the prediction records (so the hard numbers can't be hallucinated); the AI
    writes only the interpretation prose below it.
  * GUARDED. Decision-support-only, defensive-only, synthetic-mode banner, safety footer
    on every page.

Export is a real one-click PDF via fpdf2 (pure-Python, no system libraries).
"""
from __future__ import annotations

import datetime as dt

import streamlit as st
from fpdf import FPDF
from fpdf.enums import XPos, YPos
from fpdf.fonts import FontFace

from chat_assistant import build_context  # reuse grounding -> same facts as the chat

DEFAULT_MODEL = "gpt-4o-mini"

REPORT_SYSTEM_PROMPT = """\
You write the INTERPRETATION prose of a formal ANTIBIOTIC-RESPONSE REPORT for one bacterial
isolate, for "Genome Firewall", a DEFENSIVE research prototype. A structured results table
(drug, verdict, confidence, P(fail), evidence) is rendered separately ABOVE your text, so do
NOT re-tabulate the numbers and do NOT output any Markdown table. You are given a CONTEXT block
with the isolate's per-drug records, coverage scope, held-out metrics, and method/provenance.
Write GitHub-flavored Markdown using headings and short paragraphs/bullets only.

Hard rules — never break these:
1. GROUND every statement in the CONTEXT. Never invent numbers, genes, markers, verdicts, or
   metrics not present in the CONTEXT. If something is not present, omit it — do not guess.
2. DECISION SUPPORT ONLY. Give NO clinical, treatment, dosing, or prescribing advice and never
   recommend an antibiotic for a patient. State that every result must be confirmed by standard
   laboratory testing and that a qualified clinician makes treatment decisions.
3. DEFENSIVE tool: never help make an organism more resistant, more dangerous, or evade detection.
4. If the CONTEXT indicates SYNTHETIC / demonstration mode, say near the top that the verdicts are
   a pipeline demonstration, not a biological result.
5. Explain a no-call as a deliberate abstention, not a failure. Keep KNOWN resistance markers
   separate from statistical-only associations.

Required structure (use exactly these Markdown headings, no Markdown tables):
## Summary
One short paragraph: the isolate, the overall resistance picture, and the single most important caveat.
## Interpretation by drug
One `### <drug>` subsection per drug — one or two sentences on WHY that verdict, grounded in the
fields (evidence tier, supporting markers, target gate, no-call reason).
## Assembly quality
Note the QC status if present and whether it affected any verdict.
## Method & limitations
Two or three sentences on the homology-grouped evaluation, calibration, and the no-call philosophy.
## Safety
Research prototype, decision support only, confirm with standard laboratory testing, human oversight.

Keep it concise, sober, and professional — a clinical-style document, not marketing copy.
"""


# --------------------------------------------------------------------------- #
# LLM call — single seam (mirrors chat_assistant), but one-shot / temperature 0.
# --------------------------------------------------------------------------- #
def _stream_report(context: str):
    """Yield report tokens for a single grounded generation. No history; temperature=0."""
    from openai import OpenAI

    api_key = st.secrets.get("OPENAI_API_KEY")
    if not api_key:
        raise RuntimeError("OPENAI_API_KEY is not set in .streamlit/secrets.toml")
    model = st.secrets.get("OPENAI_MODEL", DEFAULT_MODEL)

    messages = [
        {"role": "system", "content": REPORT_SYSTEM_PROMPT},
        {"role": "system", "content": f"CONTEXT for this isolate:\n{context}"},
        {"role": "user", "content": "Write the report interpretation now, following the required structure."},
    ]
    stream = OpenAI(api_key=api_key).chat.completions.create(
        model=model, messages=messages, temperature=0, stream=True,
    )
    for chunk in stream:
        delta = chunk.choices[0].delta.content
        if delta:
            yield delta


# --------------------------------------------------------------------------- #
# PDF export (fpdf2) — deterministic header + results table, then the narrative.
# --------------------------------------------------------------------------- #
ACCENT = (0, 82, 163)
INK = (26, 43, 60)
MUTED = (120, 133, 148)

VERDICT_LABEL = {"likely_to_fail": "Likely to FAIL", "likely_to_work": "Likely to work",
                 "no_call": "No-call"}
VERDICT_FILL = {"likely_to_fail": (250, 233, 233), "likely_to_work": (233, 246, 238),
                "no_call": (239, 242, 246)}
VERDICT_TEXT = {"likely_to_fail": (150, 40, 40), "likely_to_work": (28, 108, 60),
                "no_call": (90, 100, 112)}
EVIDENCE_LABEL = {"known_marker": "Known marker", "statistical_only": "Statistical only",
                  "no_signal": "No signal", "no_call": "-"}

_SUBS = {"—": "-", "–": "-", "‑": "-", "’": "'", "‘": "'",
         "“": '"', "”": '"', "…": "...", "≥": ">=", "≤": "<=",
         "×": "x", "·": "-", "→": "->", "•": "-"}


def _latin1(s) -> str:
    """Core PDF fonts are latin-1 only; map common unicode and drop the rest safely."""
    s = str(s)
    for k, v in _SUBS.items():
        s = s.replace(k, v)
    return s.encode("latin-1", "replace").decode("latin-1")


def _strip_md(s: str) -> str:
    return s.replace("**", "").replace("`", "").replace("__", "")


class _ReportPDF(FPDF):
    def footer(self):
        self.set_y(-14)
        self.set_font("Helvetica", "I", 7)
        self.set_text_color(*MUTED)
        self.multi_cell(
            0, 3.3,
            _latin1("Research prototype - decision support only, not clinical advice. Every result "
                    "must be confirmed by standard laboratory testing.   Page " + str(self.page_no())),
            align="C",
        )


def _results_table(pdf: FPDF, recs: list[dict]) -> None:
    headings = FontFace(emphasis="BOLD", color=(255, 255, 255), fill_color=ACCENT)
    pdf.set_font("Helvetica", "", 8.5)
    pdf.set_text_color(*INK)
    with pdf.table(col_widths=(28, 28, 24, 16, 30), text_align=("LEFT", "LEFT", "CENTER", "CENTER", "LEFT"),
                   headings_style=headings, line_height=6, first_row_as_headings=True) as table:
        head = table.row()
        for h in ("Drug", "Verdict", "Calib. conf.", "P(fail)", "Evidence"):
            head.cell(h)
        for rec in recs:
            v = rec.get("verdict")
            r = table.row()
            r.cell(_latin1(rec.get("drug", "")))
            r.cell(_latin1(VERDICT_LABEL.get(v, v or "-")),
                   style=FontFace(emphasis="BOLD", color=VERDICT_TEXT.get(v, INK),
                                  fill_color=VERDICT_FILL.get(v, (255, 255, 255))))
            r.cell(f"{round(rec.get('confidence', 0) * 100)}%")
            r.cell(f"{rec.get('p_fail', 0):.2f}")
            r.cell(_latin1(EVIDENCE_LABEL.get(rec.get("evidence_tier"), rec.get("evidence_tier") or "-")))


def _render_narrative(pdf: FPDF, md_text: str) -> None:
    for raw in md_text.splitlines():
        line = raw.strip()
        if not line:
            pdf.ln(1.6)
            continue
        s = _latin1(_strip_md(line))
        if s.startswith("### "):
            pdf.ln(0.6); pdf.set_font("Helvetica", "B", 10); pdf.set_text_color(*INK)
            pdf.multi_cell(0, 5, s[4:], new_x=XPos.LMARGIN, new_y=YPos.NEXT)
        elif s.startswith("## "):
            pdf.ln(1.4); pdf.set_font("Helvetica", "B", 12); pdf.set_text_color(*ACCENT)
            pdf.multi_cell(0, 6, s[3:], new_x=XPos.LMARGIN, new_y=YPos.NEXT)
        elif s.startswith("# "):
            pdf.set_font("Helvetica", "B", 13); pdf.set_text_color(*INK)
            pdf.multi_cell(0, 6, s[2:], new_x=XPos.LMARGIN, new_y=YPos.NEXT)
        elif s[:2] in ("- ", "* "):
            pdf.set_font("Helvetica", "", 9.5); pdf.set_text_color(40, 50, 62)
            pdf.multi_cell(0, 4.6, "   -  " + s[2:], new_x=XPos.LMARGIN, new_y=YPos.NEXT)
        else:
            pdf.set_font("Helvetica", "", 9.5); pdf.set_text_color(40, 50, 62)
            pdf.multi_cell(0, 4.6, s, new_x=XPos.LMARGIN, new_y=YPos.NEXT)


def _model_meta(bundle: dict | None, spec: dict) -> dict | None:
    """Model version + training-cohort size for the provenance line (deterministic)."""
    if not bundle:
        return None
    n_train = 0
    for dm in (bundle.get("drug_models") or {}).values():
        n = ((dm.get("split_counts") or {}).get("train") or {}).get("n", 0)
        n_train = max(n_train, int(n or 0))
    synthetic = bool(spec.get("synthetic"))
    return {
        "version": bundle.get("model_version") or "unversioned",
        "n_train": n_train,
        "label": "synthetic fixture" if synthetic else "BV-BRC lab-measured",
        "synthetic": synthetic,
    }


def report_pdf(md_text: str, recs: list[dict], gid: str, spec: dict,
               model_meta: dict | None = None) -> bytes:
    """Build the one-click PDF: safety banner, header, deterministic results table, narrative."""
    species = spec.get("species")
    species_name = species["name"] if isinstance(species, dict) else (species or "-")
    synthetic = bool(spec.get("synthetic"))
    ts = dt.datetime.now().strftime("%Y-%m-%d %H:%M")

    pdf = _ReportPDF(format="A4")
    pdf.set_auto_page_break(auto=True, margin=18)
    pdf.set_margins(16, 14, 16)
    pdf.add_page()

    # safety banner
    if synthetic:
        pdf.set_fill_color(255, 248, 232); pdf.set_draw_color(239, 217, 163); pdf.set_text_color(138, 109, 31)
        banner = ("SYNTHETIC / DEMONSTRATION MODE - the verdicts below are a pipeline demonstration, "
                  "not a biological result.")
    else:
        pdf.set_fill_color(255, 244, 244); pdf.set_draw_color(242, 200, 200); pdf.set_text_color(154, 43, 43)
        banner = "RESEARCH PROTOTYPE - decision support only. Not for clinical use."
    pdf.set_font("Helvetica", "B", 8.5)
    pdf.multi_cell(0, 5, _latin1(banner), border=1, fill=True, align="L",
                   new_x=XPos.LMARGIN, new_y=YPos.NEXT)
    pdf.ln(3.5)

    # title + metadata
    pdf.set_text_color(*INK); pdf.set_font("Helvetica", "B", 15)
    pdf.multi_cell(0, 8, _latin1("Genome Firewall - antibiotic-response report"),
                   new_x=XPos.LMARGIN, new_y=YPos.NEXT)
    pdf.set_font("Helvetica", "", 9); pdf.set_text_color(*MUTED)
    pdf.multi_cell(0, 5, _latin1(f"Isolate {gid}    species {species_name}    generated {ts}"),
                   new_x=XPos.LMARGIN, new_y=YPos.NEXT)
    if model_meta:
        prov = (f"Model {model_meta['version']}  -  trained on {model_meta['n_train']} "
                f"{model_meta['label']} genomes (homology-grouped split)")
        pdf.set_font("Helvetica", "", 8.5); pdf.set_text_color(*MUTED)
        pdf.multi_cell(0, 4.6, _latin1(prov), new_x=XPos.LMARGIN, new_y=YPos.NEXT)
    pdf.ln(3.5)

    # deterministic results grid
    _results_table(pdf, recs)
    pdf.set_font("Helvetica", "I", 7); pdf.set_text_color(*MUTED)
    pdf.multi_cell(0, 3.6, _latin1(
        "Calib. conf. = calibrated confidence, designed to match observed frequencies "
        "(isotonic/sigmoid-calibrated on a held-out split). P(fail) = calibrated probability the drug fails."),
        new_x=XPos.LMARGIN, new_y=YPos.NEXT)
    pdf.ln(4)

    # AI interpretation
    _render_narrative(pdf, md_text)

    return bytes(pdf.output())


# --------------------------------------------------------------------------- #
# UI — button, cached generation, preview, one-click PDF download.
# --------------------------------------------------------------------------- #
# Max report generations per browser session (a courtesy cap against runaway
# spend; the real cost ceiling is a spending limit on the OpenAI key itself).
SESSION_GENERATION_CAP = 15


def render_report_section(recs: list[dict], spec: dict, gid: str, bundle: dict | None = None) -> None:
    """Render the 'Detailed report' section.

    Cost controls: the report is a deterministic (temperature=0) artifact per genome,
    so we generate it AT MOST ONCE per genome per session (cached; the Generate button
    disappears once a report exists — no Regenerate, no re-billing the same report).
    A session-wide cap bounds total generations regardless of how the user clicks.
    """
    ss = st.session_state
    cache_key = f"gf_report::{gid}"
    ss.setdefault("gf_gen_count", 0)

    st.divider()
    st.subheader("📄 Detailed report")
    st.caption("An AI-drafted, plain-language report **grounded in the predictions above** — the "
               "model phrases the pipeline's output, it never invents data. Download as a one-click PDF.")

    md = ss.get(cache_key)

    # Not yet generated for this genome: offer Generate (unless the session cap is hit).
    if md is None:
        if ss.gf_gen_count >= SESSION_GENERATION_CAP:
            st.info(f"Report-generation limit reached for this session "
                    f"({SESSION_GENERATION_CAP}). Refresh the page to reset.")
            return
        # Render the button into a placeholder so we can clear it — but only on
        # SUCCESS. On failure the button stays put, so the user can retry without
        # refreshing (and nothing was cached / counted).
        slot = st.empty()
        with slot.container():
            st.caption("_Generating uses your OpenAI credits (one call per genome)._")
            clicked = st.button("Generate report", key=f"gf_gen_{gid}", type="primary")
        if not clicked:
            return
        try:
            context = build_context(recs, spec, gid)
            with st.container(border=True):
                md = st.write_stream(_stream_report(context))  # stream it live
            ss[cache_key] = md
            ss.gf_gen_count += 1
            slot.empty()  # success -> remove the button now (kept on failure for retry)
        except Exception as e:  # missing key / network / quota — button stays, retry available
            st.warning(
                f"Couldn't generate the report: {e}\n\n"
                "Check that `OPENAI_API_KEY` is set in `.streamlit/secrets.toml`. "
                "The **Generate** button above is still available — try again."
            )
            return
    else:
        # Already generated this session: show the cached report (no re-generation).
        with st.container(border=True):
            st.markdown(md)

    try:
        pdf_bytes = report_pdf(md, recs, gid, spec, _model_meta(bundle, spec))
        st.download_button(
            "⬇ Download report (PDF)", data=pdf_bytes,
            file_name=f"genome_firewall_report_{gid}.pdf", mime="application/pdf",
            key=f"gf_dl_{gid}", type="primary",
        )
    except Exception as e:
        st.warning(f"Couldn't build the PDF: {e}")
