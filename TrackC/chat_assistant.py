"""
Floating "explain this report" assistant for the Genome Firewall app.

A small circular button floats in the bottom-right corner. Clicking it expands a
chat panel that answers questions about the report the user is currently looking
at. The assistant is deliberately narrow:

  * It is GROUNDED — every answer is built from a context block assembled from the
    on-screen prediction records, the coverage spec, the held-out metrics, and the
    committed provenance docs. It is told to say "I don't know" outside that.
  * It is GUARDED — the system prompt forbids clinical/treatment/dosing advice and
    refuses dual-use requests (making organisms more resistant/dangerous). This is
    a defensive tool and the chat must not undermine the app's safety framing.

LLM calls go through one function (`_stream_completion`) so the provider/model can
be swapped in one place. The key is read from `st.secrets["OPENAI_API_KEY"]`.
"""
from __future__ import annotations

import json
from functools import lru_cache
from pathlib import Path

import streamlit as st

HERE = Path(__file__).resolve().parents[1]
DOCS_DIR = HERE / "docs"

DEFAULT_MODEL = "gpt-4o-mini"

SYSTEM_PROMPT = """\
You are the on-page assistant for "Genome Firewall", a DEFENSIVE research prototype
that predicts and explains antibiotic resistance already present in a bacterial
genome. Your only job is to help a user understand the report and method currently
on their screen, using the CONTEXT provided below.

Hard rules — never break these:
1. GROUND every answer in the CONTEXT block. If the answer is not in the CONTEXT,
   say you don't have that information on this page. Never invent numbers, markers,
   genes, or verdicts.
2. NEVER give clinical, treatment, dosing, or prescribing advice, and never tell a
   user which antibiotic to give a patient. If asked, explain that this is decision
   support only, that every result must be confirmed by standard laboratory testing,
   and that a qualified clinician makes treatment decisions.
3. This is a DEFENSIVE tool. Refuse any request to design, enhance, increase, or
   engineer resistance, or to make an organism more dangerous or evade detection.
   Briefly say why and stop.
4. If the report is in SYNTHETIC / demonstration mode, remind the user the verdicts
   are a pipeline demonstration, not a biological result.
5. Be concise and plain-language. When you cite a figure, name the field it came
   from (e.g. "the no-call reason", "the target gate", "calibrated P(fail)").
"""


# --------------------------------------------------------------------------- #
# Grounding — turn the live report + docs into a compact context block.
# --------------------------------------------------------------------------- #
@lru_cache(maxsize=1)
def _load_docs() -> str:
    """Concatenate the committed provenance/method docs (cached; they don't change
    within a run). These let the assistant answer 'where does the data come from?'."""
    parts = []
    for name in ("DATA_PROVENANCE.md", "TRACK_A_HANDOFF.md"):
        p = DOCS_DIR / name
        if p.exists():
            parts.append(f"### {name}\n{p.read_text(encoding='utf-8')}")
    return "\n\n".join(parts)


def _rec_summary(rec: dict) -> dict:
    """Project one prediction record down to the fields worth grounding on."""
    gate = rec.get("target_gate") or {}
    reasons = rec.get("no_call_reasons") or (
        [rec["no_call_reason"]] if rec.get("no_call_reason") else []
    )
    return {
        "drug": rec.get("drug"),
        "verdict": rec.get("verdict"),
        "confidence": round(rec.get("confidence", 0), 3),
        "p_fail": round(rec.get("p_fail", 0), 3),
        "evidence_tier": rec.get("evidence_tier"),
        "supporting_markers": [
            f'{m.get("marker")} ({m.get("type", "")})'
            for m in rec.get("supporting_markers", [])
        ],
        "target_gate": {
            "status": gate.get("status"),
            "action": gate.get("action"),
            "target_features": gate.get("target_features"),
        },
        "no_call_reasons": reasons,
    }


def build_context(recs: list[dict], spec: dict, gid: str) -> str:
    """Assemble the CONTEXT block injected into every LLM call for this report."""
    species = spec.get("species")
    species_name = species["name"] if isinstance(species, dict) else species
    header = {
        "genome_id": gid,
        "synthetic_mode": bool(spec.get("synthetic")),
        "species_in_scope": species_name,
        "drugs_in_scope": spec.get("drugs"),
        "note": "Outside this species/drug set the system returns no-call.",
    }

    metrics_csv = ""
    mpath = HERE / "eval" / "overall_metrics.csv"
    if mpath.exists():
        metrics_csv = mpath.read_text(encoding="utf-8")

    return (
        "=== SCOPE ===\n"
        + json.dumps(header, ensure_ascii=False, indent=2)
        + "\n\n=== CURRENT REPORT (per-drug prediction records on screen) ===\n"
        + json.dumps([_rec_summary(r) for r in recs], ensure_ascii=False, indent=2)
        + "\n\n=== HELD-OUT PERFORMANCE (eval/overall_metrics.csv) ===\n"
        + (metrics_csv or "(metrics not generated yet)")
        + "\n\n=== METHOD & DATA PROVENANCE (project docs) ===\n"
        + _load_docs()
    )


def suggested_prompts(recs: list[dict]) -> list[str]:
    """Dynamic quick-prompt chips derived from THIS report, so the shortcuts always
    point at something actually on screen."""
    prompts = ["What is this report telling me?"]
    fails = [r["drug"] for r in recs if r.get("verdict") == "likely_to_fail"]
    nocalls = [r["drug"] for r in recs if r.get("verdict") == "no_call"]
    if fails:
        prompts.append(f"Why is {fails[0]} likely to fail?")
    if nocalls:
        prompts.append(f"Why is {nocalls[0]} a no-call?")
    prompts.append("What does 'calibrated confidence' mean here?")
    prompts.append("Where does the data come from?")
    return prompts[:4]


# --------------------------------------------------------------------------- #
# LLM call — single seam, swap provider/model here.
# --------------------------------------------------------------------------- #
def _stream_completion(context: str, history: list[dict]):
    """Yield reply tokens for the latest turn. Reads the key from st.secrets."""
    from openai import OpenAI

    api_key = st.secrets.get("OPENAI_API_KEY")
    if not api_key:
        raise RuntimeError("OPENAI_API_KEY is not set in .streamlit/secrets.toml")
    model = st.secrets.get("OPENAI_MODEL", DEFAULT_MODEL)

    messages = [
        {"role": "system", "content": SYSTEM_PROMPT},
        {"role": "system", "content": f"CONTEXT for the report on screen:\n{context}"},
        *history,
    ]
    stream = OpenAI(api_key=api_key).chat.completions.create(
        model=model, messages=messages, temperature=0.2, stream=True,
    )
    for chunk in stream:
        delta = chunk.choices[0].delta.content
        if delta:
            yield delta


# --------------------------------------------------------------------------- #
# Floating widget.
#
# Positioned with native `st.container(key=...)` + CSS `position:fixed` — no
# custom bidirectional component (streamlit-float segfaulted on rerun on macOS).
# Streamlit stamps `class="st-key-<key>"` on the container div, which the CSS
# below targets to pin the button and panel to the bottom-right corner.
# --------------------------------------------------------------------------- #
@lru_cache(maxsize=1)
def _fab_icon_uri() -> str:
    """Base64 data URI for the floating-button icon (committed under assets/).
    Empty string if the asset is missing — the CSS then falls back to a glyph."""
    import base64
    p = Path(__file__).resolve().parent / "assets" / "chatgpt-logo.png"
    if not p.exists():
        return ""
    return "data:image/png;base64," + base64.b64encode(p.read_bytes()).decode()


# Bulk CSS for the panel — plain string (single braces), no f-string interpolation.
_PANEL_CSS = """
/* ---- expanded chat panel ---- */
.st-key-gf_panel { position:fixed; right:1.9rem; bottom:6.6rem; z-index:1000;
    width:min(30rem,94vw); max-height:78vh; overflow-y:auto;
    background:#fff; border:1px solid #e5e9ee; border-radius:1rem;
    padding:1.1rem 1.2rem .3rem; box-shadow:0 18px 50px rgba(30,50,80,.24);
    font-family:'IBM Plex Sans',-apple-system,sans-serif; }
.st-key-gf_panel::-webkit-scrollbar { width:8px; }
.st-key-gf_panel::-webkit-scrollbar-thumb { background:#dbe2ea; border-radius:4px; }

/* header */
.gf-assistant-hd { font-weight:700; color:#1a2b3c; font-size:1.06rem; letter-spacing:-.01em;
    display:flex; align-items:center; gap:.4rem; }
.gf-assistant-sub { color:#8794a1; font-size:.76rem; line-height:1.35; margin-top:.15rem;
    padding-bottom:.7rem; border-bottom:1px solid #eef1f4; }

/* close button (top-right) */
.st-key-gf_closewrap button { width:2.1rem; height:2.1rem; min-height:0; padding:0;
    border-radius:.55rem; background:transparent; border:1px solid #e5e9ee; color:#8794a1;
    font-size:.95rem; box-shadow:none; }
.st-key-gf_closewrap button:hover { background:#f4f6f8; border-color:#d3dae2; color:#33404d; }

/* "quick questions" label */
.gf-ql { font-size:.68rem; font-weight:700; letter-spacing:.08em; text-transform:uppercase;
    color:#9aa6b2; margin:.7rem 0 .45rem; }

/* quick-question chips — left-aligned, stacked, calm clinical-blue cards */
.st-key-gf_chips button { width:100%; text-align:left; justify-content:flex-start;
    background:#f5f8fc; border:1px solid #e1e9f2; color:#33506e; border-radius:.6rem;
    padding:.55rem .75rem; margin-bottom:.4rem; min-height:0; font-size:.83rem; font-weight:500;
    line-height:1.35; white-space:normal; box-shadow:none; transition:all .12s ease; }
.st-key-gf_chips button:hover { background:#e9f1fb; border-color:#b9d3ee; color:#0052a3; }
.st-key-gf_chips button p { text-align:left; margin:0; }

/* "quick questions" collapse toggle — used once a conversation has started */
.st-key-gf_qtoggle button { width:100%; text-align:left; justify-content:flex-start;
    background:transparent; border:1px solid #eef1f4; border-radius:.6rem; color:#8794a1;
    padding:.38rem .7rem; min-height:0; box-shadow:none; }
.st-key-gf_qtoggle button:hover { color:#0052a3; border-color:#dbe6f2; background:#f7fafd; }
.st-key-gf_qtoggle button p { margin:0; font-size:.72rem; font-weight:600; letter-spacing:.04em;
    text-transform:uppercase; }

/* transcript — a native fixed-height scroll area so the input below it stays put */
.st-key-gf_scroll { border:none !important; padding:0 .1rem 0 0 !important; }
.st-key-gf_scroll [data-testid="stChatMessage"] { padding:.35rem .2rem; gap:.55rem; }
.st-key-gf_scroll [data-testid="stChatMessage"] p { font-size:.87rem; line-height:1.5; }

/* input box */
.st-key-gf_panel [data-testid="stChatInput"] { border-radius:.7rem; margin-top:.15rem; }
"""


def _fixed_css() -> str:
    """Positioning + look for the FAB and panel. The FAB rule is built by string
    concatenation (not an f-string) so the committed PNG can be embedded as a data
    URI without having to brace-escape the rest of the CSS."""
    icon = _fab_icon_uri()
    if icon:
        fab_face = f"background:#fff url('{icon}') center/54% no-repeat; color:transparent; font-size:0;"
    else:  # graceful fallback to a glyph on a white bubble
        fab_face = "background:#fff; color:#111; font-size:1.5rem;"
    fab_css = (
        "/* ---- floating action button (ChatGPT-logo bubble) ---- */\n"
        ".st-key-gf_fab { position:fixed; right:1.9rem; bottom:1.9rem; z-index:1001; width:3.6rem; }\n"
        ".st-key-gf_fab button { width:3.6rem; height:3.6rem; border-radius:50%; padding:0; line-height:1;"
        " border:1px solid #e3e8ee; box-shadow:0 6px 18px rgba(30,50,80,.22); " + fab_face +
        " transition:transform .12s ease, background-color .12s ease; }\n"
        ".st-key-gf_fab button:hover { background-color:#f6f8fb; transform:translateY(-2px); }\n"
        ".st-key-gf_fab button:active { transform:translateY(0); }\n"
        # Kill the label text in every state (incl. hover, where Streamlit recolors it).
        ".st-key-gf_fab button * { font-size:0 !important; color:transparent !important; }\n"
    )
    return "<style>\n" + fab_css + _PANEL_CSS + "\n</style>"


def _handle_pending(context: str) -> None:
    """If the last turn is an unanswered user message, stream the assistant reply."""
    hist = st.session_state.gf_chat_history
    if not hist or hist[-1]["role"] != "user":
        return
    with st.chat_message("assistant", avatar=":material/shield:"):
        try:
            reply = st.write_stream(_stream_completion(context, hist))
        except Exception as e:  # missing key, network, quota — fail gracefully
            reply = (
                f"I couldn't reach the assistant service: {e}\n\n"
                "Check that `OPENAI_API_KEY` is set in `.streamlit/secrets.toml`."
            )
            st.warning(reply)
    hist.append({"role": "assistant", "content": reply})


def render_floating_assistant(recs: list[dict], spec: dict, gid: str) -> None:
    """Mount the floating button + expandable chat panel. Call once, near the end
    of the page render."""
    ss = st.session_state
    ss.setdefault("gf_chat_open", False)
    ss.setdefault("gf_chat_history", [])
    ss.setdefault("gf_pending", None)

    context = build_context(recs, spec, gid)
    st.markdown(_fixed_css(), unsafe_allow_html=True)

    # A quick-prompt chip / input queued a question last run — make it a user turn.
    if ss.gf_pending:
        ss.gf_chat_history.append({"role": "user", "content": ss.gf_pending})
        ss.gf_pending = None

    # ---- expanded panel (only when open); CSS pins it via .st-key-gf_panel ----
    if ss.gf_chat_open:
        with st.container(key="gf_panel"):
            top = st.columns([0.82, 0.18], vertical_alignment="center")
            top[0].markdown(
                '<div class="gf-assistant-hd">'
                '<svg viewBox="0 0 24 24" fill="none" xmlns="http://www.w3.org/2000/svg" '
                'style="width:1.05rem;height:1.05rem;vertical-align:-2px;margin-right:.35rem">'
                '<path d="M12 2 4 5v6c0 5 3.4 8.5 8 11 4.6-2.5 8-6 8-11V5l-8-3Z" fill="#e7f0fb" '
                'stroke="#0052a3" stroke-width="1.5" stroke-linejoin="round"/>'
                '<path d="M8 12h8M12 8.5v7" stroke="#0052a3" stroke-width="1.6" stroke-linecap="round"/>'
                '</svg>Report assistant</div>'
                '<div class="gf-assistant-sub">Explains <b>this</b> report only · '
                'research prototype, not clinical advice.</div>',
                unsafe_allow_html=True,
            )
            with top[1]:
                with st.container(key="gf_closewrap"):
                    if st.button("✕", key="gf_close", help="Close"):
                        ss.gf_chat_open = False
                        st.rerun()

            # transcript — inside a native fixed-height scroll area so the chips and
            # input rendered after it always stay visible, no matter how long the
            # conversation gets. Only shown once there is something to show.
            if ss.gf_chat_history:
                with st.container(key="gf_scroll", height=280):
                    for msg in ss.gf_chat_history:
                        avatar = ":material/shield:" if msg["role"] == "assistant" else ":material/person:"
                        with st.chat_message(msg["role"], avatar=avatar):
                            st.markdown(msg["content"])
                    _handle_pending(context)  # stream the pending answer into the scroll area
            else:
                _handle_pending(context)

            # dynamic quick-prompt chips above the input. Before any conversation
            # they show as full, left-aligned cards. Once chatting has started they
            # collapse behind a one-line toggle whose open/closed state we own — so
            # sending any message (chip OR typed) always snaps it shut again. (A
            # native st.expander can't do this: `expanded=` is only the initial
            # value and the user's manual toggle sticks across reruns.)
            def _render_chips() -> None:
                with st.container(key="gf_chips"):
                    for i, q in enumerate(suggested_prompts(recs)):
                        if st.button(q, key=f"gf_chip_{i}", use_container_width=True):
                            ss.gf_pending = q
                            ss.gf_chips_open = False   # collapse on send
                            st.rerun()

            if ss.gf_chat_history:
                ss.setdefault("gf_chips_open", False)
                with st.container(key="gf_qtoggle"):
                    caret = "▾" if ss.gf_chips_open else "▸"
                    if st.button(f"{caret}  Quick questions", key="gf_qbtn",
                                 use_container_width=True):
                        ss.gf_chips_open = not ss.gf_chips_open
                        st.rerun()
                if ss.gf_chips_open:
                    _render_chips()
            else:
                st.markdown('<div class="gf-ql">Quick questions</div>', unsafe_allow_html=True)
                _render_chips()

            # free-text input — always visible: the transcript above scrolls inside
            # its own fixed-height box, so this never scrolls out of view.
            if prompt := st.chat_input("Ask about this report…", key="gf_input"):
                ss.gf_pending = prompt
                ss.gf_chips_open = False   # collapse on send
                st.rerun()

    # ---- floating toggle button (always present); CSS pins .st-key-gf_fab ----
    # Label is a zero-width space + hidden by CSS; the ChatGPT-logo is the button face.
    with st.container(key="gf_fab"):
        if st.button("​", key="gf_toggle", help="Ask about this report"):
            ss.gf_chat_open = not ss.gf_chat_open
            st.rerun()
