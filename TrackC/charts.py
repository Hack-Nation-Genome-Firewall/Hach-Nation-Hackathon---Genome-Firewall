"""
Interactive Plotly charts for the Genome Firewall report (Track C).

Palette + marks follow the project data-viz method: validated categorical slots
1-4 (blue/green/magenta/yellow), single-hue blue for the single-series
calibration facets, hairline grid, recessive axes, hover tooltips on every mark.
Rendered on the app's light clinical-blue surface.
"""
from __future__ import annotations

import numpy as np
import pandas as pd
import plotly.graph_objects as go
from plotly.subplots import make_subplots
from sklearn.calibration import calibration_curve
from sklearn.metrics import brier_score_loss

# --- cohesive clinical-blue system (validated monochrome ramp; ties to the app theme) ---
BLUE = "#2a78d6"
# monochrome blue ordinal ramp (deep -> light) — passes ordinal gate on white surface;
# reads as one intentional system instead of a rainbow of unrelated hues.
CAT = ["#103a6b", "#1c5cab", "#2a78d6", "#5598e7"]
GRID = "#e6eaef"
AXIS = "#9aa6b2"
INK = "#12233a"     # titles — near-black slate
TEXT = "#33404d"    # ticks / legend / body — dark, high-contrast
FONT = "IBM Plex Sans, system-ui, -apple-system, sans-serif"


def _style(fig: go.Figure, height: int, title: str | None = None) -> go.Figure:
    fig.update_layout(
        font=dict(family=FONT, color=TEXT, size=13),
        title=dict(text=title, font=dict(color=INK, size=16, weight=600), x=0, xanchor="left") if title else None,
        paper_bgcolor="rgba(0,0,0,0)",
        plot_bgcolor="#ffffff",
        margin=dict(l=58, r=18, t=56 if title else 42, b=48),
        height=height,
        hoverlabel=dict(font_family=FONT, font_size=13, font_color="#12233a", bgcolor="#ffffff", bordercolor=AXIS),
        legend=dict(orientation="h", yanchor="bottom", y=-0.26, x=0,
                    font=dict(size=12, color=TEXT)),
    )
    fig.update_xaxes(showgrid=False, zeroline=False, linecolor=AXIS, ticks="outside",
                     tickcolor=AXIS, tickfont=dict(color=TEXT, size=12.5))
    fig.update_yaxes(showgrid=True, gridcolor=GRID, zeroline=False, linecolor=AXIS,
                     tickfont=dict(color=TEXT, size=12.5))
    return fig


def reliability_figure(pred_df: pd.DataFrame, drugs: list[str]) -> go.Figure:
    """Small-multiples calibration curve (predicted vs observed) per antibiotic."""
    titles, data = [], {}
    for d in drugs:
        sub = pred_df[pred_df["drug"] == d]
        y = sub["y_true"].to_numpy(dtype=int)
        p = sub["p_fail"].to_numpy(dtype=float)
        brier = brier_score_loss(y, p) if len(y) else float("nan")
        titles.append(f"<b>{d}</b><br><span style='font-size:11px;color:{TEXT}'>Brier {brier:.2f} · n={len(y)}</span>")
        data[d] = (y, p)

    fig = make_subplots(rows=1, cols=len(drugs), subplot_titles=titles,
                        shared_yaxes=True, horizontal_spacing=0.035)
    fig.update_annotations(font=dict(color=INK, size=13))
    for i, d in enumerate(drugs, start=1):
        y, p = data[d]
        fig.add_trace(go.Scatter(x=[0, 1], y=[0, 1], mode="lines",
                                 line=dict(color=AXIS, dash="dot", width=1.5),
                                 hoverinfo="skip", showlegend=False), row=1, col=i)
        if len(np.unique(y)) == 2:
            n_bins = min(8, max(3, len(y) // 4))
            try:
                frac, mean_pred = calibration_curve(y, p, n_bins=n_bins, strategy="quantile")
            except ValueError:
                frac, mean_pred = calibration_curve(y, p, n_bins=n_bins, strategy="uniform")
            fig.add_trace(go.Scatter(
                x=mean_pred, y=frac, mode="lines+markers",
                line=dict(color=BLUE, width=2.5),
                marker=dict(size=10, color=BLUE, line=dict(width=2, color="#ffffff")),
                hovertemplate="predicted P(fail) %{x:.0%}<br>observed resistant %{y:.0%}<extra></extra>",
                showlegend=False), row=1, col=i)
        else:
            fig.add_annotation(text="single-class<br>held-out set", showarrow=False,
                               x=0.5, y=0.5, font=dict(color=TEXT, size=11), row=1, col=i)
        fig.update_xaxes(range=[-0.02, 1.02], row=1, col=i,
                         title_text="Predicted P(fail)", title_font=dict(size=11))
        fig.update_yaxes(range=[-0.02, 1.02], row=1, col=i)
    fig.update_yaxes(title_text="Observed fraction resistant", title_font=dict(size=11), row=1, col=1)
    return _style(fig, 330)


def performance_figure(overall_df: pd.DataFrame) -> go.Figure:
    """Grouped bars: key rubric metrics per antibiotic on the held-out split."""
    metrics = [("auroc", "AUROC"), ("balanced_accuracy", "Balanced acc"),
               ("recall_resistant", "Recall (R)"), ("recall_susceptible", "Recall (S)")]
    fig = go.Figure()
    for (col, label), color in zip(metrics, CAT):
        if col not in overall_df:
            continue
        fig.add_trace(go.Bar(
            name=label, x=overall_df["drug"], y=overall_df[col],
            marker=dict(color=color, line=dict(color="#ffffff", width=1), cornerradius=4),
            hovertemplate=f"<b>{label}</b><br>%{{x}}: %{{y:.2f}}<extra></extra>"))
    fig.update_layout(barmode="group", bargap=0.30, bargroupgap=0.06)
    fig.update_yaxes(range=[0, 1.04], tickformat=".0%")
    return _style(fig, 360, "Held-out performance by antibiotic")
