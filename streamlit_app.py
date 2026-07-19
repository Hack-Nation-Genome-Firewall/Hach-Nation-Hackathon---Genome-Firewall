"""Root entry point for Streamlit Community Cloud.

Deploy target: point the app at THIS file (streamlit_app.py). It sets the four
GENOME_FIREWALL_* env vars to the REAL trained artifacts, puts the repo on the
import path, and hands off to TrackC/app.py — so the deployed app runs on the
real 3,000-genome K. pneumoniae model (no synthetic banner), with no edits to
TrackC/app.py.
"""
from __future__ import annotations
import os
import runpy
from pathlib import Path

ROOT = Path(__file__).resolve().parent

os.environ.setdefault("GENOME_FIREWALL_SPEC", str(ROOT / "data/manifests/feature_spec.json"))
os.environ.setdefault("GENOME_FIREWALL_FEATURES", str(ROOT / "data/manifests/features.csv"))
os.environ.setdefault("GENOME_FIREWALL_SPLITS", str(ROOT / "data/manifests/split_manifest_aligned.csv"))
os.environ.setdefault("GENOME_FIREWALL_BUNDLE", str(ROOT / "models/kp_real_grouped.joblib"))
# Point the held-out performance/calibration panel at the REAL grouped eval outputs
# (the local `eval/` dir is empty in real mode, which is why the interactive chart +
# metrics table were falling back to the "run evaluate" placeholder).
os.environ.setdefault("GENOME_FIREWALL_EVAL", str(ROOT / "results/eval_grouped"))

# Make repo packages (module1_reader, module2_predictor, TrackC) importable.
import sys
for p in (str(ROOT), str(ROOT / "TrackC")):
    if p not in sys.path:
        sys.path.insert(0, p)

# Execute the actual Streamlit app in-process.
runpy.run_path(str(ROOT / "TrackC" / "app.py"), run_name="__main__")
