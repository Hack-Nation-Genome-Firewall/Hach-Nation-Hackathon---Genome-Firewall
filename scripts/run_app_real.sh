#!/usr/bin/env bash
# =============================================================================
# Launch the Streamlit decision app on the REAL trained bundle (not synthetic).
#
#   bash scripts/run_app_real.sh
#
# Points the app's four env-var paths at the real 3,000-genome K. pneumoniae
# artifacts. No edits to TrackC/app.py — the app already reads these vars and
# falls back to the synthetic fixture when they are unset.
# =============================================================================
set -euo pipefail
cd "$(dirname "$0")/.."

export GENOME_FIREWALL_SPEC="$PWD/data/manifests/feature_spec.json"
export GENOME_FIREWALL_FEATURES="$PWD/data/manifests/features.csv"
export GENOME_FIREWALL_SPLITS="$PWD/data/manifests/split_manifest_aligned.csv"
export GENOME_FIREWALL_BUNDLE="$PWD/models/kp_real_grouped.joblib"

echo "[app] REAL mode:"
echo "  spec     = $GENOME_FIREWALL_SPEC"
echo "  features = $GENOME_FIREWALL_FEATURES  (2997 genomes)"
echo "  bundle   = $GENOME_FIREWALL_BUNDLE   (kp-real-grouped-v1)"
echo

export PYTHONPATH="$PWD:${PYTHONPATH:-}"
exec streamlit run TrackC/app.py
