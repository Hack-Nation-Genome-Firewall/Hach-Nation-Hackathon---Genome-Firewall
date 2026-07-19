# Running the decision app on the REAL trained model

The Streamlit app (`TrackC/app.py`) defaults to the synthetic fixture (with a
prominent disclosure banner). To run it on the **real** 3,000-genome
K. pneumoniae bundle — the one a clinician would use — launch it with:

```bash
bash scripts/run_app_real.sh
```

That script sets the app's four env-var paths (no edits to `app.py`) to:

| Env var | Real file |
|---|---|
| `GENOME_FIREWALL_SPEC` | `data/manifests/feature_spec.json` (status `real_pending_domain_review`) |
| `GENOME_FIREWALL_FEATURES` | `data/manifests/features.csv` (2,997 genomes) |
| `GENOME_FIREWALL_SPLITS` | `data/manifests/split_manifest_aligned.csv` |
| `GENOME_FIREWALL_BUNDLE` | `models/kp_real_grouped.joblib` (`kp-real-grouped-v1`) |

Because the real `feature_spec.json` has **no `synthetic` flag**, the app's
`IS_SYNTHETIC` evaluates to `False` and the synthetic disclosure banner does not
render — the app is in genuine real-prediction mode.

## What the demo shows in real mode
- **Demo tab:** pick a real held-out test genome → per-drug verdicts
  (`likely_to_work` / `likely_to_fail` / `no_call`) with calibrated confidence,
  mechanism markers (e.g. blaKPC, blaNDM, gyrA_S83), and the OOD abstention.
- **Upload tab:** drop a FASTA → Track A's reader (`module1_reader`) runs
  AMRFinderPlus (if on PATH) → same verdicts. On a machine without AMRFinderPlus
  installed, the upload path falls back to the bundled sample so the wiring is
  demonstrable.

## Bundle provenance
`models/kp_real_grouped.joblib` = `results/models/kp_grouped.joblib`, trained by
`module2_predictor.train` on the homology-grouped split (2,098 train / 360
calibration / 539 test genomes), isotonic/sigmoid-calibrated. Held-out metrics
in `results/pitch_metrics.csv`.

> Decision support only. Every prediction must be confirmed by standard
> laboratory antimicrobial susceptibility testing before any treatment decision.
