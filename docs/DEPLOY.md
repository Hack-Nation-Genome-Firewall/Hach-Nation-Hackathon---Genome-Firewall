# Deploy the live app URL (Streamlit Community Cloud — free)

The submission needs a **live project URL**. Streamlit Community Cloud gives you
a public `https://<something>.streamlit.app` link, built straight from this repo.

## One-time setup (~5 minutes)

1. Go to **https://share.streamlit.io** and sign in with the GitHub account that
   can see this repo (`liiandy/Hach-Nation-Hackathon---Genome-Firewall`).
2. Click **Create app → Deploy a public app from GitHub**.
3. Fill in:
   - **Repository:** `liiandy/Hach-Nation-Hackathon---Genome-Firewall`
   - **Branch:** `track-a-real-features`  *(or `main` once the PRs are merged)*
   - **Main file path:** `streamlit_app.py`   ← the root entry point (real mode)
4. Click **Deploy**. First build takes ~3–5 min (installs `requirements.txt`).
5. You get a URL like `https://genome-firewall.streamlit.app` — that's your
   live project URL for the submission.

## Why `streamlit_app.py` (not `TrackC/app.py`)
The root `streamlit_app.py` sets the four `GENOME_FIREWALL_*` env vars to the
**real** trained artifacts and puts the repo packages on the import path, then
runs `TrackC/app.py`. Because the real `feature_spec.json` has no `synthetic`
flag, the app boots in **real-prediction mode** (no synthetic banner). No edits
to `TrackC/app.py`.

## What works / what's limited on the hosted app
- **Works:** the demo tab (real held-out genomes → calibrated verdicts,
  mechanism markers, OOD abstention) and the reliability/metrics views. This is
  the whole demo-video flow.
- **Limited:** the FASTA-upload tab needs the AMRFinderPlus binary, which is not
  installed on Streamlit Cloud. The app detects this (`shutil.which`) and
  defaults the upload tab to the bundled sample annotation, so the wiring is
  still demonstrable. For a *live* FASTA→annotation demo, run locally with
  AMRFinderPlus on PATH (`bash scripts/run_app_real.sh`).

## Version pins (important)
`requirements.txt` pins `scikit-learn==1.9.0`, `numpy==2.4.6`, `pandas==2.3.3`
to match the runtime the bundle was trained under — so `kp_real_grouped.joblib`
unpickles cleanly on the cloud. Do not loosen these without retraining.

## If the build fails
- **Bundle unpickle error** → a pin drifted from the training runtime; check
  `models/kp_real_grouped.joblib`'s `runtime` dict and match it in requirements.
- **ModuleNotFoundError (module1_reader / module2_predictor)** → the entry point
  must be `streamlit_app.py` at the repo root (it fixes `sys.path`), not
  `TrackC/app.py` directly.
