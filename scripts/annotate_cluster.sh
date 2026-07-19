#!/usr/bin/env bash
# =============================================================================
# Genome Firewall — Track A annotation on a cluster (Swift / any SLURM or plain box)
#
# Runs AMRFinderPlus over the 3,000 selected K. pneumoniae genomes, streaming
# each assembly (download -> annotate -> delete) so peak disk stays small.
# Resumable: re-running skips genomes already annotated.
#
# Self-contained: if conda is not on PATH (e.g. a Nix cluster like Swift) it
# bootstraps a private Miniconda into $HOME/miniconda3 automatically.
#
# TWO WAYS TO RUN
#   Plain (interactive node / login shell with many cores):
#       bash scripts/annotate_cluster.sh
#   SLURM batch:
#       sbatch scripts/annotate_cluster.sh
#
# The SBATCH lines below are read only by `sbatch`; `bash` ignores them.
# Partition/cpus are preset for Swift's `newcpu`; change if your cluster differs.
#
# NOTE ON INTERNET: this streams FASTA downloads from BV-BRC while it runs, so
# the executing node needs outbound internet. If your COMPUTE nodes are
# firewalled, run this on an internet-facing login/interactive node (bash …)
# instead of via sbatch, or pre-download the FASTAs first.
# =============================================================================
#SBATCH --job-name=gf-amrfinder
#SBATCH --output=gf-amrfinder-%j.log
#SBATCH --cpus-per-task=32
#SBATCH --mem=32G
#SBATCH --time=12:00:00
#SBATCH --partition=newcpu

set -euo pipefail

# ---- config (override via env vars) ----------------------------------------
REPO_DIR="${REPO_DIR:-$PWD}"                    # run from the repo root
IDS="${IDS:-$REPO_DIR/data/manifests/selected_genomes.csv}"
OUTDIR="${OUTDIR:-$REPO_DIR/cohort_tsv}"
FASTA_DIR="${FASTA_DIR:-$REPO_DIR/fasta_tmp}"
WORKERS="${WORKERS:-${SLURM_CPUS_PER_TASK:-8}}" # parallel genomes
AMR_THREADS="${AMR_THREADS:-2}"                 # BLAST threads per genome
ENV_NAME="${ENV_NAME:-amrfinder}"

echo "[gf] repo=$REPO_DIR  ids=$IDS  workers=$WORKERS x ${AMR_THREADS} threads"

# ---- 1. conda (bootstrap Miniconda if absent) ------------------------------
if ! command -v conda >/dev/null 2>&1; then
  CONDA_HOME="${CONDA_HOME:-$HOME/miniconda3}"
  if [ ! -x "$CONDA_HOME/bin/conda" ]; then
    echo "[gf] conda not found — bootstrapping Miniconda into $CONDA_HOME…"
    OS="$(uname -s)"; ARCH="$(uname -m)"
    case "$OS-$ARCH" in
      Linux-x86_64)  MC="Miniconda3-latest-Linux-x86_64.sh" ;;
      Linux-aarch64) MC="Miniconda3-latest-Linux-aarch64.sh" ;;
      *) echo "[gf] ERROR: unsupported platform $OS-$ARCH; install conda manually." >&2; exit 1 ;;
    esac
    curl -fsSL "https://repo.anaconda.com/miniconda/$MC" -o /tmp/miniconda_$$.sh
    bash /tmp/miniconda_$$.sh -b -p "$CONDA_HOME"
    rm -f /tmp/miniconda_$$.sh
  fi
  export PATH="$CONDA_HOME/bin:$PATH"
fi
source "$(conda info --base)/etc/profile.d/conda.sh"

# ---- 2. AMRFinderPlus env --------------------------------------------------
if ! conda env list | grep -qw "$ENV_NAME"; then
  echo "[gf] creating conda env '$ENV_NAME' (ncbi-amrfinderplus)…"
  conda create -y -n "$ENV_NAME" -c conda-forge -c bioconda ncbi-amrfinderplus pandas
fi
conda activate "$ENV_NAME"
python -m pip install --quiet pandas 2>/dev/null || true
amrfinder --version

# ---- 2. AMRFinderPlus database ---------------------------------------------
DB_DIR="${DB_DIR:-$REPO_DIR/amrfinder_db}"
if [ ! -f "$DB_DIR"/*/version.txt ] 2>/dev/null && [ -z "$(ls "$DB_DIR" 2>/dev/null)" ]; then
  echo "[gf] downloading AMRFinderPlus database to $DB_DIR…"
  mkdir -p "$DB_DIR"
  amrfinder_update -d "$DB_DIR"
fi
DB_PATH="$(ls -d "$DB_DIR"/*/ 2>/dev/null | head -1)"
DB_PATH="${DB_PATH:-$DB_DIR}"
echo "[gf] using DB: $DB_PATH"

# ---- 3. run the cohort (resumable, streaming) ------------------------------
mkdir -p "$OUTDIR" "$FASTA_DIR"
cd "$REPO_DIR"
python module1_reader/run_cohort_parallel.py \
    --ids "$IDS" \
    --db "$DB_PATH" \
    --amrfinder "$(command -v amrfinder)" \
    --outdir "$OUTDIR" \
    --fasta-dir "$FASTA_DIR" \
    --workers "$WORKERS" \
    --amr-threads "$AMR_THREADS"

# ---- 4. assemble the feature matrix + spec ---------------------------------
echo "[gf] assembling feature matrix…"
python -m module1_reader.assemble_features \
    --tsv-dir "$OUTDIR" \
    --selected "$IDS" \
    --out-features data/manifests/features.csv \
    --out-spec     data/manifests/feature_spec.json

N_TSV=$(ls "$OUTDIR"/*.tsv 2>/dev/null | wc -l | tr -d ' ')
echo "[gf] DONE. annotated $N_TSV genomes."
echo "[gf] outputs: data/manifests/features.csv, data/manifests/feature_spec.json"
echo "[gf] next (locally or here): python -m module2_predictor.train --features data/manifests/features.csv \\"
echo "         --labels data/manifests/labels.csv --splits data/manifests/split_manifest.csv \\"
echo "         --spec data/manifests/feature_spec.json --output models/kp_bundle.joblib --model-version kp-real-v1"
