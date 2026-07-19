#!/usr/bin/env bash
# =============================================================================
# GPU A/B stretch experiment: DNABERT-2 embeddings vs. AMRFinderPlus features
# on ciprofloxacin. Fire-and-forget; safe to run overnight.
#
#   sbatch scripts/gpu_ab.sbatch
#
# Targets an L40S (node9). To use a different GPU, change --gres below:
#   L40S : --gres=gpu:nvidia_l40s:1        (node9,  48GB)
#   L4   : --gres=gpu:nvidia_l4:1          (node10-12, 24GB)
#   A6000: --gres=gpu:nvidia_rtx_a6000:1   (node5,  48GB, partition oldcpu)
#   V100 : --gres=gpu:tesla_v100-pcie-16gb:1 (node7, 16GB, partition oldcpu)
# =============================================================================
#SBATCH --job-name=gf-gpu-ab
#SBATCH --output=gf-gpu-ab-%j.log
#SBATCH --partition=newcpu
#SBATCH --gres=gpu:nvidia_l40s:1
#SBATCH --cpus-per-task=4
#SBATCH --mem=32G
#SBATCH --time=08:00:00

set -euo pipefail
cd "$SLURM_SUBMIT_DIR"

# ---- conda + a torch/transformers env (reuses the bootstrapped Miniconda) ----
if ! command -v conda >/dev/null 2>&1; then
  export PATH="$HOME/miniconda3/bin:$PATH"
fi
source "$(conda info --base)/etc/profile.d/conda.sh"

ENV_NAME="gf-dnabert"
if ! conda env list | grep -qw "$ENV_NAME"; then
  echo "[ab] creating env $ENV_NAME (python 3.11)…"
  conda create -y -n "$ENV_NAME" --override-channels -c conda-forge python=3.11
fi
conda activate "$ENV_NAME"

# torch with CUDA + the deps DNABERT-2 needs. cu121 wheels work on L40S/L4/A6000/V100.
python - <<'PY' 2>/dev/null || NEED_INSTALL=1
import torch, transformers, sklearn, einops  # noqa
PY
if [ "${NEED_INSTALL:-0}" = "1" ]; then
  echo "[ab] installing torch + transformers…"
  pip install --quiet torch --index-url https://download.pytorch.org/whl/cu121
  pip install --quiet "transformers==4.44.2" einops sklearn scikit-learn pandas numpy
fi

echo "[ab] GPU check:"; python -c "import torch; print('cuda', torch.cuda.is_available(), torch.cuda.get_device_name(0) if torch.cuda.is_available() else '-')"

# ---- run the A/B --------------------------------------------------------------
python scripts/gpu_ab_dnabert.py \
    --labels data/manifests/labels.csv \
    --splits data/manifests/split_manifest.csv \
    --drug ciprofloxacin \
    --out results/gpu_ab_ciprofloxacin.txt

echo "[ab] DONE. see results/gpu_ab_ciprofloxacin.txt"
