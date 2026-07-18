"""Path setup for Track A tests.

- module1_reader/ on sys.path so tests can import `feature_annotator` / `build_features`
  by bare name.
- repo root on sys.path so the compatibility tests can `import module2_predictor.contracts`.
"""
import sys
from pathlib import Path

MODULE_DIR = Path(__file__).resolve().parents[1]   # module1_reader/
REPO_ROOT = Path(__file__).resolve().parents[2]    # repo root
for p in (str(MODULE_DIR), str(REPO_ROOT)):
    if p not in sys.path:
        sys.path.insert(0, p)
