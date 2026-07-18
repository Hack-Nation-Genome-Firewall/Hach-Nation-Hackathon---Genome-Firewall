"""Compatibility entry point for the shared Track B evaluation harness."""

from pathlib import Path
import sys


ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT))

from module2_predictor.evaluate import main  # noqa: E402


if __name__ == "__main__":
    main()
