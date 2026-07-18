"""Put the module1_reader folder on sys.path so tests can import its modules directly.

The tests import `feature_annotator` / `build_features` by bare name (not as
`module1_reader.<mod>`), so we add the module folder to the path.
"""
import sys
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parents[1]))  # -> module1_reader/
