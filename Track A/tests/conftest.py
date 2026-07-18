"""Put the Track A folder on sys.path so tests can import its modules directly.

The folder is named "Track A" (with a space), so it isn't importable as a package
by name; we add it to the path the same way the app does at runtime.
"""
import sys
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parents[1]))  # -> "Track A/"
