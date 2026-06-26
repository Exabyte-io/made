"""Resolve the shared data/ directory at the repository root."""

from pathlib import Path

# src/py/mat3ra/made/data_helper.py -> parents[4] -> repo root
DATA_DIR = Path(__file__).parents[4] / "data"
