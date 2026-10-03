"""
Reference protein sequence provider for Mini-NeoScorer.
Loads curated reference sequences (UniProt canonical isoforms)
covering target somatic loci.
"""

import json
from pathlib import Path
from typing import Dict

BASE_DIR = Path(__file__).resolve().parent.parent
JSON_PATH = BASE_DIR / "data" / "reference_proteins.json"


def load_reference_proteins() -> Dict[str, str]:
    """
    Load canonical reference protein sequences from cached JSON.
    """
    if JSON_PATH.exists():
        with open(JSON_PATH, "r", encoding="utf-8") as f:
            return json.load(f)
    raise FileNotFoundError(f"Reference protein cache not found at {JSON_PATH}")


REFERENCE_PROTEINS: Dict[str, str] = load_reference_proteins()
