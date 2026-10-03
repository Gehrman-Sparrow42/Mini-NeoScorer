"""
Data models for Mini-NeoScorer.
"""

from dataclasses import dataclass
from typing import Optional


@dataclass
class SomaticVariant:
    chrom: str
    pos_1based: int
    ref: str
    alt: str
    gene: str
    hgvsp: str
    hgvsc: str
    af: float
    depth: int
    tier: str

    @property
    def pos_0based(self) -> int:
        """
        Convert 1-based VCF coordinate to 0-based indexing (e.g. for Python slices or BED).
        """
        return self.pos_1based - 1


@dataclass
class PeptideContext:
    gene: str
    variant: SomaticVariant
    tpm: float
    mutation_aa_pos: int  # 1-based amino acid position
    wt_residue: str
    mt_residue: str
    wt_21mer: str
    mt_21mer: str
