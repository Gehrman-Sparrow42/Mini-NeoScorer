"""
Mini-NeoScorer: Phase 3 HLA Class I Binding Scoring Engine.
Implements k-mer sliding window extraction and a position weight matrix (PWM)
affinity predictor for HLA-A*02:01.
"""

import math
from dataclasses import dataclass
from typing import Dict, List, Optional, Tuple
from src.data_models import PeptideContext


@dataclass
class EpitopeCandidate:
    gene: str
    variant_hgvsp: str
    hla_allele: str
    kmer_length: int
    mut_pos_in_kmer: int  # 1-based position of the mutation inside the 9-mer
    wt_peptide: str
    mt_peptide: str
    wt_ic50_nm: float
    mt_ic50_nm: float
    wt_binder_class: str  # "Strong Binder", "Weak Binder", "Non-Binder"
    mt_binder_class: str
    agretopicity_index: float  # WT IC50 / MT IC50 (Higher = MT binds more tightly than WT)
    percentile_rank: float = 50.0
    presentation_score: float = 0.0
    mhcflurry_ic50: Optional[float] = None
    netmhcpan_ic50: Optional[float] = None
    predictor_source: str = "PWM_Baseline"


# Position Weight Matrix for HLA-A*02:01 (9-mers)
# Derived from canonical experimental binding preferences:
# - Primary anchor: Position 2 (prefers L, M, V, I; penalizes charged/proline)
# - Primary anchor: Position 9 (prefers V, L, I, A; penalizes charged/bulky)
# - Secondary anchors: Positions 1, 3, 7
HLA_A0201_PWM: Dict[int, Dict[str, float]] = {
    # Position 1
    1: {"A": 1.0, "G": 0.5, "I": 1.0, "L": 1.2, "M": 0.8, "F": 1.5, "W": 1.0, "Y": 1.5,
        "S": 0.8, "T": 0.8, "V": 1.0, "K": 0.5, "R": 0.2, "D": -0.5, "E": -0.5, "P": -1.0},
    # Position 2 (PRIMARY ANCHOR)
    2: {"L": 4.5, "M": 4.0, "V": 2.8, "I": 2.6, "A": 1.2, "T": 0.5, "Q": 0.2,
        "P": -3.5, "D": -4.0, "E": -4.0, "K": -4.0, "R": -4.0, "G": -2.0, "W": -1.5},
    # Position 3
    3: {"F": 2.0, "Y": 2.0, "W": 1.8, "L": 1.5, "I": 1.2, "V": 1.0, "M": 1.0,
        "D": -0.8, "E": -0.8, "P": -1.0, "G": 0.0},
    # Position 4
    4: {"P": 1.0, "E": 0.5, "G": 0.5, "D": 0.2, "K": 0.2, "R": 0.2, "L": 0.0, "A": 0.0},
    # Position 5
    5: {"F": 1.2, "Y": 1.2, "L": 0.8, "I": 0.8, "V": 0.5, "A": 0.5, "G": 0.2, "P": -1.0},
    # Position 6
    6: {"V": 1.2, "I": 1.2, "L": 1.0, "A": 0.8, "T": 0.5, "P": -0.5, "K": -0.5},
    # Position 7
    7: {"L": 1.5, "V": 1.2, "I": 1.2, "F": 1.0, "Y": 1.0, "M": 0.8, "K": 0.0, "E": -0.5},
    # Position 8
    8: {"K": 1.2, "R": 1.0, "T": 0.8, "S": 0.5, "A": 0.5, "L": 0.2, "D": 0.0, "P": -0.8},
    # Position 9 (PRIMARY C-TERMINAL ANCHOR)
    9: {"V": 4.5, "L": 4.2, "I": 3.2, "A": 2.0, "T": 1.0, "M": 1.2,
        "K": -4.0, "R": -4.0, "D": -4.5, "E": -4.5, "P": -4.0, "G": -2.5, "W": -2.0},
}


def score_peptide_a0201(peptide: str) -> float:
    """
    Calculate raw binding energy score for a 9-mer peptide against HLA-A*02:01.
    Higher raw score = stronger binding affinity.
    """
    if len(peptide) != 9:
        raise ValueError(f"HLA-A*02:01 model currently expects 9-mers, got len {len(peptide)}: {peptide}")

    score = 0.0
    for pos, aa in enumerate(peptide, 1):
        pos_dict = HLA_A0201_PWM.get(pos, {})
        # Default neutral background score for unlisted amino acid at this position
        score += pos_dict.get(aa, 0.0)
    return score


def raw_score_to_ic50(raw_score: str) -> float:
    """
    Convert raw matrix score to estimated IC50 (nM).
    In biological binding:
    - High positive score -> Low IC50 (< 50 nM, Strong Binder)
    - Medium score -> 50 - 500 nM (Weak Binder)
    - Low/Negative score -> High IC50 (> 500 nM, Non-Binder)
    Calibration curve: IC50 = 50000 / (1 + e^(0.8 * score - 1.5))
    """
    # Sigmoidal mapping scaled to standard 50,000 nM ceiling
    k = 0.85
    offset = 1.8
    # Clamp exponent to prevent numerical overflow
    exp_val = max(-15.0, min(15.0, -(k * raw_score - offset)))
    ic50 = 50000.0 / (1.0 + math.exp(-exp_val))
    return round(ic50, 2)


def classify_ic50(ic50_nm: float) -> str:
    """
    Standard immuno-oncology affinity classification thresholds.
    """
    if ic50_nm <= 50.0:
        return "Strong Binder"
    elif ic50_nm <= 500.0:
        return "Weak Binder"
    else:
        return "Non-Binder"


def generate_9mers_for_context(
    ctx: PeptideContext, hla_allele: str = "HLA-A*02:01"
) -> List[EpitopeCandidate]:
    """
    Slide a 9-amino acid window across the 21-mer context window.
    Only keeps 9-mers that actually contain the mutated residue.
    """
    wt_seq = ctx.wt_21mer
    mt_seq = ctx.mt_21mer

    # Locate the mutation site in the 21-mer window
    mut_idx = next(
        (i for i, (c1, c2) in enumerate(zip(wt_seq, mt_seq)) if c1 != c2), None
    )
    if mut_idx is None:
        raise ValueError(f"No sequence difference found between WT and MT for {ctx.gene}")

    candidates: List[EpitopeCandidate] = []
    k = 9

    # Slide window of size k across the 21-mer
    for start in range(0, len(mt_seq) - k + 1):
        end = start + k
        # Check if the mutation index falls inside this [start, end) window
        if start <= mut_idx < end:
            wt_kmer = wt_seq[start:end]
            mt_kmer = mt_seq[start:end]
            mut_pos_in_kmer = (mut_idx - start) + 1  # 1-based position in 9-mer

            # Score both WT and MT against HLA
            wt_raw = score_peptide_a0201(wt_kmer)
            mt_raw = score_peptide_a0201(mt_kmer)

            wt_ic50 = raw_score_to_ic50(wt_raw)
            mt_ic50 = raw_score_to_ic50(mt_raw)

            # Agretopicity Index: WT_IC50 / MT_IC50
            # Higher AI means MT binds more tightly than WT!
            agretopicity = round(wt_ic50 / max(0.01, mt_ic50), 2)

            candidate = EpitopeCandidate(
                gene=ctx.gene,
                variant_hgvsp=ctx.variant.hgvsp,
                hla_allele=hla_allele,
                kmer_length=k,
                mut_pos_in_kmer=mut_pos_in_kmer,
                wt_peptide=wt_kmer,
                mt_peptide=mt_kmer,
                wt_ic50_nm=wt_ic50,
                mt_ic50_nm=mt_ic50,
                wt_binder_class=classify_ic50(wt_ic50),
                mt_binder_class=classify_ic50(mt_ic50),
                agretopicity_index=agretopicity,
            )
            candidates.append(candidate)

    # Sort candidate 9-mers by lowest MT IC50 (strongest binding first)
    candidates.sort(key=lambda c: c.mt_ic50_nm)
    return candidates
