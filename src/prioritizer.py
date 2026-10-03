"""
Mini-NeoScorer: Phase 4 Prioritizer & Scoring Engine.
Calculates the unified Immunogenicity Prioritization Score (S_neo)
integrating clonality (AF), transcriptional abundance (TPM),
MHC presentation probability (IC50), and agretopicity (AI).
"""

import math
from dataclasses import dataclass
from typing import List, Tuple
from src.data_models import PeptideContext
from src.hla_scorer import EpitopeCandidate


@dataclass
class ScoredNeoantigen:
    gene: str
    variant_hgvsp: str
    tier: str
    hla_allele: str
    mt_peptide: str
    wt_peptide: str
    mut_pos_in_kmer: int
    af: float
    clonality_score: float      # C_AF (0.0 to 1.0)
    tpm: float
    expression_score: float     # E_TPM (log2 scale)
    mt_ic50_nm: float
    presentation_score: float   # P_MHC (0.0 to 1.0)
    agretopicity_index: float   # AI
    agretopicity_bonus: float   # A_bonus
    composite_score: float      # Final score scaled 0 - 100
    clinical_tier: str          # Tier 1 (Actionable), Tier 2 (Secondary), Tier 3 (Borderline/Deprioritized)
    recommendation: str


def compute_clonality_score(af: float) -> float:
    """
    In diploid cancer genomes, a heterozygous somatic mutation present in 100%
    of tumor cells yields an AF of ~0.50 (50%).
    Clonality score evaluates the fraction of tumor cells carrying the variant:
    min(1.0, 2.0 * AF).
    """
    return min(1.0, max(0.0, 2.0 * af))


def compute_expression_score(tpm: float) -> float:
    """
    Translates raw TPM to log2 abundance scale: log2(1 + TPM).
    Captures transcriptional volume while moderating extreme high-expression outliers.
    """
    return math.log2(1.0 + max(0.0, tpm))


def compute_presentation_score(mt_ic50_nm: float) -> float:
    """
    Maps IC50 in nM to normalized presentation probability (0.0 to 1.0):
    P_MHC = 1.0 - log10(min(50000, max(1, IC50))) / log10(50000)
    - IC50 <= 50 nM yields P_MHC >= 0.64
    - IC50 ~ 500 nM yields P_MHC ~ 0.43
    - IC50 > 10,000 nM yields P_MHC < 0.15
    """
    clamped_ic50 = min(50000.0, max(1.0, mt_ic50_nm))
    score = 1.0 - (math.log10(clamped_ic50) / math.log10(50000.0))
    return max(0.0, min(1.0, score))


def compute_agretopicity_bonus(ai: float) -> float:
    """
    Calculates bonus for gain-of-affinity mutations (AI > 1.0).
    A mutation that binds significantly more tightly than WT receives a modest boost.
    """
    if ai >= 1.0:
        return min(1.5, 1.0 + 0.15 * math.log2(ai + 1.0))
    else:
        # Slight penalty if mutation damaged binding compared to WT
        return max(0.75, 0.90 + 0.10 * ai)


def score_and_rank_candidates(
    epitope_list: List[Tuple[PeptideContext, EpitopeCandidate]]
) -> List[ScoredNeoantigen]:
    """
    Applies the Mini-NeoScorer scoring equation:
    S_neo = C_AF * E_TPM * P_MHC * A_bonus * scale_factor
    """
    scored = []
    scale_factor = 20.0  # Scales standard top scores into 0 - 100 range

    for ctx, epi in epitope_list:
        v = ctx.variant
        c_af = compute_clonality_score(v.af)
        e_tpm = compute_expression_score(ctx.tpm)
        p_mhc = compute_presentation_score(epi.mt_ic50_nm)
        a_bonus = compute_agretopicity_bonus(epi.agretopicity_index)

        raw_score = c_af * e_tpm * p_mhc * a_bonus
        final_score = round(min(100.0, raw_score * scale_factor), 1)

        # Assign clinical actionability tier
        if final_score >= 50.0 and epi.mt_ic50_nm <= 50.0:
            clinical_tier = "Tier 1: High Priority (Vaccine Candidate)"
            recommendation = "SELECT: High clonality, robust expression, strong HLA-A*02:01 presentation."
        elif final_score >= 25.0 or epi.mt_ic50_nm <= 500.0:
            clinical_tier = "Tier 2: Secondary / Backup"
            recommendation = "HOLD: Moderate presentation or limited expression. Secondary pool."
        else:
            clinical_tier = "Tier 3: Deprioritized"
            if epi.mt_ic50_nm > 500.0:
                recommendation = "REJECT: Failed HLA-A*02:01 presentation (IC50 > 500 nM)."
            elif ctx.tpm < 5.0:
                recommendation = "REJECT: Insufficient transcript abundance (TPM < 5.0)."
            else:
                recommendation = "REJECT: Subclonal variant with low tumor coverage."

        scored.append(
            ScoredNeoantigen(
                gene=ctx.gene,
                variant_hgvsp=v.hgvsp,
                tier=v.tier,
                hla_allele=epi.hla_allele,
                mt_peptide=epi.mt_peptide,
                wt_peptide=epi.wt_peptide,
                mut_pos_in_kmer=epi.mut_pos_in_kmer,
                af=v.af,
                clonality_score=round(c_af, 2),
                tpm=ctx.tpm,
                expression_score=round(e_tpm, 2),
                mt_ic50_nm=epi.mt_ic50_nm,
                presentation_score=round(p_mhc, 2),
                agretopicity_index=epi.agretopicity_index,
                agretopicity_bonus=round(a_bonus, 2),
                composite_score=final_score,
                clinical_tier=clinical_tier,
                recommendation=recommendation,
            )
        )

    # Sort descending by final composite score
    scored.sort(key=lambda s: s.composite_score, reverse=True)
    return scored
