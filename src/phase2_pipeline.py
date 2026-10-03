"""
Mini-NeoScorer: Phase 2 Pipeline
Coordinate Systems (1-based vs 0-based), RNA-seq TPM Filtering,
and Peptide Context Extraction.
"""

import re
from pathlib import Path
from typing import Dict, List, Tuple

from src.data_models import SomaticVariant, PeptideContext
from src.reference_data import REFERENCE_PROTEINS

# Standard 3-letter to 1-letter amino acid code mapping
AA_3_TO_1 = {
    "Ala": "A", "Arg": "R", "Asn": "N", "Asp": "D", "Cys": "C",
    "Gln": "Q", "Glu": "E", "Gly": "G", "His": "H", "Ile": "I",
    "Leu": "L", "Lys": "K", "Met": "M", "Phe": "F", "Pro": "P",
    "Ser": "S", "Thr": "T", "Trp": "W", "Tyr": "Y", "Val": "V",
}


def parse_vcf(vcf_path: Path) -> List[SomaticVariant]:
    """
    Parse a somatic VCF file and extract variants with their INFO tags.
    Notice POS in VCF is 1-based.
    """
    variants = []
    with open(vcf_path, "r", encoding="utf-8") as f:
        for line in f:
            line = line.strip()
            if not line or line.startswith("#"):
                continue

            parts = line.split("\t")
            chrom = parts[0]
            pos_1based = int(parts[1])  # 1-based coordinate!
            ref = parts[3]
            alt = parts[4]
            info_str = parts[7]

            # Parse key-value pairs from INFO string
            info_dict = {}
            for item in info_str.split(";"):
                if "=" in item:
                    k, v = item.split("=", 1)
                    info_dict[k] = v

            variants.append(
                SomaticVariant(
                    chrom=chrom,
                    pos_1based=pos_1based,
                    ref=ref,
                    alt=alt,
                    gene=info_dict.get("GENE", "UNKNOWN"),
                    hgvsp=info_dict.get("HGVSp", ""),
                    hgvsc=info_dict.get("HGVSc", ""),
                    af=float(info_dict.get("AF", 0.0)),
                    depth=int(info_dict.get("DP", 0)),
                    tier=info_dict.get("TIER", "UNKNOWN"),
                )
            )
    return variants


def parse_expression_tpm(tpm_path: Path) -> Dict[str, float]:
    """
    Parse gene expression TSV file into a dictionary of {gene_symbol: TPM}.
    """
    expression_map = {}
    with open(tpm_path, "r", encoding="utf-8") as f:
        header = f.readline()
        for line in f:
            line = line.strip()
            if not line:
                continue
            parts = line.split("\t")
            if len(parts) >= 2:
                gene = parts[0].strip()
                tpm = float(parts[1].strip())
                expression_map[gene] = tpm
    return expression_map


def parse_hgvsp(hgvsp_str: str) -> Tuple[str, int, str]:
    """
    Parse HGVS protein notation (e.g., 'p.Val600Glu' or 'p.Gly12Asp').
    Returns: (wt_1letter, 1_based_pos, mt_1letter)
    """
    # Regex matching 3-letter codes: p.Val600Glu
    match_3 = re.match(r"p\.([A-Z][a-z]{2})(\d+)([A-Z][a-z]{2})", hgvsp_str)
    if match_3:
        wt_3, pos_str, mt_3 = match_3.groups()
        return AA_3_TO_1[wt_3], int(pos_str), AA_3_TO_1[mt_3]

    # Regex matching 1-letter codes: p.V600E
    match_1 = re.match(r"p\.([A-Z])(\d+)([A-Z])", hgvsp_str)
    if match_1:
        wt_1, pos_str, mt_1 = match_1.groups()
        return wt_1, int(pos_str), mt_1

    raise ValueError(f"Unsupported HGVSp format: {hgvsp_str}")


def extract_peptide_context(
    variant: SomaticVariant, tpm: float, flank: int = 10
) -> PeptideContext:
    """
    Converts 1-based amino acid position to 0-based Python indexing and
    extracts flanking wild-type and mutant peptides (default 21-mer window).
    """
    wt_aa, pos_1based, mt_aa = parse_hgvsp(variant.hgvsp)
    protein_seq = REFERENCE_PROTEINS.get(variant.gene)
    if not protein_seq:
        raise KeyError(f"Reference protein sequence not found for gene {variant.gene}")

    # CRITICAL: 1-based biological coordinate to 0-based Python index
    pos_0based = pos_1based - 1

    # Quality check: ensure wild-type residue in reference matches VCF annotation
    actual_wt = protein_seq[pos_0based]
    if actual_wt != wt_aa:
        raise ValueError(
            f"Reference mismatch for {variant.gene} at pos {pos_1based}: "
            f"expected {wt_aa}, found {actual_wt}"
        )

    # Extract upstream and downstream flanking residues
    start_idx = max(0, pos_0based - flank)
    end_idx = min(len(protein_seq), pos_0based + flank + 1)

    upstream = protein_seq[start_idx:pos_0based]
    downstream = protein_seq[pos_0based + 1:end_idx]

    wt_window = upstream + wt_aa + downstream
    mt_window = upstream + mt_aa + downstream

    return PeptideContext(
        gene=variant.gene,
        variant=variant,
        tpm=tpm,
        mutation_aa_pos=pos_1based,
        wt_residue=wt_aa,
        mt_residue=mt_aa,
        wt_21mer=wt_window,
        mt_21mer=mt_window,
    )


def run_phase2(
    vcf_path: Path, tpm_path: Path, tpm_threshold: float = 1.0
) -> Tuple[List[PeptideContext], List[Tuple[SomaticVariant, float, str]]]:
    """
    Run Phase 2 filtering and peptide extraction:
    1. Parse VCF and TPM
    2. Filter out unexpressed mutations (TPM < threshold)
    3. Extract WT and MT peptide context windows
    """
    variants = parse_vcf(vcf_path)
    expression = parse_expression_tpm(tpm_path)

    passed_candidates: List[PeptideContext] = []
    filtered_out: List[Tuple[SomaticVariant, float, str]] = []

    for var in variants:
        tpm = expression.get(var.gene, 0.0)

        # RNA Expression Gate
        if tpm < tpm_threshold:
            reason = f"RNA Expression Gate: TPM={tpm:.2f} < threshold ({tpm_threshold:.2f})"
            filtered_out.append((var, tpm, reason))
            continue

        # Peptide extraction
        ctx = extract_peptide_context(var, tpm, flank=10)
        passed_candidates.append(ctx)

    return passed_candidates, filtered_out
