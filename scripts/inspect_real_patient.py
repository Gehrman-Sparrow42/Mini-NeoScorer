"""
Mini-NeoScorer: Real Patient Data Profiler
Inspects authentic clinical patient data: TCGA-D3-A1Q1 (Melanoma).
"""

import json
import sys
from pathlib import Path

# Ensure UTF-8 output
if hasattr(sys.stdout, "reconfigure"):
    sys.stdout.reconfigure(encoding="utf-8")

BASE_DIR = Path(__file__).resolve().parent.parent
DATA_DIR = BASE_DIR / "data"

MUT_FILE = DATA_DIR / "real_patient_TCGA-D3-A1Q1-06_mutations.json"
EXP_FILE = DATA_DIR / "real_patient_TCGA-D3-A1Q1-06_expression.json"


def main():
    with open(MUT_FILE, "r", encoding="utf-8") as f:
        mutations = json.load(f)

    with open(EXP_FILE, "r", encoding="utf-8") as f:
        expression = json.load(f)

    missense = [m for m in mutations if m.get("mutationType") == "Missense_Mutation"]

    print("=" * 95)
    print(" REAL CLINICAL PATIENT PROFILE: TCGA-D3-A1Q1 (Metastatic Skin Cutaneous Melanoma)")
    print(" Data Source: TCGA PanCancer Atlas (Cell, 2018) via cBioPortal API")
    print("=" * 95)
    print(f"Total Somatic Mutations Detected (WES) : {len(mutations)}")
    print(f"Missense Single-Nucleotide Variants   : {len(missense)}")
    print(f"Genes with Matched RNA-seq Expression : {len(expression)}")
    print("-" * 95)
    print(f"{'#':<4} {'Gene':<10} {'Protein Change':<16} {'Genomic Locus':<22} {'AF':<6} {'Reads (Alt/Total)':<18} {'Matched RNA (TPM)'}")
    print("-" * 95)

    for i, m in enumerate(missense[:12], 1):
        gene = m.get("gene", {}).get("hugoGeneSymbol", "UNKNOWN")
        pc = m.get("proteinChange", "UNKNOWN")
        chrom = m.get("chr", "")
        pos = m.get("startPosition", 0)
        ref = m.get("referenceAllele", "")
        alt = m.get("variantAllele", "")
        alt_reads = m.get("tumorAltCount", 0)
        ref_reads = m.get("tumorRefCount", 0)
        tot_reads = alt_reads + ref_reads
        af = alt_reads / max(1, tot_reads)
        entrez = m.get("entrezGeneId")
        tpm = expression.get(str(entrez), 0.0)

        print(
            f"{i:<4} "
            f"{gene:<10} "
            f"{pc:<16} "
            f"chr{chrom}:{pos} ({ref}>{alt}) "
            f"{af:<6.2f} "
            f"{alt_reads}/{tot_reads:<15} "
            f"{tpm:<10.1f}"
        )

    print("-" * 95)
    print(f"... and {len(missense) - 12} additional authentic patient missense mutations.")
    print("=" * 95)


if __name__ == "__main__":
    main()
