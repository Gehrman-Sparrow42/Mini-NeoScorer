"""
Mini-NeoScorer: Scale Step 3 Execution & Verification Script
Runs high-throughput 9-mer sliding window extraction and HLA-A*02:01 affinity
screening across the entire exome of patient TCGA-D3-A1Q1 (Melanoma).
"""

import json
import sys
from pathlib import Path

# Ensure UTF-8 output
if hasattr(sys.stdout, "reconfigure"):
    sys.stdout.reconfigure(encoding="utf-8")

BASE_DIR = Path(__file__).resolve().parent.parent
DATA_DIR = BASE_DIR / "data"

if str(BASE_DIR) not in sys.path:
    sys.path.insert(0, str(BASE_DIR))

MUT_FILE = DATA_DIR / "real_patient_TCGA-D3-A1Q1-06_mutations.json"
EXP_FILE = DATA_DIR / "real_patient_TCGA-D3-A1Q1-06_expression.json"

from src.scale_hla_screener import screen_patient_mutations


def main():
    print("=" * 105)
    print(" MINI-NEOSCORER SCALE STEP 3: EXOME-WIDE 9-MER SLIDING WINDOW & HLA-A*02:01 SCREENING")
    print(" Patient: TCGA-D3-A1Q1 (Metastatic Melanoma) | Target Allele: HLA-A*02:01")
    print("=" * 105)

    with open(MUT_FILE, "r", encoding="utf-8") as f:
        mutations = json.load(f)

    with open(EXP_FILE, "r", encoding="utf-8") as f:
        expression = json.load(f)

    # Execute high-throughput exome screen
    results, stats = screen_patient_mutations(
        mutations, expression, hla_allele="HLA-A*02:01", tpm_threshold=1.0
    )

    print("\n>>> HIGH-THROUGHPUT EXOME SCREENING FUNNEL:")
    print(f"  * Total somatic alterations in biopsy      : {stats['total_somatic_mutations']}")
    print(f"  * Missense protein-altering SNVs            : {stats['missense_mutations']}")
    print(f"  * Rejected by RNA Expression Gate (TPM < 1) : {stats['gated_out_low_rna']} (Silent / unexpressed)")
    print(f"  * Rejected by Isoform / Sequence QC checks  : {stats['isoform_mismatch_or_unresolved']}")
    print(f"  * Successfully screened patient variants    : {stats['successfully_screened_variants']}")
    print(f"  * Total candidate 9-mers evaluated on HLA   : {stats['total_9mers_evaluated']}")
    print(f"    |-- Strong Binders (IC50 <= 50 nM)        : {stats['strong_binders_found']}")
    print(f"    |-- Weak Binders (50 - 500 nM)            : {stats['weak_binders_found']}")
    print(f"    \\-- Non-Binders (IC50 > 500 nM)           : {stats['non_binders_found']} (Fell off HLA tray)")

    print("\n" + "=" * 105)
    print(" TOP 10 STRONGEST BINDING REAL NEOPEPTIDES DISCOVERED (HLA-A*02:01)")
    print("=" * 105)
    header = (
        f"{'#':<3} {'Gene':<10} {'Mutation':<10} {'Best 9-mer':<12} {'Pos':<5} "
        f"{'MT IC50 (nM)':<14} {'WT IC50':<12} {'Agret.':<8} {'RNA TPM':<10} {'Reads (Alt/Tot)'}"
    )
    print(header)
    print("-" * 105)

    for i, res in enumerate(results[:10], 1):
        epi = res.best_epitope
        print(
            f"{i:<3} "
            f"{res.gene:<10} "
            f"{res.protein_change:<10} "
            f"{epi.mt_peptide:<12} "
            f"P{epi.mut_pos_in_kmer:<4} "
            f"{epi.mt_ic50_nm:<14.1f} "
            f"{epi.wt_ic50_nm:<12.1f} "
            f"{epi.agretopicity_index:<8.2f} "
            f"{res.tpm:<10.1f} "
            f"{res.alt_reads}/{res.total_reads} (AF={res.af:.2f})"
        )

    print("-" * 105)
    print("Step 3 complete. Real HLA-binding neopeptides identified and mapped to patient transcriptome.")
    print("=" * 105)


if __name__ == "__main__":
    main()
