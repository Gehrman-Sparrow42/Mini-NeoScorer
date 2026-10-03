"""
Execution script for Mini-NeoScorer Phase 2:
Coordinate Translation & RNA-seq TPM Filtering.
"""

import sys
from pathlib import Path

# Ensure UTF-8 output on Windows consoles
if hasattr(sys.stdout, "reconfigure"):
    sys.stdout.reconfigure(encoding="utf-8")

BASE_DIR = Path(__file__).resolve().parent
if str(BASE_DIR) not in sys.path:
    sys.path.insert(0, str(BASE_DIR))

from src.phase2_pipeline import run_phase2


def main():
    vcf_path = BASE_DIR / "data" / "mock_mutations.vcf"
    tpm_path = BASE_DIR / "data" / "mock_expression_tpm.tsv"

    print("=" * 80)
    print(" MINI-NEOSCORER: PHASE 2 PIPELINE EXECUTION")
    print(" Coordinate Systems, RNA Expression Gate, and Peptide Context Extraction")
    print("=" * 80)

    passed, filtered = run_phase2(vcf_path, tpm_path, tpm_threshold=1.0)

    print("\n[1] RNA-seq TPM Filtering Results (Threshold >= 1.0 TPM):")
    print("-" * 80)
    for var, tpm, reason in filtered:
        print(f" [FILTERED OUT] {var.gene} ({var.hgvsp})")
        print(f"    Genomic POS: 1-based={var.pos_1based} -> 0-based={var.pos_0based}")
        print(f"    Allele Frequency (AF): {var.af:.2f}")
        print(f"    Reason: {reason}\n")

    print(f"[2] Passed Variants Generating Expressed Neopeptide Windows ({len(passed)} candidates):")
    print("-" * 80)

    for i, ctx in enumerate(passed, 1):
        v = ctx.variant
        print(f" Candidate {i}: {ctx.gene} {ctx.wt_residue}{ctx.mutation_aa_pos}{ctx.mt_residue} ({v.tier})")
        print(f"   |-- Coordinates  : Chr {v.chrom}:{v.pos_1based} (1-based VCF) -> {v.pos_0based} (0-based Python index)")
        print(f"   |-- Clonality    : AF = {v.af:.2f} ({v.af * 100:.1f}%) | Read Depth = {v.depth}x")
        print(f"   |-- Expression   : TPM = {ctx.tpm:.2f}")
        print(f"   |-- WT Context   : {ctx.wt_21mer}")
        print(f"   \\-- MT Context   : {ctx.mt_21mer}")

        # Show exactly where the single amino acid mismatch occurs
        # The mutation is centered: 10 upstream + 1 mutated + 10 downstream
        # For KRAS G12D at position 12, upstream has 11 amino acids (pos 1 to 11)
        mut_idx_in_window = len(ctx.wt_21mer) - len(ctx.variant.gene) # placeholder
        # Find exact mismatch index:
        mismatch_idx = next(
            (idx for idx, (c1, c2) in enumerate(zip(ctx.wt_21mer, ctx.mt_21mer)) if c1 != c2),
            None
        )
        if mismatch_idx is not None:
            pointer = " " * (22 + mismatch_idx) + f"^ ({ctx.wt_residue} -> {ctx.mt_residue})"
            print(pointer)
        print()

    print("=" * 80)
    print("Phase 2 complete. Ready for Phase 3 (k-mer sliding window & HLA Class I scoring).")
    print("=" * 80)


if __name__ == "__main__":
    main()
