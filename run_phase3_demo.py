"""
Mini-NeoScorer: Phase 3 Execution & Verification Script
k-mer Sliding Window Epitope Generation and HLA-A*02:01 Affinity Scoring.
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
from src.hla_scorer import generate_9mers_for_context


def main():
    vcf_path = BASE_DIR / "data" / "mock_mutations.vcf"
    tpm_path = BASE_DIR / "data" / "mock_expression_tpm.tsv"

    print("=" * 90)
    print(" MINI-NEOSCORER: PHASE 3 PIPELINE EXECUTION")
    print(" 9-mer Sliding Window Generation & HLA-A*02:01 Affinity Prediction")
    print(" Target Patient HLA Allele: HLA-A*02:01 (k = 9)")
    print("=" * 90)

    # 1. Run Phase 2 to get RNA-expressed peptide contexts
    passed_contexts, _ = run_phase2(vcf_path, tpm_path, tpm_threshold=1.0)

    all_top_epitopes = []

    for ctx in passed_contexts:
        print(f"\nVariant Target: {ctx.gene} {ctx.variant.hgvsp} (TPM = {ctx.tpm:.1f}, AF = {ctx.variant.af:.2f})")
        print(f"21-mer Context: {ctx.mt_21mer}")
        print("-" * 90)
        print(f"{'#':<3} {'Pos':<5} {'MT 9-mer':<12} {'WT 9-mer':<12} {'MT IC50 (nM)':<14} {'MT Class':<15} {'WT IC50':<12} {'Agret. Index':<12}")
        print("-" * 90)

        # Generate and score all 9-mers spanning the mutation
        candidates = generate_9mers_for_context(ctx, hla_allele="HLA-A*02:01")

        for idx, c in enumerate(candidates, 1):
            is_top = " *" if idx == 1 else ""
            print(
                f"{idx:<3} "
                f"P{c.mut_pos_in_kmer:<4} "
                f"{c.mt_peptide:<12} "
                f"{c.wt_peptide:<12} "
                f"{c.mt_ic50_nm:<14.1f} "
                f"{c.mt_binder_class + is_top:<15} "
                f"{c.wt_ic50_nm:<12.1f} "
                f"{c.agretopicity_index:<12.2f}"
            )

        top = candidates[0]
        all_top_epitopes.append((ctx, top))
        print(f"  -> Top Predicted Neoantigen: {top.mt_peptide} (IC50 = {top.mt_ic50_nm:.1f} nM | {top.mt_binder_class})")

    print("\n" + "=" * 90)
    print(" SUMMARY: TOP CANDIDATE NEOPROTEINS FOR HLA-A*02:01")
    print("=" * 90)
    print(f"{'Gene':<8} {'Mutation':<14} {'Best MT 9-mer':<14} {'MT IC50':<12} {'MT Class':<14} {'Agret. Index':<14} {'Biological Assessment'}")
    print("-" * 90)

    for ctx, top in all_top_epitopes:
        assessment = ""
        if top.mt_binder_class == "Strong Binder":
            if top.agretopicity_index > 1.5:
                assessment = "High Immunogenicity: Strong binder with gain-of-affinity over WT"
            else:
                assessment = "High Presentation: Strong binder (similar affinity to WT)"
        elif top.mt_binder_class == "Weak Binder":
            assessment = "Moderate Presentation: Weak HLA binding"
        else:
            assessment = "Low Presentation: Poor HLA-A*02:01 fit, falls off tray"

        print(
            f"{ctx.gene:<8} "
            f"{ctx.variant.hgvsp:<14} "
            f"{top.mt_peptide:<14} "
            f"{top.mt_ic50_nm:<12.1f} "
            f"{top.mt_binder_class:<14} "
            f"{top.agretopicity_index:<14.2f} "
            f"{assessment}"
        )

    print("=" * 90)
    print("Phase 3 complete. Ready for Phase 4 (Multi-parametric Neoantigen Prioritization Scoring).")
    print("=" * 90)


if __name__ == "__main__":
    main()
