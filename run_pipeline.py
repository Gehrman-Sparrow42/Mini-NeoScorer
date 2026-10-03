"""
Mini-NeoScorer: End-to-End Computational Neoantigen Prioritization Pipeline.
Orchestrates Phases 1 through 4 to produce a ranked clinical recommendation report.
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
from src.prioritizer import score_and_rank_candidates


def main():
    vcf_path = BASE_DIR / "data" / "mock_mutations.vcf"
    tpm_path = BASE_DIR / "data" / "mock_expression_tpm.tsv"
    hla_target = "HLA-A*02:01"

    print("=" * 105)
    print("       MINI-NEOSCORER: CLINICAL IN SILICO NEOANTIGEN PRIORITIZATION PIPELINE")
    print(f"       Patient Somatic Input: {vcf_path.name} | RNA Matrix: {tpm_path.name}")
    print(f"       Patient Target Class I Allele: {hla_target}")
    print("=" * 105)

    # Step 1 & 2: Parse VCF/TPM, apply RNA Expression Gate, extract 21-mer contexts
    print("\n>>> STEP 1: Somatic Mutation & RNA Expression Gating (Threshold >= 1.0 TPM)")
    passed_contexts, filtered_out = run_phase2(vcf_path, tpm_path, tpm_threshold=1.0)
    for var, tpm, reason in filtered_out:
        print(f"  [GATE REJECTED] {var.gene} ({var.hgvsp}): {reason}")
    print(f"  -> {len(passed_contexts)} variants passed the transcriptional expression gate.")

    # Step 3: 9-mer Sliding Window & HLA Class I Affinity Scoring
    print(f"\n>>> STEP 2: 9-mer Slicing & {hla_target} Presentation Prediction")
    top_epitope_pairs = []
    for ctx in passed_contexts:
        candidates = generate_9mers_for_context(ctx, hla_allele=hla_target)
        best_epitope = candidates[0]
        top_epitope_pairs.append((ctx, best_epitope))
        print(f"  * {ctx.gene:<7} {ctx.variant.hgvsp:<14} -> Best 9-mer: {best_epitope.mt_peptide} "
              f"(IC50 = {best_epitope.mt_ic50_nm:>8.1f} nM | {best_epitope.mt_binder_class})")

    # Step 4: Multi-parametric Composite Prioritization Scoring
    print("\n>>> STEP 3: Multi-Parametric Immuno-Oncology Prioritization Scoring")
    ranked_results = score_and_rank_candidates(top_epitope_pairs)

    print("\n" + "=" * 105)
    print(" FINAL RANKED CLINICAL VACCINE RECOMMENDATIONS")
    print("=" * 105)
    header = (
        f"{'Rank':<5} {'Gene':<8} {'Mutation':<14} {'Neopeptide':<12} "
        f"{'AF':<6} {'TPM':<8} {'MT IC50':<12} {'Agret.':<8} {'Score':<8} {'Tier'}"
    )
    print(header)
    print("-" * 105)

    for rank, res in enumerate(ranked_results, 1):
        print(
            f"#{rank:<4} "
            f"{res.gene:<8} "
            f"{res.variant_hgvsp:<14} "
            f"{res.mt_peptide:<12} "
            f"{res.af:<6.2f} "
            f"{res.tpm:<8.1f} "
            f"{res.mt_ic50_nm:<12.1f} "
            f"{res.agretopicity_index:<8.2f} "
            f"{res.composite_score:<8.1f} "
            f"{res.clinical_tier}"
        )

    print("-" * 105)
    print("\n>>> CLINICAL ACTIONABILITY & RATIONALE REPORT:")
    for rank, res in enumerate(ranked_results, 1):
        print(f"\n [RANK #{rank}] {res.gene} {res.variant_hgvsp} (Epitope: {res.mt_peptide}) | Score: {res.composite_score}/100")
        print(f"   |-- Clonality Component    : AF={res.af:.2f} -> Fraction of Tumor Cells Covered: {res.clonality_score * 100:.0f}%")
        print(f"   |-- Transcriptional Volume : TPM={res.tpm:.1f} -> Abundance Score: {res.expression_score:.2f}")
        print(f"   |-- HLA Presentation       : IC50={res.mt_ic50_nm:.1f} nM -> Presentation Probability: {res.presentation_score:.2f}")
        print(f"   |-- Differential Binding   : Agretopicity Index={res.agretopicity_index:.2f} (Bonus: x{res.agretopicity_bonus:.2f})")
        print(f"   \\-- Action Recommendation  : {res.recommendation}")

    print("\n" + "=" * 105)
    print("Mini-NeoScorer complete. Target ranking finalized.")
    print("=" * 105)


if __name__ == "__main__":
    main()
