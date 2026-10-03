"""
Mini-NeoScorer: Scale Step 4 - Clinical Vaccine Formulation & Prioritization.
Takes the exome-wide screening results of patient TCGA-D3-A1Q1 (Melanoma),
applies the multi-parametric prioritization math, and generates the final
Top-10 Clinical Vaccine Prescription Report.
"""

import json
import sys
from pathlib import Path

# Ensure UTF-8 output
if hasattr(sys.stdout, "reconfigure"):
    sys.stdout.reconfigure(encoding="utf-8")

BASE_DIR = Path(__file__).resolve().parent.parent
DATA_DIR = BASE_DIR / "data"
REPORT_DIR = BASE_DIR / "reports"

if str(BASE_DIR) not in sys.path:
    sys.path.insert(0, str(BASE_DIR))

MUT_FILE = DATA_DIR / "real_patient_TCGA-D3-A1Q1-06_mutations.json"
EXP_FILE = DATA_DIR / "real_patient_TCGA-D3-A1Q1-06_expression.json"

from src.scale_hla_screener import screen_patient_mutations
from src.prioritizer import (
    compute_clonality_score,
    compute_expression_score,
    compute_presentation_score,
    compute_agretopicity_bonus,
)


def main():
    print("=" * 115)
    print(" MINI-NEOSCORER SCALE STEP 4: FINAL CLINICAL VACCINE FORMULATION & MULTI-PARAMETRIC RANKING")
    print(" Patient: TCGA-D3-A1Q1 (Metastatic Cutaneous Melanoma) | Target Allele: HLA-A*02:01")
    print(" Data Source: TCGA PanCancer Atlas (Cell, 2018) via cBioPortal API")
    print("=" * 115)

    with open(MUT_FILE, "r", encoding="utf-8") as f:
        mutations = json.load(f)

    with open(EXP_FILE, "r", encoding="utf-8") as f:
        expression = json.load(f)

    # 1. Screen all candidate mutations
    print("\n>>> Screening patient exome (Gate: RNA TPM >= 1.0, 9-mer sliding window)...")
    screened_results, stats = screen_patient_mutations(
        mutations, expression, hla_allele="HLA-A*02:01", tpm_threshold=1.0
    )
    print(f"  * Total candidate 9-mers evaluated : {stats['total_9mers_evaluated']}")
    print(f"  * Strong HLA-A*02:01 Binders found : {stats['strong_binders_found']}")

    # 2. Multi-parametric Prioritization Scoring
    scored_candidates = []
    scale_factor = 20.0

    for res in screened_results:
        epi = res.best_epitope
        c_af = compute_clonality_score(res.af)
        e_tpm = compute_expression_score(res.tpm)
        p_mhc = compute_presentation_score(epi.mt_ic50_nm)
        a_bonus = compute_agretopicity_bonus(epi.agretopicity_index)

        raw_score = c_af * e_tpm * p_mhc * a_bonus
        final_score = round(min(100.0, raw_score * scale_factor), 1)

        # Clinical tier assignment
        if final_score >= 50.0 and epi.mt_ic50_nm <= 50.0:
            clinical_tier = "Tier 1: High Priority (Vaccine Payload)"
            action = "INCLUDE: High clonality, robust expression, strong HLA-A*02:01 fit."
        elif final_score >= 25.0 and epi.mt_ic50_nm <= 500.0:
            clinical_tier = "Tier 2: Secondary / Backup"
            action = "BACKUP: Moderate clonality or expression. Viable secondary target."
        else:
            clinical_tier = "Tier 3: Deprioritized"
            if epi.mt_ic50_nm > 500.0:
                action = "REJECT: Failed HLA-A*02:01 presentation (IC50 > 500 nM)."
            elif res.tpm < 5.0:
                action = "REJECT: Insufficient transcript expression in tumor."
            else:
                action = "REJECT: Subclonal / low presentation probability."

        scored_candidates.append({
            "gene": res.gene,
            "protein_change": res.protein_change,
            "best_9mer": epi.mt_peptide,
            "wt_9mer": epi.wt_peptide,
            "mut_pos_in_kmer": epi.mut_pos_in_kmer,
            "mt_ic50_nm": epi.mt_ic50_nm,
            "wt_ic50_nm": epi.wt_ic50_nm,
            "agretopicity": epi.agretopicity_index,
            "tpm": res.tpm,
            "af": res.af,
            "alt_reads": res.alt_reads,
            "total_reads": res.total_reads,
            "clonality_score": round(c_af, 2),
            "expression_score": round(e_tpm, 2),
            "presentation_score": round(p_mhc, 2),
            "agretopicity_bonus": round(a_bonus, 2),
            "composite_score": final_score,
            "clinical_tier": clinical_tier,
            "action": action,
        })

    # Sort descending by composite score
    scored_candidates.sort(key=lambda x: x["composite_score"], reverse=True)

    # 3. Print Top-15 Clinical Vaccine Formulation Table
    print("\n" + "=" * 115)
    print(" RECOMMENDED CLINICAL VACCINE FORMULATION: TOP TARGETS")
    print("=" * 115)
    header = (
        f"{'Rank':<5} {'Gene':<9} {'Mutation':<10} {'Neopeptide':<12} "
        f"{'AF':<6} {'Reads':<8} {'RNA TPM':<9} {'MT IC50':<11} {'Agret.':<8} {'Score':<7} {'Clinical Tier'}"
    )
    print(header)
    print("-" * 115)

    for rank, c in enumerate(scored_candidates[:12], 1):
        reads_str = f"{c['alt_reads']}/{c['total_reads']}"
        print(
            f"#{rank:<4} "
            f"{c['gene']:<9} "
            f"{c['protein_change']:<10} "
            f"{c['best_9mer']:<12} "
            f"{c['af']:<6.2f} "
            f"{reads_str:<8} "
            f"{c['tpm']:<9.1f} "
            f"{c['mt_ic50_nm']:<11.1f} "
            f"{c['agretopicity']:<8.2f} "
            f"{c['composite_score']:<7.1f} "
            f"{c['clinical_tier']}"
        )

    print("-" * 115)

    # 4. Detailed Biological Rationale for Top 3 Vaccine Hits
    print("\n>>> DETAILED ONCOLOGY RATIONALE FOR TOP 3 VACCINE CANDIDATES:")
    for rank, c in enumerate(scored_candidates[:3], 1):
        print(f"\n [RANK #{rank}] {c['gene']} {c['protein_change']} -> Neopeptide: '{c['best_9mer']}' (Score: {c['composite_score']}/100)")
        print(f"   |-- Clonality Profile      : AF = {c['af']:.2f} ({c['alt_reads']}/{c['total_reads']} reads) -> {c['clonality_score']*100:.0f}% tumor cell coverage")
        print(f"   |-- Transcriptional Output : TPM = {c['tpm']:.1f} (Abundance Score: {c['expression_score']:.2f})")
        print(f"   |-- HLA-A*02:01 Affinity   : MT IC50 = {c['mt_ic50_nm']:.1f} nM (Presentation Probability: {c['presentation_score']:.2f})")
        print(f"   |-- Differential Binding   : WT IC50 = {c['wt_ic50_nm']:.1f} nM -> Agretopicity Index = {c['agretopicity']:.2f}x (Bonus: x{c['agretopicity_bonus']:.2f})")
        print(f"   \\-- Clinical Assessment    : {c['action']}")

    # 5. Save persistent clinical formulation TSV report
    REPORT_DIR.mkdir(parents=True, exist_ok=True)
    report_file = REPORT_DIR / "TCGA_D3_A1Q1_clinical_vaccine_prescription.tsv"
    with open(report_file, "w", encoding="utf-8") as f:
        f.write("Rank\tGene\tProteinChange\tNeopeptide_9mer\tWT_9mer\tMutPos\tAF\tReads\tTPM\tMT_IC50_nM\tWT_IC50_nM\tAgretopicity\tScore\tClinicalTier\tAction\n")
        for rank, c in enumerate(scored_candidates, 1):
            reads_str = f"{c['alt_reads']}/{c['total_reads']}"
            f.write(
                f"{rank}\t{c['gene']}\t{c['protein_change']}\t{c['best_9mer']}\t{c['wt_9mer']}\t"
                f"{c['mut_pos_in_kmer']}\t{c['af']:.3f}\t{reads_str}\t{c['tpm']:.2f}\t{c['mt_ic50_nm']:.2f}\t"
                f"{c['wt_ic50_nm']:.2f}\t{c['agretopicity']:.2f}\t{c['composite_score']:.1f}\t"
                f"{c['clinical_tier']}\t{c['action']}\n"
            )

    print("\n" + "=" * 115)
    print(f"[OK] Full clinical vaccine prescription saved to: {report_file}")
    print("     Mini-NeoScorer Scale Step 4 complete. Personalized cancer vaccine design finalized!")
    print("=" * 115)


if __name__ == "__main__":
    main()
