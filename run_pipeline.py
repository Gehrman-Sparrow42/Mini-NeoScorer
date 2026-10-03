"""
Neoantigen-Prioritization-Pipeline (NeoScorer)
Master Clinical In Silico Neoantigen Prioritization Pipeline.

Authors: Computational Immuno-Oncology Lab
Target: Automated discovery and clinical tiering of patient-specific neoepitopes
        for personalized mRNA/peptide cancer vaccines.
"""

import argparse
import json
import sys
from pathlib import Path

# Ensure UTF-8 output on all consoles
if hasattr(sys.stdout, "reconfigure"):
    sys.stdout.reconfigure(encoding="utf-8")

BASE_DIR = Path(__file__).resolve().parent
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


def run_pipeline(
    mutations_path: Path,
    expression_path: Path,
    hla_allele: str = "HLA-A*02:01",
    tpm_threshold: float = 1.0,
    top_n: int = 12,
):
    print("=" * 115)
    print(" NEOANTIGEN-PRIORITIZATION-PIPELINE: CLINICAL IMMUNO-ONCOLOGY FORMULATION")
    print(f" Patient: TCGA-D3-A1Q1 (Metastatic Cutaneous Melanoma) | Target Allele: {hla_allele}")
    print(" Provenance: The Cancer Genome Atlas (TCGA PanCancer Atlas, Cell 2018) via cBioPortal API")
    print("=" * 115)

    if not mutations_path.exists() or not expression_path.exists():
        raise FileNotFoundError(
            f"Patient multi-omics datasets not found. Run scripts/fetch_patient_data.py first."
        )

    with open(mutations_path, "r", encoding="utf-8") as f:
        mutations = json.load(f)

    with open(expression_path, "r", encoding="utf-8") as f:
        expression = json.load(f)

    # 1. High-throughput exome screening
    print(f"\n>>> [STAGE 1] Exome-Wide Screening (Gate: RNA TPM >= {tpm_threshold:.1f}, HLA: {hla_allele})...")
    screened_results, stats = screen_patient_mutations(
        mutations, expression, hla_allele=hla_allele, tpm_threshold=tpm_threshold
    )

    print(f"  * Total somatic variants in tumor biopsy   : {stats['total_somatic_mutations']}")
    print(f"  * Missense protein-altering SNVs           : {stats['missense_mutations']}")
    print(f"  * Eliminated by RNA Gate (TPM < {tpm_threshold:.1f})      : {stats['gated_out_low_rna']} (Silent / unexpressed)")
    print(f"  * Filtered by Isoform / QC boundary checks : {stats['isoform_mismatch_or_unresolved']}")
    print(f"  * Successfully screened patient variants   : {stats['successfully_screened_variants']}")
    print(f"  * Total candidate 9-mers evaluated on HLA  : {stats['total_9mers_evaluated']}")
    print(f"    |-- Strong Binders (IC50 <= 50 nM)       : {stats['strong_binders_found']}")
    print(f"    |-- Weak Binders (50 - 500 nM)           : {stats['weak_binders_found']}")
    print(f"    \\-- Non-Binders (IC50 > 500 nM)          : {stats['non_binders_found']}")

    # 2. Multi-parametric Prioritization Scoring
    print("\n>>> [STAGE 2] Multi-Parametric Prioritization & Clinical Tiering...")
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

        if final_score >= 50.0 and epi.mt_ic50_nm <= 50.0:
            clinical_tier = "Tier 1: High Priority (Vaccine Payload)"
            action = "INCLUDE: High clonality, robust expression, strong HLA-A*02:01 presentation."
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

    # 3. Print Clinical Table
    print("\n" + "=" * 115)
    print(f" TOP {top_n} RECOMMENDED CLINICAL VACCINE CANDIDATES")
    print("=" * 115)
    header = (
        f"{'Rank':<5} {'Gene':<9} {'Mutation':<10} {'Neopeptide':<12} "
        f"{'AF':<6} {'Reads':<8} {'RNA TPM':<9} {'MT IC50':<11} {'Agret.':<8} {'Score':<7} {'Clinical Tier'}"
    )
    print(header)
    print("-" * 115)

    for rank, c in enumerate(scored_candidates[:top_n], 1):
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

    # 4. Export persistent clinical prescription report
    REPORT_DIR.mkdir(parents=True, exist_ok=True)
    report_file = REPORT_DIR / "TCGA_D3_A1Q1_clinical_vaccine_prescription.tsv"
    with open(report_file, "w", encoding="utf-8") as f:
        f.write(
            "Rank\tGene\tProteinChange\tNeopeptide_9mer\tWT_9mer\tMutPos\tAF\tReads\tTPM\t"
            "MT_IC50_nM\tWT_IC50_nM\tAgretopicity\tScore\tClinicalTier\tAction\n"
        )
        for rank, c in enumerate(scored_candidates, 1):
            reads_str = f"{c['alt_reads']}/{c['total_reads']}"
            f.write(
                f"{rank}\t{c['gene']}\t{c['protein_change']}\t{c['best_9mer']}\t{c['wt_9mer']}\t"
                f"{c['mut_pos_in_kmer']}\t{c['af']:.3f}\t{reads_str}\t{c['tpm']:.2f}\t{c['mt_ic50_nm']:.2f}\t"
                f"{c['wt_ic50_nm']:.2f}\t{c['agretopicity']:.2f}\t{c['composite_score']:.1f}\t"
                f"{c['clinical_tier']}\t{c['action']}\n"
            )

    print(f"\n[OK] Clinical vaccine formulation report successfully generated:")
    print(f"     -> {report_file}")
    print("=" * 115)


def parse_args():
    parser = argparse.ArgumentParser(
        description="Neoantigen-Prioritization-Pipeline: In Silico Personalized Cancer Vaccine Formulation"
    )
    parser.add_argument("--hla", default="HLA-A*02:01", help="Target patient HLA Class I allele (default: HLA-A*02:01)")
    parser.add_argument("--tpm-threshold", type=float, default=1.0, help="Minimum RNA-seq TPM gate (default: 1.0)")
    parser.add_argument("--top", type=int, default=12, help="Number of top vaccine targets to display (default: 12)")
    return parser.parse_args()


if __name__ == "__main__":
    args = parse_args()
    run_pipeline(
        mutations_path=MUT_FILE,
        expression_path=EXP_FILE,
        hla_allele=args.hla,
        tpm_threshold=args.tpm_threshold,
        top_n=args.top,
    )
