"""
Neoantigen-Prioritization-Pipeline (NeoScorer)
Master Clinical In Silico Neoantigen Prioritization Pipeline.

Authors: Computational Immuno-Oncology Lab
Target: Automated discovery and clinical tiering of patient-specific neoepitopes
        for personalized mRNA/peptide cancer vaccines using pre-trained Deep Neural Networks.
"""

import argparse
import json
from datetime import datetime, timezone
from patient_io import validate_inputs, report_name
from pathlib import Path
import sys

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
    predictor_name: str = "ensemble",
    top_n: int = 12,
    patient_id: str | None = None,
    output_dir: Path = REPORT_DIR,
):
    with mutations_path.open(encoding="utf-8-sig") as stream:
        mutations = json.load(stream)
    with expression_path.open(encoding="utf-8-sig") as stream:
        expression = json.load(stream)
    patient_id = validate_inputs(mutations, expression, patient_id)
    if predictor_name == "pwm" and hla_allele != "HLA-A*02:01":
        raise ValueError("The PWM baseline supports only HLA-A*02:01.")
    print("=" * 125)
    print(" NEOANTIGEN-PRIORITIZATION-PIPELINE: CLINICAL IMMUNO-ONCOLOGY FORMULATION")
    print(f" Patient: {patient_id} | Target Allele: {hla_allele}")
    print(
        f" Scoring Engine: {predictor_name.upper()} (Deep Neural Network Consensus / Presentation)"
    )
    print(
        f" Inputs: {mutations_path.name} | {expression_path.name} (cBioPortal-compatible JSON)"
    )
    print("=" * 125)

    # 1. High-throughput exome screening with selected predictor
    print(
        f"\n>>> [STAGE 1] Exome-Wide Screening (Gate: RNA TPM >= {tpm_threshold:.1f}, Predictor: {predictor_name.upper()})..."
    )
    screened_results, stats = screen_patient_mutations(
        mutations,
        expression,
        hla_allele=hla_allele,
        tpm_threshold=tpm_threshold,
        predictor_name=predictor_name,
    )

    print(
        f"  * Total somatic variants in tumor biopsy   : {stats['total_somatic_mutations']}"
    )
    print(
        f"  * Missense protein-altering SNVs           : {stats['missense_mutations']}"
    )
    print(
        f"  * Eliminated by RNA Gate (TPM < {tpm_threshold:.1f})      : {stats['gated_out_low_rna']} (Silent / unexpressed)"
    )
    print(
        f"  * Filtered by Isoform / QC boundary checks : {stats['isoform_mismatch_or_unresolved']}"
    )
    print(
        f"  * Successfully screened patient variants   : {stats['successfully_screened_variants']}"
    )
    print(
        f"  * Total candidate 9-mers evaluated on HLA  : {stats['total_9mers_evaluated']}"
    )
    print(
        f"    |-- Strong Binders (IC50 <= 50 nM)       : {stats['strong_binders_found']}"
    )
    print(
        f"    |-- Weak Binders (50 - 500 nM)           : {stats['weak_binders_found']}"
    )
    print(
        f"    \\-- Non-Binders (IC50 > 500 nM)          : {stats['non_binders_found']}"
    )

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
            action = f"INCLUDE: High clonality, robust expression, strong {hla_allele} presentation."
        elif final_score >= 25.0 and epi.mt_ic50_nm <= 500.0:
            clinical_tier = "Tier 2: Secondary / Backup"
            action = (
                "BACKUP: Moderate clonality or expression. Viable secondary target."
            )
        else:
            clinical_tier = "Tier 3: Deprioritized"
            if epi.mt_ic50_nm > 500.0:
                action = f"REJECT: Failed {hla_allele} presentation (IC50 > 500 nM)."
            elif res.tpm < 5.0:
                action = "REJECT: Insufficient transcript expression in tumor."
            else:
                action = "REJECT: Subclonal / low presentation probability."

        scored_candidates.append(
            {
                "gene": res.gene,
                "protein_change": res.protein_change,
                "best_9mer": epi.mt_peptide,
                "wt_9mer": epi.wt_peptide,
                "mut_pos_in_kmer": epi.mut_pos_in_kmer,
                "mt_ic50_nm": epi.mt_ic50_nm,
                "wt_ic50_nm": epi.wt_ic50_nm,
                "percentile_rank": epi.percentile_rank,
                "presentation_score": epi.presentation_score,
                "mhcflurry_ic50": epi.mhcflurry_ic50,
                "netmhcpan_ic50": epi.netmhcpan_ic50,
                "predictor_source": epi.predictor_source,
                "agretopicity": epi.agretopicity_index,
                "tpm": res.tpm,
                "af": res.af,
                "alt_reads": res.alt_reads,
                "total_reads": res.total_reads,
                "clonality_score": round(c_af, 2),
                "expression_score": round(e_tpm, 2),
                "presentation_prob": round(p_mhc, 2),
                "agretopicity_bonus": round(a_bonus, 2),
                "composite_score": final_score,
                "clinical_tier": clinical_tier,
                "action": action,
            }
        )

    # Sort descending by composite score
    scored_candidates.sort(key=lambda x: x["composite_score"], reverse=True)

    # 3. Print Clinical Table
    print("\n" + "=" * 125)
    print(
        f" TOP {top_n} RECOMMENDED CLINICAL VACCINE CANDIDATES (Predictor: {predictor_name.upper()})"
    )
    print("=" * 125)
    header = (
        f"{'Rank':<5} {'Gene':<9} {'Mutation':<10} {'Neopeptide':<11} "
        f"{'AF':<5} {'RNA TPM':<9} {'IC50 (nM)':<11} {'Rank %':<8} {'Agret.':<7} {'Score':<6} {'Clinical Tier'}"
    )
    print(header)
    print("-" * 125)

    for rank, c in enumerate(scored_candidates[:top_n], 1):
        rank_str = (
            f"{c['percentile_rank']:.2f}%"
            if c["percentile_rank"] is not None
            else "N/A"
        )
        print(
            f"#{rank:<4} "
            f"{c['gene']:<9} "
            f"{c['protein_change']:<10} "
            f"{c['best_9mer']:<11} "
            f"{c['af']:<5.2f} "
            f"{c['tpm']:<9.1f} "
            f"{c['mt_ic50_nm']:<11.1f} "
            f"{rank_str:<8} "
            f"{c['agretopicity']:<7.2f} "
            f"{c['composite_score']:<6.1f} "
            f"{c['clinical_tier']}"
        )
    print("-" * 125)

    # 4. Export persistent clinical prescription report
    output_dir.mkdir(parents=True, exist_ok=True)
    report_file = output_dir / report_name(patient_id)
    temporary_report = report_file.with_suffix(".tsv.tmp")
    with open(temporary_report, "w", encoding="utf-8") as f:
        f.write(
            "Rank\tGene\tProteinChange\tNeopeptide_9mer\tWT_9mer\tMutPos\tAF\tReads\tTPM\t"
            "MT_IC50_nM\tPercentile_Rank\tPresentation_Score\tMHCflurry_IC50\tNetMHCpan_IC50\t"
            "WT_IC50_nM\tAgretopicity\tScore\tClinicalTier\tAction\tPredictor\tHLA_Allele\n"
        )
        for rank, c in enumerate(scored_candidates, 1):
            reads_str = f"{c['alt_reads']}/{c['total_reads']}"
            mf_str = (
                f"{c['mhcflurry_ic50']:.1f}"
                if c["mhcflurry_ic50"] is not None
                else "N/A"
            )
            net_str = (
                f"{c['netmhcpan_ic50']:.1f}"
                if c["netmhcpan_ic50"] is not None
                else "N/A"
            )
            f.write(
                f"{rank}\t{c['gene']}\t{c['protein_change']}\t{c['best_9mer']}\t{c['wt_9mer']}\t"
                f"{c['mut_pos_in_kmer']}\t{c['af']:.3f}\t{reads_str}\t{c['tpm']:.2f}\t"
                f"{c['mt_ic50_nm']:.2f}\t{c['percentile_rank']:.2f}\t{c['presentation_score']:.4f}\t"
                f"{mf_str}\t{net_str}\t{c['wt_ic50_nm']:.2f}\t{c['agretopicity']:.2f}\t"
                f"{c['composite_score']:.1f}\t{c['clinical_tier']}\t{c['action']}\t{c['predictor_source']}\t{hla_allele}\n"
            )

    metadata = {
        "patient_id": patient_id,
        "hla_allele": hla_allele,
        "predictor": predictor_name,
        "tpm_threshold": tpm_threshold,
        "requested_top": top_n,
        "stats": stats,
        "created_at": datetime.now(timezone.utc).isoformat(),
        "mutations_file": mutations_path.name,
        "expression_file": expression_path.name,
    }
    metadata_file = report_file.with_suffix(".json")
    temporary_metadata = metadata_file.with_suffix(".json.tmp")
    temporary_metadata.write_text(json.dumps(metadata, indent=2), encoding="utf-8")
    temporary_metadata.replace(metadata_file)
    temporary_report.replace(report_file)

    print("\n[OK] Clinical vaccine formulation report successfully generated:")
    print(f"     -> {report_file}")
    print("=" * 125)


def parse_args():
    parser = argparse.ArgumentParser(
        description="Neoantigen-Prioritization-Pipeline: In Silico Personalized Cancer Vaccine Formulation"
    )
    parser.add_argument(
        "--predictor",
        choices=["ensemble", "mhcflurry", "netmhcpan", "pwm"],
        default="ensemble",
        help="HLA binding affinity & presentation engine (default: ensemble)",
    )
    parser.add_argument(
        "--hla",
        default="HLA-A*02:01",
        help="Target patient HLA Class I allele (default: HLA-A*02:01)",
    )
    parser.add_argument(
        "--tpm-threshold",
        type=float,
        default=1.0,
        help="Minimum RNA-seq TPM gate (default: 1.0)",
    )
    parser.add_argument(
        "--top",
        type=int,
        default=12,
        help="Number of top vaccine targets to display (default: 12)",
    )
    parser.add_argument("--mutations", type=Path, default=MUT_FILE)
    parser.add_argument("--expression", type=Path, default=EXP_FILE)
    parser.add_argument(
        "--patient-id", default=None, help="Inferred from mutation patientId if omitted"
    )
    parser.add_argument("--output-dir", type=Path, default=REPORT_DIR)
    return parser.parse_args()


if __name__ == "__main__":
    args = parse_args()
    run_pipeline(
        mutations_path=args.mutations,
        expression_path=args.expression,
        hla_allele=args.hla,
        tpm_threshold=args.tpm_threshold,
        predictor_name=args.predictor,
        top_n=args.top,
        patient_id=args.patient_id,
        output_dir=args.output_dir,
    )
