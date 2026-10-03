"""
Mini-NeoScorer: High-Throughput HLA Screening Module (Scale Step 3).
Executes high-throughput k-mer sliding window extraction and HLA-A*02:01
affinity scoring across hundreds of real patient somatic mutations.
"""

from dataclasses import dataclass
from typing import Dict, List, Optional, Tuple
from src.hla_scorer import EpitopeCandidate, generate_9mers_for_context
from src.protein_resolver import ProteinResolver


@dataclass
class PatientScreenResult:
    gene: str
    protein_change: str
    chrom: str
    pos: int
    ref: str
    alt: str
    alt_reads: int
    total_reads: int
    af: float
    tpm: float
    best_epitope: EpitopeCandidate
    all_9mers_evaluated: int


def screen_patient_mutations(
    mutations: List[dict],
    expression_map: Dict[str, float],
    hla_allele: str = "HLA-A*02:01",
    tpm_threshold: float = 1.0,
    resolver: Optional[ProteinResolver] = None,
) -> Tuple[List[PatientScreenResult], Dict[str, int]]:
    """
    Screens an entire patient's somatic mutations:
    1. Filters to missense mutations
    2. Applies RNA Expression Gate (TPM >= threshold)
    3. Resolves reference protein sequences
    4. Generates and scores all 9-mers against target HLA
    """
    if resolver is None:
        resolver = ProteinResolver()

    stats = {
        "total_somatic_mutations": len(mutations),
        "missense_mutations": 0,
        "gated_out_low_rna": 0,
        "isoform_mismatch_or_unresolved": 0,
        "successfully_screened_variants": 0,
        "total_9mers_evaluated": 0,
        "strong_binders_found": 0,
        "weak_binders_found": 0,
        "non_binders_found": 0,
    }

    screened_results: List[PatientScreenResult] = []

    # Prefetch uncached sequences in parallel for all expressed candidate genes
    genes_to_prefetch = []
    for m in mutations:
        if m.get("mutationType") == "Missense_Mutation":
            entrez = str(m.get("entrezGeneId", ""))
            tpm = expression_map.get(entrez, 0.0)
            if tpm >= tpm_threshold:
                g = m.get("gene", {}).get("hugoGeneSymbol")
                if g:
                    genes_to_prefetch.append(g)

    resolver.prefetch_proteins(genes_to_prefetch, max_workers=16)

    for m in mutations:
        if m.get("mutationType") != "Missense_Mutation":
            continue

        stats["missense_mutations"] += 1

        gene = m.get("gene", {}).get("hugoGeneSymbol", "")
        pc = m.get("proteinChange", "")
        entrez = str(m.get("entrezGeneId", ""))
        tpm = expression_map.get(entrez, 0.0)

        # 1. RNA Expression Gate
        if tpm < tpm_threshold:
            stats["gated_out_low_rna"] += 1
            continue

        # 2. Automated Reference Protein Resolution & Context Window
        ctx = resolver.extract_context_window(gene, pc, flank=10)
        if not ctx:
            stats["isoform_mismatch_or_unresolved"] += 1
            continue

        # Calculate Allele Frequency (AF)
        alt_reads = m.get("tumorAltCount", 0)
        ref_reads = m.get("tumorRefCount", 0)
        tot_reads = alt_reads + ref_reads
        af = alt_reads / max(1, tot_reads)

        # 3. 9-mer Sliding Window across 21-mer
        wt_seq = ctx["wt_window"]
        mt_seq = ctx["mt_window"]

        # Mock PeptideContext-like wrapper for generate_9mers_for_context
        from src.data_models import PeptideContext, SomaticVariant

        dummy_var = SomaticVariant(
            chrom=str(m.get("chr", "")),
            pos_1based=int(m.get("startPosition", 0)),
            ref=m.get("referenceAllele", ""),
            alt=m.get("variantAllele", ""),
            gene=gene,
            hgvsp=pc,
            hgvsc="",
            af=af,
            depth=tot_reads,
            tier="Clinical_Melanoma",
        )
        pep_ctx = PeptideContext(
            gene=gene,
            variant=dummy_var,
            tpm=tpm,
            mutation_aa_pos=ctx["pos_1based"],
            wt_residue=ctx["wt_aa"],
            mt_residue=ctx["mt_aa"],
            wt_21mer=wt_seq,
            mt_21mer=mt_seq,
        )

        candidates = generate_9mers_for_context(pep_ctx, hla_allele=hla_allele)
        if not candidates:
            continue

        stats["successfully_screened_variants"] += 1
        stats["total_9mers_evaluated"] += len(candidates)

        best_candidate = candidates[0]  # sorted by lowest MT IC50

        if best_candidate.mt_binder_class == "Strong Binder":
            stats["strong_binders_found"] += 1
        elif best_candidate.mt_binder_class == "Weak Binder":
            stats["weak_binders_found"] += 1
        else:
            stats["non_binders_found"] += 1

        screened_results.append(
            PatientScreenResult(
                gene=gene,
                protein_change=pc,
                chrom=str(m.get("chr", "")),
                pos=int(m.get("startPosition", 0)),
                ref=m.get("referenceAllele", ""),
                alt=m.get("variantAllele", ""),
                alt_reads=alt_reads,
                total_reads=tot_reads,
                af=af,
                tpm=tpm,
                best_epitope=best_candidate,
                all_9mers_evaluated=len(candidates),
            )
        )

    # Sort screened variants by lowest MT IC50
    screened_results.sort(key=lambda r: r.best_epitope.mt_ic50_nm)
    return screened_results, stats
