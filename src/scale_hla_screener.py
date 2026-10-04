"""
Mini-NeoScorer: High-Throughput HLA Screening Module (Scale Step 3).
Executes high-throughput k-mer sliding window extraction and HLA affinity scoring
across hundreds of patient somatic mutations using modular scoring engines:
- Ensemble (Consensus: MHCflurry 2.0 + NetMHCpan-4.1)
- MHCflurry 2.0 (Local PyTorch Deep Neural Network)
- NetMHCpan-4.1 (Official NIH IEDB REST API)
- PWM (Position Weight Matrix baseline)
"""

from dataclasses import dataclass
from typing import Dict, List, Optional, Tuple
from src.data_models import SomaticVariant
from src.dnn_predictor import (
    BaseHLAPredictor,
    EnsemblePredictor,
    MHCflurryPredictor,
    NetMHCpanIEDBPredictor,
    get_predictor,
)
from src.hla_scorer import EpitopeCandidate, classify_ic50
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
    predictor_name: str = "ensemble",
) -> Tuple[List[PatientScreenResult], Dict[str, int]]:
    """
    Screens an entire patient's somatic mutations:
    1. Filters to missense mutations
    2. Applies RNA Expression Gate (TPM >= threshold)
    3. Resolves reference protein sequences
    4. Generates and scores all 9-mers against target HLA using selected predictor
    5. Cleans up memory allocations
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
        "predictor_used": predictor_name,
    }

    screened_results: List[PatientScreenResult] = []

    # 1. Prefetch uncached sequences in parallel for all expressed candidate genes
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

    # 2. Extract context windows and candidate 9-mers
    extracted_variants = []
    all_peptides_to_score = set()

    for m in mutations:
        if m.get("mutationType") != "Missense_Mutation":
            continue

        stats["missense_mutations"] += 1

        gene = m.get("gene", {}).get("hugoGeneSymbol", "")
        pc = m.get("proteinChange", "")
        entrez = str(m.get("entrezGeneId", ""))
        tpm = expression_map.get(entrez, 0.0)

        # RNA Expression Gate
        if tpm < tpm_threshold:
            stats["gated_out_low_rna"] += 1
            continue

        # Automated Reference Protein Resolution & Context Window
        ctx = resolver.extract_context_window(gene, pc, flank=10)
        if not ctx:
            stats["isoform_mismatch_or_unresolved"] += 1
            continue

        # Calculate Allele Frequency (AF)
        alt_reads = m.get("tumorAltCount", 0)
        ref_reads = m.get("tumorRefCount", 0)
        tot_reads = alt_reads + ref_reads
        af = alt_reads / max(1, tot_reads)

        wt_seq = ctx["wt_window"]
        mt_seq = ctx["mt_window"]

        # Find mutation index in 21-mer
        mut_idx = next(
            (i for i, (c1, c2) in enumerate(zip(wt_seq, mt_seq)) if c1 != c2), None
        )
        if mut_idx is None:
            continue

        k = 9
        kmer_pairs = []
        for start in range(0, len(mt_seq) - k + 1):
            end = start + k
            if start <= mut_idx < end:
                wt_kmer = wt_seq[start:end]
                mt_kmer = mt_seq[start:end]
                mut_pos_in_kmer = (mut_idx - start) + 1
                kmer_pairs.append((wt_kmer, mt_kmer, mut_pos_in_kmer))
                all_peptides_to_score.add(wt_kmer)
                all_peptides_to_score.add(mt_kmer)

        if not kmer_pairs:
            continue

        extracted_variants.append(
            {
                "mutation": m,
                "gene": gene,
                "protein_change": pc,
                "alt_reads": alt_reads,
                "total_reads": tot_reads,
                "af": af,
                "tpm": tpm,
                "kmer_pairs": kmer_pairs,
            }
        )

    # 3. Score all peptides using selected predictor
    predictor: BaseHLAPredictor = get_predictor(predictor_name)
    all_peptides_list = list(all_peptides_to_score)
    scored_peptides: Dict[str, Dict[str, float]] = {}

    if predictor_name.lower() in ("ensemble", "consensus"):
        # Optimized two-stage ensemble:
        # Step A: High-throughput local MHCflurry pre-screening across all ~2,900 kmers (fast vectorized PyTorch)
        mf_pred = MHCflurryPredictor()
        mf_scores = mf_pred.predict_peptides(all_peptides_list, allele=hla_allele)

        # Step B: Identify candidate MT binders (IC50 < 2500 nM or top candidate per variant) and their paired WT
        priority_peptides = set()
        for item in extracted_variants:
            pairs = item["kmer_pairs"]
            # Find best MT candidate by MHCflurry
            best_pair = min(pairs, key=lambda p: mf_scores.get(p[1], {}).get("ic50", 50000.0))
            priority_peptides.add(best_pair[0])
            priority_peptides.add(best_pair[1])
            for wt_p, mt_p, _ in pairs:
                if mf_scores.get(mt_p, {}).get("ic50", 50000.0) <= 2500.0:
                    priority_peptides.add(wt_p)
                    priority_peptides.add(mt_p)

        # Step C: Query NetMHCpan-4.1 (NIH IEDB) for priority peptides with local cache
        net_pred = NetMHCpanIEDBPredictor()
        net_scores = net_pred.predict_peptides(list(priority_peptides), allele=hla_allele)

        # Step D: Compute Consensus Ensemble metrics
        import math
        for p in all_peptides_list:
            mf = mf_scores.get(p, {"ic50": 50000.0, "percentile_rank": 50.0, "presentation_score": 0.0})
            if p in net_scores:
                net = net_scores[p]
                geom_ic50 = math.sqrt(max(0.01, mf["ic50"]) * max(0.01, net["ic50"]))
                mean_rank = (mf["percentile_rank"] + net["percentile_rank"]) / 2.0
                scored_peptides[p] = {
                    "ic50": round(geom_ic50, 2),
                    "percentile_rank": round(mean_rank, 2),
                    "presentation_score": round(mf.get("presentation_score", 0.0), 4),
                    "mhcflurry_ic50": mf["ic50"],
                    "netmhcpan_ic50": net["ic50"],
                    "source": "Consensus_Ensemble(MHCflurry+NetMHCpan)",
                }
            else:
                # Obvious non-binder filtered before IEDB
                scored_peptides[p] = {
                    "ic50": mf["ic50"],
                    "percentile_rank": mf["percentile_rank"],
                    "presentation_score": mf.get("presentation_score", 0.0),
                    "mhcflurry_ic50": mf["ic50"],
                    "netmhcpan_ic50": 50000.0,
                    "source": "MHCflurry_Filter_NonBinder",
                }

        mf_pred.release_memory()
        net_pred.release_memory()
    else:
        # Direct scoring via single model (pwm, mhcflurry, or netmhcpan)
        scored_peptides = predictor.predict_peptides(all_peptides_list, allele=hla_allele)
        predictor.release_memory()

    # 4. Assemble EpitopeCandidate objects and find best epitope per variant
    for item in extracted_variants:
        m = item["mutation"]
        gene = item["gene"]
        pc = item["protein_change"]
        alt_reads = item["alt_reads"]
        tot_reads = item["total_reads"]
        af = item["af"]
        tpm = item["tpm"]

        candidates: List[EpitopeCandidate] = []
        for wt_kmer, mt_kmer, mut_pos in item["kmer_pairs"]:
            wt_info = scored_peptides.get(wt_kmer, {"ic50": 50000.0, "percentile_rank": 50.0, "presentation_score": 0.0})
            mt_info = scored_peptides.get(mt_kmer, {"ic50": 50000.0, "percentile_rank": 50.0, "presentation_score": 0.0})

            wt_ic50 = float(wt_info["ic50"])
            mt_ic50 = float(mt_info["ic50"])
            agretopicity = round(wt_ic50 / max(0.01, mt_ic50), 2)

            cand = EpitopeCandidate(
                gene=gene,
                variant_hgvsp=pc,
                hla_allele=hla_allele,
                kmer_length=9,
                mut_pos_in_kmer=mut_pos,
                wt_peptide=wt_kmer,
                mt_peptide=mt_kmer,
                wt_ic50_nm=wt_ic50,
                mt_ic50_nm=mt_ic50,
                wt_binder_class=classify_ic50(wt_ic50),
                mt_binder_class=classify_ic50(mt_ic50),
                agretopicity_index=agretopicity,
                percentile_rank=float(mt_info.get("percentile_rank", 50.0)),
                presentation_score=float(mt_info.get("presentation_score", 0.0)),
                mhcflurry_ic50=mt_info.get("mhcflurry_ic50"),
                netmhcpan_ic50=mt_info.get("netmhcpan_ic50"),
                predictor_source=str(mt_info.get("source", predictor_name)),
            )
            candidates.append(cand)

        if not candidates:
            continue

        candidates.sort(key=lambda c: c.mt_ic50_nm)
        best_candidate = candidates[0]

        stats["successfully_screened_variants"] += 1
        stats["total_9mers_evaluated"] += len(candidates)

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

    # Sort screened variants by lowest MT IC50 (strongest binding first)
    screened_results.sort(key=lambda r: r.best_epitope.mt_ic50_nm)
    return screened_results, stats
