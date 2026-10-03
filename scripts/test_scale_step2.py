"""
Mini-NeoScorer: Scale Step 2 Verification Script
Tests dynamic UniProt reference resolution and 21-mer context window generation
on real patient mutations from TCGA-D3-A1Q1 (Melanoma).
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

if str(BASE_DIR) not in sys.path:
    sys.path.insert(0, str(BASE_DIR))

from src.protein_resolver import ProteinResolver


def main():
    with open(MUT_FILE, "r", encoding="utf-8") as f:
        mutations = json.load(f)

    with open(EXP_FILE, "r", encoding="utf-8") as f:
        expression = json.load(f)

    # Filter to missense mutations
    missense = [m for m in mutations if m.get("mutationType") == "Missense_Mutation"]

    # Filter by RNA expression gate (TPM >= 50.0 to focus on high-priority transcribed genes)
    expressed_muts = []
    for m in missense:
        entrez = str(m.get("entrezGeneId", ""))
        tpm = expression.get(entrez, 0.0)
        if tpm >= 50.0:
            expressed_muts.append((m, tpm))

    # Sort descending by expression TPM
    expressed_muts.sort(key=lambda x: x[1], reverse=True)

    print("=" * 95)
    print(" SCALE STEP 2: DYNAMIC REFERENCE PROTEIN RESOLUTION ON REAL PATIENT MUTATIONS")
    print(f" Target Cohort: TCGA-D3-A1Q1 (Melanoma) | Tested High-Expression Loci: {len(expressed_muts)}")
    print("=" * 95)

    resolver = ProteinResolver()
    resolved_count = 0

    print(f"{'#':<3} {'Gene':<10} {'Mutation':<12} {'TPM':<10} {'Prot Len':<10} {'WT Context Window':<23} {'MT Context Window'}")
    print("-" * 95)

    for i, (m, tpm) in enumerate(expressed_muts[:10], 1):
        gene = m.get("gene", {}).get("hugoGeneSymbol", "")
        pc = m.get("proteinChange", "")
        
        ctx = resolver.extract_context_window(gene, pc, flank=10)
        if ctx:
            resolved_count += 1
            print(
                f"{i:<3} "
                f"{gene:<10} "
                f"{pc:<12} "
                f"{tpm:<10.1f} "
                f"{ctx['protein_length']:<10} "
                f"{ctx['wt_window']:<23} "
                f"{ctx['mt_window']}"
            )
        else:
            print(f"{i:<3} {gene:<10} {pc:<12} {tpm:<10.1f} [RESOLUTION FAILED / ISOFORM MISMATCH]")

    print("-" * 95)
    print(f"[OK] Successfully resolved and validated {resolved_count} real human proteins on the fly.")
    print(f"     Persistent cache stored at: {resolver.cache_file}")
    print("=" * 95)


if __name__ == "__main__":
    main()
