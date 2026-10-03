"""
Mini-NeoScorer: Real Clinical Patient Data Ingestion
Phase 2.1: Automated ingestion of authentic metastatic melanoma somatic mutations
and matched RNA-seq transcript expression from the TCGA PanCancer Atlas (Cell, 2018).
"""

import json
import ssl
import sys
import urllib.request
from pathlib import Path
from typing import Dict, List, Tuple

# Ensure UTF-8 output
if hasattr(sys.stdout, "reconfigure"):
    sys.stdout.reconfigure(encoding="utf-8")

BASE_DIR = Path(__file__).resolve().parent.parent
DATA_DIR = BASE_DIR / "data"

# Clinical Patient ID from TCGA-SKCM PanCancer Atlas
PATIENT_SAMPLE_ID = "TCGA-D3-A1Q1-06"
STUDY_ID = "skcm_tcga_pan_can_atlas_2018"
MUTATION_PROFILE = f"{STUDY_ID}_mutations"
RNA_PROFILE = f"{STUDY_ID}_rna_seq_v2_mrna"

CBIOPORTAL_API_BASE = "https://www.cbioportal.org/api"


def _create_ssl_context() -> ssl.SSLContext:
    ctx = ssl.create_default_context()
    ctx.check_hostname = False
    ctx.verify_mode = ssl.CERT_NONE
    return ctx


def fetch_somatic_mutations(sample_id: str) -> List[dict]:
    """
    Fetch all verified somatic mutation calls (WES) for this patient.
    """
    url = f"{CBIOPORTAL_API_BASE}/molecular-profiles/{MUTATION_PROFILE}/mutations/fetch?projection=DETAILED"
    payload = json.dumps({"sampleIds": [sample_id]}).encode("utf-8")

    req = urllib.request.Request(
        url,
        data=payload,
        headers={"Content-Type": "application/json", "User-Agent": "MiniNeoScorer/2.0"},
    )
    with urllib.request.urlopen(req, context=_create_ssl_context()) as resp:
        return json.loads(resp.read().decode("utf-8"))


def fetch_matched_rna_expression(sample_id: str, entrez_gene_ids: List[int]) -> Dict[int, float]:
    """
    Fetch matched RNA-seq normalized expression for target mutated genes.
    """
    url = f"{CBIOPORTAL_API_BASE}/molecular-profiles/{RNA_PROFILE}/molecular-data/fetch"
    payload = json.dumps({"sampleIds": [sample_id], "entrezGeneIds": entrez_gene_ids}).encode("utf-8")

    req = urllib.request.Request(
        url,
        data=payload,
        headers={"Content-Type": "application/json", "User-Agent": "MiniNeoScorer/2.0"},
    )
    with urllib.request.urlopen(req, context=_create_ssl_context()) as resp:
        records = json.loads(resp.read().decode("utf-8"))
        return {r["entrezGeneId"]: float(r.get("value", 0.0)) for r in records if "value" in r}


def ingest_real_patient_data():
    DATA_DIR.mkdir(parents=True, exist_ok=True)
    print(f">>> Fetching somatic mutations for patient {PATIENT_SAMPLE_ID} from {STUDY_ID}...")
    raw_mutations = fetch_somatic_mutations(PATIENT_SAMPLE_ID)

    # Filter to missense mutations (the primary source of single-amino-acid neoantigens)
    missense_muts = [m for m in raw_mutations if m.get("mutationType") == "Missense_Mutation"]
    print(f"  * Total somatic variants: {len(raw_mutations)}")
    print(f"  * Missense SNVs: {len(missense_muts)}")

    # Collect Entrez IDs for RNA expression lookup
    entrez_ids = list({m.get("entrezGeneId") for m in missense_muts if m.get("entrezGeneId")})
    print(f">>> Fetching matched RNA-seq expression for {len(entrez_ids)} mutated genes...")
    rna_map = fetch_matched_rna_expression(PATIENT_SAMPLE_ID, entrez_ids)
    print(f"  * Matched RNA-seq values retrieved for {len(rna_map)} genes.")

    # Save real dataset locally
    mutations_file = DATA_DIR / f"real_patient_{PATIENT_SAMPLE_ID}_mutations.json"
    expression_file = DATA_DIR / f"real_patient_{PATIENT_SAMPLE_ID}_expression.json"

    with open(mutations_file, "w", encoding="utf-8") as f:
        json.dump(raw_mutations, f, indent=2)

    with open(expression_file, "w", encoding="utf-8") as f:
        json.dump(rna_map, f, indent=2)

    print(f"\n[OK] Real patient data successfully ingested and cached:")
    print(f"  -> Mutations:  {mutations_file} ({mutations_file.stat().st_size / 1024:.1f} KB)")
    print(f"  -> Expression: {expression_file} ({expression_file.stat().st_size / 1024:.1f} KB)")


if __name__ == "__main__":
    ingest_real_patient_data()
