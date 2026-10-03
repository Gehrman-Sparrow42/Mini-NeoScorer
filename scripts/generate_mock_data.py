"""
Mini-NeoScorer: Synthetic Dataset Generator
Phase 1: Generation of mock somatic mutations (VCF) and gene expression (TPM).
"""

from pathlib import Path

# Paths
BASE_DIR = Path(__file__).resolve().parent.parent
DATA_DIR = BASE_DIR / "data"

VCF_CONTENT = """##fileformat=VCFv4.2
##fileDate=20261003
##source=MiniNeoScorerSyntheticVCF
##reference=GRCh38
##INFO=<ID=GENE,Number=1,Type=String,Description="HUGO Gene Symbol">
##INFO=<ID=HGVSp,Number=1,Type=String,Description="HGVS protein sequence alteration">
##INFO=<ID=HGVSc,Number=1,Type=String,Description="HGVS coding sequence alteration">
##INFO=<ID=AF,Number=A,Type=Float,Description="Tumor Variant Allele Frequency (0.0 to 1.0)">
##INFO=<ID=DP,Number=1,Type=Integer,Description="Total read depth at variant locus">
##INFO=<ID=TIER,Number=1,Type=String,Description="Somatic classification tier (Clonal_Driver, Subclonal, Passenger, Passenger_Silenced)">
##FORMAT=<ID=GT,Number=1,Type=String,Description="Genotype call">
##FORMAT=<ID=AD,Number=R,Type=Integer,Description="Allelic depths for REF and ALT alleles">
##FORMAT=<ID=DP,Number=1,Type=Integer,Description="Approximate read depth">
#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\tFORMAT\tTUMOR
chr7\t140753336\t.\tT\tA\t1500\tPASS\tGENE=BRAF;HGVSp=p.Val600Glu;HGVSc=c.1799T>A;AF=0.48;DP=110;TIER=Clonal_Driver\tGT:AD:DP\t0/1:57,53:110
chr12\t25245350\t.\tC\tT\t1200\tPASS\tGENE=KRAS;HGVSp=p.Gly12Asp;HGVSc=c.35G>A;AF=0.42;DP=95;TIER=Clonal_Driver\tGT:AD:DP\t0/1:55,40:95
chr3\t179218303\t.\tG\tA\t350\tPASS\tGENE=PIK3CA;HGVSp=p.Glu545Lys;HGVSc=c.1633G>A;AF=0.07;DP=120;TIER=Subclonal\tGT:AD:DP\t0/1:112,8:120
chr2\t178550100\t.\tC\tT\t850\tPASS\tGENE=TTN;HGVSp=p.Val1234Leu;HGVSc=c.3700G>T;AF=0.36;DP=80;TIER=Passenger\tGT:AD:DP\t0/1:51,29:80
chr9\t36993350\t.\tC\tT\t920\tPASS\tGENE=PAX5;HGVSp=p.Pro80Leu;HGVSc=c.239C>T;AF=0.45;DP=88;TIER=Passenger_Silenced\tGT:AD:DP\t0/1:48,40:88
"""

TPM_CONTENT = """gene_symbol\tTPM
BRAF\t54.30
KRAS\t78.10
PIK3CA\t28.60
TTN\t1.40
PAX5\t0.00
"""


def generate_mock_datasets() -> None:
    DATA_DIR.mkdir(parents=True, exist_ok=True)

    vcf_path = DATA_DIR / "mock_mutations.vcf"
    tpm_path = DATA_DIR / "mock_expression_tpm.tsv"

    vcf_path.write_text(VCF_CONTENT, encoding="utf-8")
    tpm_path.write_text(TPM_CONTENT, encoding="utf-8")

    print(f"Generated synthetic VCF at: {vcf_path}")
    print(f"Generated synthetic TPM table at: {tpm_path}")


if __name__ == "__main__":
    generate_mock_datasets()
