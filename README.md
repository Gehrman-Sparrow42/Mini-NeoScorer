# Mini-NeoScorer

A lightweight, modular, clinical-grade in silico neoantigen scoring and prioritization pipeline built from scratch in Python.

## Overview
Mini-NeoScorer implements the complete computational immuno-oncology workflow for identifying, evaluating, and prioritizing somatic neoantigens for personalized cancer vaccines:

### Module 1: Educational Prototype (Phases 1–4)
* **Phase 1**: Somatic mutation (VCF) and gene expression (TPM) ground truth.
* **Phase 2**: 1-based biological to 0-based sequence indexing, RNA expression gate ($TPM \ge 1.0$), and 21-mer flanking peptide context extraction.
* **Phase 3**: 9-mer sliding window generation, HLA-A*02:01 position weight matrix affinity prediction, and agretopicity index calculation.
* **Phase 4**: Multi-parametric prioritization scoring unifying clonality ($AF$), expression ($TPM$), presentation ($IC_{50}$), and foreignness.

### Module 2: Real-World Clinical Scale (Scale Steps 1–4)
* **Scale Step 1**: Ingestion of authentic multi-omics data from a metastatic melanoma patient (**TCGA-D3-A1Q1**, TCGA PanCancer Atlas, *Cell* 2018) containing 486 somatic variants and matched RNA-seq expression for 395 genes via cBioPortal API.
* **Scale Step 2**: Automated reference protein sequence resolution via the UniProtKB REST API with multi-threaded prefetching and persistent local disk caching.
* **Scale Step 3**: High-throughput screening across 2,912 candidate 9-mers in the patient's tumor exome.
* **Scale Step 4**: Clinical vaccine payload formulation, identifying Tier 1 targets (including clonal oncogenic driver **NRAS G12R**, **AP4B1 G525V**, and **TRIM28 R492C**) and exporting a clinical TSV prescription report.

## Repository Structure
```
Mini-NeoScorer/
├── data/
│   ├── mock_mutations.vcf
│   ├── mock_expression_tpm.tsv
│   ├── real_patient_TCGA-D3-A1Q1-06_mutations.json
│   ├── real_patient_TCGA-D3-A1Q1-06_expression.json
│   └── protein_cache.json
├── reports/
│   └── TCGA_D3_A1Q1_clinical_vaccine_prescription.tsv
├── scripts/
│   ├── generate_mock_data.py
│   ├── fetch_real_melanoma_patient.py
│   ├── inspect_real_patient.py
│   ├── test_scale_step2.py
│   ├── run_scale_step3.py
│   └── run_scale_step4_clinical_report.py
├── src/
│   ├── __init__.py
│   ├── data_models.py
│   ├── reference_data.py
│   ├── phase2_pipeline.py
│   ├── hla_scorer.py
│   ├── prioritizer.py
│   ├── protein_resolver.py
│   └── scale_hla_screener.py
├── run_phase2_demo.py
├── run_phase3_demo.py
├── run_pipeline.py
├── .gitignore
└── README.md
```

## Running the Pipeline
Run the real clinical melanoma pipeline:
```powershell
& "C:\Tools\.venv\Scripts\python.exe" scripts/run_scale_step4_clinical_report.py
```
Or run the educational prototype:
```powershell
& "C:\Tools\.venv\Scripts\python.exe" run_pipeline.py
```
Or double-click the unified desktop launcher at:
`C:\Tools\LAUNCHERS\Launch_Mini-NeoScorer.bat`
