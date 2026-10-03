# Mini-NeoScorer

A lightweight, modular, educational in silico neoantigen scoring and prioritization pipeline built from scratch in Python.

## Overview
Mini-NeoScorer demonstrates the core computational immuno-oncology workflow for identifying and prioritizing immunogenic somatic mutations:
1. **Phase 1: Somatic Mutation & Transcriptome Groundwork** (VCF structure, annotation parsing, and TPM expression profiling).
2. **Phase 2: Coordinate & Sequence Translation** (Handling genomic vs coding coordinates, 0-based vs 1-based indexing, RNA-seq TPM expression gating, and wild-type vs mutant peptide extraction).
3. **Phase 3: Epitope Generation & HLA Class I Affinity Scoring** (9-mer sliding window, HLA-A*02:01 position weight matrix energy calculations, and agretopicity index determination).
4. **Phase 4: Integrated Immunogenicity Prioritization** (Multivariate scoring integrating clonality/VAF, TPM, agretopicity index, and presentation probability into a single clinical ranking).

## Repository Structure
```
Mini-NeoScorer/
├── data/
│   ├── mock_mutations.vcf
│   ├── mock_expression_tpm.tsv
│   └── reference_proteins.json
├── scripts/
│   └── generate_mock_data.py
├── src/
│   ├── __init__.py
│   ├── data_models.py
│   ├── reference_data.py
│   ├── phase2_pipeline.py
│   ├── hla_scorer.py
│   └── prioritizer.py
├── run_phase2_demo.py
├── run_phase3_demo.py
├── run_pipeline.py
├── .gitignore
└── README.md
```

## Running the Pipeline
Run with the shared Python virtual environment:
```powershell
& "C:\Tools\.venv\Scripts\python.exe" run_pipeline.py
```
Or launch via the unified launcher at:
`C:\Tools\LAUNCHERS\Launch_Mini-NeoScorer.bat`
