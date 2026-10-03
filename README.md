# Mini-NeoScorer

A lightweight, modular, educational in silico neoantigen scoring and prioritization pipeline built from scratch in Python.

## Overview
Mini-NeoScorer demonstrates the core computational immuno-oncology workflow for identifying and prioritizing immunogenic somatic mutations:
1. **Phase 1: Somatic Mutation & Transcriptome Groundwork** (VCF structure, annotation parsing, and TPM expression profiling).
2. **Phase 2: Coordinate & Sequence Translation** (Handling genomic vs coding coordinates, 0-based vs 1-based indexing, wild-type vs mutant peptide extraction).
3. **Phase 3: Epitope Generation & HLA Class I Affinity Scoring** (k-mer sliding window, MHC-I binding estimation, mutant vs self-antigen differential binding).
4. **Phase 4: Integrated Immunogenicity Prioritization** (Multivariate scoring integrating clonality/VAF, TPM, agretopicity index, and presentation probability).

## Repository Structure
```
Mini-NeoScorer/
├── data/
│   ├── mock_mutations.vcf
│   └── mock_expression_tpm.tsv
├── scripts/
│   └── generate_mock_data.py
├── .gitignore
└── README.md
```
