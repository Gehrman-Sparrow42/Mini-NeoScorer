# Neoantigen-Prioritization-Pipeline

[![Python 3.10+](https://img.shields.io/badge/python-3.10%2B-blue.svg)](https://www.python.org/)
[![License: MIT](https://img.shields.io/badge/License-MIT-green.svg)](https://opensource.org/licenses/MIT)
[![Field: Computational Immuno-Oncology](https://img.shields.io/badge/Field-Computational%20Immuno--Oncology-purple.svg)]()
[![Data Source: TCGA PanCancer](https://img.shields.io/badge/Data%20Source-TCGA%20PanCancer%20Atlas-orange.svg)](https://gdc.cancer.gov/about-data/publications/pancanatlas)
[![Tests: 9 Passed](https://img.shields.io/badge/tests-9%20passed-success.svg)]()

A clinical-grade, modular, high-throughput computational pipeline and **interactive clinical dashboard for in silico neoantigen discovery, Deep Neural Network (DNN) HLA Class I presentation modeling, and personalized cancer vaccine payload prioritization**.

Bridges raw patient multi-omics sequencing (Whole-Exome Sequencing + RNA-seq) directly to clinically actionable, ranked mRNA/peptide vaccine formulations.

---

## Table of Contents
- [Abstract](#abstract)
- [Biological Principles](#biological-principles)
- [Pipeline Architecture & Workflow](#pipeline-architecture--workflow)
- [Mathematical Scoring Model](#mathematical-scoring-model)
- [Clinical Case Study: Metastatic Melanoma (TCGA-D3-A1Q1)](#clinical-case-study-metastatic-melanoma-tcga-d3-a1q1)
- [Interactive Multi-Patient Dashboard](#interactive-multi-patient-dashboard)
- [Repository Structure](#repository-structure)
- [Quickstart & CLI Usage](#quickstart--cli-usage)
- [Automated Testing](#automated-testing)
- [Data Provenance & Citations](#data-provenance--citations)

---

## Abstract

Personalized cancer vaccines (such as patient-specific mRNA or synthetic long peptide formulations) require selecting a restricted payload (typically 10 to 20 candidate epitopes) from hundreds of somatic mutations identified via next-generation sequencing. 

The **Neoantigen-Prioritization-Pipeline** automates this selection by systematically integrating four orthogonal biological dimensions:
1. **Somatic Clonality**: Variant Allele Frequency (VAF) from tumor-normal Whole-Exome Sequencing (WES).
2. **Transcriptional Availability**: Normalized Transcripts Per Million (TPM) from matched patient RNA-seq.
3. **MHC-I Presentation Stability**: Dual-model consensus deep learning (**MHCflurry 2.0** + **NetMHCpan-4.1**) modeling stereochemical groove binding affinity ($IC_{50}\text{ nM}$), antigen processing/cleavage, and Monte Carlo Percentile Ranks ($Rank\%$).
4. **Agretopicity (Differential Foreignness)**: Mutant versus wild-type differential binding affinity ($AI = IC_{50}^{WT} / IC_{50}^{MT}$) to maximize immunogenicity while minimizing self-tolerance and autoimmune risk.

---

## Biological Principles

```
  TUMOR CELL (Intracellular)               CELL SURFACE                    IMMUNE SYSTEM
┌────────────────────────────┐        ┌──────────────────────┐        ┌──────────────────────┐
│ DNA: Somatic Point Typo    │        │  HLA Class I Groove  │        │  CD8+ Cytotoxic      │
│  └─> mRNA Transcription    │ ──────>│  (Pockets B & F)     │ ──────>│  T-Cell Receptor     │
│  └─> Proteasome Cleavage   │        │  Cradles 9-mer       │        │  (TCR Recognition)   │
│  └─> Neopeptide Generation │        │  Mutant Epitope      │        │  Triggers Cell Death │
└────────────────────────────┘        └──────────────────────┘        └──────────────────────┘
```

1. **The Transcriptional Gate**: Somatic mutations frequently occur in epigenetically silenced or non-expressed genomic regions. A mutation cannot be translated, cleaved, or loaded onto HLA if its source transcript has zero expression ($TPM < 1.0$).
2. **Anchor Pockets on HLA-A\*02:01**: The Class I closed peptide-binding groove predominantly accommodates 9-amino-acid peptides (9-mers). Primary anchor residues at **Position 2** (Pocket B, favoring Leucine/Methionine/Valine) and **Position 9** (Pocket F, favoring Valine/Leucine) dictate physical presentation stability.
3. **Clonality and Tumor Heterogeneity**: Vaccines targeting subclonal mutations ($AF \ll 0.20$) allow unmutated tumor clones to escape. Prioritizing clonal driver mutations ($AF \to 0.50$ in diploid tumors) ensures tumor-wide immune eradication.
4. **Consensus Deep Learning**: Single neural networks are prone to model-specific false positives. Combining **MHCflurry 2.0** (antigen processing and TAP transport) with **NetMHCpan-4.1** (pan-allele surface presentation) ensures robust dual-validated candidates.

---

## Pipeline Architecture & Workflow

```text
┌──────────────────────────────────────────────────────────────────────────────────┐
│ STAGE 1: Multi-Omics Patient Data Ingestion                                      │
│  * WES Somatic Calls: Mutation types, read counts (Alt/Total), VAF calculation   │
│  * Matched RNA-seq: Gene-level RSEM / TPM abundance quantification                │
└────────────────────────────────────────┬─────────────────────────────────────────┘
                                         │
                                         ▼
┌──────────────────────────────────────────────────────────────────────────────────┐
│ STAGE 2: RNA-seq Expression Gating                                               │
│  * Filter: Drop silent/unexpressed loci (Threshold: TPM >= 1.0)                  │
│  * Reduces downstream search space by 60% - 80%                                  │
└────────────────────────────────────────┬─────────────────────────────────────────┘
                                         │
                                         ▼
┌──────────────────────────────────────────────────────────────────────────────────┐
│ STAGE 3: Dynamic UniProt Proteome Resolution & Context Extraction                 │
│  * 16-worker multi-threaded prefetching from UniProtKB Swiss-Prot                │
│  * Translates 1-based clinical coordinates to 0-based sequence indices           │
│  * Strict reference residue verification (rejects alternative splice mismatches) │
│  * Extracts 21-mer flanking windows (10 aa upstream + MUT + 10 aa downstream)    │
└────────────────────────────────────────┬─────────────────────────────────────────┘
                                         │
                                         ▼
┌──────────────────────────────────────────────────────────────────────────────────┐
│ STAGE 4: High-Throughput 9-mer Slicing & Deep Neural Network Ensemble            │
│  * Slices all overlapping 9-mers spanning the somatic alteration                 │
│  * Local PyTorch MHCflurry 2.0 (groove affinity + antigen processing/cleavage)    │
│  * Official NIH IEDB NetMHCpan-4.1 (pan-allele neural network + Percentile Rank) │
│  * Geometric consensus IC50 & Arithmetic consensus Percentile Rank calculation   │
└────────────────────────────────────────┬─────────────────────────────────────────┘
                                         │
                                         ▼
┌──────────────────────────────────────────────────────────────────────────────────┐
│ STAGE 5: Multi-Parametric Prioritization & Clinical Payload Formulation          │
│  * Multiplies Clonality * Expression * Presentation * Agretopicity Bonus         │
│  * Assigns Clinical Actionability Tiers (Tier 1: Vaccine, Tier 2: Backup)        │
│  * Exports verified Clinical Prescription Table (TSV) & JSON metadata sidecar    │
│  * Interactive Clinical Web Dashboard & Visual Analytics (FastAPI + Chart.js)   │
└──────────────────────────────────────────────────────────────────────────────────┘
```

---

## Mathematical Scoring Model

The final **Immunogenicity Prioritization Score** ($S_{neo}$, scaled $0 - 100$) integrates four multiplicative biological terms:

$$S_{neo} = C_{AF} \times E_{TPM} \times P_{MHC} \times A_{bonus} \times 20.0$$

### 1. Clonality Score ($C_{AF}$)
In diploid human tumor genomes, a heterozygous somatic alteration present in 100% of malignant cells yields an allele frequency $AF \approx 0.50$:
$$C_{AF} = \min(1.0, \max(0.0, 2.0 \times AF))$$

### 2. Transcriptional Abundance Score ($E_{TPM}$)
Logarithmic transformation reflecting saturating biological antigen availability:
$$E_{TPM} = \log_2(1.0 + \text{TPM})$$

### 3. MHC-I Presentation Probability ($P_{MHC}$)
Calibrated transformation of consensus $IC_{50}$ dissociation constant (nM) normalized to the standard 50,000 nM ceiling:
$$P_{MHC} = 1.0 - \frac{\log_{10}(\text{clamped } IC_{50})}{\log_{10}(50000)}$$
* $IC_{50} \le 50\text{ nM}$ (Strong Binder) $\to P_{MHC} \ge 0.65$
* $IC_{50} \approx 500\text{ nM}$ (Weak Binder threshold) $\to P_{MHC} \approx 0.43$
* $IC_{50} > 10,000\text{ nM}$ (Non-Binder) $\to P_{MHC} < 0.15$

### 4. Agretopicity Foreignness Factor ($A_{bonus}$)
Evaluates the ratio of wild-type to mutant binding affinity ($AI = IC_{50}^{WT} / IC_{50}^{MT}$):
* If $AI \ge 1.0$ (Gain-of-Binding mutation):
  $$A_{bonus} = \min(1.50, 1.0 + 0.15 \times \log_2(AI + 1.0))$$
* If $AI < 1.0$ (Loss of affinity compared to self):
  $$A_{bonus} = \max(0.75, 0.90 + 0.10 \times AI)$$

### 5. Consensus Ensemble Equations
$$\text{Ensemble IC}_{50} = \sqrt{\text{IC}_{50}^{\text{MHCflurry}} \times \text{IC}_{50}^{\text{NetMHCpan}}}$$
$$\text{Consensus Rank \%} = \frac{\text{Rank}_{\text{MHCflurry}} + \text{Rank}_{\text{NetMHCpan}}}{2}$$

---

## Clinical Case Study: Metastatic Melanoma (TCGA-D3-A1Q1)

The pipeline was validated against authentic clinical patient data from **The Cancer Genome Atlas (TCGA) Skin Cutaneous Melanoma cohort** (*Cell*, 2018):

### Exome Screening Funnel
| Metric | Count | Biological Significance |
| :--- | :--- | :--- |
| **Total Somatic Mutations (WES)** | **486** | High tumor mutation burden typical of UV-induced melanoma |
| **Missense Alterations** | **410** | Non-synonymous single-amino-acid changes |
| **Eliminated by RNA Gate ($TPM < 1.0$)** | **66** | Transcriptionally silent genes (e.g. *ALPP*, *CXCR5*) safely dropped |
| **Filtered by Isoform / QC boundary** | **18** | Alternative splice discrepancies safely flagged |
| **Candidate 9-mers Screened** | **2,912** | Overlapping peptide windows evaluated against HLA-A\*02:01 |
| **Strong Binders ($IC_{50} \le 50\text{ nM}$ or $Rank \le 0.5\%$)** | **28** | Dual-consensus high-affinity surface presentation candidates |
| **Weak Binders ($50 - 500\text{ nM}$)** | **51** | Viable secondary/backup presentation pool |

### Top Recommended Vaccine Payload (Consensus Ensemble)
```text
Rank  Gene       Mutation   Neopeptide   AF     RNA TPM   IC50 (nM)   Rank %   Agret.    Score   Clinical Tier
-----------------------------------------------------------------------------------------------------------------------------
#1    AP4B1      G525V      LLLVVIDEV    0.65   325.3     17.8        0.14%    0.81x     100.0   Tier 1: Vaccine Payload
#2    TNFRSF10B  R145L      FLEEDSPEM    0.20   879.7     24.0        0.10%    340.78x   82.9    Tier 1: Vaccine Payload
#3    TMEM245    G754W      FLWTYWAAV    0.15   1634.3    4.3         0.04%    1.25x     64.3    Tier 1: Vaccine Payload
#4    PCNX3      L1607I     SLEPFIYGL    0.18   1194.8    22.8        0.16%    1.07x     61.2    Tier 1: Vaccine Payload
#5    GPR158     L535F      RMLAVILFV    0.13   2328.6    5.1         0.02%    1.29x     59.8    Tier 1: Vaccine Payload
#6    PML        S388Y      RLQDLSYCI    0.12   1176.6    18.0        0.12%    7.28x     52.3    Tier 1: Vaccine Payload
```

> **Key Biological Insights from Deep Learning**:
> - **`TNFRSF10B R145L` (`FLEEDSPEM`)**: Introduces a canonical Leucine anchor at Pocket B (Position 2), converting an unpresentable wild-type ($IC_{50} = 8,185\text{ nM}$) into a high-affinity neoepitope ($IC_{50} = 24.0\text{ nM}$)—a remarkable **340.78× Agretopicity gain-of-affinity**.
> - **Elimination of Matrix Artifacts (`NRAS G12R`)**: Simple position weight matrices falsely favor Arginine at position 7. The Deep Neural Network correctly recognized that the bulky, positively charged Arginine disrupts hydrophobic groove packing ($IC_{50} = 4,599\text{ nM}$), safely preventing an unpresentable decoy from entering the vaccine payload.

---

## Interactive Multi-Patient Dashboard

A production-grade, local web dashboard (FastAPI + Vanilla JS + Chart.js) is included for interactive exploration, multi-patient cohort management, and clinical formulation review.

```text
┌────────────────────────────────────────────────────────────────────────────────────────┐
│  NeoScorer · PRECISION IMMUNO-ONCOLOGY                   [ Active Patient: TCGA-D3-A1Q1 ▼ ] ↻ │
├────────────────────────────────────────────────────────────────────────────────────────┤
│  [Alterations: 326]    [Strong Binders: 28]    [Tier 1 Payload: 6]    [Secondary: 35]  │
├────────────────────────────────────────────────────┬───────────────────────────────────┤
│  EXPRESSION × PRESENTATION LANDSCAPE               │  MUTANT / WILD-TYPE AGRETOPICITY  │
│  (Interactive Bubble Scatter: TPM vs IC50 vs VAF)  │  (Log-scale bar chart: MT vs WT)  │
├────────────────────────────────────────────────────┴───────────────────────────────────┤
│  CLINICAL VACCINE PRESCRIPTION EXPLORER                  [⌕ Search] [Tier: All ▼] [↓ Export] │
│  Rank | Gene | Mutation | Neopeptide | VAF | TPM | IC50 (nM) | Rank% | Agret. | Tier     │
├────────────────────────────────────────────────────────────────────────────────────────┤
│  PREDICTION CONSENSUS EVIDENCE CARD (Selected Candidate)                                │
│  MHCflurry 2.0 (Groove + Cleavage) ↔ NetMHCpan-4.1 (IEDB Pan-Allele Presentation)       │
└────────────────────────────────────────────────────────────────────────────────────────┘
```

### Key Features:
- **Dynamic Patient Discovery**: Automatically scans `reports/` for all patient prescription reports; seamlessly switch between patients via the navbar dropdown.
- **Interactive Visualizations**: Four-quadrant immuno-oncology scatter plots and logarithmic agretopicity bar charts (bundled locally with zero CDN dependencies).
- **Candidate Filtering & Multi-Level Sorting**: Instant search by gene, mutation, or peptide sequence; filter by clinical tier; Shift-click table headers for multi-level sorting.
- **Model Consensus Audit**: Inspect individual MHCflurry and NetMHCpan predictions for any selected candidate.
- **"Run New Patient" Interface**: Upload or select patient WES/RNA-seq JSONs and execute the background pipeline with streaming console logs.
- **One-Click Export**: Download filtered prescriptions as clinical CSV or TSV files.

---

## Repository Structure

```text
Neoantigen-Prioritization-Pipeline/
├── data/                                    # Multi-omics inputs & reference caches
│   ├── real_patient_TCGA-D3-A1Q1-06_mutations.json   # 486 WES somatic variant calls
│   ├── real_patient_TCGA-D3-A1Q1-06_expression.json  # Matched RNA-seq TPM values
│   ├── protein_cache.json                   # Local cache of 305+ UniProt protein sequences
│   └── netmhcpan_cache.json                 # Persistent cache of NIH IEDB NetMHCpan queries
├── reports/                                 # Clinical outputs
│   ├── TCGA_D3_A1Q1_clinical_vaccine_prescription.tsv # Formatted prescription table
│   └── TCGA_D3_A1Q1_clinical_vaccine_prescription.json # Patient sidecar metadata
├── dashboard_backend/                       # Local dashboard server modules
│   ├── __init__.py
│   ├── jobs.py                              # Asynchronous pipeline execution & queueing
│   └── reports.py                           # Dynamic patient discovery, parsing & filtering
├── dashboard_static/                        # Self-contained frontend assets
│   ├── index.html                           # Accessible semantic dashboard UI
│   ├── styles.css                           # Modern slate/emerald dark theme tokens
│   ├── app.js                               # Patient state management & grid controller
│   ├── charts.js                            # Interactive Chart.js analytics
│   ├── data.js                              # Client API bindings
│   └── vendor/chart.umd.js                  # Bundled offline Chart.js v4.4.8
├── scripts/                                 # Auxiliary data utilities
│   ├── fetch_patient_data.py                # Ingestion client for cBioPortal REST API
│   └── inspect_patient.py                   # Multi-omics data profiler
├── src/                                     # Core scientific engine
│   ├── __init__.py
│   ├── data_models.py                       # Dataclasses and coordinate conversions
│   ├── dnn_predictor.py                     # Deep Neural Network engine (MHCflurry + NetMHCpan)
│   ├── protein_resolver.py                  # Multi-threaded UniProt REST API sequence resolver
│   ├── hla_scorer.py                        # 9-mer sliding window & HLA-A*02:01 PWM model
│   ├── prioritizer.py                       # Multi-parametric scoring mathematics
│   └── scale_hla_screener.py                # Exome screening funnel orchestrator
├── tests/                                   # Automated regression test suite
│   └── test_dashboard.py                    # Pytest suite for API, jobs, and exports
├── dashboard.py                             # FastAPI dashboard entrypoint (loopback-only)
├── patient_io.py                            # Shared input validation & file naming logic
├── run_pipeline.py                          # Master CLI pipeline entrypoint
├── requirements.txt                         # Runtime dependencies (torch, mhcflurry, fastapi)
├── requirements-dev.txt                     # Test dependencies (pytest)
├── DASHBOARD.md                             # Technical architecture documentation for dashboard
├── .gitignore
└── README.md
```

---

## Quickstart & CLI Usage

### Prerequisites
* Python 3.10+ (Tested on Python 3.13 Windows 64-bit AMD64).
* Install dependencies:
```bash
pip install -r requirements.txt
```

### 1. Launch the Interactive Web Dashboard
Run the dashboard and automatically open it in your browser (`http://localhost:8520`):
```bash
python dashboard.py --open-browser
```
*(Alternatively, on Windows, double-click `C:\Tools\LAUNCHERS\Launch_NeoScorer_Dashboard.bat`)*

### 2. Run the Command-Line Pipeline Directly
Execute the full clinical exome prioritization pipeline with the Consensus Ensemble:
```bash
python run_pipeline.py --predictor ensemble
```

Customize clinical parameters via CLI arguments:
```bash
python run_pipeline.py \
  --mutations data/real_patient_TCGA-D3-A1Q1-06_mutations.json \
  --expression data/real_patient_TCGA-D3-A1Q1-06_expression.json \
  --hla HLA-A*02:01 \
  --predictor ensemble \
  --tpm-threshold 1.0 \
  --top 12
```

### CLI Options:
* `--predictor`: HLA presentation engine:
  * `ensemble` (Default): Consensus combining local PyTorch **MHCflurry 2.0** and official NIH IEDB **NetMHCpan-4.1** via geometric mean IC50 and mean percentile rank.
  * `mhcflurry`: Local deep neural network modeling groove affinity and antigen processing (proteasome cleavage + TAP transport).
  * `netmhcpan`: Official NIH IEDB REST API for NetMHCpan-4.1 pan-allele artificial neural network.
  * `pwm`: Ultra-fast position weight matrix baseline.
* `--hla`: Target patient HLA Class I allele (default: `HLA-A*02:01`).
* `--tpm-threshold`: Minimum RNA expression cutoff gate in TPM (default: `1.0`).
* `--top`: Number of top-ranked vaccine candidates to display in console log (default: `12`).
* `--mutations`: Path to input somatic mutations JSON.
* `--expression`: Path to input RNA-seq TPM JSON.
* `--patient-id`: Optional patient identifier (inferred from mutation file if omitted).
* `--output-dir`: Output directory for generated prescriptions (default: `reports/`).

---

## Automated Testing

Execute the comprehensive test suite to verify patient isolation, data validation, exports, and background workers:

```bash
pytest tests/
```

---

## Data Provenance & Citations

1. **TCGA PanCancer Atlas**:
   * Hoadley, K. A. et al. *Cell-of-Origin Patterns Dominate the Molecular Classification of 10,000 Tumors from 33 Types of Cancer*. **Cell** 173, 291–304.e6 (2018). [doi:10.1016/j.cell.2018.03.022](https://doi.org/10.1016/j.cell.2018.03.022).
2. **MHCflurry 2.0**:
   * O'Donnell, T. J. et al. *MHCflurry 2.0: Improved Pan-Allele Prediction of MHC Class I-Presented Peptides by Incorporating Antigen Processing*. **Cell Systems** 11, 42–48.e7 (2020). [doi:10.1016/j.cels.2020.06.010](https://doi.org/10.1016/j.cels.2020.06.010).
3. **NetMHCpan-4.1 (IEDB)**:
   * Reynisson, B. et al. *NetMHCpan-4.1 and NetMHCIIpan-4.0: Improved predictions of MHC antigen presentation by concurrent motif deconvolution and integration of MS eluted ligand data*. **Nucleic Acids Research** 48, W449–W454 (2020). [doi:10.1093/nar/gkaa379](https://doi.org/10.1093/nar/gkaa379).
4. **cBioPortal for Cancer Genomics**:
   * Cerami, E. et al. *The cBio Cancer Genomics Portal: An Open Platform for Exploring Multidimensional Cancer Genomics Data*. **Cancer Discovery** 2, 401–404 (2012). [doi:10.1158/2159-8290.CD-12-0095](https://doi.org/10.1158/2159-8290.CD-12-0095).
5. **UniProtKB / Swiss-Prot**:
   * The UniProt Consortium. *UniProt: the Universal Protein Knowledgebase in 2023*. **Nucleic Acids Research** 51, D523–D531 (2023). [doi:10.1093/nar/gkac1052](https://doi.org/10.1093/nar/gkac1052).
