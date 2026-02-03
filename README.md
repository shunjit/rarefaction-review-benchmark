# Rarefaction Review Benchmark

[![License: MIT](https://img.shields.io/badge/License-MIT-yellow.svg)](https://opensource.org/licenses/MIT)

Reproducibility package for:

**"Rarefaction in QIIME 2 workflows: when it is justified for diversity, when it is harmful for differential abundance, and how to report it reproducibly"**

Submitted to *ISME Communications*

---

## Overview

This repository contains all scripts, metadata, and documentation required to reproduce the benchmark analyses presented in the manuscript. The analyses evaluate rarefaction practices for microbiome diversity analysis using three public 16S rRNA amplicon datasets.

### Key Findings

1. **Depth–group confounding diagnosis** using Kruskal–Wallis test with η² effect size provides an actionable framework for assessing whether rarefaction-based analyses may be compromised.

2. **Sensitivity analysis across multiple rarefaction depths** demonstrates that diversity conclusions can remain stable even under potential confounding, provided depth selection is conservative.

3. **Repeated rarefaction (q2-boots)** reveals that stochastic variability from single rarefaction is negligible (CV < 1%), supporting the use of single rarefaction for routine diversity analyses.

---

## Repository Structure

```
rarefaction-review-benchmark/
├── README.md                      # This file
├── LICENSE                        # MIT License
├── CITATION.cff                   # Citation information
├── environment.yml                # Conda environment specification
│
├── data/
│   ├── README.md                  # Data acquisition instructions
│   └── metadata/
│       ├── metadata_soil.tsv      # Soil dataset metadata
│       ├── metadata_bioethanol.tsv # Bioethanol dataset metadata
│       └── metadata_rice.tsv      # Rice dataset metadata
│
├── scripts/
│   ├── run_all.sh                 # Master script (runs all steps)
│   ├── 01_download_data.sh        # Download FASTQ from ENA (wrapper)
│   ├── 02_dada2_processing.sh     # DADA2 denoising pipeline
│   ├── 03_build_phylogeny.sh      # Phylogenetic tree construction
│   ├── 04_confounding_diagnosis.sh # Depth–group confounding diagnosis
│   ├── 05_alpha_rarefaction.sh    # Alpha-rarefaction curve generation
│   ├── 06_sensitivity_analysis.sh # Multi-depth sensitivity analysis
│   ├── 07_q2boots_analysis.sh     # Repeated rarefaction with q2-boots
│   ├── depth_group_confounding_diagnosis.py  # Python diagnosis script
│   └── utils/
│       └── ena_fastq_download.sh  # ENA/NCBI/DDBJ Universal FASTQ Downloader
│
├── results/
│   ├── tables/
│   │   ├── Table_S1_sensitivity_analysis.md
│   │   └── Table_2_benchmark_summary.md
│   └── figures/
│       └── (generated during analysis)
│
└── provenance/
    └── README.md                  # QIIME 2 provenance documentation
```

---

## Requirements

### Software Dependencies

| Software | Version | Purpose |
|----------|---------|---------|
| QIIME 2 | 2025.10.1 (amplicon distribution) | Core analysis platform |
| q2-boots | 2025.10.1 | Repeated rarefaction/bootstrapping |
| Python | 3.10.14 | Scripting and visualization |

### System Requirements

- **Operating System**: macOS (Apple Silicon recommended) or Linux
- **RAM**: 16 GB minimum, 32 GB recommended for Rice dataset
- **Storage**: ~50 GB for raw FASTQ files and intermediate outputs

---

## Quick Start

### 1. Clone this repository

```bash
git clone https://github.com/shunjit/rarefaction-review-benchmark.git
cd rarefaction-review-benchmark
```

### 2. Install QIIME 2 environment

The analysis requires QIIME 2 amplicon distribution 2025.10. Install using the official specification file:

```bash
# For Apple Silicon (M1/M2/M3/M4) Macs:
wget https://data.qiime2.org/distro/amplicon/qiime2-amplicon-2025.10-py310-osx-conda-arm64.yml
conda env create -n qiime2-amplicon-2025.10 --file qiime2-amplicon-2025.10-py310-osx-conda-arm64.yml

# Activate the environment
conda activate qiime2-amplicon-2025.10

# Verify installation
qiime --version        # Expected: q2cli version 2025.10.1
qiime boots --help     # Should display q2-boots help
```

### 3. Run the complete pipeline

```bash
# Option A: Run all steps sequentially (8-12 hours total)
bash scripts/run_all.sh

# Option B: Run in background (recommended for full pipeline)
nohup caffeinate -ims bash scripts/run_all.sh > pipeline_log.txt 2>&1 &
```

### 4. Run individual steps (alternative)

If you prefer to run steps individually or resume from a specific point:

```bash
# Step 1: Download FASTQ files (~2-4 hours)
bash scripts/01_download_data.sh

# Step 2: DADA2 denoising (~2-3 hours)
bash scripts/02_dada2_processing.sh

# Step 3: Build phylogenetic trees (~30-60 minutes)
bash scripts/03_build_phylogeny.sh

# Step 4: Diagnose depth-group confounding (~10 minutes)
bash scripts/04_confounding_diagnosis.sh

# Step 5: Generate alpha-rarefaction curves (~30-60 minutes)
bash scripts/05_alpha_rarefaction.sh

# Step 6: Run multi-depth sensitivity analysis (~2-3 hours)
bash scripts/06_sensitivity_analysis.sh

# Step 7: Run q2-boots repeated rarefaction (~2-3 hours)
bash scripts/07_q2boots_analysis.sh
```

---

## Datasets

Three public 16S rRNA amplicon datasets were selected from Schloss (2024) to represent distinct regimes of sequencing depth, sample size, and confounding severity.

| Dataset | Accession | Samples | Primary Group | Confounding Status |
|---------|-----------|---------|---------------|-------------------|
| Soil | PRJEB10725 | 18 | treatment (3 groups) | POTENTIAL (η²=0.28) |
| Bioethanol | PRJNA276052 | 95 | batch (19) / time_point (12) | NONE / POTENTIAL |
| Rice | PRJNA255789 | 488* | compartment (4 groups) | MINOR (η²=0.02) |

*Five samples excluded due to incomplete paired-end data in ENA.

See `data/README.md` for detailed data acquisition instructions.

---

## Key Parameters

### DADA2 Processing

| Dataset | --trunc-len-f | --trunc-len-r | Rationale |
|---------|---------------|---------------|-----------|
| Soil | 133 | 133 | Quality drop at position 180+; conservative trimming |
| Bioethanol | 200 | 200 | V4 region; adequate overlap maintained |
| Rice | 133 | 133 | Quality degradation after position 180 |

### Primary Rarefaction Depths

| Dataset | Depth | Selection Criteria |
|---------|-------|-------------------|
| Soil | 20,000 | Plateau reached; 100% sample retention |
| Bioethanol | 35,000 | Near-plateau; 100% sample retention |
| Rice | 20,000 | Plateau reached; 99% sample retention |

### Confounding Diagnosis Thresholds

| Status | Criteria |
|--------|----------|
| NONE | Kruskal–Wallis p ≥ 0.05 |
| MINOR | p < 0.05 and η² < 0.06 |
| POTENTIAL | p < 0.05 and η² ≥ 0.06 |

---

## Reproducibility Notes

### QIIME 2 Provenance

All QIIME 2 artifacts (.qza) and visualizations (.qzv) contain embedded provenance tracking. Upload any .qzv file to [https://view.qiime2.org](https://view.qiime2.org) to inspect:

- Software versions and plugin information
- Input file checksums
- Complete parameter specifications
- Execution timestamps

### Random Seeds

Where applicable, random seeds are fixed for reproducibility:

- DADA2: default (deterministic given same input)
- Rarefaction: controlled via `--p-seed` parameter
- q2-boots: `--p-random-state` parameter

---

## Citation

If you use this code or data, please cite:

```
Tokoro S. Rarefaction in QIIME 2 workflows: when it is justified for diversity, 
when it is harmful for differential abundance, and how to report it reproducibly. 
ISME Communications. [Year];[Volume]:[Pages]. doi:[DOI]
```

See `CITATION.cff` for machine-readable citation information.

---

## References

1. Schloss PD. Rarefaction is currently the best approach to control for uneven sequencing effort in amplicon sequence analyses. *mSphere*. 2024;9:e00354-23.

2. McMurdie PJ, Holmes S. Waste not, want not: why rarefying microbiome data is inadmissible. *PLoS Comput Biol*. 2014;10:e1003531.

3. Raspet I, et al. Facilitating bootstrapped and rarefaction-based microbiome diversity analysis with q2-boots. *F1000Research*. 2025;14:87.

4. Bolyen E, et al. Reproducible, interactive, scalable and extensible microbiome data science using QIIME 2. *Nat Biotechnol*. 2019;37:852–857.

---

## License

This project is licensed under the MIT License - see the [LICENSE](LICENSE) file for details.

---

## Contact

For questions or issues, please open a GitHub issue or contact:

**Shunji Tokoro**  
Gifu University of Medical Science  
ORCID: [0009-0003-8042-9869](https://orcid.org/0009-0003-8042-9869)
