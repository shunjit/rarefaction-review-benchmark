# Rarefaction Review Benchmark

[![License: MIT](https://img.shields.io/badge/License-MIT-yellow.svg)](https://opensource.org/licenses/MIT)

Reproducibility package for:

**"Rarefaction validity is task-dependent: justified for microbiome diversity, harmful for differential abundance, reportable via a decision framework"**

Submitted to *ISME Communications*

---

## Overview

This repository contains all scripts, metadata, and documentation required to reproduce the benchmark analyses presented in the manuscript. The analyses evaluate rarefaction practices for microbiome diversity analysis using three public 16S rRNA amplicon datasets.

### Key Findings

1. **Depth–group confounding diagnosis** using Kruskal–Wallis test with η² effect size provides an actionable framework for assessing whether rarefaction-based analyses may be compromised.

2. **Sensitivity analysis across multiple rarefaction depths** demonstrates that diversity conclusions can remain stable even under potential confounding, provided depth selection is conservative.

3. **Repeated rarefaction (q2-boots)** reveals that stochastic variability from single rarefaction is negligible (CV < 1%), supporting the use of single rarefaction for routine diversity analyses.

---

## Scope

This repository reproduces the **diversity-side benchmark** presented in the manuscript (Steps 0–5 of the decision framework). The differential abundance (DA) recommendations in the manuscript are based on published literature and do not involve novel benchmarking; accordingly, no DA analysis scripts are included here. For DA method implementations, see the original R packages: [ANCOM-BC2](https://github.com/FrederickHuangLin/ANCOMBC), [LinDA](https://github.com/zhouhj1994/LinDA), [MaAsLin2](https://github.com/biobakery/Maaslin2), and [ALDEx2](https://bioconductor.org/packages/ALDEx2/).

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
│       └── ena_fastq_download.sh  # ENA/NCBI/DDBJ universal FASTQ downloader
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
| QIIME 2 amplicon distribution | 2025.10 | Core analysis platform |
| q2cli | 2025.10.1 | QIIME 2 command-line interface |
| q2-boots | 2025.10.1 | Repeated rarefaction/bootstrapping |
| q2-composition | 2025.10.1 | ANCOM-BC2 (via QIIME 2 plugin) |
| Python | 3.10.14 | Scripting and visualization |

### System Requirements

- **Operating System**: macOS (Intel or Apple Silicon via Rosetta 2) or Linux
- **RAM**: 16 GB minimum, 32 GB recommended for Rice dataset (n=488)
- **Storage**: ~50 GB for raw FASTQ files and intermediate outputs

---

## Quick Start

### 1. Clone this repository

```bash
git clone https://github.com/shunjit/rarefaction-review-benchmark.git
cd rarefaction-review-benchmark
```

### 2. Install QIIME 2 environment

The analysis requires QIIME 2 amplicon distribution 2025.10. Install using the official specification file for your platform.

**Apple Silicon (M1/M2/M3/M4) Macs — via Rosetta 2 emulation:**

```bash
CONDA_SUBDIR=osx-64 conda env create \
  --name qiime2-amplicon-2025.10 \
  --file https://raw.githubusercontent.com/qiime2/distributions/refs/heads/dev/2025.10/amplicon/released/qiime2-amplicon-macos-latest-conda.yml
conda activate qiime2-amplicon-2025.10
conda config --env --set subdir osx-64
```

**Intel Macs:**

```bash
conda env create \
  --name qiime2-amplicon-2025.10 \
  --file https://raw.githubusercontent.com/qiime2/distributions/refs/heads/dev/2025.10/amplicon/released/qiime2-amplicon-macos-latest-conda.yml
conda activate qiime2-amplicon-2025.10
```

**Linux / WSL:**

```bash
conda env create \
  --name qiime2-amplicon-2025.10 \
  --file https://raw.githubusercontent.com/qiime2/distributions/refs/heads/dev/2025.10/amplicon/released/qiime2-amplicon-ubuntu-latest-conda.yml
conda activate qiime2-amplicon-2025.10
```

**Verify installation:**

```bash
qiime --version        # Expected: q2cli version 2025.10.1
qiime boots --help     # Should display q2-boots help
python --version       # Expected: Python 3.10.14
```

> **Note:** As of 2025.10, QIIME 2 does not provide native ARM64 builds. Apple Silicon Macs run QIIME 2 via Rosetta 2 emulation. The `CONDA_SUBDIR=osx-64` prefix and the `conda config --env --set subdir osx-64` command are both required.

### 3. Run the complete pipeline

```bash
# Option A: Run all steps sequentially (estimated 8–12 hours total)
bash scripts/run_all.sh

# Option B: Run in background with sleep prevention (recommended)
nohup caffeinate -ims bash scripts/run_all.sh > pipeline_log.txt 2>&1 &
```

### 4. Run individual steps

```bash
bash scripts/01_download_data.sh         # Download FASTQ files (~2–4 hours)
bash scripts/02_dada2_processing.sh      # DADA2 denoising (~2–3 hours)
bash scripts/03_build_phylogeny.sh       # Phylogenetic trees (~30–60 min)
bash scripts/04_confounding_diagnosis.sh # Depth–group confounding (~10 min)
bash scripts/05_alpha_rarefaction.sh     # Alpha-rarefaction curves (~30–60 min)
bash scripts/06_sensitivity_analysis.sh  # Multi-depth sensitivity (~2–3 hours)
bash scripts/07_q2boots_analysis.sh      # q2-boots repeated rarefaction (~2–3 hours)
```

---

## Datasets

Three public 16S rRNA amplicon datasets were selected from Schloss (2024) [1] to represent distinct regimes of sequencing depth, sample size, and confounding severity.

| Dataset | Accession | Samples | Primary Group | Confounding Status |
|---------|-----------|---------|---------------|-------------------|
| Soil | PRJEB10725 | 18 | treatment (k=3) | POTENTIAL (η²=0.28) |
| Bioethanol | PRJNA276052 | 95 | batch (k=19) / time_point (k=12) | NONE / POTENTIAL (η²=0.12) |
| Rice | PRJNA255789 | 488* | compartment (k=4) | MINOR (η²=0.02) |

\*Five samples excluded due to ENA paired-end file registration issues.

See `data/README.md` for detailed data acquisition instructions.

---

## Key Parameters

### Rarefaction Depths

| Dataset    | Primary Depth | Selection Criteria                         | Sample Retention |
| ---------- | ------------- | ------------------------------------------ | ---------------- |
| Soil       | 20,000        | Shannon plateau reached                    | 100% (18/18)     |
| Bioethanol | 35,000        | Shannon plateau behaviour; minimal dropout | 100% (95/95)     |
| Rice       | 20,000        | Shannon plateau reached                    | 99% (483/488)    |

### Confounding Diagnosis Thresholds

| Status | Criteria | Interpretation |
|--------|----------|---------------|
| NONE | Kruskal–Wallis p ≥ 0.05 | No evidence of depth–group association |
| MINOR | p < 0.05 and η² < 0.06 | Statistically significant but small effect |
| POTENTIAL | p < 0.05 and η² ≥ 0.06 | Meaningful confounding; sensitivity analysis required |

η² is calculated as: η² = (H − k + 1) / (N − k), where H is the Kruskal–Wallis statistic, k is the number of groups, and N is the total sample size.

---

## Reproducibility Notes

### QIIME 2 Provenance

All QIIME 2 artifacts (.qza) and visualizations (.qzv) contain embedded provenance tracking. Upload any .qzv file to [https://view.qiime2.org](https://view.qiime2.org) to inspect software versions, input file checksums, complete parameter specifications, and execution timestamps. See the `provenance/` directory for archived provenance from this study.

### Random Seeds

Where applicable, random seeds are fixed for reproducibility:

- DADA2: deterministic given identical input
- Rarefaction: controlled via `--p-seed` parameter
- q2-boots: controlled via `--p-random-state` parameter

---

## Citation

If you use this code or data, please cite:

```
Tokoro S. Rarefaction validity is task-dependent: justified for microbiome diversity, harmful for differential abundance, reportable via a decision framework. 
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

This project is licensed under the MIT License — see the [LICENSE](LICENSE) file for details.

---

## Contact

For questions or issues, please open a GitHub issue or contact:

**Shunji Tokoro**  
Gifu University of Medical Science  
ORCID: [0009-0003-8042-9869](https://orcid.org/0009-0003-8042-9869)
