# Rarefaction Review Benchmark

[![License: MIT](https://img.shields.io/badge/License-MIT-yellow.svg)](https://opensource.org/licenses/MIT)

Reproducibility package for:

**"Rarefaction validity is task-dependent: justified for microbiome diversity, harmful for differential abundance, reportable via a decision framework"**

Submitted to *ISME Communications*

---

## Overview

This repository contains the scripts, frozen run configurations, summary outputs, metadata, and documentation required to reproduce the benchmark analyses presented in the manuscript:

1. the **diversity-side benchmark** (Steps 0–5 of the decision framework) on three public 16S rRNA amplicon datasets, and
2. the **depth-confounding differential abundance (DA) benchmark** introduced in revision (manuscript Supplementary Text S2), comprising a calibrated simulation study and a real-data Type I error study (`da_benchmark/`).

### Key Findings

1. **Depth–group confounding diagnosis** using the Kruskal–Wallis test with η² effect size provides an actionable framework for assessing whether rarefaction-based analyses may be compromised.

2. **Sensitivity analysis across multiple rarefaction depths** demonstrates that diversity conclusions can remain stable even under potential confounding, provided depth selection is conservative.

3. **Repeated rarefaction (q2-boots)** reveals that stochastic variability from single rarefaction is negligible (CV < 1%), supporting the use of single rarefaction for routine diversity analyses.

4. **Depth-confounding DA benchmark**: when sequencing depth is associated with the comparison groups, rarefaction-based DA testing inflates false positives relative to unrarefied analysis of the same data; the effect is method-dependent (ANCOM-BC2, MaAsLin 3, Wilcoxon). Code, frozen run configurations, and summary outputs are archived in `da_benchmark/`.

---

## Scope

This repository reproduces **both benchmarks presented in the manuscript**:

- **Diversity-side benchmark** (Steps 0–5 of the decision framework): pipeline scripts `script/01`–`07` plus the Step 0 diagnostic `script/depth_group_confounding_diagnosis.py`.
- **Depth-confounding DA benchmark** (Supplementary Text S2): simulation and real-data harness code, frozen run configurations, per-run environment records, and summary outputs under `da_benchmark/`. Full per-replicate run bodies (~8.8 GB) are not distributable in a repository; they are regenerable from the archived code, configurations, and seeds (see `da_benchmark/README.md`).

The DA **method implementations themselves** are not re-implemented here. The benchmark calls the original R packages: [ANCOM-BC2](https://github.com/FrederickHuangLin/ANCOMBC), [MaAsLin 3](https://github.com/biobakery/maaslin3), and base-R Wilcoxon; [LinDA](https://github.com/zhouhj1994/LinDA) and [ALDEx2](https://bioconductor.org/packages/ALDEx2/) are discussed in the manuscript.

---

## Repository Structure

```
rarefaction-review-benchmark/
├── README.md                      # This file
├── LICENSE                        # MIT License
├── CITATION.cff                   # Citation information
├── environment.yml                # Conda environment specification (QIIME 2)
│
├── data/
│   ├── README.md                  # Data acquisition instructions (ENA accessions)
│   └── metadata/
│       ├── metadata_soil.tsv      # Soil dataset metadata
│       ├── metadata_bioethanol.tsv # Bioethanol dataset metadata
│       └── metadata_rice.tsv      # Rice dataset metadata
│
├── script/                        # Diversity-side pipeline (frozen as executed)
│   ├── run_all.sh                 # Master script (runs all steps)
│   ├── 01_download_data.sh        # Download FASTQ from ENA (wrapper)
│   ├── 02_dada2_processing.sh     # DADA2 denoising pipeline
│   ├── 03_build_phylogeny.sh      # Phylogenetic tree construction
│   ├── 04_confounding_diagnosis.sh # Step 0: depth–group confounding diagnosis
│   ├── 05_alpha_rarefaction.sh    # Step 2: alpha-rarefaction depth selection
│   ├── 06_sensitivity_analysis.sh # Multi-depth sensitivity analysis
│   ├── 07_q2boots_analysis.sh     # Step 3: repeated rarefaction with q2-boots
│   ├── depth_group_confounding_diagnosis.py  # Step 0 diagnostic (Python)
│   ├── run_with_base_dir.sh       # Helper: run a frozen script with your data root
│   └── utils/
│       └── ena_fastq_download.sh  # ENA/NCBI/DDBJ universal FASTQ downloader
│
├── da_benchmark/                  # Depth-confounding DA benchmark (Suppl. Text S2)
│   ├── README.md                  # Layout, environments, how to run
│   ├── simulation/                # Calibrated simulation study (D-sim)
│   ├── realdata_typeI/            # Real-data Type I error study (D-real)
│   ├── common_set/                # Common-sample-set re-analysis (n=486)
│   └── reaggregate_b/             # Post-hoc re-aggregations of frozen outputs
│
├── figures/                       # Manuscript figure generation code + input tables
│   ├── README.md
│   ├── make_figure1_r1.py         # Figure 1 (decision framework)
│   ├── make_figure2_r1.py         # Figure 2 (alpha diversity vs depth)
│   ├── figure3_alpha_stability_v2_1.py  # Figure 3 (q2-boots CV)
│   ├── make_figure4_r1.py         # Figure 4 (DA rarefaction schematic)
│   └── data/                      # Input summary tables (CSV) for Figures 2–3
│
└── provenance/
    ├── README.md                  # QIIME 2 provenance documentation + archive scope
    └── rice_sensitivity_permanova/  # Archived primary PERMANOVA .qzv (4 depths)
```

---

## Requirements

### Software Dependencies (diversity-side benchmark)

| Software | Version | Purpose |
|----------|---------|---------|
| QIIME 2 amplicon distribution | 2025.10 | Core analysis platform |
| q2cli | 2025.10.1 | QIIME 2 command-line interface |
| q2-boots | 2025.10.1 | Repeated rarefaction/bootstrapping |
| q2-composition | 2025.10.1 | ANCOM-BC2 (via QIIME 2 plugin) |
| Python | 3.10.14 | Scripting and visualization |

### Software Dependencies (DA benchmark)

The DA benchmark harness runs in R. Two separate R environments were used, and their exact package versions are archived per run in `da_benchmark/*/runs/<run_id>/sessionInfo/`:

| Environment | Key packages | Recorded in |
|-------------|--------------|-------------|
| ANCOM-BC2 environment | ANCOMBC (R) | `sessionInfo/ancombc_env_R.txt` |
| MaAsLin 3 environment | maaslin3 1.2.0 (R 4.5.3) | `sessionInfo/maaslin3_env_R.txt` |

See `da_benchmark/README.md` for how the harness locates these environments (`ANCOMBC_RSCRIPT` / `MAASLIN3_RSCRIPT` environment variables).

### System Requirements

- **Operating System**: macOS (Intel or Apple Silicon via Rosetta 2) or Linux
- **RAM**: 16 GB minimum, 32 GB recommended for the Rice dataset
- **Storage**: ~50 GB for raw FASTQ files and intermediate outputs (the repository itself is ~55 MB)

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

### 3. Point the pipeline at your data root

The numbered pipeline scripts are archived **exactly as executed** for the manuscript, including the absolute data root of the original analysis machine (`BASE_DIR="/Volumes/PS3000/benchmark_data"`). To run them on your machine without modifying the frozen scripts, use the provided helper, which patches only that single assignment in a temporary copy and refuses to run if anything else would change:

```bash
BASE_DIR=/path/to/your/data_root bash script/run_with_base_dir.sh 04_confounding_diagnosis.sh
```

### 4. Run the pipeline

```bash
# Individual steps via the helper (recommended)
BASE_DIR=... bash script/run_with_base_dir.sh 01_download_data.sh   # FASTQ download (~2–4 h)
BASE_DIR=... bash script/run_with_base_dir.sh 02_dada2_processing.sh # DADA2 (~2–3 h)
BASE_DIR=... bash script/run_with_base_dir.sh 03_build_phylogeny.sh  # Trees (~30–60 min)
BASE_DIR=... bash script/run_with_base_dir.sh 04_confounding_diagnosis.sh # Step 0 (~10 min)
BASE_DIR=... bash script/run_with_base_dir.sh 05_alpha_rarefaction.sh # Step 2 (~30–60 min)
BASE_DIR=... bash script/run_with_base_dir.sh 06_sensitivity_analysis.sh # (~2–3 h)
BASE_DIR=... bash script/run_with_base_dir.sh 07_q2boots_analysis.sh # Step 3 (~2–3 h)
```

(`script/run_all.sh` is the frozen master script from the original run; it invokes the numbered scripts directly and therefore uses the original `BASE_DIR`. For a full pipeline run on another machine, loop the steps through `run_with_base_dir.sh` as shown in its header.)

The Step 0 diagnostic can also be called directly on any feature table:

```bash
python script/depth_group_confounding_diagnosis.py \
  --table-qza table.qza --metadata metadata.tsv \
  --group-column treatment --dataset-name mydata --output-dir out
```

---

## Datasets

Three public 16S rRNA amplicon datasets were selected from Schloss (2024) [1] to represent distinct regimes of sequencing depth, sample size, and confounding severity.

| Dataset | Accession | Samples imported | Primary Group | Confounding Status |
|---------|-----------|------------------|---------------|-------------------|
| Soil | PRJEB10725 | 18 | treatment (k=3) | POTENTIAL (η²=0.28) |
| Bioethanol | PRJNA276052 | 95 | batch (k=19) / time_point (k=12) | INDETERMINATE (p=0.110) / POTENTIAL (η²=0.12) |
| Rice | PRJNA255789 | 493* | compartment (k=4) | MINOR (η²=0.02) |

\*Of 498 registered samples, 5 were excluded due to ENA paired-end file registration issues (see `data/README.md`); 493 samples were imported.

See `data/README.md` for detailed data acquisition instructions.

---

## Key Parameters

### Rarefaction Depths

| Dataset    | Primary Depth | Selection Criteria                         | Sample Retention |
| ---------- | ------------- | ------------------------------------------ | ---------------- |
| Soil       | 20,000        | Shannon plateau reached                    | 100% (18/18)     |
| Bioethanol | 35,000        | Shannon plateau behaviour; minimal dropout | 100% (95/95)     |
| Rice       | 20,000        | Shannon plateau reached                    | 99% (488/493)    |

### Confounding Diagnosis Thresholds

The manuscript's revised decision framework (Table 1; Supplementary Text S1 "Rule versions") uses a four-category rule:

| Status | Criteria | Interpretation |
|--------|----------|---------------|
| NONE | p ≥ 0.20 | No evidence of depth–group association |
| INDETERMINATE | 0.05 ≤ p < 0.20 | Evidence inconclusive; report and interpret with context |
| MINOR | p < 0.05 and η² < 0.06 | Statistically significant but small effect |
| POTENTIAL | p < 0.05 and η² ≥ 0.06 | Meaningful confounding; sensitivity analysis required |

η² is calculated as: η² = (H − k + 1) / (N − k), where H is the Kruskal–Wallis statistic, k is the number of groups, and N is the total sample size.

> **Note:** `depth_group_confounding_diagnosis.py` is frozen as executed and prints the original three-category labels (NONE / MINOR / POTENTIAL with NONE at p ≥ 0.05). The four-category classification above is derivable from the p-value and η² it reports; both rule versions are stated in Supplementary Text S1.

---

## Depth-Confounding DA Benchmark (`da_benchmark/`)

The revision adds a two-part DA benchmark (manuscript Supplementary Text S2):

- **Simulation study (D-sim)** — Dirichlet-multinomial simulation calibrated to the Rice rhizoplane compartment (`simulation/harness/calibrate/`), with depth–group confounding imposed by design; null and spike-in cells; ANCOM-BC2, MaAsLin 3, and Wilcoxon fits on rarefied vs unrarefied tables.
- **Real-data Type I error study (D-real)** — artificial two-group splits of the Rice dataset (depth-ordered vs restricted permutation), where any detection is a false positive.

What is archived per run: the frozen cell configuration (`cells.tsv`, `config.yaml`), the execution record (`run.env`), exact package environments (`sessionInfo/`), and all summary outputs (`summary/`). The harness code, calibration outputs, seeds, and configurations regenerate the full runs; per-replicate bodies (~8.8 GB) are not stored in the repository. See `da_benchmark/README.md` for layout, environment setup, and a ~3-minute end-to-end smoke run.

---

## Reproducibility Notes

### Frozen scripts

All analysis code in `script/`, `da_benchmark/`, and `figures/` is archived byte-identical to the versions that produced the manuscript results (verified against the project's freeze ledger). Machine-specific absolute paths visible in `script/` and in `da_benchmark` execution records (`run.env`) are part of that frozen record; use `script/run_with_base_dir.sh` (diversity pipeline) or the `REVISION_ROOT` environment variable (DA harness) to run on another machine.

### QIIME 2 Provenance

All QIIME 2 artifacts (.qza) and visualizations (.qzv) contain embedded provenance tracking. Upload any .qzv file to [https://view.qiime2.org](https://view.qiime2.org) to inspect software versions, input file checksums, complete parameter specifications, and execution timestamps. The `provenance/` directory documents this system and archives the primary PERMANOVA .qzv files (4 rarefaction depths, weighted UniFrac, Rice) that underpin the manuscript's beta-diversity stability values.

### Random Seeds

Where applicable, random seeds are fixed for reproducibility:

- DADA2: deterministic given identical input
- Rarefaction: controlled via `--p-seed` parameter
- q2-boots: controlled via `--p-random-state` parameter
- DA benchmark: hierarchical seeding from a single `BASE_SEED` per run (recorded in `run.env`), making every replicate deterministically reproducible

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
5. Lin H, Peddada SD. Multigroup analysis of compositions of microbiomes with covariate adjustments and repeated measures. *Nat Methods*. 2024;21:83–91.
6. Nickols WA, et al. MaAsLin 3: refining and extending generalized multivariable linear models for meta-omic association discovery. (maaslin3 R package.)

---

## License

This project is licensed under the MIT License — see the [LICENSE](LICENSE) file for details.

---

## Contact

For questions or issues, please open a GitHub issue or contact:

**Shunji Tokoro**  
Gifu University of Medical Science  
ORCID: [0009-0003-8042-9869](https://orcid.org/0009-0003-8042-9869)
