# Data Acquisition

This document describes how to obtain the raw sequencing data used in the benchmark analyses.

---

## Overview

Three public 16S rRNA amplicon datasets were selected from Schloss (2024) to represent distinct regimes of sequencing depth, sample size, and depth–group confounding severity.

| Dataset | Study Accession | Samples | Region | Platform | Source |
|---------|-----------------|---------|--------|----------|--------|
| Soil | PRJEB10725 | 18 | V4 | Illumina MiSeq | ENA |
| Bioethanol | PRJNA276052 | 95 | V4 | Illumina MiSeq | ENA |
| Rice | PRJNA255789 | 493* | V4 | Illumina MiSeq | ENA |

*Five samples excluded due to incomplete paired-end data (see Known Issues below).

---

## Data Source

All data are downloaded from the **European Nucleotide Archive (ENA)**, which mirrors NCBI SRA data in a more accessible format (direct FASTQ downloads without requiring SRA Toolkit).

ENA URLs follow this pattern:
```
https://www.ebi.ac.uk/ena/browser/view/{ACCESSION}
```

---

## Automated Download

The recommended method is to use the provided download script:

```bash
# From repository root directory (see script/run_with_base_dir.sh to
# redirect the download target to your own data root)
BASE_DIR=/path/to/your/data_root bash script/run_with_base_dir.sh 01_download_data.sh
```

This script will:
1. Create the necessary directory structure
2. Download run information tables from ENA
3. Download paired-end FASTQ files for all samples
4. Verify file integrity using MD5 checksums (when available)

**Estimated download time**: 2-4 hours (depending on connection speed)  
**Storage required**: ~30 GB

---

## Manual Download Instructions

If the automated script fails or you prefer manual download, follow these steps for each dataset.

### Soil Dataset (PRJEB10725)

1. Navigate to: https://www.ebi.ac.uk/ena/browser/view/PRJEB10725

2. Click "Download All" or use wget:
   ```bash
   # Create directory
   mkdir -p data/raw/soil
   
   # Download run info
   wget -O data/raw/soil/filereport_PRJEB10725.tsv \
     "https://www.ebi.ac.uk/ena/portal/api/filereport?accession=PRJEB10725&result=read_run&fields=run_accession,fastq_ftp,fastq_md5&format=tsv"
   
   # Download FASTQ files (example for one sample)
   wget -P data/raw/soil/ ftp://ftp.sra.ebi.ac.uk/vol1/fastq/ERR102/ERR1021???/ERR1021???_1.fastq.gz
   wget -P data/raw/soil/ ftp://ftp.sra.ebi.ac.uk/vol1/fastq/ERR102/ERR1021???/ERR1021???_2.fastq.gz
   ```

### Bioethanol Dataset (PRJNA276052)

1. Navigate to: https://www.ebi.ac.uk/ena/browser/view/PRJNA276052

2. Download procedure:
   ```bash
   mkdir -p data/raw/bioethanol
   
   wget -O data/raw/bioethanol/filereport_PRJNA276052.tsv \
     "https://www.ebi.ac.uk/ena/portal/api/filereport?accession=PRJNA276052&result=read_run&fields=run_accession,fastq_ftp,fastq_md5&format=tsv"
   ```

### Rice Dataset (PRJNA255789)

1. Navigate to: https://www.ebi.ac.uk/ena/browser/view/PRJNA255789

2. Download procedure:
   ```bash
   mkdir -p data/raw/rice
   
   wget -O data/raw/rice/filereport_PRJNA255789.tsv \
     "https://www.ebi.ac.uk/ena/portal/api/filereport?accession=PRJNA255789&result=read_run&fields=run_accession,fastq_ftp,fastq_md5&format=tsv"
   ```

---

## Known Issues

### Rice Dataset: Missing Paired-End Files

Five samples in the Rice dataset are registered as paired-end sequencing but only have single merged FASTQ files available in ENA:

| Run Accession | Compartment | Issue |
|---------------|-------------|-------|
| SRR1524332 | Pre_planting | Single file only |
| SRR1524667 | Pre_planting | Single file only |
| SRR1524669 | Pre_planting | Single file only |
| SRR1524671 | Pre_planting | Single file only |
| SRR1524672 | Bulk Soil | Single file only |

**Impact**: The entire Pre_planting compartment (4 samples) is excluded from the analysis. This limitation is documented in the manuscript Methods section.

**Identification**: These samples can be identified by their `AvgSpotLen` values (245-247 bp), which is approximately half the expected length for paired 250 bp reads (~494 bp).

---

## Metadata Files

Pre-formatted metadata files compatible with QIIME 2 are provided in `data/metadata/`:

| File | Dataset | Key Variables |
|------|---------|---------------|
| `metadata_soil.tsv` | Soil | treatment (C/T/H) |
| `metadata_bioethanol.tsv` | Bioethanol | batch, time_point |
| `metadata_rice.tsv` | Rice | compartment |

These metadata files were generated from ENA sample information and formatted to meet QIIME 2 requirements (tab-separated, `sample-id` as first column, `#q2:types` directive row).

---

## File Organization After Download

After successful download, your data directory should have this structure:

```
data/
├── README.md                      # This file
├── raw/
│   ├── soil/
│   │   ├── filereport_PRJEB10725.tsv
│   │   ├── ERR1021XXX_1.fastq.gz
│   │   ├── ERR1021XXX_2.fastq.gz
│   │   └── ...
│   ├── bioethanol/
│   │   ├── filereport_PRJNA276052.tsv
│   │   ├── SRR18XXXXX_1.fastq.gz
│   │   ├── SRR18XXXXX_2.fastq.gz
│   │   └── ...
│   └── rice/
│       ├── filereport_PRJNA255789.tsv
│       ├── SRR15XXXXX_1.fastq.gz
│       ├── SRR15XXXXX_2.fastq.gz
│       └── ...
└── metadata/
    ├── metadata_soil.tsv
    ├── metadata_bioethanol.tsv
    └── metadata_rice.tsv
```

---

## Verification

After download, verify file counts:

```bash
# Expected file counts (paired-end, so 2 files per sample)
echo "Soil: $(ls data/raw/soil/*.fastq.gz 2>/dev/null | wc -l) files (expected: 36)"
echo "Bioethanol: $(ls data/raw/bioethanol/*.fastq.gz 2>/dev/null | wc -l) files (expected: 190)"
echo "Rice: $(ls data/raw/rice/*.fastq.gz 2>/dev/null | wc -l) files (expected: ~986)"
```

---

## References

Schloss PD. Rarefaction is currently the best approach to control for uneven sequencing effort in amplicon sequence analyses. *mSphere*. 2024;9:e00354-23. doi:10.1128/msphere.00354-23

Original dataset publications are cited in the manuscript.
