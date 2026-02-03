#!/bin/bash
# ============================================================================
# Script: 02_dada2_processing.sh
# Purpose: Run DADA2 denoising pipeline for all benchmark datasets
# 
# This script performs DADA2 paired-end denoising on the downloaded FASTQ
# files, generating ASV feature tables and representative sequences.
#
# Prerequisites:
#   - QIIME 2 environment activated (qiime2-amplicon-2025.10)
#   - FASTQ files downloaded (01_download_data.sh completed)
#   - Sufficient disk space (~10 GB for intermediate files)
#
# Usage:
#   conda activate qiime2-amplicon-2025.10
#   bash scripts/02_dada2_processing.sh
#
# Estimated time: 2-3 hours (Rice dataset takes longest)
#
# DADA2 Parameters (determined from quality profiles):
#   Soil:      --trunc-len-f 133 --trunc-len-r 133
#   Bioethanol: --trunc-len-f 200 --trunc-len-r 200
#   Rice:      --trunc-len-f 133 --trunc-len-r 133
#
# Output per dataset:
#   02_denoised/table.qza        - Feature table (ASV counts)
#   02_denoised/rep-seqs.qza     - Representative sequences
#   02_denoised/stats.qza        - Denoising statistics
#   02_denoised/table.qzv        - Feature table visualization
#   02_denoised/stats.qzv        - Statistics visualization
#
# Reference: Manuscript Methods - Sequence Processing
# ============================================================================

set -euo pipefail

# ============================================================================
# Configuration
# ============================================================================

# Base directory - MODIFY THIS to match your environment
BASE_DIR="/Volumes/PS3000/benchmark_data"

# Metadata directory (from Schloss 2024 repository)
METADATA_DIR="${BASE_DIR}/00_github/Schloss_Rarefaction_mSphere_2024/data"

# Number of threads for DADA2
N_THREADS=8

# Color output
RED='\033[0;31m'
GREEN='\033[0;32m'
YELLOW='\033[0;33m'
BLUE='\033[0;34m'
BOLD='\033[1m'
NC='\033[0m'

log_info() { echo -e "${BLUE}[INFO]${NC} $(date '+%Y-%m-%d %H:%M:%S') $*"; }
log_success() { echo -e "${GREEN}[OK]${NC} $(date '+%Y-%m-%d %H:%M:%S') $*"; }
log_warn() { echo -e "${YELLOW}[WARN]${NC} $(date '+%Y-%m-%d %H:%M:%S') $*"; }
log_error() { echo -e "${RED}[ERROR]${NC} $(date '+%Y-%m-%d %H:%M:%S') $*" >&2; }
log_header() {
    echo ""
    echo -e "${BOLD}========================================${NC}"
    echo -e "${BOLD}  $*${NC}"
    echo -e "${BOLD}========================================${NC}"
}

# ============================================================================
# Pre-flight Checks
# ============================================================================

log_header "Pre-flight Checks"

# Check QIIME 2 environment
if ! command -v qiime &> /dev/null; then
    log_error "QIIME 2 not found. Please activate the conda environment:"
    log_error "  conda activate qiime2-amplicon-2025.10"
    exit 1
fi

QIIME_VERSION=$(qiime --version | head -1)
log_success "QIIME 2 detected: ${QIIME_VERSION}"

# ============================================================================
# 1. Soil Dataset
# ============================================================================

log_header "1. Soil Dataset (PRJEB10725)"

SOIL_DIR="${BASE_DIR}/Soil_v2"
SOIL_OUT="${SOIL_DIR}/02_denoised"

log_info "Input: ${SOIL_DIR}/01_imported/demux.qza"
log_info "Parameters: --trunc-len-f 133 --trunc-len-r 133"
log_info "Rationale: Quality drop observed at position 180+; conservative trimming"

mkdir -p "${SOIL_OUT}"

# DADA2 denoising
qiime dada2 denoise-paired \
    --i-demultiplexed-seqs "${SOIL_DIR}/01_imported/demux.qza" \
    --p-trunc-len-f 133 \
    --p-trunc-len-r 133 \
    --p-n-threads "${N_THREADS}" \
    --o-table "${SOIL_OUT}/table.qza" \
    --o-representative-sequences "${SOIL_OUT}/rep-seqs.qza" \
    --o-denoising-stats "${SOIL_OUT}/stats.qza" \
    --verbose

# Generate visualizations
qiime metadata tabulate \
    --m-input-file "${SOIL_OUT}/stats.qza" \
    --o-visualization "${SOIL_OUT}/stats.qzv"

qiime feature-table summarize \
    --i-table "${SOIL_OUT}/table.qza" \
    --m-sample-metadata-file "${METADATA_DIR}/soil/metadata_soil.tsv" \
    --o-visualization "${SOIL_OUT}/table_with_metadata.qzv"

qiime feature-table tabulate-seqs \
    --i-data "${SOIL_OUT}/rep-seqs.qza" \
    --o-visualization "${SOIL_OUT}/rep-seqs.qzv"

log_success "Soil DADA2 processing complete"

# ============================================================================
# 2. Bioethanol Dataset
# ============================================================================

log_header "2. Bioethanol Dataset (PRJNA276052)"

BIO_DIR="${BASE_DIR}/Bioethanol_v2"
BIO_OUT="${BIO_DIR}/02_denoised"

log_info "Input: ${BIO_DIR}/01_imported/demux.qza"
log_info "Parameters: --trunc-len-f 200 --trunc-len-r 200"
log_info "Rationale: V4 region; adequate overlap maintained"

mkdir -p "${BIO_OUT}"

# DADA2 denoising
qiime dada2 denoise-paired \
    --i-demultiplexed-seqs "${BIO_DIR}/01_imported/demux.qza" \
    --p-trunc-len-f 200 \
    --p-trunc-len-r 200 \
    --p-n-threads "${N_THREADS}" \
    --o-table "${BIO_OUT}/table.qza" \
    --o-representative-sequences "${BIO_OUT}/rep-seqs.qza" \
    --o-denoising-stats "${BIO_OUT}/stats.qza" \
    --verbose

# Generate visualizations
qiime metadata tabulate \
    --m-input-file "${BIO_OUT}/stats.qza" \
    --o-visualization "${BIO_OUT}/stats.qzv"

qiime feature-table summarize \
    --i-table "${BIO_OUT}/table.qza" \
    --m-sample-metadata-file "${METADATA_DIR}/bioethanol/metadata_bioethanol.tsv" \
    --o-visualization "${BIO_OUT}/table_with_metadata.qzv"

qiime feature-table tabulate-seqs \
    --i-data "${BIO_OUT}/rep-seqs.qza" \
    --o-visualization "${BIO_OUT}/rep-seqs.qzv"

log_success "Bioethanol DADA2 processing complete"

# ============================================================================
# 3. Rice Dataset
# ============================================================================

log_header "3. Rice Dataset (PRJNA255789)"

RICE_DIR="${BASE_DIR}/Rice_v2"
RICE_OUT="${RICE_DIR}/02_denoised"

log_info "Input: ${RICE_DIR}/01_imported/demux.qza"
log_info "Parameters: --trunc-len-f 133 --trunc-len-r 133"
log_info "Rationale: Quality degradation after position 180"
log_warn "Note: This dataset is large (n=493) and will take longer to process"

mkdir -p "${RICE_OUT}"

# DADA2 denoising
qiime dada2 denoise-paired \
    --i-demultiplexed-seqs "${RICE_DIR}/01_imported/demux.qza" \
    --p-trunc-len-f 133 \
    --p-trunc-len-r 133 \
    --p-n-threads "${N_THREADS}" \
    --o-table "${RICE_OUT}/table.qza" \
    --o-representative-sequences "${RICE_OUT}/rep-seqs.qza" \
    --o-denoising-stats "${RICE_OUT}/stats.qza" \
    --verbose

# Generate visualizations
qiime metadata tabulate \
    --m-input-file "${RICE_OUT}/stats.qza" \
    --o-visualization "${RICE_OUT}/stats.qzv"

qiime feature-table summarize \
    --i-table "${RICE_OUT}/table.qza" \
    --m-sample-metadata-file "${METADATA_DIR}/rice/metadata_rice.tsv" \
    --o-visualization "${RICE_OUT}/table_with_metadata.qzv"

qiime feature-table tabulate-seqs \
    --i-data "${RICE_OUT}/rep-seqs.qza" \
    --o-visualization "${RICE_OUT}/rep-seqs.qzv"

log_success "Rice DADA2 processing complete"

# ============================================================================
# Summary
# ============================================================================

log_header "DADA2 Processing Complete"

echo ""
log_info "Output locations:"
echo "  Soil:      ${SOIL_OUT}/"
echo "  Bioethanol: ${BIO_OUT}/"
echo "  Rice:      ${RICE_OUT}/"
echo ""
log_info "Key output files per dataset:"
echo "  - table.qza              (Feature table)"
echo "  - rep-seqs.qza           (Representative sequences)"
echo "  - stats.qza/.qzv         (Denoising statistics)"
echo "  - table_with_metadata.qzv (Feature table with group info)"
echo ""
log_info "Next step: Build phylogenetic trees"
echo "  bash scripts/03_build_phylogeny.sh"
