#!/bin/bash
# ============================================================================
# Script: 01_download_data.sh
# Purpose: Download FASTQ files from ENA for all benchmark datasets
# 
# This wrapper script uses the ENA Universal FASTQ Downloader (utils/ena_fastq_download.sh)
# to download raw sequencing data for three 16S rRNA amplicon datasets:
#   - Soil (PRJEB10725): 18 samples
#   - Bioethanol (PRJNA276052): 95 samples
#   - Rice (PRJNA255789): 493 samples (498 - 5 excluded)
#
# The utility script provides:
#   - Parallel downloads with configurable connections
#   - MD5 checksum verification
#   - Resume capability for interrupted downloads
#   - Filtering by library strategy or platform
#
# Prerequisites:
#   - wget and curl installed
#   - Approximately 30 GB free disk space
#   - Stable internet connection
#
# Usage:
#   bash scripts/01_download_data.sh
#
# Options (passed to underlying downloader):
#   -p N    : Number of parallel downloads (default: 4)
#   --no-md5: Skip MD5 verification for faster downloads
#
# Estimated time: 2-4 hours (depending on connection speed)
#
# Output structure:
#   {OUTPUT_DIR}/{dataset}/01_fastq/  - FASTQ files
#   {OUTPUT_DIR}/{dataset}/00_meta/   - Metadata from ENA
#
# Reference: Manuscript Methods - Data Acquisition
# ============================================================================

set -euo pipefail

# ============================================================================
# Configuration
# ============================================================================

# Script locations
SCRIPT_DIR="$(cd "$(dirname "$0")" && pwd)"
DOWNLOADER="${SCRIPT_DIR}/utils/ena_fastq_download.sh"

# Output directory - MODIFY THIS to match your environment
OUTPUT_BASE="/Volumes/PS3000/benchmark_data"

# Dataset accessions (from Schloss 2024 mSphere study)
SOIL_ACCESSION="PRJEB10725"
BIOETHANOL_ACCESSION="PRJNA276052"
RICE_ACCESSION="PRJNA255789"

# Download settings
PARALLEL_DOWNLOADS=4    # Number of concurrent downloads
VERIFY_MD5=true         # Set to false for faster downloads (--no-md5)

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

# Check downloader script exists
if [[ ! -f "${DOWNLOADER}" ]]; then
    log_error "ENA downloader not found: ${DOWNLOADER}"
    log_error "Please ensure utils/ena_fastq_download.sh is in place."
    exit 1
fi
log_success "ENA downloader found"

# Make downloader executable
chmod +x "${DOWNLOADER}"

# Check required commands
for cmd in wget curl awk; do
    if ! command -v "$cmd" &> /dev/null; then
        log_error "Required command not found: $cmd"
        exit 1
    fi
done
log_success "All dependencies available"

# Build common options
COMMON_OPTS="-p ${PARALLEL_DOWNLOADS}"
if [[ "${VERIFY_MD5}" != true ]]; then
    COMMON_OPTS="${COMMON_OPTS} --no-md5"
fi

log_info "Output base: ${OUTPUT_BASE}"
log_info "Parallel downloads: ${PARALLEL_DOWNLOADS}"
log_info "MD5 verification: ${VERIFY_MD5}"

# ============================================================================
# Download Datasets
# ============================================================================

log_header "Downloading Benchmark Datasets"

# ----------------------------------------------------------------------------
# 1. Soil Dataset (PRJEB10725)
# ----------------------------------------------------------------------------
log_header "1/3: Soil Dataset (${SOIL_ACCESSION})"

log_info "Description: Agricultural soil microbiome study"
log_info "Samples: 18 (treatment groups: C, T, H)"
log_info "Region: V4 (515F-806R)"

"${DOWNLOADER}" "${SOIL_ACCESSION}" \
    -o "${OUTPUT_BASE}/Soil_v2" \
    ${COMMON_OPTS}

log_success "Soil dataset download complete"

# ----------------------------------------------------------------------------
# 2. Bioethanol Dataset (PRJNA276052)
# ----------------------------------------------------------------------------
log_header "2/3: Bioethanol Dataset (${BIOETHANOL_ACCESSION})"

log_info "Description: Bioethanol fermentation microbiome time series"
log_info "Samples: 95 (19 batches × ~5 time points)"
log_info "Region: V4"

"${DOWNLOADER}" "${BIOETHANOL_ACCESSION}" \
    -o "${OUTPUT_BASE}/Bioethanol_v2" \
    ${COMMON_OPTS}

log_success "Bioethanol dataset download complete"

# ----------------------------------------------------------------------------
# 3. Rice Dataset (PRJNA255789)
# ----------------------------------------------------------------------------
log_header "3/3: Rice Dataset (${RICE_ACCESSION})"

log_info "Description: Rice rhizosphere microbiome compartment study"
log_info "Samples: 493 (4 compartments: Endosphere, Rhizosphere, Rhizoplane, Bulk Soil)"
log_info "Region: V4"
log_warn "Note: 5 samples excluded due to incomplete paired-end data in ENA"
log_warn "This is the largest dataset and will take the longest to download"

"${DOWNLOADER}" "${RICE_ACCESSION}" \
    -o "${OUTPUT_BASE}/Rice_v2" \
    ${COMMON_OPTS}

log_success "Rice dataset download complete"

# ============================================================================
# Summary
# ============================================================================

log_header "All Downloads Complete"

echo ""
log_info "Downloaded data locations:"
echo "  Soil:       ${OUTPUT_BASE}/Soil_v2/01_fastq/"
echo "  Bioethanol: ${OUTPUT_BASE}/Bioethanol_v2/01_fastq/"
echo "  Rice:       ${OUTPUT_BASE}/Rice_v2/01_fastq/"
echo ""
log_info "ENA metadata saved to:"
echo "  ${OUTPUT_BASE}/{dataset}/00_meta/"
echo ""
log_info "Next step: Import FASTQ files into QIIME 2 and run DADA2"
echo "  bash scripts/02_dada2_processing.sh"
