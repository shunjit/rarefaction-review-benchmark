#!/bin/bash
# ============================================================================
# Script: 05_alpha_rarefaction.sh
# Purpose: Generate alpha-rarefaction curves for all benchmark datasets
# 
# This script generates alpha-rarefaction curves to:
#   1. Determine if sequencing depth is sufficient (plateau detection)
#   2. Identify appropriate rarefaction depth for diversity analysis
#   3. Visualize diversity accumulation patterns across groups
#
# Prerequisites:
#   - QIIME 2 environment activated (qiime2-amplicon-2025.10)
#   - DADA2 processing and phylogeny construction completed
#   - Confounding diagnosis completed (to inform depth selection)
#
# Usage:
#   conda activate qiime2-amplicon-2025.10
#   bash scripts/05_alpha_rarefaction.sh
#
# Estimated time: 30-60 minutes
#
# Output per dataset:
#   03_diversity/alpha_rarefaction.qzv - Interactive rarefaction curves
#
# Reference: Manuscript Section 5.1 - Step 1: Evaluate alpha-rarefaction plateau
# ============================================================================

set -euo pipefail

# ============================================================================
# Configuration
# ============================================================================

# Base directory - MODIFY THIS to match your environment
BASE_DIR="/Volumes/PS3000/benchmark_data"
METADATA_DIR="${BASE_DIR}/00_github/Schloss_Rarefaction_mSphere_2024/data"

# Color output
RED='\033[0;31m'
GREEN='\033[0;32m'
YELLOW='\033[0;33m'
BLUE='\033[0;34m'
BOLD='\033[1m'
NC='\033[0m'

log_info() { echo -e "${BLUE}[INFO]${NC} $*"; }
log_success() { echo -e "${GREEN}[OK]${NC} $*"; }
log_warn() { echo -e "${YELLOW}[WARN]${NC} $*"; }
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

if ! command -v qiime &> /dev/null; then
    echo "Error: QIIME 2 not found. Please activate the conda environment."
    exit 1
fi
log_success "QIIME 2 environment detected"

# ============================================================================
# 1. Soil Dataset
# ============================================================================

log_header "1. Soil Dataset"

log_info "Depth range: 36,251-56,114 reads"
log_info "max-depth: 35,000 (conservative, below H group minimum of 36,834)"
log_info "Confounding status: POTENTIAL (η²=0.28)"

SOIL_OUT="${BASE_DIR}/Soil_v2/03_diversity"
mkdir -p "${SOIL_OUT}"

qiime diversity alpha-rarefaction \
    --i-table "${BASE_DIR}/Soil_v2/02_denoised/table.qza" \
    --i-phylogeny "${BASE_DIR}/Soil_v2/02_denoised/rooted-tree.qza" \
    --p-max-depth 35000 \
    --p-min-depth 5000 \
    --p-steps 15 \
    --p-iterations 10 \
    --m-metadata-file "${METADATA_DIR}/soil/metadata_soil.tsv" \
    --o-visualization "${SOIL_OUT}/alpha_rarefaction.qzv"

log_success "Output: ${SOIL_OUT}/alpha_rarefaction.qzv"

# ============================================================================
# 2. Bioethanol Dataset
# ============================================================================

log_header "2. Bioethanol Dataset"

log_info "Depth range: 36,110-383,223 reads"
log_info "max-depth: 100,000 (exploring up to median range)"
log_info "Confounding status: batch=NONE, time_point=POTENTIAL (η²=0.12)"

BIO_OUT="${BASE_DIR}/Bioethanol_v2/03_diversity"
mkdir -p "${BIO_OUT}"

qiime diversity alpha-rarefaction \
    --i-table "${BASE_DIR}/Bioethanol_v2/02_denoised/table.qza" \
    --i-phylogeny "${BASE_DIR}/Bioethanol_v2/02_denoised/rooted-tree.qza" \
    --p-max-depth 100000 \
    --p-min-depth 10000 \
    --p-steps 15 \
    --p-iterations 10 \
    --m-metadata-file "${METADATA_DIR}/bioethanol/metadata_bioethanol.tsv" \
    --o-visualization "${BIO_OUT}/alpha_rarefaction.qzv"

log_success "Output: ${BIO_OUT}/alpha_rarefaction.qzv"

# ============================================================================
# 3. Rice Dataset
# ============================================================================

log_header "3. Rice Dataset"

log_info "Depth range: 2,112-201,037 reads"
log_info "max-depth: 60,000 (exploring up to median range)"
log_info "Confounding status: MINOR (η²=0.02)"
log_warn "Note: Large dataset (n=493), this step may take longer"

RICE_OUT="${BASE_DIR}/Rice_v2/03_diversity"
mkdir -p "${RICE_OUT}"

qiime diversity alpha-rarefaction \
    --i-table "${BASE_DIR}/Rice_v2/02_denoised/table.qza" \
    --i-phylogeny "${BASE_DIR}/Rice_v2/02_denoised/rooted-tree.qza" \
    --p-max-depth 60000 \
    --p-min-depth 2000 \
    --p-steps 15 \
    --p-iterations 10 \
    --m-metadata-file "${METADATA_DIR}/rice/metadata_rice.tsv" \
    --o-visualization "${RICE_OUT}/alpha_rarefaction.qzv"

log_success "Output: ${RICE_OUT}/alpha_rarefaction.qzv"

# ============================================================================
# Summary
# ============================================================================

log_header "Alpha-Rarefaction Complete"

echo ""
log_info "Generated files:"
echo "  1. ${BASE_DIR}/Soil_v2/03_diversity/alpha_rarefaction.qzv"
echo "  2. ${BASE_DIR}/Bioethanol_v2/03_diversity/alpha_rarefaction.qzv"
echo "  3. ${BASE_DIR}/Rice_v2/03_diversity/alpha_rarefaction.qzv"
echo ""
log_info "Next steps:"
echo "  1. Open each .qzv file at https://view.qiime2.org/"
echo "  2. Identify plateau depth for each dataset"
echo "  3. Proceed to sensitivity analysis"
echo ""
echo "  bash scripts/06_sensitivity_analysis.sh"
