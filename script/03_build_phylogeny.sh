#!/bin/bash
# ============================================================================
# Script: 03_build_phylogeny.sh
# Purpose: Build phylogenetic trees for all benchmark datasets
# 
# This script constructs rooted phylogenetic trees from representative
# sequences using MAFFT alignment and FastTree. The trees are required
# for phylogenetic diversity metrics (UniFrac, Faith's PD).
#
# Prerequisites:
#   - QIIME 2 environment activated (qiime2-amplicon-2025.10)
#   - DADA2 processing completed (02_dada2_processing.sh)
#
# Usage:
#   conda activate qiime2-amplicon-2025.10
#   bash scripts/03_build_phylogeny.sh
#
# Estimated time: 30-60 minutes (Rice dataset takes longest due to ~80k ASVs)
#
# Pipeline:
#   1. MAFFT alignment of representative sequences
#   2. Masking of highly variable positions
#   3. FastTree construction of unrooted tree
#   4. Midpoint rooting
#
# Output per dataset:
#   02_denoised/aligned-rep-seqs.qza        - MAFFT alignment
#   02_denoised/masked-aligned-rep-seqs.qza - Masked alignment
#   02_denoised/unrooted-tree.qza           - Unrooted tree
#   02_denoised/rooted-tree.qza             - Rooted tree (midpoint)
#
# Reference: Manuscript Methods - Phylogenetic Analysis
# ============================================================================

set -euo pipefail

# ============================================================================
# Configuration
# ============================================================================

# Base directory - MODIFY THIS to match your environment
BASE_DIR="/Volumes/PS3000/benchmark_data"

# Number of threads
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

log_success "QIIME 2 environment detected"

# ============================================================================
# 1. Soil Dataset
# ============================================================================

log_header "1. Soil Dataset Phylogeny"

SOIL_DIR="${BASE_DIR}/Soil_v2/02_denoised"

log_info "Input: ${SOIL_DIR}/rep-seqs.qza"
log_info "ASV count: ~6,485"

qiime phylogeny align-to-tree-mafft-fasttree \
    --i-sequences "${SOIL_DIR}/rep-seqs.qza" \
    --o-alignment "${SOIL_DIR}/aligned-rep-seqs.qza" \
    --o-masked-alignment "${SOIL_DIR}/masked-aligned-rep-seqs.qza" \
    --o-tree "${SOIL_DIR}/unrooted-tree.qza" \
    --o-rooted-tree "${SOIL_DIR}/rooted-tree.qza" \
    --p-n-threads "${N_THREADS}"

log_success "Soil phylogeny complete"

# ============================================================================
# 2. Bioethanol Dataset
# ============================================================================

log_header "2. Bioethanol Dataset Phylogeny"

BIO_DIR="${BASE_DIR}/Bioethanol_v2/02_denoised"

log_info "Input: ${BIO_DIR}/rep-seqs.qza"
log_info "ASV count: ~2,194"

qiime phylogeny align-to-tree-mafft-fasttree \
    --i-sequences "${BIO_DIR}/rep-seqs.qza" \
    --o-alignment "${BIO_DIR}/aligned-rep-seqs.qza" \
    --o-masked-alignment "${BIO_DIR}/masked-aligned-rep-seqs.qza" \
    --o-tree "${BIO_DIR}/unrooted-tree.qza" \
    --o-rooted-tree "${BIO_DIR}/rooted-tree.qza" \
    --p-n-threads "${N_THREADS}"

log_success "Bioethanol phylogeny complete"

# ============================================================================
# 3. Rice Dataset
# ============================================================================

log_header "3. Rice Dataset Phylogeny"

RICE_DIR="${BASE_DIR}/Rice_v2/02_denoised"

log_info "Input: ${RICE_DIR}/rep-seqs.qza"
log_info "ASV count: ~79,807"
log_warn "This dataset is large and will take longer to process"

qiime phylogeny align-to-tree-mafft-fasttree \
    --i-sequences "${RICE_DIR}/rep-seqs.qza" \
    --o-alignment "${RICE_DIR}/aligned-rep-seqs.qza" \
    --o-masked-alignment "${RICE_DIR}/masked-aligned-rep-seqs.qza" \
    --o-tree "${RICE_DIR}/unrooted-tree.qza" \
    --o-rooted-tree "${RICE_DIR}/rooted-tree.qza" \
    --p-n-threads "${N_THREADS}"

log_success "Rice phylogeny complete"

# ============================================================================
# Summary
# ============================================================================

log_header "Phylogeny Construction Complete"

echo ""
log_info "Output files created:"
echo "  Soil:      ${SOIL_DIR}/rooted-tree.qza"
echo "  Bioethanol: ${BIO_DIR}/rooted-tree.qza"
echo "  Rice:      ${RICE_DIR}/rooted-tree.qza"
echo ""
log_info "Next step: Run confounding diagnosis"
echo "  bash scripts/04_confounding_diagnosis.sh"
