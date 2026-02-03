#!/bin/bash
# ============================================================================
# Script: 07_q2boots_analysis.sh
# Purpose: Run repeated rarefaction analysis using q2-boots for all datasets
# 
# This script quantifies the stochastic uncertainty introduced by single
# rarefaction using q2-boots with n=100 iterations. Results demonstrate that
# stochastic variability from rarefaction is negligible (CV < 1%).
#
# Key outputs:
#   - Alpha diversity collections (100 Shannon vectors per dataset)
#   - Beta diversity collections (100 distance matrices per dataset)
#   - Average distance matrices for downstream visualization
#
# Prerequisites:
#   - QIIME 2 environment activated (qiime2-amplicon-2025.10)
#   - q2-boots plugin installed (verify: qiime boots --help)
#   - Sensitivity analysis completed (to determine primary depths)
#
# Usage:
#   conda activate qiime2-amplicon-2025.10
#   bash scripts/07_q2boots_analysis.sh
#
# Background execution (recommended):
#   nohup caffeinate -ims ./scripts/07_q2boots_analysis.sh > q2boots_log.txt 2>&1 &
#
# Estimated time: 2-3 hours
#
# Output structure per dataset:
#   05_q2boots/
#     - boots-alpha-shannon-collection.qza  (n=100 alpha vectors)
#     - boots-beta-wunifrac-collection.qza  (n=100 distance matrices)
#     - boots-beta-wunifrac-avg.qza         (averaged distance matrix)
#
# Reference: Manuscript Section 5.2 - Repeated rarefaction (q2-boots)
# ============================================================================

set -euo pipefail

# ============================================================================
# Configuration
# ============================================================================

# Base directory - MODIFY THIS to match your environment
BASE_DIR="/Volumes/PS3000/benchmark_data"
METADATA_DIR="${BASE_DIR}/00_github/Schloss_Rarefaction_mSphere_2024/data"

# Number of iterations (standard for uncertainty quantification)
N_ITERATIONS=100

# Beta averaging method
# Options: non-metric-mean, non-metric-median, medoid
# non-metric-mean: arithmetic mean of each distance matrix element (most common)
AVERAGE_METHOD="non-metric-mean"

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
log_error() { echo -e "${RED}[ERROR]${NC} $(date '+%Y-%m-%d %H:%M:%S') $*"; }
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
    log_error "QIIME 2 not found. Please run: conda activate qiime2-amplicon-2025.10"
    exit 1
fi

# Check q2-boots installation
if ! qiime boots --help &> /dev/null; then
    log_error "q2-boots not found. Please install the plugin."
    exit 1
fi

log_success "QIIME 2 and q2-boots detected"

# Verify input files exist
SOIL_TABLE="${BASE_DIR}/Soil_v2/02_denoised/table.qza"
SOIL_TREE="${BASE_DIR}/Soil_v2/02_denoised/rooted-tree.qza"
BIO_TABLE="${BASE_DIR}/Bioethanol_v2/02_denoised/table.qza"
BIO_TREE="${BASE_DIR}/Bioethanol_v2/02_denoised/rooted-tree.qza"
RICE_TABLE="${BASE_DIR}/Rice_v2/02_denoised/table.qza"
RICE_TREE="${BASE_DIR}/Rice_v2/02_denoised/rooted-tree.qza"

for f in "${SOIL_TABLE}" "${SOIL_TREE}" "${BIO_TABLE}" "${BIO_TREE}" "${RICE_TABLE}" "${RICE_TREE}"; do
    if [[ ! -f "$f" ]]; then
        log_error "File not found: $f"
        exit 1
    fi
done

log_success "All input files verified"

echo ""
log_info "=== Script started: $(date) ==="
log_info "Iterations: ${N_ITERATIONS}"
log_info "Averaging method: ${AVERAGE_METHOD}"
echo ""

# ============================================================================
# 1. Soil Dataset
# ============================================================================

log_header "1. Soil Dataset: q2-boots (n=${N_ITERATIONS})"

SOIL_DIR="${BASE_DIR}/Soil_v2"
SOIL_OUT="${SOIL_DIR}/05_q2boots"
SOIL_DEPTH=20000  # Primary depth from sensitivity analysis

log_info "Primary depth: ${SOIL_DEPTH} reads"
log_info "Iterations: ${N_ITERATIONS}"
log_info "Resampling: rarefaction (without replacement)"

mkdir -p "${SOIL_OUT}"

# Alpha diversity collection (Shannon)
# --p-no-replacement: specifies rarefaction (without replacement sampling)
# Shannon does not require phylogeny
log_info "Running alpha diversity (Shannon) collection..."
qiime boots alpha-collection \
    --i-table "${SOIL_TABLE}" \
    --p-sampling-depth "${SOIL_DEPTH}" \
    --p-metric shannon \
    --p-n "${N_ITERATIONS}" \
    --p-no-replacement \
    --o-alpha-diversities "${SOIL_OUT}/boots-alpha-shannon-collection.qza"

log_success "Alpha collection complete"

# Beta diversity collection (Weighted UniFrac)
# --i-phylogeny: required for UniFrac metrics
log_info "Running beta diversity (Weighted UniFrac) collection..."
qiime boots beta-collection \
    --i-table "${SOIL_TABLE}" \
    --i-phylogeny "${SOIL_TREE}" \
    --p-sampling-depth "${SOIL_DEPTH}" \
    --p-metric weighted_unifrac \
    --p-n "${N_ITERATIONS}" \
    --p-no-replacement \
    --o-distance-matrices "${SOIL_OUT}/boots-beta-wunifrac-collection.qza"

log_success "Beta collection complete"

# Compute average distance matrix
log_info "Computing average distance matrix (method: ${AVERAGE_METHOD})..."
qiime boots beta-average \
    --i-data "${SOIL_OUT}/boots-beta-wunifrac-collection.qza" \
    --p-average-method "${AVERAGE_METHOD}" \
    --o-average-distance-matrix "${SOIL_OUT}/boots-beta-wunifrac-avg.qza"

log_success "Soil q2-boots complete"

# ============================================================================
# 2. Bioethanol Dataset
# ============================================================================

log_header "2. Bioethanol Dataset: q2-boots (n=${N_ITERATIONS})"

BIO_DIR="${BASE_DIR}/Bioethanol_v2"
BIO_OUT="${BIO_DIR}/05_q2boots"
BIO_DEPTH=35000  # Primary depth from sensitivity analysis

log_info "Primary depth: ${BIO_DEPTH} reads"
log_info "Iterations: ${N_ITERATIONS}"
log_info "Resampling: rarefaction (without replacement)"

mkdir -p "${BIO_OUT}"

# Alpha diversity collection (Shannon)
log_info "Running alpha diversity (Shannon) collection..."
qiime boots alpha-collection \
    --i-table "${BIO_TABLE}" \
    --p-sampling-depth "${BIO_DEPTH}" \
    --p-metric shannon \
    --p-n "${N_ITERATIONS}" \
    --p-no-replacement \
    --o-alpha-diversities "${BIO_OUT}/boots-alpha-shannon-collection.qza"

log_success "Alpha collection complete"

# Beta diversity collection (Weighted UniFrac)
log_info "Running beta diversity (Weighted UniFrac) collection..."
qiime boots beta-collection \
    --i-table "${BIO_TABLE}" \
    --i-phylogeny "${BIO_TREE}" \
    --p-sampling-depth "${BIO_DEPTH}" \
    --p-metric weighted_unifrac \
    --p-n "${N_ITERATIONS}" \
    --p-no-replacement \
    --o-distance-matrices "${BIO_OUT}/boots-beta-wunifrac-collection.qza"

log_success "Beta collection complete"

# Compute average distance matrix
log_info "Computing average distance matrix (method: ${AVERAGE_METHOD})..."
qiime boots beta-average \
    --i-data "${BIO_OUT}/boots-beta-wunifrac-collection.qza" \
    --p-average-method "${AVERAGE_METHOD}" \
    --o-average-distance-matrix "${BIO_OUT}/boots-beta-wunifrac-avg.qza"

log_success "Bioethanol q2-boots complete"

# ============================================================================
# 3. Rice Dataset
# ============================================================================

log_header "3. Rice Dataset: q2-boots (n=${N_ITERATIONS})"

RICE_DIR="${BASE_DIR}/Rice_v2"
RICE_OUT="${RICE_DIR}/05_q2boots"
RICE_DEPTH=20000  # Primary depth from sensitivity analysis

log_info "Primary depth: ${RICE_DEPTH} reads"
log_info "Iterations: ${N_ITERATIONS}"
log_info "Resampling: rarefaction (without replacement)"
log_warn "Rice is a large dataset (n=493, ~80k ASVs) - this will take time"

mkdir -p "${RICE_OUT}"

# Alpha diversity collection (Shannon)
log_info "Running alpha diversity (Shannon) collection..."
qiime boots alpha-collection \
    --i-table "${RICE_TABLE}" \
    --p-sampling-depth "${RICE_DEPTH}" \
    --p-metric shannon \
    --p-n "${N_ITERATIONS}" \
    --p-no-replacement \
    --o-alpha-diversities "${RICE_OUT}/boots-alpha-shannon-collection.qza"

log_success "Alpha collection complete"

# Beta diversity collection (Weighted UniFrac)
log_info "Running beta diversity (Weighted UniFrac) collection..."
qiime boots beta-collection \
    --i-table "${RICE_TABLE}" \
    --i-phylogeny "${RICE_TREE}" \
    --p-sampling-depth "${RICE_DEPTH}" \
    --p-metric weighted_unifrac \
    --p-n "${N_ITERATIONS}" \
    --p-no-replacement \
    --o-distance-matrices "${RICE_OUT}/boots-beta-wunifrac-collection.qza"

log_success "Beta collection complete"

# Compute average distance matrix
log_info "Computing average distance matrix (method: ${AVERAGE_METHOD})..."
qiime boots beta-average \
    --i-data "${RICE_OUT}/boots-beta-wunifrac-collection.qza" \
    --p-average-method "${AVERAGE_METHOD}" \
    --o-average-distance-matrix "${RICE_OUT}/boots-beta-wunifrac-avg.qza"

log_success "Rice q2-boots complete"

# ============================================================================
# Summary
# ============================================================================

log_header "All q2-boots Analyses Complete"

echo ""
log_info "Execution parameters:"
echo "  Iterations: ${N_ITERATIONS}"
echo "  Averaging method: ${AVERAGE_METHOD}"
echo "  Resampling: rarefaction (without replacement)"
echo ""
log_info "Output files:"
echo ""
echo "Soil (depth=${SOIL_DEPTH}):"
echo "  ${SOIL_OUT}/"
echo "    - boots-alpha-shannon-collection.qza  (n=${N_ITERATIONS} alpha vectors)"
echo "    - boots-beta-wunifrac-collection.qza  (n=${N_ITERATIONS} distance matrices)"
echo "    - boots-beta-wunifrac-avg.qza         (averaged distance matrix)"
echo ""
echo "Bioethanol (depth=${BIO_DEPTH}):"
echo "  ${BIO_OUT}/"
echo "    - boots-alpha-shannon-collection.qza"
echo "    - boots-beta-wunifrac-collection.qza"
echo "    - boots-beta-wunifrac-avg.qza"
echo ""
echo "Rice (depth=${RICE_DEPTH}):"
echo "  ${RICE_OUT}/"
echo "    - boots-alpha-shannon-collection.qza"
echo "    - boots-beta-wunifrac-collection.qza"
echo "    - boots-beta-wunifrac-avg.qza"
echo ""
log_info "Next steps:"
echo "  1. Export alpha collections and compute per-sample CV:"
echo "     - qiime tools export --input-path boots-alpha-*.qza --output-path exported/"
echo "     - Calculate: CV = (std / mean) * 100 for each sample"
echo ""
echo "  2. Use averaged distance matrices for PCoA visualization:"
echo "     - qiime emperor plot --i-pcoa ... --m-metadata-file ..."
echo ""
echo "  3. Update Table 2 with Alpha CV values"
echo ""
echo "  4. Generate Figure 3 (repeated rarefaction stability)"
echo ""
echo "=== Script completed: $(date) ==="
