#!/bin/bash
# ============================================================================
# Script: run_all.sh
# Purpose: Execute the complete analysis pipeline for the rarefaction benchmark
# 
# This master script runs all analysis steps sequentially:
#   1. Download FASTQ files from ENA
#   2. DADA2 denoising
#   3. Phylogenetic tree construction
#   4. Depth-group confounding diagnosis
#   5. Alpha-rarefaction curve generation
#   6. Multi-depth sensitivity analysis
#   7. q2-boots repeated rarefaction
#
# Prerequisites:
#   - QIIME 2 environment activated (qiime2-amplicon-2025.10)
#   - q2-boots plugin installed
#   - ~50 GB disk space for raw data and outputs
#   - Stable internet connection (for data download)
#
# Usage:
#   conda activate qiime2-amplicon-2025.10
#   bash scripts/run_all.sh
#
# Background execution (recommended for full pipeline):
#   nohup caffeinate -ims bash scripts/run_all.sh > pipeline_log.txt 2>&1 &
#
# Estimated total time: 8-12 hours
#   - Download: 2-4 hours
#   - DADA2: 2-3 hours
#   - Phylogeny: 30-60 minutes
#   - Confounding diagnosis: 10 minutes
#   - Alpha-rarefaction: 30-60 minutes
#   - Sensitivity analysis: 2-3 hours
#   - q2-boots: 2-3 hours
#
# Note: Individual scripts can be run separately if needed.
# ============================================================================

set -euo pipefail

# ============================================================================
# Configuration
# ============================================================================

SCRIPT_DIR="$(cd "$(dirname "$0")" && pwd)"
LOG_DIR="${SCRIPT_DIR}/../logs"
mkdir -p "${LOG_DIR}"

TIMESTAMP=$(date '+%Y%m%d_%H%M%S')
MAIN_LOG="${LOG_DIR}/pipeline_${TIMESTAMP}.log"

# Color output
RED='\033[0;31m'
GREEN='\033[0;32m'
YELLOW='\033[0;33m'
BLUE='\033[0;34m'
BOLD='\033[1m'
NC='\033[0m'

log_info() { echo -e "${BLUE}[INFO]${NC} $(date '+%Y-%m-%d %H:%M:%S') $*" | tee -a "${MAIN_LOG}"; }
log_success() { echo -e "${GREEN}[OK]${NC} $(date '+%Y-%m-%d %H:%M:%S') $*" | tee -a "${MAIN_LOG}"; }
log_warn() { echo -e "${YELLOW}[WARN]${NC} $(date '+%Y-%m-%d %H:%M:%S') $*" | tee -a "${MAIN_LOG}"; }
log_error() { echo -e "${RED}[ERROR]${NC} $(date '+%Y-%m-%d %H:%M:%S') $*" | tee -a "${MAIN_LOG}"; }
log_header() {
    echo "" | tee -a "${MAIN_LOG}"
    echo -e "${BOLD}============================================================${NC}" | tee -a "${MAIN_LOG}"
    echo -e "${BOLD}  $*${NC}" | tee -a "${MAIN_LOG}"
    echo -e "${BOLD}============================================================${NC}" | tee -a "${MAIN_LOG}"
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
log_success "QIIME 2 environment detected"

# Check q2-boots
if ! qiime boots --help &> /dev/null; then
    log_error "q2-boots not found. Please install the plugin."
    exit 1
fi
log_success "q2-boots plugin detected"

# Check required commands
for cmd in wget curl awk; do
    if ! command -v "$cmd" &> /dev/null; then
        log_error "Required command not found: $cmd"
        exit 1
    fi
done
log_success "All required commands available"

log_info "Log file: ${MAIN_LOG}"
echo ""

# ============================================================================
# Pipeline Execution
# ============================================================================

PIPELINE_START=$(date +%s)

# Step 1: Download FASTQ files
log_header "Step 1/7: Download FASTQ Files"
log_info "Estimated time: 2-4 hours"
bash "${SCRIPT_DIR}/01_download_data.sh" 2>&1 | tee -a "${MAIN_LOG}"
log_success "Step 1 complete"

# Step 2: DADA2 processing
log_header "Step 2/7: DADA2 Denoising"
log_info "Estimated time: 2-3 hours"
bash "${SCRIPT_DIR}/02_dada2_processing.sh" 2>&1 | tee -a "${MAIN_LOG}"
log_success "Step 2 complete"

# Step 3: Phylogenetic tree construction
log_header "Step 3/7: Phylogenetic Tree Construction"
log_info "Estimated time: 30-60 minutes"
bash "${SCRIPT_DIR}/03_build_phylogeny.sh" 2>&1 | tee -a "${MAIN_LOG}"
log_success "Step 3 complete"

# Step 4: Confounding diagnosis
log_header "Step 4/7: Depth-Group Confounding Diagnosis"
log_info "Estimated time: 10 minutes"
bash "${SCRIPT_DIR}/04_confounding_diagnosis.sh" 2>&1 | tee -a "${MAIN_LOG}"
log_success "Step 4 complete"

# Step 5: Alpha-rarefaction curves
log_header "Step 5/7: Alpha-Rarefaction Curves"
log_info "Estimated time: 30-60 minutes"
bash "${SCRIPT_DIR}/05_alpha_rarefaction.sh" 2>&1 | tee -a "${MAIN_LOG}"
log_success "Step 5 complete"

# Step 6: Sensitivity analysis
log_header "Step 6/7: Multi-Depth Sensitivity Analysis"
log_info "Estimated time: 2-3 hours"
bash "${SCRIPT_DIR}/06_sensitivity_analysis.sh" 2>&1 | tee -a "${MAIN_LOG}"
log_success "Step 6 complete"

# Step 7: q2-boots analysis
log_header "Step 7/7: q2-boots Repeated Rarefaction"
log_info "Estimated time: 2-3 hours"
bash "${SCRIPT_DIR}/07_q2boots_analysis.sh" 2>&1 | tee -a "${MAIN_LOG}"
log_success "Step 7 complete"

# ============================================================================
# Pipeline Summary
# ============================================================================

PIPELINE_END=$(date +%s)
PIPELINE_DURATION=$((PIPELINE_END - PIPELINE_START))
HOURS=$((PIPELINE_DURATION / 3600))
MINUTES=$(((PIPELINE_DURATION % 3600) / 60))

log_header "Pipeline Complete"

echo "" | tee -a "${MAIN_LOG}"
log_success "All 7 steps completed successfully!"
echo "" | tee -a "${MAIN_LOG}"
log_info "Total execution time: ${HOURS}h ${MINUTES}m"
log_info "Log file: ${MAIN_LOG}"
echo "" | tee -a "${MAIN_LOG}"
log_info "Next steps:"
echo "  1. Review confounding diagnosis results" | tee -a "${MAIN_LOG}"
echo "  2. Examine alpha-rarefaction curves at https://view.qiime2.org/" | tee -a "${MAIN_LOG}"
echo "  3. Compare sensitivity analysis results across depths" | tee -a "${MAIN_LOG}"
echo "  4. Calculate Alpha CV from q2-boots collections" | tee -a "${MAIN_LOG}"
echo "  5. Generate figures and tables for manuscript" | tee -a "${MAIN_LOG}"
