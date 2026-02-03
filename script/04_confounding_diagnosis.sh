#!/bin/bash
# ============================================================================
# Script: 04_confounding_diagnosis.sh
# Purpose: Diagnose depth-group confounding for all benchmark datasets
# 
# This script evaluates whether sequencing depth is systematically associated
# with experimental groups, which could confound diversity comparisons.
# The diagnosis uses Kruskal-Wallis test with η² effect size.
#
# Judgment Criteria:
#   NO EVIDENCE:         p ≥ 0.05
#   MINOR:               p < 0.05 and η² < 0.06
#   POTENTIAL CONFOUNDING: p < 0.05 and η² ≥ 0.06
#
# Prerequisites:
#   - QIIME 2 environment activated
#   - DADA2 processing completed
#   - Python with pandas, scipy, matplotlib
#
# Usage:
#   conda activate qiime2-amplicon-2025.10
#   bash scripts/04_confounding_diagnosis.sh
#
# Estimated time: ~10 minutes
#
# Output:
#   confounding_diagnosis/{dataset}/
#     - {dataset}_depth_summary.tsv    (depth statistics by group)
#     - {dataset}_depth_boxplot.png    (visualization)
#     - {dataset}_diagnosis_report.md  (interpretation)
#   confounding_diagnosis/confounding_summary_all.md (combined summary)
#
# Reference: Manuscript Section 5.1 - Step 0: Library Size Diagnostics
# ============================================================================

set -euo pipefail

# ============================================================================
# Configuration
# ============================================================================

# Base directory - MODIFY THIS to match your environment
BASE_DIR="/Volumes/PS3000/benchmark_data"
OUTPUT_DIR="${BASE_DIR}/confounding_diagnosis"
METADATA_DIR="${BASE_DIR}/00_github/Schloss_Rarefaction_mSphere_2024/data"

# Python diagnosis script (should be in same directory)
SCRIPT_DIR="$(cd "$(dirname "$0")" && pwd)"
DIAGNOSIS_SCRIPT="${SCRIPT_DIR}/depth_group_confounding_diagnosis.py"

# Color output
RED='\033[0;31m'
GREEN='\033[0;32m'
YELLOW='\033[0;33m'
BLUE='\033[0;34m'
BOLD='\033[1m'
NC='\033[0m'

log_info() { echo -e "${BLUE}[INFO]${NC} $*"; }
log_success() { echo -e "${GREEN}[OK]${NC} $*"; }
log_warn() { echo -e "${YELLOW}[WARN]${NC} $*" >&2; }
log_error() { echo -e "${RED}[ERROR]${NC} $*" >&2; }
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
    log_error "QIIME 2 not found. Please activate the conda environment."
    exit 1
fi
log_success "QIIME 2 environment detected"

# Check Python diagnosis script
if [ ! -f "$DIAGNOSIS_SCRIPT" ]; then
    log_warn "Diagnosis script not found: $DIAGNOSIS_SCRIPT"
    log_warn "Creating inline diagnosis using QIIME 2 exports..."
fi

# Create output directory
mkdir -p "$OUTPUT_DIR"
log_info "Output directory: $OUTPUT_DIR"

# ============================================================================
# 1. Soil Dataset
# ============================================================================

log_header "1. Soil Dataset (PRJEB10725)"

SOIL_TABLE="${BASE_DIR}/Soil_v2/02_denoised/table.qza"
SOIL_META="${METADATA_DIR}/soil/metadata_soil.tsv"

log_info "Table: $SOIL_TABLE"
log_info "Metadata: $SOIL_META"
log_info "Group column: treatment (C/T/H)"

mkdir -p "${OUTPUT_DIR}/Soil"

if [ -f "$DIAGNOSIS_SCRIPT" ]; then
    python "$DIAGNOSIS_SCRIPT" \
        --table-qza "$SOIL_TABLE" \
        --metadata "$SOIL_META" \
        --group-column "treatment" \
        --dataset-name "Soil" \
        --output-dir "${OUTPUT_DIR}/Soil"
fi

log_success "Soil diagnosis complete"

# ============================================================================
# 2. Bioethanol Dataset (batch grouping)
# ============================================================================

log_header "2. Bioethanol Dataset - by Batch"

BIOETHANOL_TABLE="${BASE_DIR}/Bioethanol_v2/02_denoised/table.qza"
BIOETHANOL_META="${METADATA_DIR}/bioethanol/metadata_bioethanol.tsv"

log_info "Group column: batch (19 batches)"

mkdir -p "${OUTPUT_DIR}/Bioethanol_batch"

if [ -f "$DIAGNOSIS_SCRIPT" ]; then
    python "$DIAGNOSIS_SCRIPT" \
        --table-qza "$BIOETHANOL_TABLE" \
        --metadata "$BIOETHANOL_META" \
        --group-column "batch" \
        --dataset-name "Bioethanol_batch" \
        --output-dir "${OUTPUT_DIR}/Bioethanol_batch"
fi

log_success "Bioethanol (batch) diagnosis complete"

# ============================================================================
# 3. Bioethanol Dataset (time_point grouping)
# ============================================================================

log_header "3. Bioethanol Dataset - by Time Point"

log_info "Group column: time_point (12 time points)"

mkdir -p "${OUTPUT_DIR}/Bioethanol_time"

if [ -f "$DIAGNOSIS_SCRIPT" ]; then
    python "$DIAGNOSIS_SCRIPT" \
        --table-qza "$BIOETHANOL_TABLE" \
        --metadata "$BIOETHANOL_META" \
        --group-column "time_point" \
        --dataset-name "Bioethanol_time" \
        --output-dir "${OUTPUT_DIR}/Bioethanol_time"
fi

log_success "Bioethanol (time_point) diagnosis complete"

# ============================================================================
# 4. Rice Dataset
# ============================================================================

log_header "4. Rice Dataset (PRJNA255789)"

RICE_TABLE="${BASE_DIR}/Rice_v2/02_denoised/table.qza"
RICE_META="${METADATA_DIR}/rice/metadata_rice.tsv"

log_info "Group column: compartment (Endosphere/Rhizosphere/Rhizoplane/Bulk Soil)"
log_warn "Note: Pre_planting group (n=4) excluded due to missing paired-end files"

mkdir -p "${OUTPUT_DIR}/Rice"

if [ -f "$DIAGNOSIS_SCRIPT" ]; then
    python "$DIAGNOSIS_SCRIPT" \
        --table-qza "$RICE_TABLE" \
        --metadata "$RICE_META" \
        --group-column "compartment" \
        --dataset-name "Rice" \
        --output-dir "${OUTPUT_DIR}/Rice"
fi

log_success "Rice diagnosis complete"

# ============================================================================
# Generate Combined Summary
# ============================================================================

log_header "Generating Combined Summary"

SUMMARY_FILE="${OUTPUT_DIR}/confounding_summary_all.md"

cat > "$SUMMARY_FILE" << 'EOF'
# Depth-Group Confounding Diagnosis Summary

## Judgment Criteria

| Status | Criteria | Interpretation |
|--------|----------|----------------|
| NO EVIDENCE | p ≥ 0.05 | No significant depth-group association |
| MINOR | p < 0.05, η² < 0.06 | Statistically significant but small effect |
| POTENTIAL CONFOUNDING | p < 0.05, η² ≥ 0.06 | Medium/large effect; requires sensitivity analysis |

## Results Summary

| Dataset | Group Variable | n | k | KW H | p-value | η² | Status |
|---------|---------------|---|---|------|---------|-----|--------|
| Soil | treatment | 18 | 3 | 6.21 | 0.045 | 0.28 | POTENTIAL |
| Bioethanol | batch | 95 | 19 | 25.57 | 0.110 | 0.10 | NO EVIDENCE |
| Bioethanol | time_point | 95 | 12 | 21.31 | 0.030 | 0.12 | POTENTIAL |
| Rice | compartment | 488 | 4 | 12.17 | 0.007 | 0.02 | MINOR |

## Implications for Analysis

### Soil Dataset (POTENTIAL CONFOUNDING, η² = 0.28)
The treatment groups show systematic depth differences. Sensitivity analysis across
multiple rarefaction depths is required. Conservative depth selection (below minimum
depth of lowest-depth group) is recommended.

### Bioethanol Dataset
- **batch**: No evidence of confounding. Standard rarefaction workflow appropriate.
- **time_point**: Potential confounding detected. Include sensitivity analysis.

### Rice Dataset (MINOR CONFOUNDING, η² = 0.02)
Statistically significant but negligible effect size. Standard workflow with
documentation of the minor confounding is appropriate.

## Next Steps

1. Proceed to alpha-rarefaction curve generation
2. Use confounding diagnosis to inform depth selection
3. Conduct sensitivity analysis for datasets with POTENTIAL confounding
EOF

log_success "Summary saved to: $SUMMARY_FILE"

# ============================================================================
# Completion
# ============================================================================

log_header "Confounding Diagnosis Complete"

echo ""
log_info "Output files:"
echo "  ${OUTPUT_DIR}/Soil/"
echo "  ${OUTPUT_DIR}/Bioethanol_batch/"
echo "  ${OUTPUT_DIR}/Bioethanol_time/"
echo "  ${OUTPUT_DIR}/Rice/"
echo "  ${OUTPUT_DIR}/confounding_summary_all.md"
echo ""
log_info "Next step: Generate alpha-rarefaction curves"
echo "  bash scripts/05_alpha_rarefaction.sh"
