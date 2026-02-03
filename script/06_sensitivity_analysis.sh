#!/bin/bash
# ============================================================================
# Script: 06_sensitivity_analysis.sh
# Purpose: Run multi-depth sensitivity analysis for all benchmark datasets
# 
# This script evaluates the stability of diversity conclusions across multiple
# rarefaction depths, which is critical for datasets with potential depth-group
# confounding.
#
# For each depth, the script computes:
#   - Alpha diversity metrics (Shannon, Faith's PD, Observed Features)
#   - Beta diversity distance matrices (Weighted/Unweighted UniFrac, Bray-Curtis)
#   - Group significance tests (Kruskal-Wallis for alpha, PERMANOVA for beta)
#
# Prerequisites:
#   - QIIME 2 environment activated (qiime2-amplicon-2025.10)
#   - Alpha-rarefaction curves completed (to determine plateau depths)
#
# Usage:
#   conda activate qiime2-amplicon-2025.10
#   bash scripts/06_sensitivity_analysis.sh
#
# Background execution (recommended for long runs):
#   nohup ./scripts/06_sensitivity_analysis.sh > sensitivity_log.txt 2>&1 &
#
# Estimated time: 2-3 hours (Rice dataset takes longest)
#
# Output structure per dataset:
#   04_sensitivity_analysis/depth_{D}/
#     - core_metrics/           (alpha/beta diversity results)
#     - shannon_group_significance.qzv
#     - permanova_*.qzv         (beta diversity significance)
#
# Reference: Manuscript Section 5.1 - Step 2: Sensitivity analysis
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
# 1. Soil Dataset (POTENTIAL CONFOUNDING)
# ============================================================================

log_header "1. Soil Dataset: Sensitivity Analysis"

SOIL_DIR="${BASE_DIR}/Soil_v2"
SOIL_TABLE="${SOIL_DIR}/02_denoised/table.qza"
SOIL_TREE="${SOIL_DIR}/02_denoised/rooted-tree.qza"
SOIL_META="${METADATA_DIR}/soil/metadata_soil.tsv"
SOIL_OUT="${SOIL_DIR}/04_sensitivity_analysis"

# Depth parameters (based on plateau depth of 20,000)
# All depths retain all 18 samples (H group min = 36,834)
SOIL_DEPTHS=(15000 20000 25000 30000)

log_info "Diagnosis: POTENTIAL CONFOUNDING (η² = 0.28)"
log_info "Plateau depth: 20,000 reads"
log_info "Test depths: ${SOIL_DEPTHS[*]}"
log_info "Expected retention: 18/18 samples at all depths"

mkdir -p "${SOIL_OUT}"

for D in "${SOIL_DEPTHS[@]}"; do
    log_info "Processing Soil at depth: ${D}"
    
    OUTDIR="${SOIL_OUT}/depth_${D}"
    mkdir -p "${OUTDIR}"
    
    # Core metrics (alpha/beta diversity suite)
    qiime diversity core-metrics-phylogenetic \
        --i-table "${SOIL_TABLE}" \
        --i-phylogeny "${SOIL_TREE}" \
        --p-sampling-depth "${D}" \
        --m-metadata-file "${SOIL_META}" \
        --output-dir "${OUTDIR}/core_metrics"
    
    # Alpha diversity group significance (Shannon)
    qiime diversity alpha-group-significance \
        --i-alpha-diversity "${OUTDIR}/core_metrics/shannon_vector.qza" \
        --m-metadata-file "${SOIL_META}" \
        --o-visualization "${OUTDIR}/shannon_group_significance.qzv"
    
    # Alpha diversity group significance (Faith's PD)
    qiime diversity alpha-group-significance \
        --i-alpha-diversity "${OUTDIR}/core_metrics/faith_pd_vector.qza" \
        --m-metadata-file "${SOIL_META}" \
        --o-visualization "${OUTDIR}/faith_pd_group_significance.qzv"
    
    # Alpha diversity group significance (Observed Features)
    qiime diversity alpha-group-significance \
        --i-alpha-diversity "${OUTDIR}/core_metrics/observed_features_vector.qza" \
        --m-metadata-file "${SOIL_META}" \
        --o-visualization "${OUTDIR}/observed_features_group_significance.qzv"
    
    # PERMANOVA (Weighted UniFrac)
    qiime diversity beta-group-significance \
        --i-distance-matrix "${OUTDIR}/core_metrics/weighted_unifrac_distance_matrix.qza" \
        --m-metadata-file "${SOIL_META}" \
        --m-metadata-column treatment \
        --p-method permanova \
        --p-pairwise \
        --p-permutations 999 \
        --o-visualization "${OUTDIR}/permanova_weighted_unifrac.qzv"
    
    # PERMANOVA (Unweighted UniFrac)
    qiime diversity beta-group-significance \
        --i-distance-matrix "${OUTDIR}/core_metrics/unweighted_unifrac_distance_matrix.qza" \
        --m-metadata-file "${SOIL_META}" \
        --m-metadata-column treatment \
        --p-method permanova \
        --p-pairwise \
        --p-permutations 999 \
        --o-visualization "${OUTDIR}/permanova_unweighted_unifrac.qzv"
    
    # PERMANOVA (Bray-Curtis)
    qiime diversity beta-group-significance \
        --i-distance-matrix "${OUTDIR}/core_metrics/bray_curtis_distance_matrix.qza" \
        --m-metadata-file "${SOIL_META}" \
        --m-metadata-column treatment \
        --p-method permanova \
        --p-pairwise \
        --p-permutations 999 \
        --o-visualization "${OUTDIR}/permanova_bray_curtis.qzv"
    
    log_success "Soil depth ${D} complete"
done

log_success "Soil sensitivity analysis complete"

# ============================================================================
# 2. Bioethanol Dataset (time_point: POTENTIAL CONFOUNDING)
# ============================================================================

log_header "2. Bioethanol Dataset: Sensitivity Analysis"

BIO_DIR="${BASE_DIR}/Bioethanol_v2"
BIO_TABLE="${BIO_DIR}/02_denoised/table.qza"
BIO_TREE="${BIO_DIR}/02_denoised/rooted-tree.qza"
BIO_META="${METADATA_DIR}/bioethanol/metadata_bioethanol.tsv"
BIO_OUT="${BIO_DIR}/04_sensitivity_analysis"

# Depth parameters (based on plateau depth of ~35,714)
# Min depth 36,110, so 35,000 and below retains all samples
BIO_DEPTHS=(20000 25000 30000 35000)

log_info "Diagnosis: batch=NONE, time_point=POTENTIAL (η² = 0.12)"
log_info "Plateau depth: ~35,714 reads"
log_info "Test depths: ${BIO_DEPTHS[*]}"
log_info "Expected retention: 95/95 samples at all depths"

mkdir -p "${BIO_OUT}"

for D in "${BIO_DEPTHS[@]}"; do
    log_info "Processing Bioethanol at depth: ${D}"
    
    OUTDIR="${BIO_OUT}/depth_${D}"
    mkdir -p "${OUTDIR}"
    
    # Core metrics
    qiime diversity core-metrics-phylogenetic \
        --i-table "${BIO_TABLE}" \
        --i-phylogeny "${BIO_TREE}" \
        --p-sampling-depth "${D}" \
        --m-metadata-file "${BIO_META}" \
        --output-dir "${OUTDIR}/core_metrics"
    
    # Alpha diversity group significance (Shannon)
    qiime diversity alpha-group-significance \
        --i-alpha-diversity "${OUTDIR}/core_metrics/shannon_vector.qza" \
        --m-metadata-file "${BIO_META}" \
        --o-visualization "${OUTDIR}/shannon_group_significance.qzv"
    
    # PERMANOVA by batch (NO EVIDENCE, included for comparison)
    qiime diversity beta-group-significance \
        --i-distance-matrix "${OUTDIR}/core_metrics/weighted_unifrac_distance_matrix.qza" \
        --m-metadata-file "${BIO_META}" \
        --m-metadata-column batch \
        --p-method permanova \
        --p-permutations 999 \
        --o-visualization "${OUTDIR}/permanova_wunifrac_batch.qzv"
    
    # PERMANOVA by time_point (POTENTIAL CONFOUNDING - primary focus)
    qiime diversity beta-group-significance \
        --i-distance-matrix "${OUTDIR}/core_metrics/weighted_unifrac_distance_matrix.qza" \
        --m-metadata-file "${BIO_META}" \
        --m-metadata-column time_point \
        --p-method permanova \
        --p-permutations 999 \
        --o-visualization "${OUTDIR}/permanova_wunifrac_timepoint.qzv"
    
    # PERMANOVA (Bray-Curtis by time_point)
    qiime diversity beta-group-significance \
        --i-distance-matrix "${OUTDIR}/core_metrics/bray_curtis_distance_matrix.qza" \
        --m-metadata-file "${BIO_META}" \
        --m-metadata-column time_point \
        --p-method permanova \
        --p-permutations 999 \
        --o-visualization "${OUTDIR}/permanova_bc_timepoint.qzv"
    
    log_success "Bioethanol depth ${D} complete"
done

log_success "Bioethanol sensitivity analysis complete"

# ============================================================================
# 3. Rice Dataset (MINOR CONFOUNDING)
# ============================================================================

log_header "3. Rice Dataset: Sensitivity Analysis"

RICE_DIR="${BASE_DIR}/Rice_v2"
RICE_TABLE="${RICE_DIR}/02_denoised/table.qza"
RICE_TREE="${RICE_DIR}/02_denoised/rooted-tree.qza"
RICE_META="${METADATA_DIR}/rice/metadata_rice.tsv"
RICE_OUT="${RICE_DIR}/04_sensitivity_analysis"

# Depth parameters (based on plateau depth of ~22,714)
# At 20,000: 488/493 retained; at 22,714: 486/493 retained
RICE_DEPTHS=(10000 15000 20000 22714)

log_info "Diagnosis: MINOR CONFOUNDING (η² = 0.02)"
log_info "Plateau depth: ~22,714 reads"
log_info "Test depths: ${RICE_DEPTHS[*]}"
log_info "Retention at 20,000: 488/493 (5 dropout)"
log_warn "Rice is a large dataset (n=493, ~80k ASVs) - this will take time"

mkdir -p "${RICE_OUT}"

for D in "${RICE_DEPTHS[@]}"; do
    log_info "Processing Rice at depth: ${D}"
    
    OUTDIR="${RICE_OUT}/depth_${D}"
    mkdir -p "${OUTDIR}"
    
    # Core metrics
    qiime diversity core-metrics-phylogenetic \
        --i-table "${RICE_TABLE}" \
        --i-phylogeny "${RICE_TREE}" \
        --p-sampling-depth "${D}" \
        --m-metadata-file "${RICE_META}" \
        --output-dir "${OUTDIR}/core_metrics"
    
    # Alpha diversity group significance (Shannon)
    qiime diversity alpha-group-significance \
        --i-alpha-diversity "${OUTDIR}/core_metrics/shannon_vector.qza" \
        --m-metadata-file "${RICE_META}" \
        --o-visualization "${OUTDIR}/shannon_group_significance.qzv"
    
    # Alpha diversity group significance (Faith's PD)
    qiime diversity alpha-group-significance \
        --i-alpha-diversity "${OUTDIR}/core_metrics/faith_pd_vector.qza" \
        --m-metadata-file "${RICE_META}" \
        --o-visualization "${OUTDIR}/faith_pd_group_significance.qzv"
    
    # PERMANOVA by compartment (Weighted UniFrac)
    qiime diversity beta-group-significance \
        --i-distance-matrix "${OUTDIR}/core_metrics/weighted_unifrac_distance_matrix.qza" \
        --m-metadata-file "${RICE_META}" \
        --m-metadata-column compartment \
        --p-method permanova \
        --p-pairwise \
        --p-permutations 999 \
        --o-visualization "${OUTDIR}/permanova_wunifrac_compartment.qzv"
    
    # PERMANOVA by compartment (Unweighted UniFrac)
    qiime diversity beta-group-significance \
        --i-distance-matrix "${OUTDIR}/core_metrics/unweighted_unifrac_distance_matrix.qza" \
        --m-metadata-file "${RICE_META}" \
        --m-metadata-column compartment \
        --p-method permanova \
        --p-pairwise \
        --p-permutations 999 \
        --o-visualization "${OUTDIR}/permanova_uunifrac_compartment.qzv"
    
    log_success "Rice depth ${D} complete"
done

log_success "Rice sensitivity analysis complete"

# ============================================================================
# Summary
# ============================================================================

log_header "All Sensitivity Analyses Complete"

echo ""
log_info "Output structure:"
echo ""
echo "Soil (POTENTIAL CONFOUNDING):"
echo "  ${SOIL_OUT}/"
for D in "${SOIL_DEPTHS[@]}"; do
    echo "    depth_${D}/"
done
echo ""
echo "Bioethanol (time_point: POTENTIAL CONFOUNDING):"
echo "  ${BIO_OUT}/"
for D in "${BIO_DEPTHS[@]}"; do
    echo "    depth_${D}/"
done
echo ""
echo "Rice (MINOR CONFOUNDING):"
echo "  ${RICE_OUT}/"
for D in "${RICE_DEPTHS[@]}"; do
    echo "    depth_${D}/"
done
echo ""
log_info "Next steps:"
echo "  1. Open .qzv files at https://view.qiime2.org/"
echo "  2. Compare metrics across depths:"
echo "     - Alpha: Kruskal-Wallis H, p-value, group medians"
echo "     - Beta: PERMANOVA pseudo-F, R², p-value"
echo "  3. Evaluate conclusion stability"
echo "  4. Proceed to q2-boots analysis"
echo ""
echo "  bash scripts/07_q2boots_analysis.sh"
