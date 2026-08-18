# QIIME 2 Provenance Documentation

This document explains how QIIME 2's built-in provenance tracking system supports the reproducibility of this analysis, and states exactly which primary QIIME 2 outputs are archived in this repository.

---

## Archived Primary Outputs in This Repository

The complete set of intermediate QIIME 2 artifacts for the three datasets is ~15 GB and is not distributable in a repository; every artifact is regenerable from the archived scripts, pinned environment, and fixed seeds (see the repository README). What **is** archived here is the primary evidence for the manuscript's central beta-diversity stability claim (Table 2, Rice; Supplementary Table S1 (D)): the PERMANOVA visualizations for the Rice dataset at all four sensitivity-analysis rarefaction depths.

`rice_sensitivity_permanova/depth_{10000,15000,20000,22714}/permanova_wunifrac_compartment.qzv`
(weighted UniFrac, compartment, 999 permutations):

| Depth | Samples (n) | pseudo-F | R² | Artifact UUID | sha256 |
|---|---|---|---|---|---|
| 10,000 | 492 | 165.488004 | 0.504 | `8add901a-626a-4867-81a7-e7484661a139` | `84d79258b516…c9f7` |
| 15,000 | 490 | 165.015866 | 0.505 | `270caebf-55d9-4ce2-8298-d796890ebd0c` | `1acea25bc60a…6ab1` |
| 20,000 | 488 | 164.840491 | 0.505 | `ad70edf2-a284-4e29-bca3-324052df1e03` | `95eb02b16e8b…efbe` |
| 22,714 | 486 | 165.417978 | 0.507 | `00fe4f4d-3595-4b51-a8a1-92ceb60f2058` | `b1fbeac8d19c…50f4` |

(R² = F·(k−1) / (F·(k−1) + N−k) with k = 4 compartments; full 64-character checksums are reproducible with `shasum -a 256`.)

Each file can be dropped onto [view.qiime2.org](https://view.qiime2.org) to inspect the full computational lineage — input table checksums, rarefaction seed and depth, distance metric, software versions, and timestamps — which is what makes these four files sufficient anchors for third-party verification: the embedded provenance ties them back to the public FASTQ accessions through every intermediate step.

---

## What is QIIME 2 Provenance?

QIIME 2 automatically records comprehensive metadata about every computational step in a decentralized provenance system. Each QIIME 2 artifact (.qza) and visualization (.qzv) file contains a complete record of its computational history, including all ancestor artifacts, software versions, and parameters.

This design means that reproducibility information is never separated from the data it describes, and any artifact can be fully traced back to its origins.

---

## Information Captured in Provenance

### Software Environment

Each artifact records the exact versions of all software components used in its creation, including QIIME 2 framework version, plugin name and version, Python version, and all dependency versions.

### Computational Parameters

Every parameter passed to every command is recorded, including user-specified parameters, default values that were used, and random seeds (when applicable).

### Input Lineage

The complete history of input artifacts is recorded, forming a directed acyclic graph (DAG) that traces back to the original imported data. This allows full reconstruction of the analysis pipeline.

### Execution Metadata

Timestamps for when each action was performed are recorded, along with the unique identifier (UUID) for each artifact.

---

## Viewing Provenance

### Using QIIME 2 View (Recommended)

The easiest way to inspect provenance is through the QIIME 2 View web interface.

1. Navigate to https://view.qiime2.org
2. Drag and drop any .qzv file onto the page
3. Click the "Provenance" tab to explore the computational history

The provenance viewer provides an interactive graph showing all computational steps and their relationships. Clicking on any node reveals the complete details for that step.

### Using Command Line

You can also export provenance information using the QIIME 2 command line interface:

```bash
# Export provenance as a ZIP archive
qiime tools export \
  --input-path artifact.qza \
  --output-path exported_provenance/

# The provenance directory contains:
# - action.yaml (parameters and software info)
# - citations.bib (relevant citations)
# - VERSION (QIIME 2 archive format version)
```

---

## Key Artifacts in This Analysis

The following table lists the key QIIME 2 artifacts generated during the analysis and their locations.

### Phase 1: DADA2 Processing

| Artifact | Description | Location |
|----------|-------------|----------|
| `table.qza` | Feature table (ASV counts) | `{dataset}/02_denoised/` |
| `rep-seqs.qza` | Representative sequences | `{dataset}/02_denoised/` |
| `denoising-stats.qza` | DADA2 processing statistics | `{dataset}/02_denoised/` |
| `rooted-tree.qza` | Phylogenetic tree | `{dataset}/02_denoised/` |

### Phase 2: Diversity Analysis

| Artifact | Description | Location |
|----------|-------------|----------|
| `alpha-rarefaction.qzv` | Alpha-rarefaction curves | `{dataset}/03_diversity/` |
| `core-metrics-results/` | Directory containing diversity metrics | `{dataset}/03_diversity/` |
| `shannon_vector.qza` | Shannon diversity values | `{dataset}/03_diversity/core-metrics-results/` |
| `weighted_unifrac_distance_matrix.qza` | Beta diversity matrix | `{dataset}/03_diversity/core-metrics-results/` |

### Phase 3: q2-boots Analysis

| Artifact | Description | Location |
|----------|-------------|----------|
| `boots-alpha-shannon.qza` | Repeated rarefaction alpha results | `{dataset}/03_diversity/` |
| `boots-beta-wunifrac.qza` | Repeated rarefaction beta results | `{dataset}/03_diversity/` |

---

## Reproducing the Analysis from Provenance

If you have access to only the final artifacts (without the original scripts), you can reconstruct the complete analysis pipeline from provenance.

### Step 1: Extract the Action Graph

Upload the final artifact to view.qiime2.org and navigate to the Provenance tab. The graph shows all computational steps in order.

### Step 2: Extract Parameters for Each Step

Click on each node in the provenance graph to see the exact parameters used, then reconstruct the corresponding QIIME 2 command.

### Step 3: Re-execute

Execute the reconstructed commands in order, using the same input data and parameters.

---

## Provenance Limitations

While QIIME 2 provenance is comprehensive for QIIME 2 operations, it does have some limitations to be aware of.

### External Tools

Operations performed outside QIIME 2 (such as data preprocessing with custom scripts) are not captured in artifact provenance. This is why we provide separate scripts for data download and preprocessing.

### System Environment

Hardware specifications, operating system details, and conda environment state are not fully captured. The environment.yml file in this repository addresses this gap.

### Random Number Generation

While random seeds are recorded, exact reproducibility of stochastic operations may depend on software versions and hardware. We recommend using the same QIIME 2 version specified in environment.yml.

---

## Best Practices Used in This Project

To maximize reproducibility, this project follows several best practices.

### Artifact Preservation

Key intermediate artifacts are preserved, not just final results. This allows verification of any intermediate step.

### Parameter Documentation

All non-default parameters are explicitly documented in both scripts and the README.

### Version Pinning

Exact software versions are specified in environment.yml to ensure consistent behavior.

### Seed Specification

Random seeds are explicitly set where applicable to enable exact reproduction of stochastic results.

---

## References

Bolyen E, et al. Reproducible, interactive, scalable and extensible microbiome data science using QIIME 2. *Nature Biotechnology*. 2019;37:852–857.

QIIME 2 Provenance Documentation: https://docs.qiime2.org/
