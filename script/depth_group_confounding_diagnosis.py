#!/usr/bin/env python3
"""
Depth–Group Confounding Diagnosis Script
=========================================
Purpose: Diagnose depth–group confounding for rarefaction review paper benchmark

This script extracts depth information from QIIME2 FeatureTable[Frequency] artifacts
and tests whether sequencing depth is systematically associated with experimental groups.

Usage:
    python depth_group_confounding_diagnosis.py \
        --table-qza /path/to/table.qza \
        --metadata /path/to/metadata.tsv \
        --group-column treatment \
        --output-dir /path/to/output

Outputs:
    - depth_summary.tsv: Descriptive statistics by group
    - depth_boxplot.png: Visualization of depth distribution by group
    - confounding_test.tsv: Statistical test results
    - diagnosis_report.md: Human-readable summary

Author: For rarefaction review paper benchmark (Schloss 2024 datasets)
Date: 2026-01
"""

import argparse
import os
import tempfile
import subprocess
from pathlib import Path

import pandas as pd
import numpy as np
from scipy import stats
import matplotlib.pyplot as plt
import seaborn as sns

# ==============================================================================
# Helper Functions
# ==============================================================================

def extract_depth_from_qza(table_qza_path: str) -> pd.Series:
    """
    Extract per-sample depth (total read count) from a QIIME2 FeatureTable[Frequency] artifact.
    
    Parameters
    ----------
    table_qza_path : str
        Path to table.qza file
    
    Returns
    -------
    pd.Series
        Series with sample-id as index and depth as values
    """
    with tempfile.TemporaryDirectory() as tmpdir:
        # Export BIOM table from QZA
        subprocess.run([
            "qiime", "tools", "export",
            "--input-path", table_qza_path,
            "--output-path", tmpdir
        ], check=True, capture_output=True)
        
        # Convert BIOM to TSV
        biom_path = os.path.join(tmpdir, "feature-table.biom")
        tsv_path = os.path.join(tmpdir, "feature-table.tsv")
        subprocess.run([
            "biom", "convert",
            "-i", biom_path,
            "-o", tsv_path,
            "--to-tsv"
        ], check=True, capture_output=True)
        
        # Read and compute depth
        df = pd.read_csv(tsv_path, sep="\t", skiprows=1, index_col=0)
        depth = df.sum(axis=0)
        depth.name = "depth"
        depth.index.name = "sample-id"
        
    return depth


def cliffs_delta(x: np.ndarray, y: np.ndarray) -> tuple:
    """
    Calculate Cliff's delta effect size for two groups.
    
    Interpretation:
        |d| < 0.147: negligible
        |d| < 0.33: small
        |d| < 0.474: medium
        |d| >= 0.474: large
    
    Parameters
    ----------
    x, y : np.ndarray
        Arrays of values for two groups
    
    Returns
    -------
    tuple
        (delta, interpretation)
    """
    n_x, n_y = len(x), len(y)
    
    # Count dominance
    greater = sum((xi > yj for xi in x for yj in y))
    less = sum((xi < yj for xi in x for yj in y))
    
    delta = (greater - less) / (n_x * n_y)
    
    # Interpretation
    abs_delta = abs(delta)
    if abs_delta < 0.147:
        interp = "negligible"
    elif abs_delta < 0.33:
        interp = "small"
    elif abs_delta < 0.474:
        interp = "medium"
    else:
        interp = "large"
    
    return delta, interp


def eta_squared_kruskal(h_stat: float, n: int, k: int) -> float:
    """
    Calculate eta-squared effect size for Kruskal-Wallis test.
    
    Formula: eta² = (H - k + 1) / (n - k)
    
    Interpretation:
        eta² < 0.01: small
        eta² < 0.06: medium
        eta² >= 0.14: large
    
    Parameters
    ----------
    h_stat : float
        Kruskal-Wallis H statistic
    n : int
        Total sample size
    k : int
        Number of groups
    
    Returns
    -------
    float
        eta-squared value
    """
    eta2 = (h_stat - k + 1) / (n - k)
    return max(0, eta2)  # Clip to non-negative


def diagnose_depth_group_confounding(
    depth: pd.Series,
    metadata: pd.DataFrame,
    group_column: str,
    dataset_name: str,
    output_dir: str
) -> dict:
    """
    Perform comprehensive depth–group confounding diagnosis.
    
    Parameters
    ----------
    depth : pd.Series
        Per-sample depth
    metadata : pd.DataFrame
        Sample metadata with group_column
    group_column : str
        Column name for grouping variable
    dataset_name : str
        Name of the dataset (for output filenames)
    output_dir : str
        Directory for output files
    
    Returns
    -------
    dict
        Dictionary containing all diagnosis results
    """
    # Ensure output directory exists
    os.makedirs(output_dir, exist_ok=True)
    
    # Merge depth with metadata
    merged = metadata.copy()
    merged["depth"] = depth
    merged = merged.dropna(subset=[group_column, "depth"])
    
    # Group-wise statistics
    groups = merged.groupby(group_column)["depth"]
    
    summary_stats = groups.agg([
        ("n", "count"),
        ("mean", "mean"),
        ("median", "median"),
        ("std", "std"),
        ("min", "min"),
        ("max", "max"),
        ("q25", lambda x: x.quantile(0.25)),
        ("q75", lambda x: x.quantile(0.75))
    ]).round(1)
    
    # Log-transformed statistics (more relevant for comparing distributions)
    log_depth = np.log10(merged["depth"])
    merged["log10_depth"] = log_depth
    
    log_summary = merged.groupby(group_column)["log10_depth"].agg([
        ("log10_mean", "mean"),
        ("log10_std", "std")
    ]).round(3)
    
    summary_stats = summary_stats.join(log_summary)
    
    # Save summary statistics
    summary_path = os.path.join(output_dir, f"{dataset_name}_depth_summary.tsv")
    summary_stats.to_csv(summary_path, sep="\t")
    
    # Statistical tests
    group_names = list(groups.groups.keys())
    group_values = [groups.get_group(g).values for g in group_names]
    n_groups = len(group_names)
    n_total = len(merged)
    
    results = {
        "dataset": dataset_name,
        "group_column": group_column,
        "n_groups": n_groups,
        "n_samples": n_total,
        "summary_stats": summary_stats
    }
    
    # Kruskal-Wallis test (non-parametric ANOVA)
    if n_groups >= 2:
        h_stat, p_kw = stats.kruskal(*group_values)
        eta2 = eta_squared_kruskal(h_stat, n_total, n_groups)
        
        if eta2 < 0.01:
            eta2_interp = "small (negligible)"
        elif eta2 < 0.06:
            eta2_interp = "medium"
        else:
            eta2_interp = "large"
        
        results["kruskal_wallis"] = {
            "H": h_stat,
            "p_value": p_kw,
            "eta_squared": eta2,
            "eta_squared_interpretation": eta2_interp
        }
    
    # Pairwise comparisons (Cliff's delta) for 2-3 groups
    if 2 <= n_groups <= 4:
        pairwise = {}
        for i, g1 in enumerate(group_names):
            for g2 in group_names[i+1:]:
                v1 = groups.get_group(g1).values
                v2 = groups.get_group(g2).values
                
                # Mann-Whitney U test
                u_stat, p_mw = stats.mannwhitneyu(v1, v2, alternative='two-sided')
                
                # Cliff's delta
                delta, delta_interp = cliffs_delta(v1, v2)
                
                pairwise[f"{g1}_vs_{g2}"] = {
                    "mann_whitney_U": u_stat,
                    "p_value": p_mw,
                    "cliffs_delta": delta,
                    "effect_size": delta_interp
                }
        
        results["pairwise_comparisons"] = pairwise
    
    # Visualization
    fig, axes = plt.subplots(1, 2, figsize=(12, 5))
    
    # Box plot (raw scale)
    ax1 = axes[0]
    merged.boxplot(column="depth", by=group_column, ax=ax1)
    ax1.set_title(f"{dataset_name}: Depth by {group_column}")
    ax1.set_xlabel(group_column)
    ax1.set_ylabel("Sequencing Depth (reads)")
    plt.suptitle("")  # Remove automatic title
    
    # Box plot (log scale)
    ax2 = axes[1]
    merged.boxplot(column="log10_depth", by=group_column, ax=ax2)
    ax2.set_title(f"{dataset_name}: log10(Depth) by {group_column}")
    ax2.set_xlabel(group_column)
    ax2.set_ylabel("log10(Sequencing Depth)")
    plt.suptitle("")
    
    plt.tight_layout()
    plot_path = os.path.join(output_dir, f"{dataset_name}_depth_boxplot.png")
    plt.savefig(plot_path, dpi=150, bbox_inches="tight")
    plt.close()
    
    results["plot_path"] = plot_path
    
    # Generate diagnosis report
    report = generate_diagnosis_report(results, dataset_name, group_column, output_dir)
    results["report"] = report
    
    return results


def generate_diagnosis_report(results: dict, dataset_name: str, group_column: str, output_dir: str) -> str:
    """
    Generate human-readable diagnosis report in Markdown format.
    """
    report_lines = [
        f"# Depth–Group Confounding Diagnosis Report",
        f"## Dataset: {dataset_name}",
        f"",
        f"**Grouping variable:** {group_column}",
        f"**Number of groups:** {results['n_groups']}",
        f"**Total samples:** {results['n_samples']}",
        f"",
        f"---",
        f"",
        f"## 1. Descriptive Statistics",
        f"",
        f"```",
        results["summary_stats"].to_string(),
        f"```",
        f"",
    ]
    
    # Kruskal-Wallis results
    if "kruskal_wallis" in results:
        kw = results["kruskal_wallis"]
        report_lines.extend([
            f"## 2. Kruskal-Wallis Test (Non-parametric ANOVA)",
            f"",
            f"- **H statistic:** {kw['H']:.2f}",
            f"- **p-value:** {kw['p_value']:.4e}",
            f"- **Effect size (η²):** {kw['eta_squared']:.4f} ({kw['eta_squared_interpretation']})",
            f"",
        ])
        
        # Interpretation
        if kw["p_value"] < 0.05:
            report_lines.append(f"**Interpretation:** Statistically significant depth difference between groups (p < 0.05).")
            if kw["eta_squared"] >= 0.06:
                report_lines.append(f"**⚠️ WARNING:** Effect size is medium or large (η² ≥ 0.06). "
                                   f"Depth–group confounding may affect diversity comparisons.")
            else:
                report_lines.append(f"**Note:** Effect size is small (η² < 0.06). "
                                   f"Statistical significance may be driven by large sample size; practical effect is limited.")
        else:
            report_lines.append(f"**Interpretation:** No statistically significant depth difference between groups (p ≥ 0.05).")
        report_lines.append(f"")
    
    # Pairwise comparisons
    if "pairwise_comparisons" in results:
        report_lines.extend([
            f"## 3. Pairwise Comparisons (Mann-Whitney U + Cliff's Delta)",
            f"",
            f"| Comparison | U statistic | p-value | Cliff's δ | Effect size |",
            f"|------------|-------------|---------|-----------|-------------|",
        ])
        for comp, vals in results["pairwise_comparisons"].items():
            report_lines.append(
                f"| {comp} | {vals['mann_whitney_U']:.0f} | {vals['p_value']:.4e} | "
                f"{vals['cliffs_delta']:.3f} | {vals['effect_size']} |"
            )
        report_lines.append(f"")
    
    # Summary judgment
    report_lines.extend([
        f"## 4. Summary Judgment for Benchmark",
        f"",
    ])
    
    # Determine confounding status
    kw = results.get("kruskal_wallis", {})
    p_val = kw.get("p_value", 1.0)
    eta2 = kw.get("eta_squared", 0.0)
    
    if p_val < 0.05 and eta2 >= 0.06:
        confounding_status = "POTENTIAL CONFOUNDING"
        confounding_note = ("This dataset shows statistically significant depth differences with medium/large effect size. "
                          "Rarefaction depth selection requires extra caution to avoid differential dropout.")
    elif p_val < 0.05 and eta2 < 0.06:
        confounding_status = "MINOR CONFOUNDING (statistically significant but small effect)"
        confounding_note = ("Depth differences are statistically significant but effect size is small. "
                          "Proceed with standard depth selection, but report group-stratified dropout.")
    else:
        confounding_status = "NO EVIDENCE OF CONFOUNDING"
        confounding_note = ("No statistically significant depth–group association detected. "
                          "Rarefaction depth can be selected based on plateau and overall dropout.")
    
    report_lines.extend([
        f"**Status:** {confounding_status}",
        f"",
        f"{confounding_note}",
        f"",
        f"---",
        f"",
        f"*Generated by depth_group_confounding_diagnosis.py*",
    ])
    
    report = "\n".join(report_lines)
    
    # Save report
    report_path = os.path.join(output_dir, f"{dataset_name}_diagnosis_report.md")
    with open(report_path, "w") as f:
        f.write(report)
    
    print(f"Report saved to: {report_path}")
    
    return report


# ==============================================================================
# Main Execution
# ==============================================================================

def main():
    parser = argparse.ArgumentParser(
        description="Diagnose depth–group confounding for rarefaction benchmark",
        formatter_class=argparse.RawDescriptionHelpFormatter
    )
    parser.add_argument("--table-qza", required=True, help="Path to table.qza")
    parser.add_argument("--metadata", required=True, help="Path to metadata.tsv")
    parser.add_argument("--group-column", required=True, help="Metadata column for grouping")
    parser.add_argument("--dataset-name", default="dataset", help="Name for output files")
    parser.add_argument("--output-dir", default="./confounding_diagnosis", help="Output directory")
    
    args = parser.parse_args()
    
    # Load metadata
    metadata = pd.read_csv(args.metadata, sep="\t", index_col=0)
    
    # Extract depth
    print(f"Extracting depth from {args.table_qza}...")
    depth = extract_depth_from_qza(args.table_qza)
    
    # Run diagnosis
    print(f"Running depth–group confounding diagnosis...")
    results = diagnose_depth_group_confounding(
        depth=depth,
        metadata=metadata,
        group_column=args.group_column,
        dataset_name=args.dataset_name,
        output_dir=args.output_dir
    )
    
    print(f"\n{'='*60}")
    print(f"Diagnosis complete for {args.dataset_name}")
    print(f"{'='*60}")
    print(results["report"])


if __name__ == "__main__":
    main()
