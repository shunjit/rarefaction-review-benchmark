#!/usr/bin/env python3
"""R2-1 figures from the dreal_main_v1 summary (phase 1, B=0: raw/rarefied).

D-real is a depth-associated empirical stress test on Rice (design section 5):
per rep, the SAME 80% subsample is split by median depth within compartment
(real_depth) and, as the empirical null discovery background, the SAME sample
set gets its labels randomly permuted within compartment (real_perm).  There
is no ground truth, so everything here is descriptive distributions over the
R = 100 resampled splits (median / IQR / points), never SEs or tests.

Three print figures (PNG 300 dpi + PDF), condition colors as in R2-2
(validated palette; raw=blue, rarefied=orange), scenario encoded by fill
(depth split = filled, restricted permutation = open):

  fig_dreal_strat_main      rice_strat (stratified, main analysis)
                            row 1: detected/tested per scenario x condition
                            row 2: within-rep paired difference depth - perm
  fig_dreal_rhizo_robust    the same for rice_rhizo (single compartment)
  fig_dreal_achieved_ratio  achieved depth ratios per assignment, and the
                            paired difference vs achieved ratio (rice_strat)

No nominal-alpha line is drawn: the tested universe is not a known null
(depth may covary with real biology); the permutation arm is the reference.

Usage:
    python make_figures_r2_1.py <summary_dir> <out_dir>
"""

import sys
from pathlib import Path

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
from scipy.stats import spearmanr

SUMMARY = Path(sys.argv[1]) if len(sys.argv) > 1 else Path(
    "/Volumes/PS3000/benchmark_data/revision_r1/realdata_typeI/runs/"
    "dreal_main_v1/summary")
OUT = Path(sys.argv[2]) if len(sys.argv) > 2 else Path(__file__).parent / "output"
OUT.mkdir(parents=True, exist_ok=True)

INK = "#0b0b0b"
INK2 = "#52514e"
GRID = "#d9d8d4"
COND_COLOR = {"raw": "#2a78d6", "rarefied": "#eb6834"}
CONDS = ["raw", "rarefied"]
METHODS = [("ancombc2", "ANCOM-BC2"), ("maaslin3", "MaAsLin3 (abundance)"),
           ("wilcoxon", "Wilcoxon (TSS)")]
RNG = np.random.RandomState(20260814)          # deterministic jitter

plt.rcParams.update({
    "font.size": 9, "axes.titlesize": 10, "axes.labelsize": 9.5,
    "xtick.labelsize": 8.5, "ytick.labelsize": 8.5,
    "text.color": INK, "axes.edgecolor": INK2, "axes.labelcolor": INK,
    "xtick.color": INK2, "ytick.color": INK2,
    "axes.spines.top": False, "axes.spines.right": False,
    "figure.dpi": 120, "savefig.dpi": 300,
})

pr = pd.read_csv(SUMMARY / "per_rep_detections.tsv", sep="\t")
pr = pr[pr["method"] != "maaslin3_prev"]       # prevalence model: secondary

# within-rep paired difference (same subsample, same tables, labels differ)
dep = pr[pr.scenario == "real_depth"].set_index(
    ["cell_id", "method", "condition", "rep"])
per = pr[pr.scenario == "real_perm"].set_index(
    ["cell_id", "method", "condition", "rep"])
pair = dep[["n_detected", "achieved_ratio_depth", "achieved_ratio_perm"]].copy()
pair["diff"] = dep["n_detected"] - per["n_detected"]
pair = pair.dropna(subset=["diff"]).reset_index()


def style_ax(ax):
    ax.grid(axis="y", color=GRID, linewidth=0.7, zorder=0)
    ax.set_axisbelow(True)


def box(ax, x, vals, width=0.46):
    ax.boxplot([vals], positions=[x], widths=[width], showfliers=False,
               medianprops=dict(color=INK, linewidth=1.4),
               boxprops=dict(color=INK2, linewidth=0.9),
               whiskerprops=dict(color=INK2, linewidth=0.9),
               capprops=dict(color=INK2, linewidth=0.9), zorder=2)


def dots(ax, x, vals, color, filled, width=0.15):
    jx = x + RNG.uniform(-width, width, size=len(vals))
    if filled:
        ax.scatter(jx, vals, s=7, facecolor=color, edgecolor="none",
                   alpha=0.5, zorder=3)
    else:
        ax.scatter(jx, vals, s=8, facecolor="white", edgecolor=color,
                   linewidth=0.6, alpha=0.75, zorder=3)


def sel(cell, method, cond, scen):
    r = pr[(pr.cell_id == cell) & (pr.method == method)
           & (pr.condition == cond) & (pr.scenario == scen)]
    if len(r) == 0:
        sys.exit(f"lookup failed: {cell}/{method}/{cond}/{scen}")
    return r


def legend_fig(fig, y=1.0):
    handles = [
        plt.Line2D([], [], color=COND_COLOR["raw"], marker="s", ms=7,
                   linewidth=0, label="raw"),
        plt.Line2D([], [], color=COND_COLOR["rarefied"], marker="s", ms=7,
                   linewidth=0, label="rarefied"),
        plt.Line2D([], [], color=INK2, marker="o", ms=6.5, linewidth=0,
                   label="depth split"),
        plt.Line2D([], [], color=INK2, marker="o", ms=6.5, mfc="white",
                   mew=1.2, linewidth=0, label="restricted permutation"),
    ]
    fig.legend(handles=handles, loc="upper center", ncol=4, frameon=False,
               bbox_to_anchor=(0.5, y), fontsize=9)


def save(fig, name):
    for ext in ("png", "pdf"):
        fig.savefig(OUT / f"{name}.{ext}", bbox_inches="tight",
                    facecolor="white")
    plt.close(fig)
    print(f"wrote {name}.png/.pdf")


# ------------------------------------------- figs 1-2: per-cell distributions
X4 = [1.0, 2.0, 3.5, 4.5]                      # raw d/p, rarefied d/p
PAIR_MID = [(X4[0] + X4[1]) / 2, (X4[2] + X4[3]) / 2]


def cell_figure(cell, name):
    fig, axes = plt.subplots(2, 3, figsize=(7.2, 5.0),
                             height_ratios=[1.15, 1.0])
    n_med = int(pr[pr.cell_id == cell]["n_tested"].median())
    for j, (m, mlab) in enumerate(METHODS):
        ax = axes[0][j]
        for k, cond in enumerate(CONDS):
            for s, scen in enumerate(["real_depth", "real_perm"]):
                x = X4[2 * k + s]
                v = sel(cell, m, cond, scen)["frac_detected"].values
                box(ax, x, v)
                dots(ax, x, v, COND_COLOR[cond], filled=(scen == "real_depth"))
        ax.set_title(mlab)
        style_ax(ax)
        ax.set_xticks(X4, ["depth", "perm", "depth", "perm"], fontsize=8)
        ax.set_xlim(0.35, 5.15)
        for mid, cond in zip(PAIR_MID, CONDS):
            ax.text(mid, -0.24, cond, transform=ax.get_xaxis_transform(),
                    ha="center", fontsize=9, color=INK)
        ax.margins(y=0.06)
        ax.set_ylim(bottom=min(0, ax.get_ylim()[0]) - 0.015 * ax.get_ylim()[1])
        if j == 0:
            ax.set_ylabel(f"detected / tested  (q < 0.05)\n"
                          f"~{n_med:,} tested taxa")

        ax = axes[1][j]
        for k, cond in enumerate(CONDS):
            g = pair[(pair.cell_id == cell) & (pair.method == m)
                     & (pair.condition == cond)]["diff"].values
            box(ax, [1.0, 2.0][k], g, width=0.4)
            dots(ax, [1.0, 2.0][k], g, COND_COLOR[cond], filled=True,
                 width=0.13)
        ax.axhline(0, color=INK2, linewidth=1.0, linestyle=(0, (4, 3)),
                   zorder=1)
        style_ax(ax)
        ax.set_xticks([1.0, 2.0], ["raw", "rarefied"])
        ax.set_xlim(0.45, 2.55)
        ax.margins(y=0.08)
        if j == 0:
            ax.set_ylabel("Δ detections per rep\n(depth − permutation)")
    fig.subplots_adjust(top=0.87, hspace=0.42, wspace=0.30)
    legend_fig(fig, y=0.985)
    save(fig, name)


cell_figure("rice_strat", "fig_dreal_strat_main")
cell_figure("rice_rhizo", "fig_dreal_rhizo_robust")

# ------------------------------------- fig 3: achieved depth ratio and effect
fig, axes = plt.subplots(2, 2, figsize=(7.2, 5.6))

ax = axes[0][0]
ratios = pr[["cell_id", "rep", "achieved_ratio_depth",
             "achieved_ratio_perm"]].drop_duplicates(["cell_id", "rep"])
for i, cell in enumerate(["rice_strat", "rice_rhizo"]):
    g = ratios[ratios.cell_id == cell]
    for s, col in enumerate(["achieved_ratio_depth", "achieved_ratio_perm"]):
        x = X4[2 * i + s]
        box(ax, x, g[col].values)
        dots(ax, x, g[col].values, INK2, filled=(s == 0))
ax.axhline(1.0, color=INK2, linewidth=1.0, linestyle=(0, (4, 3)), zorder=1)
style_ax(ax)
ax.set_xticks(X4, ["depth", "perm", "depth", "perm"], fontsize=8)
ax.set_xlim(0.35, 5.15)
for mid, lab in zip(PAIR_MID, ["stratified", "Rhizoplane"]):
    ax.text(mid, -0.20, lab, transform=ax.get_xaxis_transform(),
            ha="center", fontsize=9, color=INK)
ax.set_ylabel("achieved depth ratio\n(median High / median Low)")
ax.set_title("achieved ratio per assignment", fontsize=9.5)

for p, (m, mlab) in enumerate(METHODS):
    ax = axes.flat[p + 1]
    txt = []
    for cond in CONDS:
        g = pair[(pair.cell_id == "rice_strat") & (pair.method == m)
                 & (pair.condition == cond)]
        ax.scatter(g["achieved_ratio_depth"], g["diff"], s=8,
                   facecolor=COND_COLOR[cond], edgecolor="none", alpha=0.55,
                   zorder=3)
        rho = spearmanr(g["achieved_ratio_depth"], g["diff"]).statistic
        txt.append(f"{cond} {rho:+.2f}")
    ax.axhline(0, color=INK2, linewidth=1.0, linestyle=(0, (4, 3)), zorder=1)
    style_ax(ax)
    ax.set_title(mlab, fontsize=9.5)
    ax.text(0.03, 0.97, "Spearman ρ: " + " · ".join(txt),
            transform=ax.transAxes, va="top", fontsize=7.5, color=INK2)
    ax.margins(y=0.10)
    if p >= 1:
        ax.set_xlabel("achieved depth ratio (depth split)")
    if p <= 1:
        ax.set_ylabel("Δ detections per rep\n(depth − permutation)")
fig.text(0.5, 0.004, "scatter panels: rice_strat; R = 100 resampled splits; "
         "descriptive only (subsamples overlap across reps)", ha="center",
         fontsize=8, color=INK2)
legend_fig(fig, y=0.995)
fig.subplots_adjust(top=0.90, bottom=0.13, hspace=0.38, wspace=0.32)
save(fig, "fig_dreal_achieved_ratio")

print(f"\nall figures -> {OUT}")
