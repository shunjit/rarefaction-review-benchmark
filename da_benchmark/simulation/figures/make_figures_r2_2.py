#!/usr/bin/env python3
"""R2-2 figures from the dsim_main_v1 final summary (metrics_by_cell.tsv).

Four print-oriented figures (PNG 300 dpi + PDF), fixed condition colors
(validated 3-slot categorical palette; raw=blue, rarefied=orange,
boots=aqua), thin marks, recessive grid, ink-colored text/legends:

  fig_dsim_null_main        taxon-FPR (log) and FWER vs depth ratio rho,
                            3 methods, null main cells (150 reps)
  fig_dsim_spike_main       empirical FDR (+-1 MC SE) and power_itt vs rho,
                            3 methods, spike f=3 main cells (150 reps)
  fig_dsim_null_sens        taxon-FPR (log) one-factor sensitivity
                            (theta, D0, n) at rho=100, ANCOM-BC2 & MaAsLin3
  fig_dsim_spike_sens       FDR and power one-factor sensitivity
                            (f, DA fraction, D0) at rho=100, both methods

Zeros on the log FPR axis are drawn at a floor (open marker, "0" label);
the floor sits below the smallest positive value in the data.

Usage:
    python make_figures_r2_2.py <summary_dir> <out_dir>
"""

import sys
from pathlib import Path

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import pandas as pd

SUMMARY = Path(sys.argv[1]) if len(sys.argv) > 1 else Path(
    "/Volumes/PS3000/benchmark_data/revision_r1/simulation/runs/dsim_main_v1/summary")
OUT = Path(sys.argv[2]) if len(sys.argv) > 2 else Path(__file__).parent / "output"
OUT.mkdir(parents=True, exist_ok=True)

INK = "#0b0b0b"
INK2 = "#52514e"
GRID = "#d9d8d4"
COND_COLOR = {"raw": "#2a78d6", "rarefied": "#eb6834", "boots": "#1baf7a"}
COND_LABEL = {"raw": "raw", "rarefied": "rarefied", "boots": "boots (≥50%)"}
CONDS = ["raw", "rarefied", "boots"]
METHODS = [("ancombc2", "ANCOM-BC2"), ("maaslin3", "MaAsLin3 (abundance)"),
           ("wilcoxon", "Wilcoxon (TSS)")]
FLOOR = 2e-5          # log-axis stand-in for exact zeros

plt.rcParams.update({
    "font.size": 9, "axes.titlesize": 10, "axes.labelsize": 9.5,
    "xtick.labelsize": 8.5, "ytick.labelsize": 8.5,
    "text.color": INK, "axes.edgecolor": INK2, "axes.labelcolor": INK,
    "xtick.color": INK2, "ytick.color": INK2,
    "axes.spines.top": False, "axes.spines.right": False,
    "figure.dpi": 120, "savefig.dpi": 300,
})

d = pd.read_csv(SUMMARY / "metrics_by_cell.tsv", sep="\t")
d = d[d["method"] != "maaslin3_prev"]          # prevalence model: secondary


def get(cell, method, cond, col):
    r = d[(d.cell_id == cell) & (d.method == method) & (d.condition == cond)]
    if len(r) != 1:
        sys.exit(f"lookup failed: {cell}/{method}/{cond} -> {len(r)} rows")
    return float(r.iloc[0][col])


def style_ax(ax):
    ax.grid(axis="y", color=GRID, linewidth=0.7, zorder=0)
    ax.set_axisbelow(True)


def series(ax, xs, ys, cond, log=False):
    col = COND_COLOR[cond]
    plot_y = [max(y, FLOOR) if log else y for y in ys]
    ax.plot(xs, plot_y, color=col, linewidth=1.8, zorder=3)
    for x, y, py in zip(xs, ys, plot_y):
        if log and y == 0:
            ax.plot([x], [py], marker="o", ms=6.5, mfc="white", mec=col,
                    mew=1.6, zorder=4)
            ax.annotate("0", (x, py), textcoords="offset points",
                        xytext=(0, 7), ha="center", fontsize=7.5, color=INK2)
        else:
            ax.plot([x], [py], marker="o", ms=6.5, mfc=col, mec="white",
                    mew=1.1, zorder=4)


def nominal(ax, y=0.05):
    ax.axhline(y, color=INK2, linewidth=1.0, linestyle=(0, (4, 3)), zorder=2)


def legend_fig(fig, with_nominal=True, y=1.0):
    handles = [plt.Line2D([], [], color=COND_COLOR[c], marker="o", ms=6.5,
                          mec="white", mew=1.1, linewidth=1.8,
                          label=COND_LABEL[c]) for c in CONDS]
    if with_nominal:
        handles.append(plt.Line2D([], [], color=INK2, linewidth=1.0,
                                  linestyle=(0, (4, 3)), label="nominal 0.05"))
    fig.legend(handles=handles, loc="upper center", ncol=len(handles),
               frameon=False, bbox_to_anchor=(0.5, y), fontsize=9)


def save(fig, name):
    for ext in ("png", "pdf"):
        fig.savefig(OUT / f"{name}.{ext}", bbox_inches="tight",
                    facecolor="white")
    plt.close(fig)
    print(f"wrote {name}.png/.pdf")


# ---------------------------------------------------------------- fig 1: null
RHO_CELLS = {"null": ["null_r1", "null_r10", "null_r100"],
             "spike": ["spike_r1_f3", "spike_r10_f3", "spike_r100_f3"]}
RHO_X = [0, 1, 2]
RHO_LBL = ["1", "10", "100"]

fig, axes = plt.subplots(2, 3, figsize=(7.2, 4.8), sharex=True)
for j, (m, mlab) in enumerate(METHODS):
    ax = axes[0][j]
    for cond in CONDS:
        series(ax, RHO_X, [get(c, m, cond, "fpr_taxon")
                           for c in RHO_CELLS["null"]], cond, log=True)
    ax.set_yscale("log")
    ax.set_ylim(FLOOR * 0.6, 0.4)
    nominal(ax)
    ax.set_title(mlab)
    style_ax(ax)
    if j == 0:
        ax.set_ylabel("taxon-level FPR")

    ax = axes[1][j]
    for cond in CONDS:
        series(ax, RHO_X, [get(c, m, cond, "fwer")
                           for c in RHO_CELLS["null"]], cond)
    ax.set_ylim(-0.04, 1.06)
    nominal(ax)
    style_ax(ax)
    ax.set_xticks(RHO_X, RHO_LBL)
    ax.set_xlabel("depth ratio ρ")
    if j == 0:
        ax.set_ylabel("FWER  P(≥1 false positive)")
fig.subplots_adjust(top=0.86, hspace=0.14, wspace=0.28)
legend_fig(fig, y=0.985)
save(fig, "fig_dsim_null_main")

# --------------------------------------------------------------- fig 2: spike
fig, axes = plt.subplots(2, 3, figsize=(7.2, 4.8), sharex=True)
for j, (m, mlab) in enumerate(METHODS):
    ax = axes[0][j]
    for cond in CONDS:
        ys = [get(c, m, cond, "mean_fdp") for c in RHO_CELLS["spike"]]
        se = [get(c, m, cond, "fdr_mc_se") for c in RHO_CELLS["spike"]]
        ax.errorbar(RHO_X, ys, yerr=se, color=COND_COLOR[cond], linewidth=0,
                    elinewidth=1.1, capsize=2.5, capthick=1.1, zorder=2)
        series(ax, RHO_X, ys, cond)
    ax.set_ylim(-0.04, 1.06)
    nominal(ax)
    ax.set_title(mlab)
    style_ax(ax)
    if j == 0:
        ax.set_ylabel("empirical FDR (±1 MC SE)")

    ax = axes[1][j]
    for cond in CONDS:
        series(ax, RHO_X, [get(c, m, cond, "power_itt")
                           for c in RHO_CELLS["spike"]], cond)
    ax.set_ylim(-0.02, 0.42)
    style_ax(ax)
    ax.set_xticks(RHO_X, RHO_LBL)
    ax.set_xlabel("depth ratio ρ")
    if j == 0:
        ax.set_ylabel("power (intention-to-test)")
fig.subplots_adjust(top=0.86, hspace=0.14, wspace=0.28)
legend_fig(fig, y=0.985)
save(fig, "fig_dsim_spike_main")

# ------------------------------------------------- fig 3: null sensitivity
NULL_FACTORS = [
    ("θ (overdispersion)", ["null_r100_th15", "null_r100", "null_r100_th200"],
     ["15", "52.65*", "200"]),
    ("base depth D0", ["null_r100_d2500", "null_r100", "null_r100_d40000"],
     ["2,500", "10,000*", "40,000"]),
    ("group size n", ["null_r100_n10", "null_r100", "null_r100_n50"],
     ["10", "20*", "50"]),
]
SENS_METHODS = [("ancombc2", "ANCOM-BC2"), ("maaslin3", "MaAsLin3 (abundance)")]

fig, axes = plt.subplots(2, 3, figsize=(7.2, 4.6))
for i, (m, mlab) in enumerate(SENS_METHODS):
    for j, (flab, cells, ticks) in enumerate(NULL_FACTORS):
        ax = axes[i][j]
        for cond in CONDS:
            series(ax, [0, 1, 2],
                   [get(c, m, cond, "fpr_taxon") for c in cells], cond, log=True)
        ax.set_yscale("log")
        ax.set_ylim(FLOOR * 0.6, 0.6)
        nominal(ax)
        style_ax(ax)
        ax.set_xticks([0, 1, 2], ticks)
        ax.set_xlim(-0.35, 2.35)
        if i == 0:
            ax.set_title(flab, fontsize=9.5)
        if i == 1:
            ax.set_xlabel(flab, fontsize=9)
        if j == 0:
            ax.set_ylabel(f"{mlab}\ntaxon-level FPR", fontsize=9)
fig.subplots_adjust(top=0.84, bottom=0.17, hspace=0.32, wspace=0.28)
legend_fig(fig, y=0.985)
fig.text(0.5, 0.02, "one-factor sensitivity around the null ρ=100 cell "
         "(* = main-analysis value; 100 reps, B=10)",
         ha="center", fontsize=8, color=INK2)
save(fig, "fig_dsim_null_sens")

# ------------------------------------------------ fig 4: spike sensitivity
SPIKE_FACTORS = [
    ("effect size f", ["spike_r100_f15", "spike_r100_f3", "spike_r100_f5"],
     ["1.5", "3*", "5"]),
    ("DA fraction", ["spike_r100_f3_da5", "spike_r100_f3", "spike_r100_f3_da20"],
     ["5%", "10%*", "20%"]),
    ("base depth D0", ["spike_r100_f3_d2500", "spike_r100_f3", "spike_r100_f3_d40000"],
     ["2,500", "10,000*", "40,000"]),
]

fig, axes = plt.subplots(4, 3, figsize=(7.2, 8.6))
for i, (m, mlab) in enumerate(SENS_METHODS):
    for j, (flab, cells, ticks) in enumerate(SPIKE_FACTORS):
        ax = axes[2 * i][j]
        for cond in CONDS:
            ys = [get(c, m, cond, "mean_fdp") for c in cells]
            se = [get(c, m, cond, "fdr_mc_se") for c in cells]
            ax.errorbar([0, 1, 2], ys, yerr=se, color=COND_COLOR[cond],
                        linewidth=0, elinewidth=1.1, capsize=2.5,
                        capthick=1.1, zorder=2)
            series(ax, [0, 1, 2], ys, cond)
        ax.set_ylim(-0.04, 1.06)
        nominal(ax)
        style_ax(ax)
        ax.set_xticks([0, 1, 2], ticks)
        ax.set_xlim(-0.35, 2.35)
        if i == 0:
            ax.set_title(flab, fontsize=9.5)
        if j == 0:
            ax.set_ylabel(f"{mlab}\nempirical FDR", fontsize=9)

        ax = axes[2 * i + 1][j]
        for cond in CONDS:
            series(ax, [0, 1, 2],
                   [get(c, m, cond, "power_itt") for c in cells], cond)
        ax.set_ylim(-0.02, 0.42)
        style_ax(ax)
        ax.set_xticks([0, 1, 2], ticks)
        ax.set_xlim(-0.35, 2.35)
        if j == 0:
            ax.set_ylabel(f"{mlab}\npower (ITT)", fontsize=9)
fig.subplots_adjust(top=0.90, bottom=0.07, hspace=0.38, wspace=0.28)
legend_fig(fig, y=0.975)
fig.text(0.5, 0.015, "one-factor sensitivity around the spike ρ=100, f=3 cell "
         "(* = main-analysis value; 100 reps, B=10)",
         ha="center", fontsize=8, color=INK2)
save(fig, "fig_dsim_spike_sens")

print(f"\nall figures -> {OUT}")
