#!/usr/bin/env python3
"""R1-1 figures from the dsim_main_v1 final summary.

R1-1 (Reviewer #1): does rarefaction inflate the FDR of DA testing?
The spike-in scenario of D-sim measures empirical FDR (mean FDP +- MC SE)
and empirical power on Rice-calibrated synthetic data, so the raw vs
rarefied contrast answers the claim directly, per cell.

Two print figures (PNG 300 dpi + PDF), colors as in R2-2 (validated
palette; raw=blue, rarefied=orange, boots=aqua):

  fig_r11_fdr_direction   dumbbell raw -> rarefied per spike cell
                          (9 cells x 3 methods; FDR row + power row);
                          boots is omitted here (in the TSVs) to keep
                          the reviewer contrast readable
  fig_r11_stratum         abundance-stratum breakdown (low/mid/high)
                          for the main spike cell rho=100, f=3
                          (150 reps; FDR row + conditional power row)

Usage:
    python make_figures_r1_1.py <summary_dir> <out_dir>
"""

import sys
from pathlib import Path

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import pandas as pd

SUMMARY = Path(sys.argv[1]) if len(sys.argv) > 1 else Path(
    "/Volumes/PS3000/benchmark_data/revision_r1/simulation/runs/"
    "dsim_main_v1/summary")
OUT = Path(sys.argv[2]) if len(sys.argv) > 2 else Path(__file__).parent / "output"
OUT.mkdir(parents=True, exist_ok=True)

INK = "#0b0b0b"
INK2 = "#52514e"
GRID = "#d9d8d4"
COND_COLOR = {"raw": "#2a78d6", "rarefied": "#eb6834", "boots": "#1baf7a"}
COND_LABEL = {"raw": "raw", "rarefied": "rarefied", "boots": "boots (≥50%)"}
METHODS = [("ancombc2", "ANCOM-BC2"), ("maaslin3", "MaAsLin3 (abundance)"),
           ("wilcoxon", "Wilcoxon (TSS)")]

plt.rcParams.update({
    "font.size": 9, "axes.titlesize": 10, "axes.labelsize": 9.5,
    "xtick.labelsize": 8.5, "ytick.labelsize": 8.5,
    "text.color": INK, "axes.edgecolor": INK2, "axes.labelcolor": INK,
    "xtick.color": INK2, "ytick.color": INK2,
    "axes.spines.top": False, "axes.spines.right": False,
    "figure.dpi": 120, "savefig.dpi": 300,
})

d = pd.read_csv(SUMMARY / "metrics_by_cell.tsv", sep="\t")
d = d[d["method"] != "maaslin3_prev"]
ds = pd.read_csv(SUMMARY / "metrics_by_stratum.tsv", sep="\t")
ds = ds[ds["method"] != "maaslin3_prev"]


def get(frame, col, **kv):
    r = frame
    for k, v in kv.items():
        r = r[r[k] == v]
    if len(r) != 1:
        sys.exit(f"lookup failed: {kv} -> {len(r)} rows")
    return float(r.iloc[0][col])


def style_ax(ax, axis="x"):
    ax.grid(axis=axis, color=GRID, linewidth=0.7, zorder=0)
    ax.set_axisbelow(True)


def save(fig, name):
    for ext in ("png", "pdf"):
        fig.savefig(OUT / f"{name}.{ext}", bbox_inches="tight",
                    facecolor="white")
    plt.close(fig)
    print(f"wrote {name}.png/.pdf")


# ------------------------------------------------ fig 1: dumbbells per cell
CELLS = [
    ("spike_r1_f3", "ρ=1, f=3 (main)"),
    ("spike_r10_f3", "ρ=10, f=3 (main)"),
    ("spike_r100_f3", "ρ=100, f=3 (main)"),
    ("spike_r100_f15", "ρ=100, f=1.5"),
    ("spike_r100_f5", "ρ=100, f=5"),
    ("spike_r100_f3_da5", "ρ=100, f=3, DA 5%"),
    ("spike_r100_f3_da20", "ρ=100, f=3, DA 20%"),
    ("spike_r100_f3_d2500", "ρ=100, f=3, D0=2,500"),
    ("spike_r100_f3_d40000", "ρ=100, f=3, D0=40,000"),
]
YS = list(range(len(CELLS)))[::-1]              # main cells on top

fig, axes = plt.subplots(2, 3, figsize=(7.4, 6.2), sharey=True)
for j, (m, mlab) in enumerate(METHODS):
    for i, (col, se_col, xlim, xlab, nom) in enumerate([
            ("mean_fdp", "fdr_mc_se", (-0.03, 1.03),
             "empirical FDR (±1 MC SE)", 0.05),
            ("power_itt", None, (-0.02, 0.45),
             "power (intention-to-test)", None)]):
        ax = axes[i][j]
        for (cell, _), y in zip(CELLS, YS):
            xr = get(d, col, cell_id=cell, method=m, condition="raw")
            xf = get(d, col, cell_id=cell, method=m, condition="rarefied")
            ax.plot([xr, xf], [y, y], color=INK2, linewidth=1.0, zorder=2)
            for cond, x in (("raw", xr), ("rarefied", xf)):
                if se_col:
                    se = get(d, se_col, cell_id=cell, method=m,
                             condition=cond)
                    ax.errorbar([x], [y], xerr=[se], color=COND_COLOR[cond],
                                linewidth=0, elinewidth=1.0, capsize=2.0,
                                capthick=1.0, zorder=3)
                ax.plot([x], [y], marker="o", ms=6, mfc=COND_COLOR[cond],
                        mec="white", mew=1.0, linewidth=0, zorder=4)
        if nom is not None:
            ax.axvline(nom, color=INK2, linewidth=1.0,
                       linestyle=(0, (4, 3)), zorder=1)
        style_ax(ax)
        ax.set_xlim(*xlim)
        ax.set_ylim(-0.6, len(CELLS) - 0.4)
        ax.set_yticks(YS, [lab for _, lab in CELLS], fontsize=8)
        if i == 0:
            ax.set_title(mlab)
        if i == 1:
            ax.set_xlabel(xlab, fontsize=9)
        else:
            ax.set_xlabel(xlab, fontsize=9)
handles = [plt.Line2D([], [], color=COND_COLOR[c], marker="o", ms=6,
                      mec="white", mew=1.0, linewidth=0, label=COND_LABEL[c])
           for c in ("raw", "rarefied")]
handles.append(plt.Line2D([], [], color=INK2, linewidth=1.0,
                          linestyle=(0, (4, 3)), label="nominal 0.05"))
fig.legend(handles=handles, loc="upper center", ncol=3, frameon=False,
           bbox_to_anchor=(0.5, 1.0), fontsize=9)
fig.text(0.5, 0.005, "all 9 spike-in cells of dsim_main_v1 "
         "(main 150 reps, sensitivity 100 reps); boots in TSVs",
         ha="center", fontsize=8, color=INK2)
fig.subplots_adjust(top=0.90, bottom=0.09, hspace=0.30, wspace=0.10)
save(fig, "fig_r11_fdr_direction")

# ---------------------------------------- fig 2: abundance-stratum breakdown
STRATA = ["low", "mid", "high"]
SX = [0, 1, 2]
CELL = "spike_r100_f3"

fig, axes = plt.subplots(2, 3, figsize=(7.2, 4.8), sharex=True)
for j, (m, mlab) in enumerate(METHODS):
    for i, (col, ylim, ylab, nom) in enumerate([
            ("fdr", (-0.04, 1.06), "per-stratum FDR", 0.05),
            ("power_cot", (-0.02, 0.62),
             "per-stratum power\n(cond. on tested)", None)]):
        ax = axes[i][j]
        for cond in ("raw", "rarefied", "boots"):
            ys = [get(ds, col, cell_id=CELL, method=m, condition=cond,
                      stratum=s) for s in STRATA]
            ax.plot(SX, ys, color=COND_COLOR[cond], linewidth=1.8, zorder=3)
            ax.plot(SX, ys, marker="o", ms=6.5, mfc=COND_COLOR[cond],
                    mec="white", mew=1.1, linewidth=0, zorder=4)
        if nom is not None:
            ax.axhline(nom, color=INK2, linewidth=1.0,
                       linestyle=(0, (4, 3)), zorder=2)
        style_ax(ax, axis="y")
        ax.set_ylim(*ylim)
        if i == 0:
            ax.set_title(mlab)
        if i == 1:
            ax.set_xticks(SX, STRATA)
            ax.set_xlabel("baseline-abundance stratum")
        if j == 0:
            ax.set_ylabel(ylab)
handles = [plt.Line2D([], [], color=COND_COLOR[c], marker="o", ms=6.5,
                      mec="white", mew=1.1, linewidth=1.8,
                      label=COND_LABEL[c]) for c in ("raw", "rarefied",
                                                     "boots")]
handles.append(plt.Line2D([], [], color=INK2, linewidth=1.0,
                          linestyle=(0, (4, 3)), label="nominal 0.05"))
fig.legend(handles=handles, loc="upper center", ncol=4, frameon=False,
           bbox_to_anchor=(0.5, 1.0), fontsize=9)
fig.text(0.5, 0.005, "main spike cell ρ=100, f=3 (150 reps); spiked taxa "
         "drawn evenly from the three strata", ha="center", fontsize=8,
         color=INK2)
fig.subplots_adjust(top=0.85, bottom=0.13, hspace=0.24, wspace=0.28)
save(fig, "fig_r11_stratum")

print(f"\nall figures -> {OUT}")
