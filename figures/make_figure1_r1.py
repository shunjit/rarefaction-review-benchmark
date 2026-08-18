#!/usr/bin/env python3
"""
Figure 1 R1 — revision update of v9_2 (submitted version)

R1 changes (2026-08-17, C-4d; base = make_figure_v9_2.py sha256 2775d57e5a9c…,
whose output text-matches the submitted fig1.pdf):
  1. Confounding status box: rule_v2_beta020 4-category
     (NONE p>=0.20 | INDETERMINATE 0.05<=p<0.20 | MINOR | POTENTIAL)
  2. Decision box: 'Confounding POTENTIAL?' -> 'Confounding classification?'
     with branch labels 'POTENTIAL (required) / INDETERMINATE (recommended)'
     (right, sensitivity) and 'NONE / MINOR' (left, proceed)
  3. Sensitivity-methods box: MaAsLin2 [24] -> MaAsLin3 [34]
  4. Output names fig1_r1.*

Key changes from v8 (inherited notes):
  1. Font: Arial / Liberation Sans / Helvetica (sans-serif fallback)
  2. Figure size: ~180 mm output width ≈ ISME full-width figure (~170 mm),
     so font sizes in script ≈ font sizes at publication (minimal scaling).
  3. Font sizes: 7.5–9.5 pt → ≥ 6–8 pt after minor journal scaling.
  4. Colorblind-safe palette: warm/cool distinction for adjacent boxes;
     border colors adjusted for deuteranopia/protanopia contrast.
  5. Line/border widths increased for visibility after print scaling.
  6. Title & legend removed per journal policy (lines 595, Instructions to Authors).
  LEFT  → Differential Abundance (DA)
  RIGHT → Diversity / Ordination
"""

import matplotlib.pyplot as plt
from matplotlib.patches import FancyBboxPatch
import os

# ============================================================
# Global style — Arial / Liberation Sans / Helvetica (sans-serif fallback)
# ============================================================
plt.rcParams.update({
    'font.family': 'sans-serif',                                          # 汎用ファミリー名を指定
    'font.sans-serif': ['Arial', 'Liberation Sans', 'Helvetica', 'DejaVu Sans'],  # 優先順位リスト
    'font.size': 8,
    'figure.dpi': 150,
    'savefig.dpi': 300,
    'savefig.bbox': 'tight',
})

# ============================================================
# Helper functions
# ============================================================
def draw_box(ax, x, y, width, height, text,
             facecolor='white', edgecolor='black',
             fontsize=8, fontweight='normal', text_color='black',
             linewidth=1.8, alpha=1.0):
    """Draw a rounded box with centred multi-line text."""
    box = FancyBboxPatch(
        (x - width / 2, y - height / 2), width, height,
        boxstyle="round,pad=0.02,rounding_size=0.02",
        facecolor=facecolor, edgecolor=edgecolor,
        linewidth=linewidth, alpha=alpha,
    )
    ax.add_patch(box)
    ax.text(x, y, text,
            ha='center', va='center',
            fontsize=fontsize, fontweight=fontweight,
            color=text_color, multialignment='center')


def draw_arrow(ax, start, end, color='#555555', text='',
               text_pos=0.5, fontsize=7.5, linestyle='solid'):
    """Draw an arrow with an optional mid-path label."""
    ax.annotate(
        '', xy=end, xytext=start,
        arrowprops=dict(arrowstyle='->', color=color,
                        lw=1.4, linestyle=linestyle),
    )
    if text:
        mx = start[0] + (end[0] - start[0]) * text_pos
        my = start[1] + (end[1] - start[1]) * text_pos
        ax.text(mx, my, text,
                ha='center', va='center', fontsize=fontsize,
                color=color, fontweight='bold',
                bbox=dict(boxstyle='round,pad=0.2',
                          facecolor='white', edgecolor='none', alpha=0.85))


# ============================================================
# Figure setup — sized close to ISME full-width column
#   Data coords 0–14 (x) × −0.9–13.4 (y) map onto ~7.5 × 8.0 in.
#   After bbox_inches='tight', output ≈ 185 mm wide.
#   At 170 mm print width → scale ≈ 0.92;
#     8 pt in script → ~7.4 pt at print (≥ 6 pt minimum).
# ============================================================
fig, ax = plt.subplots(figsize=(7.5, 8.0))
ax.set_xlim(0, 14)
ax.set_ylim(-0.9, 13.4)
ax.axis('off')

# ============================================================
# Colour palette — colorblind-accessible
#   • Fill colours: pastel tints (high contrast with black text)
#   • Edge colours: saturated, differ in both hue AND luminance
#     so that deuteranopia/protanopia readers can still
#     distinguish adjacent boxes (e.g., Yes/No pair).
#   • Warm/cool axis reinforces the red ↔ green distinction.
# ============================================================
C_START    = '#2C5F9E'   # dark blue fill  (white text)
C_STEP0    = '#D9E2F3'   # pale blue fill
C_DECISION = '#FFF2CC'   # pale amber fill
C_DIVERSITY = '#D4EDDA'  # pale green fill
C_DA       = '#FCE4D6'   # pale peach fill
C_ACTION   = '#DEEBF7'   # pale sky-blue fill
C_WARNING  = '#FDDDD5'   # warm peach fill   (shifted from pink)
C_SUCCESS  = '#C6EFCE'   # cool mint fill
C_Q2BOOTS  = '#E6E6FA'   # lavender fill
C_REPORT   = '#F0F0F0'   # pale grey fill

# Edge/border semantic colours
E_BLUE     = '#2166AC'   # diagnostics & action steps
E_AMBER    = '#B8860B'   # decision points (dark goldenrod)
E_GREEN    = '#1A7431'   # diversity path & proceed boxes
E_ORANGE   = '#D55E00'   # DA path header (Okabe-Ito orange)
E_RED      = '#B91C1C'   # caution / warning
E_PURPLE   = '#5B4A8A'   # q2-boots (optional)
E_DARK     = '#333333'   # report box

# Default border width (thicker for print visibility)
BW = 1.8

# ============================================================
# === START BOX ===
# ============================================================
draw_box(ax, 7, 12.8, 5.5, 0.6,
         'Microbiome data with uneven sequencing depth',
         facecolor=C_START, edgecolor=C_START,
         fontsize=9, fontweight='bold', text_color='white', linewidth=BW)

# === STEP 0: Library size diagnostics ===
draw_arrow(ax, (7, 12.5), (7, 12.1))
draw_box(ax, 7, 11.6, 6.5, 0.85,
         'Step 0: Library size diagnostics\n'
         'Summarize depth distribution; test depth–group\n'
         'association (Kruskal–Wallis + η²) [2,6,26]',
         facecolor=C_STEP0, edgecolor=E_BLUE,
         fontsize=8, fontweight='bold', linewidth=BW)

# Confounding status output (R1: rule_v2_beta020, 4 categories on 3 lines)
draw_arrow(ax, (7, 11.15), (7, 10.83))
draw_box(ax, 7, 10.35, 7.8, 0.95,
         'Confounding status:\n'
         'NONE (p ≥ 0.20)  |  INDETERMINATE (0.05 ≤ p < 0.20)\n'
         'MINOR (p < 0.05, η² < 0.06)  |  POTENTIAL (p < 0.05, η² ≥ 0.06)',
         facecolor=C_STEP0, edgecolor=E_BLUE,
         fontsize=7.5, linewidth=1.2)

# === STEP 1: Task bifurcation ===
draw_arrow(ax, (7, 9.87), (7, 9.55))
draw_box(ax, 7, 9.15, 4.2, 0.65,
         'Step 1: What is your\nanalytic task?',
         facecolor=C_DECISION, edgecolor=E_AMBER,
         fontsize=9, fontweight='bold', linewidth=BW)

# === BIFURCATION ARROWS ===
draw_arrow(ax, (4.9, 9.15), (3.5, 8.4), text='DA', text_pos=0.4)
draw_arrow(ax, (9.1, 9.15), (10.5, 8.4), text='Diversity', text_pos=0.4)

# ================================================================
# LEFT BRANCH — Differential Abundance
# ================================================================
draw_box(ax, 3.5, 8.05, 4.2, 0.65,
         'Differential abundance\n(feature-wise inference)',
         facecolor=C_DA, edgecolor=E_ORANGE,
         fontsize=8, fontweight='bold', linewidth=BW)

draw_arrow(ax, (3.5, 7.72), (3.5, 7.3))
draw_box(ax, 3.5, 6.8, 4.4, 0.9,
         'Rarefaction discouraged\nDiscards reads; increases variance;\n'
         'can introduce depth-dependent bias',
         facecolor=C_WARNING, edgecolor=E_RED,
         fontsize=8, fontweight='bold', linewidth=BW)

draw_arrow(ax, (3.5, 6.35), (3.5, 5.95))
draw_box(ax, 3.5, 5.4, 4.4, 0.9,
         'Use bias-aware models:\nANCOM-BC2 (primary) [19]\n'
         '+ covariate adjustment',
         facecolor=C_SUCCESS, edgecolor=E_GREEN,
         fontsize=8, fontweight='bold', linewidth=BW)

draw_arrow(ax, (3.5, 4.95), (3.5, 4.55))
draw_box(ax, 3.5, 4.0, 4.4, 0.9,
         'Sensitivity analysis:\nLinDA [20] / MaAsLin3 [34]\nALDEx2 [8]',
         facecolor=C_ACTION, edgecolor=E_BLUE,
         fontsize=8, linewidth=BW)

draw_arrow(ax, (3.5, 3.55), (3.5, 3.15))
draw_box(ax, 3.5, 2.6, 4.4, 0.9,
         'Report concordant findings\nacross methods\n'
         '(reduce method dependence) [21]',
         facecolor=C_SUCCESS, edgecolor=E_GREEN,
         fontsize=8, fontweight='bold', linewidth=BW)

# ================================================================
# RIGHT BRANCH — Diversity / Ordination
# ================================================================
draw_box(ax, 10.5, 8.05, 4.2, 0.65,
         'Diversity / Ordination\n(α/β diversity, distance-based tests)',
         facecolor=C_DIVERSITY, edgecolor=E_GREEN,
         fontsize=8, fontweight='bold', linewidth=BW)

draw_arrow(ax, (10.5, 7.72), (10.5, 7.35))
draw_box(ax, 10.5, 6.9, 4.2, 0.85,
         'Step 2: Select sampling depth\nRun alpha-rarefaction [1];\n'
         'check plateau + dropout (Fig. 2)',
         facecolor=C_ACTION, edgecolor=E_BLUE,
         fontsize=8, linewidth=BW)

draw_arrow(ax, (10.5, 6.45), (10.5, 6.1))
draw_box(ax, 10.5, 5.65, 4.2, 0.85,
         'Step 3: Quantify stability\nq2-boots / beta-rarefaction [1,10];\n'
         'report CV at chosen depth (Fig. 3)',
         facecolor=C_ACTION, edgecolor=E_BLUE,
         fontsize=8, linewidth=BW)

draw_arrow(ax, (10.5, 5.2), (10.5, 4.85))
draw_box(ax, 10.5, 4.4, 4.2, 0.85,
         'Step 4: Hypothesis testing\nPERMANOVA [28] + PERMDISP [29];\n'
         'report metric, permutations, covariates',
         facecolor=C_ACTION, edgecolor=E_BLUE,
         fontsize=8, linewidth=BW)

draw_arrow(ax, (10.5, 3.95), (10.5, 3.7))
draw_box(ax, 10.5, 3.3, 3.6, 0.7,
         'Confounding\nclassification?',
         facecolor=C_DECISION, edgecolor=E_AMBER,
         fontsize=8.5, fontweight='bold', linewidth=BW)

# ================================================================
# Dashed line: Step 0 result → Confounding POTENTIAL?
#   Routed along the right margin to avoid overlapping panels.
# ================================================================
_dash_x_right = 13.3
_dash_y_top   = 10.35
_dash_y_bot   = 3.3
_dash_x_entry = 12.3     # right edge of classification? box

ax.plot([10.9, _dash_x_right], [_dash_y_top, _dash_y_top],
        color=E_BLUE, lw=1.4, linestyle='dashed', clip_on=False)
ax.plot([_dash_x_right, _dash_x_right], [_dash_y_top, _dash_y_bot],
        color=E_BLUE, lw=1.4, linestyle='dashed', clip_on=False)
ax.annotate('', xy=(_dash_x_entry, _dash_y_bot),
            xytext=(_dash_x_right, _dash_y_bot),
            arrowprops=dict(arrowstyle='->', color=E_BLUE,
                            lw=1.4, linestyle='dashed'))

ax.text(_dash_x_right + 0.05, 6.8, 'Step 0\nresult',
        ha='left', va='center', fontsize=7.5, color=E_BLUE,
        fontstyle='italic',
        bbox=dict(boxstyle='round,pad=0.3', facecolor='white',
                  edgecolor=E_BLUE, alpha=0.9, linewidth=0.8))

# --- Classification branches (R1: labels replace Yes/No) ---
draw_arrow(ax, (11.6, 2.95), (11.95, 2.52))
ax.text(13.9, 2.88, 'POTENTIAL /\nINDETERMINATE',
        ha='right', va='top', fontsize=7, color='#555555', fontweight='bold',
        multialignment='right')
draw_box(ax, 12.5, 2.05, 2.8, 0.85,
         'Run sensitivity\nanalysis at\nmultiple depths',
         facecolor=C_WARNING, edgecolor=E_RED,
         fontsize=8, linewidth=BW)

draw_arrow(ax, (8.7, 3.05), (8.45, 2.52))
ax.text(8.35, 2.88, 'NONE / MINOR',
        ha='right', va='top', fontsize=7, color='#555555', fontweight='bold')
draw_box(ax, 8.5, 2.05, 2.8, 0.85,
         'Proceed with\nselected depth\n(report dropout)',
         facecolor=C_SUCCESS, edgecolor=E_GREEN,
         fontsize=8, linewidth=BW)

# --- q2-boots ---
draw_arrow(ax, (12.5, 1.6), (11.2, 1.25))
draw_arrow(ax, (8.5, 1.6), (9.8, 1.25))
draw_box(ax, 10.5, 0.9, 4.6, 0.65,
         'Optional: q2-boots [10] (repeated rarefaction\n'
         'for uncertainty quantification)',
         facecolor=C_Q2BOOTS, edgecolor=E_PURPLE,
         fontsize=7.5, linewidth=BW)

# ================================================================
# STEP 5 — Converging report step
# ================================================================
draw_arrow(ax, (3.5, 2.1), (5.5, 0.2))
draw_arrow(ax, (10.5, 0.55), (8.5, 0.2))
draw_box(ax, 7, -0.3, 7.2, 0.75,
         'Step 5: Report with QIIME 2 provenance [1]\n'
         '— depth, dropout, stability, method versions (Table 1)',
         facecolor=C_REPORT, edgecolor=E_DARK,
         fontsize=8, fontweight='bold', linewidth=BW)

# ============================================================
# Save
# ============================================================
plt.tight_layout()

output_png = './fig1_r1.png'
output_pdf = './fig1_r1.pdf'

fig.savefig(output_png, dpi=300, bbox_inches='tight', facecolor='white')
fig.savefig(output_pdf, bbox_inches='tight', facecolor='white')
plt.close()

print(f"✅ Figure 1 R1: {os.getcwd()}")
print(f"   {output_png}")
print(f"   {output_pdf}")
