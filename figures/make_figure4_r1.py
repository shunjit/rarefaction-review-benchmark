#!/usr/bin/env python3
"""
Figure 4 R1 (base = figure4_da_rarefaction_problems_v6_2.py sha256 dffc50200f9f…)
R1 changes (2026-08-17, C-4d): bottom bar "Raw counts" -> "Unrarefied counts"
(terminology pass, CG-10); output fig4_r1.*

Figure 4 v6 — Why rarefaction is problematic for differential abundance
       and how bias-aware models address depth variation

ISME Communications 投稿用
Version: v6 (2026-02-28)
引用番号: 本文 v3.5.6.24 参考文献リスト準拠

変更履歴 (v5 → v6):
    - figsize: (14, 10) → (7, 5) に変更
      「作成サイズ ≈ 出版サイズ」方式への統一。
      Figure 1（make_figure_v9.py, figsize=7.5×8.0）および
      Figure 2（v2a_v2, figsize=7.5×3.0）と同一のアプローチ。
    - フォントサイズ: 全体を約半分に調整
      v5: 14–16pt（出版時50%縮小で ~7–8pt）
      v6: 7–8.5pt（出版時~7%縮小で ~6.5–7.9pt）
    - ライン幅: Figure 1 と同等の視覚的重みに調整
    - 座標系 (0–14 × 0–10) は変更なし
      → 全ボックス座標・レイアウト定数をそのまま維持
    - その他の変更なし（テキスト内容、カラーパレット、レイアウト構造は v5 と同一）

出版時サイズに関する注記:
    この図は figsize=(7, 5) インチで作成されています。
    ISME Communications の2カラム幅（約170mm ≈ 6.7インチ）に配置する場合、
    縮小率は約 6.7/7 ≈ 0.96（約4%縮小）です。

    フォントサイズの対応（スクリプト → 出版時）:
        7.0pt → ~6.7pt（本文ボックス、ワークフローバー）
        7.5pt → ~7.2pt（まとめボックス）
        8.5pt → ~8.2pt（ヘッダー）
    すべて ISME 最小基準（6pt）を十分に上回ります。

    参考: Figures 間のサイズ方式比較
        Figure 1: figsize=(7.5, 8.0), 出力幅≈185mm, 縮小率≈0.92
        Figure 2: figsize=(7.5, 3.0), 出力幅≈190mm, 縮小率≈0.89
        Figure 4: figsize=(7.0, 5.0), 出力幅≈178mm, 縮小率≈0.96

実行方法:
    python3 figure4_da_rarefaction_problems_v6.py

出力:
    figure4_da_rarefaction_problems.pdf  （投稿用）
    figure4_da_rarefaction_problems.png  （確認用）

依存ライブラリ:
    matplotlib >= 3.5
"""

import matplotlib.pyplot as plt
from matplotlib.patches import FancyBboxPatch


# =====================================================================
# ① グローバル設定
# =====================================================================

plt.rcParams.update({
    'font.family': 'sans-serif',                                          # ← 汎用ファミリー名を指定
    'font.sans-serif': ['Arial', 'Liberation Sans', 'Helvetica', 'DejaVu Sans'],  # ← 優先順位リスト
    'font.size': 7,                    # v5: 14 → v6: 7（出版時 ~6.7pt）
    'figure.dpi': 150,
    'savefig.dpi': 300,
    'savefig.bbox': 'tight',
    'pdf.fonttype': 42,
})

# --- カラーパレット（v5 から変更なし）---
COLORS = {
    # 問題カラム（vermillion 系）
    'problem_header':   '#D55E00',
    'problem_bg':       '#FDE8D0',
    'problem_border':   '#C45200',
    'problem_summary_bg':     '#FADCC4',
    'problem_summary_border': '#A04000',

    # 解決策カラム（blue 系）
    'solution_header':  '#0072B2',
    'solution_bg':      '#D6EAF8',
    'solution_border':  '#005A8E',
    'solution_summary_bg':     '#B8D8F0',
    'solution_summary_border': '#004A7A',

    # 接続要素
    'flow':             '#2C3E6B',
    'connector':        '#555555',
    'connector_label':  '#333333',

    # テキスト
    'text_dark':        '#1A1A1A',
    'text_white':       '#FFFFFF',
}

# --- レイアウト定数（v5 から変更なし）---
# 座標系は (0–14) × (0–10) を維持。
# figsize を (7, 5) に変更したことで、
# 1データ単位 = 0.5 インチ（v5 では 1.0 インチ）。
# ボックス座標はすべてそのまま使用可能。

LEFT_X = 3.5
RIGHT_X = 10.5
BOX_W = 5.2
BOX_H = 1.1
HEADER_W = 5.5
HEADER_H = 0.75

Y_HEADER = 9.2
Y_ROW1 = 7.85
Y_ROW2 = 6.50
Y_ROW3 = 5.15
Y_ROW4 = 3.80
Y_SUMMARY = 2.30
Y_WORKFLOW = 0.85

# --- ライン幅設定（v6 で新設）---
# v5 では figsize=(14,10) でライン幅 1.8pt だったが、
# figsize=(7,5) ではボックスが物理的に半分になるため、
# 同じ 1.8pt ではボックスに対して相対的に2倍太く見える。
# Figure 1（make_figure_v9.py）は figsize=(7.5,8.0) で BW=1.8 を使用しているが、
# Figure 4 はボックス内テキストが3行と密なため、やや細めに設定して
# テキストの可読性を優先する。
LW_BOX = 1.2          # 通常ボックス枠線（v5: 1.8）
LW_SUMMARY = 1.5      # まとめボックス枠線（v5: 2.0）
LW_ARROW = 1.2        # 通常矢印（v5: 1.8）
LW_CONNECTOR = 1.4    # 接続矢印（v5: 2.0）
LW_CENTER_ARROW = 2.0 # 中央 "Use instead" 矢印（v5: 3.0）


# =====================================================================
# ② ヘルパー関数
# =====================================================================

def draw_box(ax, x, y, width, height, text,
             facecolor='white', edgecolor='black',
             fontsize=7, fontweight='normal', text_color=None,
             linewidth=None):
    """
    角丸のボックスを描画し、中央にテキストを配置する。

    Parameters
    ----------
    ax : matplotlib.axes.Axes
    x, y : float — ボックスの中心座標（データ座標）
    width, height : float — ボックスのサイズ（データ座標）
    text : str — ボックス内テキスト（\\n で改行）
    facecolor, edgecolor : str — 背景色・枠線色
    fontsize : float — フォントサイズ（v6: 出版時ほぼそのまま）
    fontweight : str — 'normal' or 'bold'
    text_color : str or None — テキスト色
    linewidth : float or None — 枠線の太さ（None の場合 LW_BOX を使用）
    """
    if text_color is None:
        text_color = COLORS['text_dark']
    if linewidth is None:
        linewidth = LW_BOX

    box = FancyBboxPatch(
        (x - width / 2, y - height / 2),
        width, height,
        boxstyle="round,pad=0.02,rounding_size=0.03",
        facecolor=facecolor,
        edgecolor=edgecolor,
        linewidth=linewidth,
    )
    ax.add_patch(box)
    ax.text(
        x, y, text,
        ha='center', va='center',
        fontsize=fontsize, fontweight=fontweight,
        color=text_color,
        multialignment='center',
        linespacing=1.25,
    )


def draw_arrow(ax, start, end, color='gray', style='->', lw=None):
    """矢印を描画する。"""
    if lw is None:
        lw = LW_ARROW
    ax.annotate(
        '',
        xy=end, xytext=start,
        arrowprops=dict(arrowstyle=style, color=color, lw=lw),
    )


# =====================================================================
# ③–⑤ メイン描画
# =====================================================================

def create_figure():
    """Figure 4 を生成して保存する。"""

    # --- キャンバスの作成 ---
    # v6: figsize を (14, 10) → (7, 5) に変更。
    # 座標系 (0–14) × (0–10) は維持するため、
    # 1データ単位 = 0.5 インチ（v5 では 1.0 インチ）。
    fig, ax = plt.subplots(figsize=(7, 5))
    ax.set_xlim(0, 14)
    ax.set_ylim(0, 10)
    ax.axis('off')

    # =================================================================
    # ③ 左カラム: Rarefaction の問題点（vermillion 系）
    # =================================================================

    # --- ヘッダー ---
    draw_box(ax, LEFT_X, Y_HEADER, HEADER_W, HEADER_H,
             'Rarefaction for DA: Problems',
             facecolor=COLORS['problem_header'],
             edgecolor=COLORS['problem_header'],
             fontsize=8.5, fontweight='bold',       # v5: 16 → v6: 8.5
             text_color=COLORS['text_white'])

    # --- 問題1: Read discard ---
    draw_box(ax, LEFT_X, Y_ROW1, BOX_W, BOX_H,
             '1. Read discard\n'
             'High-depth samples lose information\n'
             '\u2192 Reduced statistical power',
             facecolor=COLORS['problem_bg'],
             edgecolor=COLORS['problem_border'])

    # --- 問題2: Feature loss ---
    draw_box(ax, LEFT_X, Y_ROW2, BOX_W, BOX_H,
             '2. Feature loss\n'
             'Low-abundance taxa become zero\n'
             '\u2192 Biased inference on rare taxa',
             facecolor=COLORS['problem_bg'],
             edgecolor=COLORS['problem_border'])

    # --- 問題3: Resampling variance ---
    draw_box(ax, LEFT_X, Y_ROW3, BOX_W, BOX_H,
             '3. Resampling variance\n'
             'Single rarefaction adds stochastic noise\n'
             '\u2192 Increased estimator variance',
             facecolor=COLORS['problem_bg'],
             edgecolor=COLORS['problem_border'])

    # --- 問題4: Sample dropout ---
    draw_box(ax, LEFT_X, Y_ROW4, BOX_W, BOX_H,
             '4. Sample dropout\n'
             'Low-depth samples excluded entirely\n'
             '\u2192 Systematic bias if depth ~ group',
             facecolor=COLORS['problem_bg'],
             edgecolor=COLORS['problem_border'])

    # --- まとめ（警告）---
    draw_box(ax, LEFT_X, Y_SUMMARY, BOX_W, BOX_H,
             'CAUTION: For DA, these effects\n'
             'directly reduce detection power\n'
             'and can introduce systematic bias [5]',
             facecolor=COLORS['problem_summary_bg'],
             edgecolor=COLORS['problem_summary_border'],
             fontsize=7.5,                          # まとめは本文より少し大きく
             fontweight='bold',
             linewidth=LW_SUMMARY)

    # =================================================================
    # ④ 右カラム: Bias-aware models の解決策（blue 系）
    # =================================================================

    # --- ヘッダー ---
    draw_box(ax, RIGHT_X, Y_HEADER, HEADER_W, HEADER_H,
             'Bias-aware Models: Solutions',
             facecolor=COLORS['solution_header'],
             edgecolor=COLORS['solution_header'],
             fontsize=8.5, fontweight='bold',       # v5: 16 → v6: 8.5
             text_color=COLORS['text_white'])

    # --- 解決策1: Use all reads ---
    draw_box(ax, RIGHT_X, Y_ROW1, BOX_W, BOX_H,
             '1. Use all reads\n'
             'No subsampling required\n'
             '\u2192 Preserves statistical power',
             facecolor=COLORS['solution_bg'],
             edgecolor=COLORS['solution_border'])

    # --- 解決策2: Estimate sampling fraction ---
    draw_box(ax, RIGHT_X, Y_ROW2, BOX_W, BOX_H,
             '2. Estimate sampling fraction\n'
             'ANCOM-BC/BC2 models depth as bias [18,19]\n'
             '\u2192 Proper compositional inference',
             facecolor=COLORS['solution_bg'],
             edgecolor=COLORS['solution_border'])

    # --- 解決策3: Covariate adjustment ---
    draw_box(ax, RIGHT_X, Y_ROW3, BOX_W, BOX_H,
             '3. Covariate adjustment\n'
             'Control for confounders (batch, etc.)\n'
             '\u2192 Reduced false positives',
             facecolor=COLORS['solution_bg'],
             edgecolor=COLORS['solution_border'])

    # --- 解決策4: Sensitivity analysis ---
    draw_box(ax, RIGHT_X, Y_ROW4, BOX_W, BOX_H,
             '4. Sensitivity analysis\n'
             'Compare with LinDA [20], ALDEx2 [8]\n'
             '\u2192 Method-robust conclusions',
             facecolor=COLORS['solution_bg'],
             edgecolor=COLORS['solution_border'])

    # --- まとめ（推奨）---
    draw_box(ax, RIGHT_X, Y_SUMMARY, BOX_W, BOX_H,
             'RECOMMENDED for DA:\n'
             'ANCOM-BC2 [19] + sensitivity analysis\n'
             'with LinDA [20] / ALDEx2 [8]',
             facecolor=COLORS['solution_summary_bg'],
             edgecolor=COLORS['solution_summary_border'],
             fontsize=7.5,                          # まとめは本文より少し大きく
             fontweight='bold',
             linewidth=LW_SUMMARY)

    # =================================================================
    # ⑤ 中央・下部の接続要素
    # =================================================================

    # --- 中央の "Use instead" 矢印 ---
    arrow_y = 5.82
    arrow_x_start = 6.2
    arrow_x_end = 7.8

    ax.annotate(
        '',
        xy=(arrow_x_end, arrow_y),
        xytext=(arrow_x_start, arrow_y),
        arrowprops=dict(
            # v6: 矢印ヘッドサイズを縮小（物理的に半分の図面に合わせる）
            arrowstyle='->,head_width=0.3,head_length=0.2',  # v5: 0.4, 0.3
            color=COLORS['connector'],
            lw=LW_CENTER_ARROW,
        ),
    )

    # ラベル（矢印の上に配置）
    ax.text(
        (arrow_x_start + arrow_x_end) / 2,
        arrow_y + 0.25,                   # v5: +0.35 → v6: +0.25（小フォントに対応）
        'Use instead',
        ha='center', va='bottom',
        fontsize=7, fontweight='bold',    # v5: 14 → v6: 7
        color=COLORS['connector_label'],
        fontstyle='italic',
    )

    # --- 下部ワークフローバー ---
    draw_box(ax, 7, Y_WORKFLOW, 12.5, 0.9,
             'Workflow: Unrarefied counts \u2192 ANCOM-BC2 (primary) '
             '\u2192 LinDA + ALDEx2 (sensitivity) '
             '\u2192 Report concordant findings',
             facecolor=COLORS['flow'],
             edgecolor=COLORS['flow'],
             fontsize=7, fontweight='bold',          # v5: 14 → v6: 7
             text_color=COLORS['text_white'],
             linewidth=LW_SUMMARY)

    # --- まとめボックスからワークフローバーへの接続矢印 ---
    workflow_top = Y_WORKFLOW + 0.45 + 0.12
    draw_arrow(ax, (LEFT_X, Y_SUMMARY - BOX_H / 2),
               (4.5, workflow_top),
               color=COLORS['connector'], lw=LW_CONNECTOR)
    draw_arrow(ax, (RIGHT_X, Y_SUMMARY - BOX_H / 2),
               (9.5, workflow_top),
               color=COLORS['connector'], lw=LW_CONNECTOR)

    # =================================================================
    # 保存
    # =================================================================
    plt.tight_layout()

    fig.savefig(
        './fig4_r1.pdf',
        bbox_inches='tight',
        facecolor='white',
    )
    print("  PDF saved: fig4_r1.pdf")

    fig.savefig(
        './fig4_r1.png',
        dpi=300,
        bbox_inches='tight',
        facecolor='white',
    )
    print("  PNG saved: fig4_r1.png")

    plt.close()


# =====================================================================
# エントリーポイント
# =====================================================================

if __name__ == '__main__':
    print("Generating Figure 4 (v6)...")
    print()
    print("  Design approach: creation size ≈ publication size")
    print("  figsize: (7, 5) inches")
    print("  Coordinate system: (0–14) × (0–10) [unchanged from v5]")
    print("  1 data unit = 0.5 inches")
    print()
    print("  Publication scaling estimate:")
    print("    ISME 2-column width ≈ 170mm (6.7in)")
    print("    Output width ≈ 178mm (7.0in)")
    print("    Scale factor ≈ 0.96 (4% reduction)")
    print()
    print("  Font sizes at publication:")
    for label, size in [('Content boxes', 7.0), ('Summary boxes', 7.5),
                         ('Headers', 8.5), ('Workflow bar', 7.0),
                         ('"Use instead"', 7.0)]:
        print(f"    {label}: {size}pt → ~{size * 6.7 / 7:.1f}pt at print")
    print()

    create_figure()
    print()
    print("Done.")


# =====================================================================
# カスタマイズ FAQ（v6 更新）
# =====================================================================
#
# Q: v5 との視覚的な違いは？
# A: 最終出版物としてはほぼ同一です。v5 は 14×10 インチで作成し
#    出版時に ~50% 縮小される前提でしたが、v6 は 7×5 インチで作成し
#    出版時に ~4% しか縮小されません。スクリプト内のフォントサイズが
#    出版時のサイズとほぼ一致するため、プレビュー時の見た目が
#    そのまま最終出版物に近くなります。
#
# Q: フォント優先順位を変更したい
# A: plt.rcParams['font.sans-serif'] のリスト順を変更してください。
#    font.family は 'sans-serif' のまま維持すること。
#
# Q: カラーパレットの根拠は？
# A: v5 から変更なし。Okabe-Ito カラーユニバーサルデザイン推奨パレットを
#    基調とし、vermillion (#D55E00) と blue (#0072B2) で対比を表現。
#
# Q: ライン幅の設計根拠は？
# A: v6 では figsize 縮小に伴い、ボックスの物理サイズが半減したため、
#    ライン幅もスケールダウン（1.8→1.2, 2.0→1.5 等）。
#    Figure 1（make_figure_v9.py）は BW=1.8 を使用しているが、
#    Figure 4 はボックス内テキストが3行と密なため、
#    やや細めに設定してテキスト可読性を優先した。
#
# Q: 引用番号のマッピング（本文 v3.5.6.24 準拠）
# A:   [5]  = McMurdie & Holmes 2014 (rarefaction critique)
#      [8]  = Fernandes 2014 (ALDEx2)
#      [18] = Lin & Peddada 2020 (ANCOM-BC)
#      [19] = Lin & Peddada 2024 (ANCOM-BC2)
#      [20] = Zhou 2022 (LinDA)
