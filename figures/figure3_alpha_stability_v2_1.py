#!/usr/bin/env python3
"""
Figure 3: Stability of alpha diversity estimates under repeated rarefaction
==========================================================================
v2: アクセシビリティ改訂版

【v1 → v2 変更点】
1. P1: 作成サイズ ≈ 出版サイズ方式に統一
   - figsize: (14, 4.5) → (7.5, 3.0) — Figure 2 v2a_v2 と同一
   - 縮小率: ~0.48 → ~0.89（ISME 2カラム幅 ~6.7in 基準）
   - 全フォントサイズを再設計（出版時 6pt 以上を保証）
   - 線幅を縮小後の視認性に合わせて調整

2. P3: Panel C にマーカー形状を追加（色以外の識別手段）
   - Soil: 丸 (o)、Bioethanol: 四角 (s)、Rice: 三角 (^)

3. Panel A 中央値アノテーション位置修正
   - v1: xy=(i, median), xytext=(12, 0) — ボックス内で点群と重なる
   - v2: 各データセットの最大値の上に配置し、重なりを完全回避

4. 微調整
   - Panel C マーカーサイズ: 20 → 25（縮小後の視認性向上）
   - stripplot ジッター幅: 0.15 → 0.12（縮小後のはみ出し防止）
   - boxplot 線幅: 1.2 → 1.0（縮小後のバランス調整）

【必要なデータファイル】
- alpha_cv_summary_soil.csv
- alpha_cv_summary_bioethanol.csv
- alpha_cv_summary_rice.csv

作成日: 2026-03-01
対象論文: ISME Communications Rarefaction Review
"""

import os
import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
import seaborn as sns

# ============================================================
# 基本設定（論文品質のMatplotlib設定）
# ============================================================
# v2: 「作成サイズ ≈ 出版サイズ」方式に合わせたフォントサイズ設計
# ISME 2カラム幅 ~6.7in に対して figsize=7.5in → 縮小率 ~0.89
# すべてのフォントサイズが出版時 6pt 以上を満たすよう設計
plt.rcParams.update({
    'font.family': 'sans-serif',                                          # 汎用ファミリー名を指定
    'font.sans-serif': ['Arial', 'Liberation Sans', 'Helvetica', 'DejaVu Sans'],  # 優先順位リスト
    'font.size': 8,                # 本文 → 出版時 ~7.1pt
    'axes.titlesize': 8,           # パネルタイトル → 出版時 ~8.0pt
    'axes.labelsize': 8,           # 軸ラベル → 出版時 ~7.1pt
    'xtick.labelsize': 7,          # X軸目盛り → 出版時 ~6.2pt
    'ytick.labelsize': 7,          # Y軸目盛り → 出版時 ~6.2pt
    'legend.fontsize': 7,          # 凡例 → 出版時 ~6.2pt
    'figure.dpi': 150,             # 画面表示用DPI
    'savefig.dpi': 300,            # 保存時DPI
    'savefig.bbox': 'tight',       # 余白を自動調整
    'pdf.fonttype': 42,            # PDFにフォントを埋め込む
})


# ============================================================
# 調整可能なパラメータ（ここを編集して微調整）
# ============================================================

# --- データディレクトリ ---
DATA_DIR = '/Volumes/PS3000/benchmark_data/data_for_figure3'

# --- 出力ディレクトリ ---
OUTPUT_DIR = '/Volumes/PS3000/benchmark_data/data_for_figure3/output'

# --- 図のサイズ（インチ）---
# v2: 出版サイズに近い値に変更
# ISME Communications 2カラム幅 ~6.7in に対して縮小率 ~0.89
FIGURE_WIDTH = 7.5     # v1: 14 → v2: 7.5
FIGURE_HEIGHT = 3.0    # v1: 4.5 → v2: 3.0

# --- CV閾値（参照線）---
CV_THRESHOLD = 1.0  # 1%

# --- カラーパレット（Okabe-Itoベース）---
COLORS = {
    'Soil': '#E69F00',       # オレンジ（Okabe-Ito公式名: orange）
    'Bioethanol': '#56B4E9', # スカイブルー（Okabe-Ito公式名: sky blue）
    'Rice': '#009E73',       # ティール（Okabe-Ito公式名: bluish green）
}

# --- P3: マーカー形状（Panel C 散布図用）---
# 色以外の識別手段を提供（ISME投稿規定 line 610 準拠）
MARKERS = {
    'Soil': 'o',        # 丸（circle）
    'Bioethanol': 's',  # 四角（square）
    'Rice': '^',        # 三角（triangle）
}

# --- データセット情報 ---
DATASET_INFO = {
    'Soil': {
        'depth': 20000,
        'n': 18,
    },
    'Bioethanol': {
        'depth': 35000,
        'n': 95,
    },
    'Rice': {
        'depth': 20000,
        'n': 488,
    },
}

# --- データセット表示順序 ---
DATASET_ORDER = ['Soil', 'Bioethanol', 'Rice']

# --- 線幅定数（v2: 出版サイズに合わせて調整）---
LW_BOX = 1.0          # v1: 1.2 → v2: 1.0（boxplotの枠線）
LW_REFLINE = 0.8      # v1: 1.0 → v2: 0.8（CV=1%参照線）

# --- Panel A アノテーション設定（v2: 重なり解消）---
# v1 の問題: xytext=(12, 0) でボックス内の中央値位置に配置
#   → stripplot の点群とテキストが重なり、数値が判読不能
# v2 の解決策: 各データセットの最大値の上にテキストを配置
#   → ボックス・ヒゲ・点群との重なりを完全に回避
ANNOT_Y_OFFSET_PT = 4   # 最大値からの上方オフセット（ポイント単位）
ANNOT_FONTSIZE = 7       # v1: 9 → v2: 7（出版時 ~6.2pt）


# ============================================================
# データ読み込み
# ============================================================

def load_cv_data():
    """
    各データセットのCV summaryファイルを読み込み、
    統合したDataFrameを返す
    """
    files = {
        'Soil': f'{DATA_DIR}/alpha_cv_summary_soil.csv',
        'Bioethanol': f'{DATA_DIR}/alpha_cv_summary_bioethanol.csv',
        'Rice': f'{DATA_DIR}/alpha_cv_summary_rice.csv',
    }
    
    all_data = []
    
    for dataset_name, filepath in files.items():
        if not os.path.exists(filepath):
            print(f"警告: {dataset_name} のファイルが見つかりません: {filepath}")
            continue
        
        df = pd.read_csv(filepath, index_col=0)
        df['dataset'] = dataset_name
        df['sample_id'] = df.index
        df = df.reset_index(drop=True)
        
        all_data.append(df)
        print(f"  {dataset_name}: {len(df)} samples loaded")
    
    if not all_data:
        raise FileNotFoundError("データファイルが見つかりません")
    
    return pd.concat(all_data, ignore_index=True)


# ============================================================
# プロット関数
# ============================================================

def plot_panel_a_boxplot(ax, data):
    """
    パネルA: 箱ひげ図 + ストリッププロット
    
    v2 変更点:
    - 中央値アノテーションを最大値の上に移動（重なり解消）
    - 線幅・フォントサイズを出版サイズに最適化
    - stripplot ジッター幅を縮小（はみ出し防止）
    """
    palette = [COLORS[d] for d in DATASET_ORDER]
    
    import warnings
    with warnings.catch_warnings():
        warnings.filterwarnings("ignore", category=FutureWarning)
        
        # 箱ひげ図
        sns.boxplot(
            data=data, 
            x='dataset', 
            y='cv_percent',
            order=DATASET_ORDER, 
            palette=palette,
            width=0.5, 
            linewidth=LW_BOX,     # v2: 1.2 → 1.0
            fliersize=0,
            ax=ax
        )
        
        # 個々のサンプルを点で表示
        sns.stripplot(
            data=data, 
            x='dataset', 
            y='cv_percent',
            order=DATASET_ORDER, 
            palette=palette,
            size=3,
            alpha=0.6,
            jitter=0.12,          # v2: 0.15 → 0.12（縮小後のはみ出し防止）
            ax=ax
        )
    
    # ----- v2: 中央値アノテーションの位置修正 -----
    # 設計方針:
    #   v1 では xy=(i, median_val), xytext=(12, 0) で配置していたが、
    #   テキストが boxplot 本体および stripplot 点群と重なっていた。
    #   v2 では各データセットの最大値の上方にテキストを配置することで、
    #   すべてのデータ要素との重なりを完全に回避する。
    #   テキストは水平中央揃えで boxplot の真上に表示される。
    for i, dataset in enumerate(DATASET_ORDER):
        subset = data[data['dataset'] == dataset]
        median_val = subset['cv_percent'].median()
        max_val = subset['cv_percent'].max()
        
        # 最大値の上方にテキストを配置
        ax.annotate(
            f'{median_val:.2f}%',
            xy=(i, max_val),               # アンカー: 最大値の位置
            xytext=(0, ANNOT_Y_OFFSET_PT), # 上方にオフセット
            textcoords='offset points',
            fontsize=ANNOT_FONTSIZE,        # v2: 9 → 7
            fontweight='bold',
            color=COLORS[dataset],
            ha='center',                    # v2: 水平中央揃え（v1は左揃え）
            va='bottom',                    # v2: 下端を基準（v1はva='center'）
        )
    
    # CV = 1% の参照線
    ax.axhline(
        y=CV_THRESHOLD, 
        color='gray', 
        linestyle='--', 
        linewidth=LW_REFLINE,   # v2: 1.0 → 0.8
        alpha=0.7
    )
    ax.text(
        2.55, CV_THRESHOLD, 
        f'CV = {CV_THRESHOLD:.0f}%',
        fontsize=7,             # v2: 8 → 7（出版時 ~6.2pt）
        color='gray', 
        va='center'
    )
    
    # 軸設定
    ax.set_xlabel('')
    ax.set_ylabel('Coefficient of Variation (%)')
    
    # パネルラベル（括弧なし、左上端に配置）
    ax.text(
        -0.02, 1.05,
        'A',
        transform=ax.transAxes,
        fontsize=10,            # v2: 12 → 10（出版時 ~8.9pt）
        fontweight='bold',
        va='bottom',
        ha='right',
    )
    ax.set_title(
        'Alpha diversity stability\n(Shannon index, 100 iterations)',
        fontweight='bold'
    )
    
    # Y軸範囲
    y_max = max(data['cv_percent'].max() * 1.15, CV_THRESHOLD * 1.2)
    ax.set_ylim(0, y_max)
    
    # X軸ラベルにn数を追加
    new_labels = [
        f'{d}\n(n={DATASET_INFO[d]["n"]})' 
        for d in DATASET_ORDER
    ]
    ax.set_xticks(range(len(DATASET_ORDER)))
    ax.set_xticklabels(new_labels)


def plot_panel_b_histogram(ax, data):
    """
    パネルB: CV分布のヒストグラム + KDE曲線
    
    v2 変更点:
    - 線幅・フォントサイズを出版サイズに最適化
    """
    for dataset in DATASET_ORDER:
        subset = data[data['dataset'] == dataset]
        
        sns.histplot(
            data=subset, 
            x='cv_percent',
            kde=True,
            color=COLORS[dataset],
            alpha=0.4,
            label=f'{dataset} (n={DATASET_INFO[dataset]["n"]})',
            bins=20,
            stat='density',
            line_kws={'linewidth': 1.2},  # v2追加: KDE曲線の線幅を明示
            ax=ax
        )
    
    # CV = 1% の参照線
    ax.axvline(
        x=CV_THRESHOLD, 
        color='gray', 
        linestyle='--', 
        linewidth=LW_REFLINE,   # v2: 1.0 → 0.8
        alpha=0.7
    )
    
    # 軸設定
    ax.set_xlabel('Coefficient of Variation (%)')
    ax.set_ylabel('Density')
    
    # パネルラベル
    ax.text(
        -0.02, 1.05,
        'B',
        transform=ax.transAxes,
        fontsize=10,            # v2: 12 → 10
        fontweight='bold',
        va='bottom',
        ha='right',
    )
    ax.set_title(
        'Distribution of Alpha CV\nacross samples',
        fontweight='bold'
    )
    
    ax.legend(loc='upper right', framealpha=0.9, fontsize=7)
    
    x_max = max(data['cv_percent'].max() * 1.1, CV_THRESHOLD * 1.2)
    ax.set_xlim(0, x_max)


def plot_panel_c_scatter(ax, data):
    """
    パネルC: Shannon多様性 vs CV の散布図
    
    v2 変更点:
    - P3: データセットごとに異なるマーカー形状を使用
      Soil: 丸(o)、Bioethanol: 四角(s)、Rice: 三角(^)
    - マーカーサイズ: 20 → 25（縮小後の視認性向上）
    """
    for dataset in DATASET_ORDER:
        subset = data[data['dataset'] == dataset]
        
        ax.scatter(
            subset['mean'],
            subset['cv_percent'],
            c=COLORS[dataset],
            alpha=0.6,
            s=25,                       # v2: 20 → 25（縮小後の視認性向上）
            label=dataset,
            marker=MARKERS[dataset],    # v2追加: P3 マーカー形状
            edgecolors='white',
            linewidths=0.3
        )
    
    # CV = 1% の参照線
    ax.axhline(
        y=CV_THRESHOLD, 
        color='gray', 
        linestyle='--', 
        linewidth=LW_REFLINE,   # v2: 1.0 → 0.8
        alpha=0.7
    )
    
    # 軸設定
    ax.set_xlabel('Mean Shannon index (100 iterations)')
    ax.set_ylabel('Coefficient of Variation (%)')
    
    # パネルラベル
    ax.text(
        -0.02, 1.05,
        'C',
        transform=ax.transAxes,
        fontsize=10,            # v2: 12 → 10
        fontweight='bold',
        va='bottom',
        ha='right',
    )
    ax.set_title(
        'Shannon diversity vs\nresampling variability',
        fontweight='bold'
    )
    
    ax.legend(loc='upper right', framealpha=0.9, fontsize=7)


# ============================================================
# メイン処理
# ============================================================

def main():
    """
    メイン処理：データ読み込み → Figure作成 → 保存
    """
    print("=" * 60)
    print("Figure 3: Stability of Alpha Diversity Estimates (v2)")
    print("=" * 60)
    print()
    
    # ----- v2 設計情報の表示 -----
    print("【v2 設計パラメータ】")
    print(f"  作成サイズ: {FIGURE_WIDTH} × {FIGURE_HEIGHT} inches")
    scale_factor = 6.7 / FIGURE_WIDTH
    print(f"  出版時縮小率: {scale_factor:.2f} (ISME 2カラム幅 ~6.7in)")
    print(f"  フォント設計:")
    for key, label in [
        ('axes.titlesize', 'パネルタイトル'),
        ('axes.labelsize', '軸ラベル'),
        ('xtick.labelsize', '目盛りラベル'),
        ('legend.fontsize', '凡例'),
    ]:
        script_pt = plt.rcParams[key]
        pub_pt = script_pt * scale_factor
        status = "✓" if pub_pt >= 6.0 else "✗"
        print(f"    {label}: {script_pt}pt → 出版時 ~{pub_pt:.1f}pt {status}")
    print(f"    パネルラベル: 10pt → 出版時 ~{10*scale_factor:.1f}pt ✓")
    print(f"    中央値アノテーション: {ANNOT_FONTSIZE}pt "
          f"→ 出版時 ~{ANNOT_FONTSIZE*scale_factor:.1f}pt "
          f"{'✓' if ANNOT_FONTSIZE*scale_factor >= 6.0 else '✗'}")
    print()
    
    # ----- データ読み込み -----
    print("データを読み込み中...")
    data = load_cv_data()
    print(f"合計サンプル数: {len(data)}")
    print()
    
    # ----- サマリー統計の表示 -----
    print("Summary Statistics (Alpha CV):")
    print("-" * 50)
    for dataset in DATASET_ORDER:
        subset = data[data['dataset'] == dataset]
        cv = subset['cv_percent']
        
        print(f"{dataset}:")
        print(f"  n = {len(subset)}, depth = {DATASET_INFO[dataset]['depth']:,}")
        print(f"  CV median = {cv.median():.3f}%")
        print(f"  CV mean   = {cv.mean():.3f}%")
        print(f"  CV range  = {cv.min():.3f}% - {cv.max():.3f}%")
        
        below_threshold = (cv < CV_THRESHOLD).sum() / len(cv) * 100
        print(f"  CV < {CV_THRESHOLD}%: {below_threshold:.1f}% of samples")
    print("-" * 50)
    print()
    
    # ----- Figure作成 -----
    print("Figure を作成中...")
    
    fig, axes = plt.subplots(1, 3, figsize=(FIGURE_WIDTH, FIGURE_HEIGHT))
    
    # パネルA: 箱ひげ図（中央値アノテーション位置修正済み）
    plot_panel_a_boxplot(axes[0], data)
    
    # パネルB: ヒストグラム
    plot_panel_b_histogram(axes[1], data)
    
    # パネルC: 散布図（P3: マーカー形状追加済み）
    plot_panel_c_scatter(axes[2], data)
    
    plt.tight_layout()
    
    # ----- 保存 -----
    os.makedirs(OUTPUT_DIR, exist_ok=True)
    
    # PDF（投稿用）
    pdf_path = f'{OUTPUT_DIR}/figure3_alpha_stability.pdf'
    fig.savefig(pdf_path, bbox_inches='tight', facecolor='white')
    print(f"PDF保存: {pdf_path}")
    
    # PNG（確認用）
    png_path = f'{OUTPUT_DIR}/figure3_alpha_stability.png'
    fig.savefig(png_path, dpi=300, bbox_inches='tight', facecolor='white')
    print(f"PNG保存: {png_path}")
    
    plt.close()
    
    print()
    print("=" * 60)
    print("Figure 3 v2 作成完了!")
    print("=" * 60)
    print()
    print("【v1 → v2 変更サマリー】")
    print("  P1: figsize (14, 4.5) → (7.5, 3.0) — 出版サイズ方式統一")
    print("  P1: 全フォントサイズ再設計 — 出版時 6pt 以上保証")
    print("  P3: Panel C マーカー形状追加 — ○/□/△")
    print("  Fix: Panel A 中央値アノテーション位置 — 最大値上方に移動")
    print()
    print("【Cross-Figure 作成サイズ統一状況】")
    print("  Figure 1: 7.5 × 8.0 in (scale ~0.92) ✓")
    print("  Figure 2: 7.5 × 3.0 in (scale ~0.89) ✓")
    print("  Figure 3: 7.5 × 3.0 in (scale ~0.89) ✓ ← 本改訂")
    print("  Figure 4: 7.0 × 5.0 in (scale ~0.96) ✓")


if __name__ == '__main__':
    main()
