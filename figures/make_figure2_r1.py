#!/usr/bin/env python3
"""
Figure 2 R1 — Alpha-rarefaction curves (base = figure2_v2a_v2_1_linestyles.py sha256 1c1c93efd363…)
R1 changes (2026-08-17, C-4d): Panel C title n=488 -> n=493 (imported-sample accounting,
C-2/R3-8; curves and retention annotation unchanged); output fig2_r1.*; local output dir.
==========================================================================

v2a からの変更点:
  [P4] アノテーション・凡例の位置をパネルごとに個別制御する仕組みを導入。
       v2a では全パネル共通の固定比率（0.88, 0.75）で位置を決定していたが、
       パネル間でデータの分布域が大きく異なるため、以下の重なりが発生していた:
         - Panel A: "Selected depth" が Heated 曲線に重なる
         - Panel A: 凡例が Thaw 曲線に重なる
         - Panel C: "Selected depth" が Bulk Soil/Rhizosphere 曲線に重なる
         - Panel C: "Retention" が Rhizoplane 曲線に重なる
         - Panel C: 凡例が Endosphere 曲線に重なる
       各パネルのDATASETS辞書に位置パラメータを追加し、個別調整可能にした。

  優先度1–3の変更（v2aから継承）:
  [P1] FIGURE_WIDTH: 14 → 7.5 インチ
  [P2] Rhizoplane 黄色: #F0E442 → #C4A000
  [P3a] 群ごとに異なる線スタイル

出力:
    figure2_v2a_v2_linestyles.pdf（投稿用）
    figure2_v2a_v2_linestyles.png（確認用）
"""

import os
import re
import pandas as pd
import numpy as np
import matplotlib.pyplot as plt

# ============================================================
# 基本設定
# ============================================================
plt.rcParams.update({
    'font.family': 'sans-serif',                                          # 汎用ファミリー名を指定
    'font.sans-serif': ['Arial', 'Liberation Sans', 'Helvetica', 'DejaVu Sans'],  # 優先順位リスト
    'font.size': 8,
    'axes.titlesize': 9,
    'axes.labelsize': 8,
    'xtick.labelsize': 7,
    'ytick.labelsize': 7,
    'legend.fontsize': 7,
    'figure.dpi': 150,
    'savefig.dpi': 300,
    'savefig.bbox': 'tight',
    'pdf.fonttype': 42,
})


# ============================================================
# 調整可能なパラメータ
# ============================================================

DATA_DIR = '/Volumes/PS3000/benchmark_data/data_for_figure2'
OUTPUT_DIR = '/Volumes/PS3000/benchmark_data/data_for_figure2/output'

FIGURE_WIDTH = 7.5
FIGURE_HEIGHT = 3.0

SELECTED_DEPTHS = {
    'Soil': 20000,
    'Bioethanol': 35000,
    'Rice': 20000,
}

RETENTION_RATES = {
    'Soil': '100%',
    'Bioethanol': '100%',
    'Rice': '99%',
}

DISPLAY_N = {
    'Soil': 18,
    'Bioethanol': 95,
    'Rice': 493,
}

# --- カラーパレット ---
SOIL_COLORS = {
    'C': '#56B4E9',
    'T': '#009E73',
    'H': '#E69F00',
}

SOIL_LABELS = {
    'C': 'Control',
    'T': 'Thaw',
    'H': 'Heated',
}

SOIL_LINESTYLES = {
    'C': '-',
    'H': '--',
    'T': ':',
}

RICE_COLORS = {
    'Endosphere': '#D55E00',
    'Rhizoplane': '#C4A000',
    'Rhizosphere': '#CC79A7',
    'Bulk Soil': '#999999',
}

RICE_LINESTYLES = {
    'Bulk Soil': '-',
    'Endosphere': '--',
    'Rhizoplane': ':',
    'Rhizosphere': '-.',
}

BIOETHANOL_COLOR = '#0072B2'


# ============================================================
# [P4] パネルごとのアノテーション・凡例位置設定
# ============================================================
# 各値は y 軸範囲に対する比率（0.0=下端, 1.0=上端）。
# 位置調整が必要な場合はこのセクションの値のみ変更すればよい。
#
# --- 設計メモ ---
# Panel A (Soil): y軸範囲 ≈ 6.8–9.6
#   Heated ≈ 9.3–9.5 (比率 ~0.89–0.96)
#   Control ≈ 7.7–8.2 (比率 ~0.32–0.50)
#   Thaw ≈ 6.9–7.0 (比率 ~0.04–0.07)
#   → "Selected depth" を Control と Heated の間（~0.62）に配置
#   → "Retention" をその下（~0.48）に配置
#   → 凡例を Control と Thaw の間に引き上げ
#
# Panel B (Bioethanol): y軸範囲 ≈ 1.82–1.98
#   単一曲線で重なりの問題なし → デフォルト位置を維持
#
# Panel C (Rice): y軸範囲 ≈ 4.8–10.2
#   上位クラスタ (Bulk Soil ≈ 9.7, Rhizosphere ≈ 9.4): 比率 ~0.85–0.91
#   下位クラスタ (Endosphere ≈ 5.1, Rhizoplane ≈ 5.0): 比率 ~0.04–0.06
#   クラスタ間ギャップ: ~5.3–8.5 (比率 ~0.09–0.69)
#   → "Selected depth" をギャップ上部（~0.58）に配置
#   → "Retention" をギャップ中央（~0.42）に配置
#   → 凡例をギャップ内に引き上げ

ANNOTATION_POSITIONS = {
    'Soil': {
        'selected_depth_y': 0.8,     # v2a: 0.88 → Control–Heated間に下降
        'retention_y': 0.7,          # v2a: 0.75 → "Selected depth"の下に下降
        # 凡例位置: Axes座標 (x, y) で指定。
        # Thaw(~0.07)より上、Control(~0.32)より下に配置。
        'legend_loc': 'lower right',
        'legend_bbox': (0.98, 0.12),  # Thaw曲線の上方に引き上げ
    },
    'Bioethanol': {
        'selected_depth_y': 0.88,     # 変更なし
        'retention_y': 0.75,          # 変更なし
        'legend_loc': None,           # 凡例なし（"(overall mean)"テキスト使用）
        'legend_bbox': None,
    },
    'Rice': {
        'selected_depth_y': 0.58,     # v2a: 0.88 → クラスタ間ギャップ上部に下降
        'retention_y': 0.3,          # v2a: 0.75 → クラスタ間ギャップ中央に下降
        'legend_loc': 'center right',
        'legend_bbox': (0.98, 0.55),  # Endosphereの上方（ギャップ内）に引き上げ
    },
}


# ============================================================
# データセット設定
# ============================================================
DATASETS = {
    'Soil': {
        'file': f'{DATA_DIR}/shannon_soil.csv',
        'group_col': 'treatment',
        'selected_depth': SELECTED_DEPTHS['Soil'],
        'colors': SOIL_COLORS,
        'labels': SOIL_LABELS,
        'linestyles': SOIL_LINESTYLES,
        'n_display': DISPLAY_N['Soil'],
        'retention': RETENTION_RATES['Soil'],
    },
    'Bioethanol': {
        'file': f'{DATA_DIR}/shannon_bioethanol.csv',
        'group_col': 'time_point',
        'selected_depth': SELECTED_DEPTHS['Bioethanol'],
        'colors': None,
        'labels': None,
        'linestyles': None,
        'n_display': DISPLAY_N['Bioethanol'],
        'retention': RETENTION_RATES['Bioethanol'],
    },
    'Rice': {
        'file': f'{DATA_DIR}/shannon_rice.csv',
        'group_col': 'compartment',
        'selected_depth': SELECTED_DEPTHS['Rice'],
        'colors': RICE_COLORS,
        'labels': None,
        'linestyles': RICE_LINESTYLES,
        'n_display': DISPLAY_N['Rice'],
        'retention': RETENTION_RATES['Rice'],
    },
}


# ============================================================
# データ処理関数
# ============================================================

def parse_rarefaction_csv(filepath):
    """QIIME2 alpha-rarefaction CSVを解析し、整形済みDataFrameを返す。"""
    df = pd.read_csv(filepath)

    depth_cols = [c for c in df.columns if c.startswith('depth-')]
    depths = sorted(set(
        int(re.search(r'depth-(\d+)_', c).group(1))
        for c in depth_cols
    ))

    meta_cols = [
        c for c in df.columns
        if not c.startswith('depth-') and c != 'sample-id'
    ]

    records = []
    for _, row in df.iterrows():
        sample_id = row['sample-id']
        meta = {mc: row[mc] for mc in meta_cols if mc in df.columns}

        for depth in depths:
            iter_cols = [c for c in depth_cols if f'depth-{depth}_' in c]
            values = row[iter_cols].values.astype(float)
            valid_values = values[~np.isnan(values)]

            if len(valid_values) > 0:
                records.append({
                    'sample_id': sample_id,
                    'depth': depth,
                    'shannon_mean': np.mean(valid_values),
                    'shannon_std': np.std(valid_values),
                    'n_iter': len(valid_values),
                    **meta
                })

    return pd.DataFrame(records)


# ============================================================
# プロット関数
# ============================================================

def plot_rarefaction_panel(ax, df, config, panel_name):
    """
    1パネル分のrarefaction曲線をプロット。

    Parameters
    ----------
    ax : matplotlib.axes.Axes
    df : pd.DataFrame
    config : dict
        データセット設定（色、ラベル、線スタイルなど）
    panel_name : str
        データセット名（'Soil', 'Bioethanol', 'Rice'）。
        ANNOTATION_POSITIONS の参照キーとして使用。
    """
    group_col = config['group_col']
    selected_depth = config['selected_depth']
    colors = config['colors']
    labels = config['labels']
    linestyles = config['linestyles']
    retention = config['retention']

    # [P4] パネル固有のアノテーション位置を取得
    ann_pos = ANNOTATION_POSITIONS[panel_name]

    # ===== 曲線の描画 =====
    if colors and group_col in df.columns:
        groups = sorted(df[group_col].dropna().unique())

        for group in groups:
            group_df = df[df[group_col] == group]

            summary = group_df.groupby('depth').agg(
                mean=('shannon_mean', 'mean'),
                std=('shannon_mean', 'std'),
                n=('shannon_mean', 'count')
            ).reset_index()
            summary['se'] = summary['std'] / np.sqrt(summary['n'])

            color = colors.get(group, '#666666')
            label = labels.get(group, group) if labels else str(group)
            ls = linestyles.get(group, '-') if linestyles else '-'

            ax.plot(
                summary['depth'],
                summary['mean'],
                color=color,
                linewidth=1.5,
                label=label,
                alpha=0.9,
                linestyle=ls,
            )

            ax.fill_between(
                summary['depth'],
                summary['mean'] - summary['se'],
                summary['mean'] + summary['se'],
                color=color,
                alpha=0.2
            )
    else:
        # Bioethanol: 全サンプル平均
        summary = df.groupby('depth').agg(
            mean=('shannon_mean', 'mean'),
            std=('shannon_mean', 'std'),
            n=('shannon_mean', 'count')
        ).reset_index()
        summary['se'] = summary['std'] / np.sqrt(summary['n'])

        ax.plot(
            summary['depth'],
            summary['mean'],
            color=BIOETHANOL_COLOR,
            linewidth=1.5,
            label='Mean ± SE'
        )
        ax.fill_between(
            summary['depth'],
            summary['mean'] - summary['se'],
            summary['mean'] + summary['se'],
            color=BIOETHANOL_COLOR,
            alpha=0.3
        )

    # ===== 選択深度の垂直線 =====
    ax.axvline(
        x=selected_depth,
        color='red',
        linestyle='--',
        linewidth=1.2,
        alpha=0.7,
        zorder=10
    )

    # ===== Selected depth アノテーション =====
    # [P4] パネル固有の y 位置比率を使用
    y_lo, y_hi = ax.get_ylim()
    y_span = y_hi - y_lo
    y_pos = y_lo + y_span * ann_pos['selected_depth_y']

    ax.annotate(
        'Selected\ndepth',
        xy=(selected_depth, y_pos),
        xytext=(selected_depth * 0.55, y_pos),  # x方向をやや左に調整（0.60→0.55）
        fontsize=7,
        color='red',
        ha='center',
        va='center',
        arrowprops=dict(arrowstyle='->', color='red', lw=1)
    )

    # ===== Retention率の表示 =====
    # [P4] パネル固有の y 位置比率を使用
    y_ret = y_lo + y_span * ann_pos['retention_y']
    ax.text(
        selected_depth * 1.03,
        y_ret,
        f'Retention\n{retention}',
        fontsize=7,              # 6.5pt → 7pt に引き上げ（出版時 ~6.2pt、6pt基準超過）
        color='red',
        va='center',
        ha='left',
        style='italic',
        bbox=dict(
            boxstyle='round,pad=0.2',
            facecolor='white',
            edgecolor='red',
            alpha=0.85,
            linewidth=0.5
        )
    )

    # ===== 軸の設定 =====
    ax.set_xlabel('Sequencing depth')
    ax.set_ylabel('Shannon index')

    ax.xaxis.set_major_formatter(
        plt.FuncFormatter(lambda x, p: f'{x/1000:.0f}k')
    )


# ============================================================
# メイン処理
# ============================================================

def main():
    """メイン処理：データ読み込み → Figure作成 → 保存"""
    print("=" * 60)
    print("Figure 2 v2a_v2: Alpha-rarefaction Curves")
    print("  (Line Styles + Per-panel Position Tuning)")
    print("=" * 60)
    print()
    print(f"Figure size: {FIGURE_WIDTH} x {FIGURE_HEIGHT} inches")
    print(f"  → Publication width ~170mm (6.7in), scale ≈ {6.7/FIGURE_WIDTH:.2f}")
    print()

    # ----- データ読み込み -----
    data_dict = {}
    for name, config in DATASETS.items():
        filepath = config['file']

        if not os.path.exists(filepath):
            print(f"警告: {name} のデータファイルが見つかりません")
            print(f"  パス: {filepath}")
            data_dict[name] = None
            continue

        print(f"{name}: 読み込み中...")
        df = parse_rarefaction_csv(filepath)
        data_dict[name] = df

        n_samples = df['sample_id'].nunique()
        depth_min = df['depth'].min()
        depth_max = df['depth'].max()
        print(f"  サンプル数: {n_samples}")
        print(f"  深度範囲: {depth_min:,} - {depth_max:,}")

        if config['group_col'] in df.columns:
            groups = df[config['group_col']].unique()
            print(f"  群: {list(groups)}")
        print()

    # ----- Figure作成 -----
    print("Figure を作成中...")

    fig, axes = plt.subplots(1, 3, figsize=(FIGURE_WIDTH, FIGURE_HEIGHT))

    panel_labels = ['A', 'B', 'C']
    panel_names = list(DATASETS.keys())  # ['Soil', 'Bioethanol', 'Rice']

    for idx, (name, config) in enumerate(DATASETS.items()):
        ax = axes[idx]
        df = data_dict.get(name)

        if df is not None and len(df) > 0:
            # [P4] panel_name を渡してパネル固有の位置制御を有効化
            plot_rarefaction_panel(ax, df, config, panel_name=name)

            # パネルラベル（ISME規定準拠: 括弧なし大文字）
            ax.text(
                -0.02, 1.05,
                panel_labels[idx],
                transform=ax.transAxes,
                fontsize=10,
                fontweight='bold',
                va='bottom',
                ha='right',
            )

            # パネルタイトル
            n_display = config['n_display']
            ax.set_title(
                f'{name} (n={n_display})',
                fontweight='bold',
                fontsize=9
            )

            # [P4] 凡例: パネル固有の位置設定を適用
            ann_pos = ANNOTATION_POSITIONS[name]

            if config['colors'] and len(config['colors']) <= 6:
                legend_kwargs = {
                    'framealpha': 0.9,
                    'fontsize': 7,
                }
                if ann_pos.get('legend_bbox'):
                    # bbox_to_anchor で精密な位置指定
                    legend_kwargs['loc'] = ann_pos['legend_loc']
                    legend_kwargs['bbox_to_anchor'] = ann_pos['legend_bbox']
                else:
                    legend_kwargs['loc'] = ann_pos.get('legend_loc', 'lower right')

                ax.legend(**legend_kwargs)

            elif name == 'Bioethanol':
                ax.text(
                    0.98, 0.02,
                    '(overall mean)',
                    transform=ax.transAxes,
                    fontsize=6.5,
                    ha='right',
                    va='bottom',
                    style='italic',
                    color='#666666'
                )
        else:
            ax.text(
                0.5, 0.5,
                f'{name}\n(データなし)',
                ha='center',
                va='center',
                transform=ax.transAxes,
                fontsize=10
            )
            ax.set_title(f'{name}', fontweight='bold')
            ax.text(
                -0.02, 1.05,
                panel_labels[idx],
                transform=ax.transAxes,
                fontsize=10,
                fontweight='bold',
                va='bottom',
                ha='right',
            )

    plt.tight_layout()

    # ----- 保存 -----
    os.makedirs(OUTPUT_DIR, exist_ok=True)

    pdf_path = './fig2_r1.pdf'
    fig.savefig(pdf_path, bbox_inches='tight', facecolor='white')
    print(f"PDF保存: {pdf_path}")

    png_path = './fig2_r1.png'
    fig.savefig(png_path, dpi=300, bbox_inches='tight', facecolor='white')
    print(f"PNG保存: {png_path}")

    plt.close()

    print()
    print("=" * 60)
    print("Figure 2 v2a_v2 作成完了!")
    print("=" * 60)
    print()
    print("【v2a_v2 の変更点（v2a からの差分）】")
    print("  [P4] パネルごとのアノテーション位置調整:")
    print()
    for pname, pos in ANNOTATION_POSITIONS.items():
        print(f"  {pname}:")
        print(f"    Selected depth y = {pos['selected_depth_y']}"
              f"  (v2a: 0.88)")
        print(f"    Retention y      = {pos['retention_y']}"
              f"  (v2a: 0.75)")
        if pos.get('legend_bbox'):
            print(f"    Legend: loc='{pos['legend_loc']}', "
                  f"bbox_to_anchor={pos['legend_bbox']}")
        else:
            print(f"    Legend: (none / overall mean text)")
        print()
    print("【位置の微調整方法】")
    print("  ANNOTATION_POSITIONS 辞書の数値を変更してください。")
    print("  各値は y 軸範囲に対する比率（0.0=下端, 1.0=上端）です。")
    print("  例: Soil の 'selected_depth_y' を 0.62 → 0.55 に変更すると、")
    print("      アノテーションがさらに下に移動します。")


if __name__ == '__main__':
    main()
