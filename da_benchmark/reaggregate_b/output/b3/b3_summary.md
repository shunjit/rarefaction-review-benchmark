# B-3: closure シフトの定量 — 理論値と検出床の対比（CG-08 部分採用分）

- 実施日: 2026-08-16／入力 = `dsim_main_v1` スパイク 9 セル（1,050 rep）の per-rep 出力＋truth
- 理論: 潜在段階注入の総質量係数 c により、**非注入 taxa の相対存在量は正確に 1/c に一律縮小**する（truth の `rel_fc_analytic`＝1/c）。estimand はこれを absolute null と固定（設計書 §6）。B-3 はこの closure シフト（log c）が、各手法が実際に検出する効果量（検出床）と比べてどの規模かを実測する
- 尺度: ancombc2 = 自然対数（bias-corrected LFC）、maaslin3 = log2（LOG 変換の abundance 係数）。**wilcoxon は effect が平均差（非対数尺度）のため対比の対象外**（メモに明記の上除外）
- 判定規則は default（q<0.05）。boots は効果量が複製中央値であり対象外

## 理論値（セル別）

c は rep 間でほぼ一定（min–max は TSV 参照）。log2 シフト = log2(c)。

| セル | c 平均 | log2 シフト |
|---|---|---|
| spike_r100_f15 | 1.009 | 0.0126 |
| spike_r100_f3_da5 | 1.034 | 0.0466 |
| spike_r100_f3_d2500 | 1.056 | 0.0774 |
| spike_r10_f3 | 1.063 | 0.0854 |
| spike_r100_f3_d40000 | 1.068 | 0.0931 |
| **spike_r1_f3 / r10 / r100 (f3)** | **1.063–1.073** | **0.085–0.098** |
| spike_r100_f5 | 1.126 | 0.1676 |
| spike_r100_f3_da20 | 1.140 | 0.1842 |

## 検出床との対比（代表 = spike_r100_f3、全 36 行は TSV）

| 手法 | 条件 | closure シフト（同尺度） | FP \|effect\| p05 | 床/シフト比 | FP のうちシフト以下 |
|---|---|---|---|---|---|
| ancombc2 | raw | 0.0683 (ln) | 1.762 | **25.8×** | 0% |
| ancombc2 | rarefied | 0.0683 (ln) | 1.135 | **16.6×** | 0% |
| maaslin3 | rarefied | 0.0985 (log2) | 2.662 | **27.0×** | 0% |
| maaslin3 | raw | 0.0985 (log2) | 0.0234 | 0.24× | **23.1%** |

## 観察

1. **ancombc2（両条件）と maaslin3 rarefied では、FP の効果量はいずれも closure シフトの 17〜27 倍以上で、シフト以下の FP は 0%。** closure が偽陽性の説明になりえない、という estimand 弁護（CG-08 回答の中核）を直接支持する。
2. **maaslin3 raw のみ**、ρ≥10 のセルで FP の 12〜29% が closure 尺度以下の微小効果（p05 ≈ 0.014–0.033）。ただし **closure は ρ に依存しない**（c ≈ 1.06–1.07 で ρ=1/10/100 同一）のに、この微小効果 FP は **ρ=1 では消失**（spike_r1_f3: FP 9 件、床 33×、シフト以下 0%）。したがってこれは closure ではなく**深度交絡（unrarefied＋TSS+LOG の深度依存ゼロ構造）の産物**であり、むしろ「ρ=1 で FDR 挙動が正常＝closure 効果は検出床未満」という裁定時の傍証（90_反映メモ CG-08 欄）を precise にした形。回答文書では maaslin3 raw の但し書きとして扱う。
3. 床/シフト比が最小になる設計は da20（シフト最大 0.184）でも ancombc2 床は桁上（TSV 参照）。

## 成果物

- `b3_closure_theory_by_cell.tsv`（9 セル）
- `b3_detection_floor_vs_closure.tsv`（9 セル × 2 手法 × 2 条件 = 36 行）
