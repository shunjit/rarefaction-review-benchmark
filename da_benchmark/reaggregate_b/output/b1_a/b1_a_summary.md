# B-1 拡張（b1_a）: diff_robust 並記の A 系適用 — n 感度＋極端 quantile（CG-05, GM-01×CG-05）

- 実施日: 2026-08-17（A 系完了処理と同朝）／**採用決定: 2026-08-17（所）— 8/20
  チェックポイントの表現強度判断の正式材料とし、同夜の freeze イベントで凍結する**
- 規則・範囲は B-1 と同一: **default** = `detected`（q_group < 0.05）／**diff_robust** =
  `detected & passed_ss`。method = ancombc2、condition ∈ {raw, rarefied}
  （boots は pseudo_sens=FALSE のため設計上対象外）。再フィットなし（凍結ローカルミラーの純再集計）
- 入力: `dsim_a2_nsens_v1`（200 rep、2 セル = spike_r100_f3_n10/_n50）＋
  `dreal_a1_x25_v1`（100 rep、1 セル = rice_strat_x25 × 2 割当）
- 自己検査: default 規則の再計算値は両 run の凍結 summary と全指標一致
  （dsim 最大差 7.1e-15、dreal 0 = 機械精度）
- スクリプト: `reaggregate_b/R/b1_diff_robust.R` を**無改修**で再利用（B-1 と同一・凍結済み）。
  出力 TSV 5 本＋本メモ

## D-sim: n 感度セル（B-1 主要セル表の n=10/50 拡張）

FWER = P(偽陽性 ≥ 1)、FDR = mean FDP、power = intention-to-test。
survival = diff_robust 検出数 / default 検出数（総数比）。
比較用の n=20 行は B-1（`output/b1/b1_summary.md`、dsim_main_v1 spike_r100_f3）。

| セル | 条件 | FWER 既定→robust | FDR 既定→robust | power_itt 既定→robust | survival |
|---|---|---|---|---|---|
| spike_r100_f3_n10 | raw | 0.72 → 0.34 | 0.615 → 0.327 | 0.035 → 0.002 | 0.126 |
| spike_r100_f3（n20、B-1） | raw | 0.980 → 0.460 | 0.813 → 0.447 | 0.084 → 0.003 | 0.066 |
| spike_r100_f3_n50 | raw | 1.000 → 0.39 | 0.789 → 0.339 | 0.197 → 0.011 | 0.040 |
| spike_r100_f3_n10 | rarefied | 0.87 → 0.74 | 0.748 → 0.653 | 0.059 → 0.023 | 0.411 |
| spike_r100_f3（n20、B-1） | rarefied | 1.000 → 1.000 | 0.842 → 0.841 | 0.240 → 0.124 | 0.522 |
| spike_r100_f3_n50 | rarefied | 1.000 → 1.000 | 0.777 → 0.792 | 0.436 → 0.228 | 0.569 |

## D-real: rice_strat_x25（達成深度比 2.53×。記述統計・検出数の中央値、tested 2,766）

| 条件 | 割当 | det 中央値 既定→robust |
|---|---|---|
| raw | depth | 1,760.5 → 303.5 |
| raw | perm | 1,173.5 → 286.5 |
| rarefied | depth | 1,497.5 → 259.5 |
| rarefied | perm | 1,249.5 → 75 |

ペア内差（depth − perm、中央値 [IQR]、深度側が多いペアの割合。比較列 = B-1 の
rice_strat main、達成深度比 1.62×）:

| 条件 | 規則 | x25（2.53×） | main（1.62×、B-1） |
|---|---|---|---|
| raw | default | +591.5 [+517, +641.75]、100% | +453.5、100% |
| raw | diff_robust | **+17.5 [−53, +75]、57%（中立化）** | **−107 [−152.25, −68.75]、2%（反転）** |
| rarefied | default | +256.5 [+219.25, +342]、100% | +211.5、100% |
| rarefied | diff_robust | **+171.5 [+134.75, +196.25]、96%（明確に残存）** | +5.5 [−17.25, +20]、55%（消失） |

## 観察（8/20 チェックポイントの判断材料 — B-1 観察 1〜4 への追補）

1. **B-1 観察 1（rarefied の inflation は diff_robust で消えない）は n に頑健。**
   n50 で FDR 0.777→0.792・FWER 1→1 と不変（n10 のみ 0.87→0.74 と微減 — 検出総数が
   小さく survival 0.41 の間引きが FWER に届く、という規模効果の範疇）。
2. **B-1 観察 2（raw の超過は sensitivity-fragile、robust 化で検出ごと消える）も n に頑健。**
   survival は n10 0.126 / n20 0.066 / n50 0.040 と n とともにさらに低下、
   power は 0.002〜0.011 でほぼ消失。「diff_robust は raw の救済ではない」が全 n で成立。
3. **D-real は増強分割で内訳が変わる（B-1 観察 4 の重要な修正）**: main の
   「diff_robust で raw 反転・rarefied 消失」のうち、**rarefied 側の消失は 1.62× に
   固有で一般化しない**。x25（2.53×）では rarefied の diff_robust 超過が
   +171.5・96% と明確に残存し、しかも置換背景が 75 まで下がるため**背景比では約 2.3 倍**
   （default の約 0.2 倍より相対的にむしろ大きい）。raw 側は反転（2%）→中立（57%）へ。
4. **表現強度への含意（確定は 8/20）**: 「既定判定では大規模超過、公式 sensitivity filter
   併記では raw は大幅減衰・rarefied の inflation は非感受」という B-1 の中間形要約は
   維持できるが、D-real の filter 後挙動は「raw 反転・rarefied 消失」と**一括では書けない**
   （分割強度依存）。回答文書・MEMO では「filter 後の残余は分割の深度信号強度に依存し、
   増強分割では rarefied 側にも depth-associated 超過が残る」と記述する（因果語は使わない）。

## 成果物（一次正本 = 本ディレクトリ。OneDrive 確定コピーは b1a_ プレフィクス）

- `b1_dsim_ancombc2_rules_by_cell.tsv` / `b1_dsim_ancombc2_rules_wide.tsv`
  （→ OneDrive R1-1/b1a_dsim_*。a2_* と同居）
- `b1_dreal_ancombc2_rules_by_cell.tsv` / `_paired.tsv` / `_per_rep.tsv`
  （→ OneDrive R2-1/b1a_dreal_*。a1_* と同居）
- `b1_provenance.txt`（→ R1-1/b1a_provenance.txt）／本メモ（→ 両ディレクトリへ
  b1a_summary.md として複製）
- **凍結**: 8/20 チェックポイントの freeze イベントで PS3000 6＋OneDrive 8 を登録予定。
  R2-1/R1-1 README への b1a_ 追記と R2-1/R2-2 MEMO の diff_robust 文言改版は
  8/20 の決定文言とセットで行い、同 freeze で新版登録する
