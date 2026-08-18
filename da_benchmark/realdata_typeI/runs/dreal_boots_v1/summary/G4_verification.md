# G4 再現性検証 — dreal_boots_v1 (第2相) vs dreal_main_v1 (第1相)

実施: 2026-08-15 朝（第2相完走・正本化直後）。
スクリプト: `realdata_typeI/harness/verify_g4_dreal.sh`（正本）。

## 目的

第2相 `dreal_boots_v1`（B=10 × R=25 × 2 セル、BASE_SEED=20260813）は、第1相
`dreal_main_v1`（B=0 × R=100、同一 BASE_SEED）の rep 1–25 と同一のサブサンプル・
同一の rarefied table を決定論再現し、raw/rarefied 条件を再計算する設計
（`D_HARNESS_DESIGN.md` v2 §14）。この再現が実際に成立していることを確認し、
第2相の boots 行を第1相の raw/rarefied 分布と rep 単位でペアリングする操作を
正当化する（G4）。

## 方法

全 50 (cell × rep) ペア（rice_strat / rice_rhizo × rep 1–25）で以下 3 点の
md5 一致を照合:

1. **結果 TSV の raw/rarefied 行**: `results/<cell>/rep_NNNN.tsv.gz` から
   condition ∈ {raw, rarefied} の行を抽出（ヘッダ保持・書き出し順の影響を
   除くため sort 後）— effect/se/p/q/detected/reason/passed_ss を含む全列
2. **truth ファイル**: `rep_NNNN_truth.tsv.gz` 全体（同一サブサンプル ⇒
   同一 tested set の確認）
3. **geninfo**: `rep_NNNN_geninfo.tsv` の B 列（第 25 列。設計上の唯一の差:
   10 vs 0）を除く全 27 列 — サブサンプル構成・置換割当・達成深度比・
   rarefied depth・seed 系列

## 結果

**PASS 50 / FAIL 0 — 全ペアでバイト同一**。

- 3 手法（ancombc2 / maaslin3(+maaslin3_prev) / wilcoxon）× 2 シナリオ
  （real_depth / real_perm）× raw/rarefied の fit 出力（数値列含む）が
  第1相と完全一致
- geninfo の一致により、サブサンプル抽出・High/Low 割当・restricted
  permutation・rarefied table（raref_depth・dropout 0）・seed 系列
  （dataset_seed）まで同一であることを確認
- 唯一の差は設計どおり B（10 vs 0）のみ

## 含意

- ハーネスの rep 単位決定論（seed 階層化）が第1相→第2相の run 跨ぎで実証された
- 第2相 boots / boots@NN 行は、第1相 rep 1–25 の raw/rarefied 行と**同一データ・
  同一 fit に対する条件追加**として扱える（R2-1 追補の前提成立）
- 第2相 summary 内の raw/rarefied 行（n_reps=25）は第1相 rep 1–25 の
  サブセット集計と等価
