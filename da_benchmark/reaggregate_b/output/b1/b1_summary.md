# B-1: ANCOM-BC2 diff_robust 再集計 — 既定判定との並記（CG-05, GM-01）

- 実施日: 2026-08-16（前倒し着手）／再フィットなし（凍結ローカルミラーの純再集計）
- 規則: **default** = `detected`（q_group < 0.05、主判定）／**diff_robust** = `detected & passed_ss`（擬似カウント感度フィルタ通過）
- 範囲: method = ancombc2、condition ∈ {raw, rarefied}。**boots は判断2 のハイブリッド規則で pseudo_sens=FALSE のため設計上対象外**
- 入力: `dsim_main_v1`（2,100 rep、18 セル）＋ `dreal_main_v1`（200 rep、2 セル × 2 割当）
- 自己検査: default 規則の再計算値は両 run の凍結 summary と全指標一致（最大差 4e-14 = 機械精度）
- スクリプト: `reagg_b/R/b1_diff_robust.R`（成果物 TSV 5 本 + 本メモ）

## D-sim 主要セル（並記表の要約）

FWER = P(偽陽性 ≥ 1)、FDR = mean FDP、power = intention-to-test。survival = diff_robust 検出数 / default 検出数（総数比）。

| セル | 条件 | FWER 既定→robust | FDR 既定→robust | power_itt 既定→robust | survival |
|---|---|---|---|---|---|
| null_r1 | raw | 1.000 → 1.000 | 1.000 → 1.000 | — | 0.563 |
| null_r10 | raw | 0.993 → 0.767 | 0.993 → 0.767 | — | 0.158 |
| **null_r100** | **raw** | **0.927 → 0.373** | 0.927 → 0.373 | — | **0.070** |
| null_r1 | rarefied | 1.000 → 1.000 | 1.000 → 1.000 | — | 0.500 |
| null_r10 | rarefied | 1.000 → 1.000 | 1.000 → 1.000 | — | 0.480 |
| **null_r100** | **rarefied** | **1.000 → 1.000** | 1.000 → 1.000 | — | 0.494 |
| spike_r1_f3 | raw | 1.000 → 1.000 | 0.841 → 0.838 | 0.185 → 0.104 | 0.549 |
| spike_r10_f3 | raw | 0.993 → 0.840 | 0.823 → 0.746 | 0.124 → 0.019 | 0.189 |
| **spike_r100_f3** | **raw** | 0.980 → 0.460 | **0.813 → 0.447** | **0.084 → 0.003** | **0.066** |
| spike_r1_f3 | rarefied | 1.000 → 1.000 | 0.844 → 0.834 | 0.250 → 0.132 | 0.508 |
| spike_r10_f3 | rarefied | 1.000 → 1.000 | 0.845 → 0.837 | 0.230 → 0.117 | 0.494 |
| **spike_r100_f3** | **rarefied** | 1.000 → 1.000 | **0.842 → 0.841** | 0.240 → 0.124 | 0.522 |

感度 12 セル分は `b1_dsim_ancombc2_rules_by_cell.tsv` / `_wide.tsv` に全収載。

## D-real（記述統計・検出数の中央値、tested ≈ 2,800）

| セル | 条件 | 割当 | det 中央値 既定→robust |
|---|---|---|---|
| rice_strat | raw | depth | 1,696.5 → 194.5 |
| rice_strat | raw | perm | 1,238 → 306.5 |
| rice_strat | rarefied | depth | 1,588 → 84.5 |
| rice_strat | rarefied | perm | 1,378 → 73.5 |
| rice_rhizo | raw | depth | 1,484 → 140.5 |
| rice_rhizo | raw | perm | 1,144.5 → 225 |
| rice_rhizo | rarefied | depth | 1,161 → 122 |
| rice_rhizo | rarefied | perm | 989 → 56 |

ペア内差（depth − perm、中央値 [IQR]、深度側が多いペアの割合）:

| セル | 条件 | default | diff_robust |
|---|---|---|---|
| rice_strat | raw | +453.5 [+398.5, +515.25]、100% | **−107 [−152.25, −68.75]、2%** |
| rice_strat | rarefied | +211.5 [+173.75, +279.25]、100% | +5.5 [−17.25, +20]、55% |
| rice_rhizo | raw | +345 [+277.5, +403.25]、100% | **−85 [−111, −44]、6%** |
| rice_rhizo | rarefied | +173 [+129.5, +229]、99% | +62 [+28, +82.25]、84% |

## 観察（8/20 チェックポイントの判断材料 — 表現強度はそこで確定）

1. **rarefied の inflation は diff_robust で消えない。** D-sim 帰無の rarefied は全 ρ で FWER 1.000 のまま（検出の約半数は個別に落ちるが、残存検出だけで毎 rep 偽陽性が出る）。スパイクの FDR も 0.84 → 0.84 と不変。
2. **raw の超過は大部分が sensitivity-fragile。** ρ=100 で検出の 93% が passed_ss 不通過。robust 化で FWER 0.927→0.373 と大幅減衰するが、**名目 0.05 の約 7.5 倍が残存**。一方 robust 化の代償で raw の検出力はほぼ消失（ITT 0.084→0.003 @ρ=100）。
3. CG-05 裁定の枠組みでは「残存なら維持、消えるなら緩和」の**中間形**: 「既定判定では大規模な超過検出。公式 sensitivity filter（diff_robust）併記では raw の超過は大幅に減衰するが名目水準への回復はなく、rarefied の inflation はフィルタ非感受」が実測に忠実な要約。
4. **D-real の深度分割超過（CG-03 の depth-associated excess）は、既定判定では全条件で明瞭（+173〜+454）だが、diff_robust では raw で逆転（perm 側が多い）、rice_strat rarefied で消失。** 超過検出の主成分が擬似カウント感度フラグ付き検出であることを示す。R2-1 MEMO の文言（GM-03 の実測 +128〜+611 引用部）に diff_robust 並記を追補する際はこの反転も記述する。

## 成果物

- `b1_dsim_ancombc2_rules_by_cell.tsv`（18 セル × 2 条件 × 2 規則、長形式）
- `b1_dsim_ancombc2_rules_wide.tsv`（並記＋survival）
- `b1_dreal_ancombc2_rules_by_cell.tsv` / `b1_dreal_ancombc2_rules_paired.tsv` / `b1_dreal_ancombc2_rules_per_rep.tsv`
- `b1_provenance.txt`
