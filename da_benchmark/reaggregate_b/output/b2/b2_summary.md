# B-2: D0 別の実現ゼロ率表（CG-07 後段）

- 実施日: 2026-08-16／入力 = `dsim_main_v1` 全 2,100 rep の per-rep geninfo（新規計算なし）
- D0 は事前登録された「サンプリングゼロのレバー」（裁定 CG-07: 呼称限定つき）。名目でなく**実現**ゼロ率を引用するための表

## D0 レバー系列（tested taxa、rep 平均）

| セル | D0 | 単一 rarefaction 深度（中央値） | ゼロ率 raw | ゼロ率 rarefied |
|---|---|---|---|---|
| null_r100_d2500 | 2,500 | 1,396 | 0.491 | 0.651 |
| null_r100 | 10,000 | 5,711 | 0.424 | 0.557 |
| null_r100_d40000 | 40,000 | 23,056 | 0.368 | 0.478 |
| spike_r100_f3_d2500 | 2,500 | 1,467 | 0.492 | 0.651 |
| spike_r100_f3 | 10,000 | 5,858 | 0.426 | 0.558 |
| spike_r100_f3_d40000 | 40,000 | 22,849 | 0.369 | 0.481 |

- D0 の 16 倍レンジで実現ゼロ率は raw 0.37→0.49、rarefied 0.48→0.65 と単調。rarefaction は各セルでゼロ率を約 +0.11〜+0.17 追加する（最小深度への縮約のため）
- 全 18 セル分（min/max・SD・achieved ρ・tested_n 含む）は `b2_realized_zero_rates_by_cell.tsv`

## 成果物

- `b2_realized_zero_rates_by_cell.tsv`（18 セル）
- `b2_realized_zero_rates_d0_lever.tsv`（D0 レバー 6 セル抜粋）
