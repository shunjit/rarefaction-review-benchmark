# Depth-Confounding Differential Abundance Benchmark

This directory archives the differential abundance (DA) benchmark introduced in
the manuscript revision (Supplementary Text S2): the code, frozen run
configurations, per-run environment records, and summary outputs for

- **D-sim** (`simulation/`) — a Dirichlet-multinomial simulation study
  calibrated to the Rice rhizoplane compartment, with depth–group confounding
  imposed by design (null and spike-in cells), and
- **D-real** (`realdata_typeI/`) — a real-data Type I error study using
  artificial two-group splits of the Rice dataset, where any detection is a
  false positive, plus
- **common_set** — a common-sample-set re-analysis (n = 486) that re-runs the
  phase-1 Rice sensitivity PERMANOVA on the identical sample set across all
  four rarefaction depths, and
- **reaggregate_b/** — post-hoc re-aggregations of the frozen run outputs
  (no re-fitting), e.g. ANCOM-BC2 detection-rule sensitivity.

Three DA procedures are compared on rarefied vs unrarefied versions of the same
tables: **ANCOM-BC2**, **MaAsLin 3**, and the **Wilcoxon rank-sum test**
(representing bias-correction, linear-model, and non-parametric baseline
method classes).

---

## What is (and is not) archived

| Archived here | Not archived (regenerable) |
|---|---|
| Harness code (`harness/R`, shell drivers) | Per-replicate run bodies (`results/`, `reps/`, logs; ~8.8 GB) |
| Calibration code and outputs (`simulation/harness/calibrate/`, incl. the extracted rhizoplane ASV table) | Intermediate QIIME 2 artifacts of the phase-1 pipeline (see `../provenance/`) |
| Frozen run configurations (`runs/<id>/cells.tsv`, `config.yaml`) | |
| Execution records (`runs/<id>/run.env`) | |
| Exact R environments (`runs/<id>/sessionInfo/`) | |
| All summary outputs (`runs/<id>/summary/`) | |
| Figure scripts for the supplementary DA figures (`*/figures/*.py`) | Figure PDFs (regenerable from the summaries) |
| Re-aggregation code and outputs (`reaggregate_b/R`, `reaggregate_b/output`) | |

Every replicate is deterministically reproducible: each run derives all random
seeds hierarchically from the single `BASE_SEED` recorded in its `run.env`
(D-sim: 20260811, D-real: 20260813), so re-running a cell/replicate with the
archived code and configuration reproduces the frozen summaries. The second
phase of D-real (`dreal_boots_v1`) deliberately re-computes replicates 1–25 of
the first phase and was verified byte-identical on all 50 cell × replicate
pairs (`runs/dreal_boots_v1/summary/G4_verification.md`).

`run.env` files are **execution records**, not editable configuration: they
document the absolute paths and environment of the original analysis machine
and are archived verbatim. To run the harness elsewhere, set the environment
variables described below instead of editing any archived file.

---

## Archived runs

| Run | Design | Scale | Role in the manuscript |
|---|---|---|---|
| `simulation/runs/dsim_main_v1` | 18 cells (null + spike-in grid) | 2,100/2,100 reps, 0 failures | Main D-sim results (FWER/FDR, detection direction) |
| `simulation/runs/dsim_conv_b100_v1` | 2 cells, B=100 | 100/100 reps | Convergence check for the boots stability score |
| `simulation/runs/dsim_a2_nsens_v1` | 2 spike cells, n=10/50 per group | 200/200 reps | Sample-size sensitivity |
| `realdata_typeI/runs/dreal_main_v1` | 2 cells (stratified, rhizoplane), B=0 | 200/200 reps | Main D-real Type I error (depth-ordered vs restricted permutation) |
| `realdata_typeI/runs/dreal_boots_v1` | same cells, B=10 | 50/50 reps | Rarefaction-stability supplement; reproducibility check |
| `realdata_typeI/runs/dreal_a1_x25_v1` | 1 cell, extreme-quantile split (2.53×) | 100/100 reps | Split-strength sensitivity |
| `common_set/runs/a3_common486_v1` | 4 depths × fixed 486-sample table | 4/4 depths | Direct demonstration that Rice beta-diversity stability is not a sample-composition artefact (wUniFrac R² 0.506–0.508, all p = 0.001) |

`reaggregate_b/output/` contains the corresponding re-aggregations
(`b1`–`b6`, plus `b1_a` extending `b1` to the two sensitivity runs). Each
`b*_summary.md` states its inputs and method; `b1` self-checks by reproducing
the frozen default-rule summaries to machine precision before re-aggregating.

---

## Environments

The harness calls two R installations through environment variables (defaults
in parentheses):

- `ANCOMBC_RSCRIPT` — Rscript with the ANCOMBC package
  (`$HOME/miniforge3/envs/ancombc_native/bin/Rscript`)
- `MAASLIN3_RSCRIPT` — Rscript with maaslin3
  (`$HOME/miniforge3/envs/maaslin3_env/bin/Rscript`)

The exact package versions used for the manuscript are archived per run in
`runs/<id>/sessionInfo/ancombc_env_R.txt` and `maaslin3_env_R.txt`
(maaslin3 1.2.0 on R 4.5.3). Wilcoxon uses base R from the ANCOM-BC2
environment. No QIIME 2 installation is required for the DA harness itself.

---

## Running the harness

Set `REVISION_ROOT` to this directory (the harness resolves everything under
it: `simulation/`, `realdata_typeI/`, shared fit scripts in
`simulation/harness/R`):

```bash
export REVISION_ROOT=/path/to/rarefaction-review-benchmark/da_benchmark
export RUNS_ROOT=$HOME/dsim_runs   # hot run state on a local disk
```

**Smoke test (~3 minutes)** — one replicate of two micro cells end-to-end
(generate → Wilcoxon → ANCOM-BC2 → MaAsLin 3 → collect):

```bash
cd "$REVISION_ROOT/simulation/harness"
./night_launch.sh dsim_micro_$(date +%Y%m%d) configs/cells_micro.tsv 2
```

Progress: `ls $RUNS_ROOT/<run_id>/markers | wc -l` (2 markers = done). The
launcher runs detached via `nohup`, writes a log under `<run_dir>/log/`, and
rsyncs the finished run to `$REVISION_ROOT/simulation/runs/<run_id>/`.

**Main runs** — same drivers with the archived configurations, e.g.:

```bash
./night_launch.sh my_dsim_main configs/cells_main_v1.tsv 8          # D-sim main
cd "$REVISION_ROOT/realdata_typeI/harness"
./night_launch_real.sh my_dreal_main configs/cells_real_v1.tsv 8    # D-real main
```

D-real additionally requires the extracted Rice input tables (regenerate with
`realdata_typeI/harness/extract_rice_full.py` from the phase-1 Rice feature
table; input checksums are archived in each run's
`sessionInfo/input_checksums.sha256`). Expect ~30 h wall time per main run on
8 cores; `BASE_SEED`, timeouts, and pseudo-sensitivity switches are recorded
in each archived `run.env`.

**Common-set re-analysis** (requires QIIME 2 and the phase-1 Rice artifacts):

```bash
cd "$REVISION_ROOT/common_set/runs/a3_common486_v1/script"
REVISION_ROOT=... RAW_ROOT=/path/to/phase1_data ./run_a3_common486.sh
```

**Re-aggregations** read finished run directories (they need the per-replicate
bodies, i.e. local mirrors of completed runs — not just this archive) and write
to `reaggregate_b/output/`.

---

## Reading the summaries

Per-run `summary/` files are tab-separated with self-describing headers; the
main ones are `metrics_by_cell.tsv` / `metrics_real_by_cell.tsv` (per-cell
FWER/FDR/power aggregates), `per_rep_detections.tsv` (per-replicate detection
counts), `paired_depth_vs_perm.tsv` (paired split-type contrasts),
`timings_summary.tsv`, and for the common set `a3_permanova_summary.tsv` /
`a3_alpha_summary.tsv`. The supplementary DA figures are produced by
`simulation/figures/make_figures_r2_2.py`, `simulation/figures/make_figures_r1_1.py`,
and `realdata_typeI/figures/make_figures_r2_1.py` from these summaries.
