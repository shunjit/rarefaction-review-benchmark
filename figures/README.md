# Manuscript Figure Generation Code

Scripts and input tables for the manuscript's main figures. The scripts are
archived byte-identical to the versions that produced the submitted/revised
figures; where a script embeds an absolute input path from the original
analysis machine, the one- or two-line change needed to run it elsewhere is
documented below (the archived files themselves are left frozen).

| Figure | Script | Data inputs | Notes |
|---|---|---|---|
| Figure 1 (decision framework) | `make_figure1_r1.py` | none (drawn programmatically) | Revised Step 0 rule (four-category), MaAsLin 3 reference |
| Figure 2 (alpha diversity vs rarefaction depth) | `make_figure2_r1.py` | `data/shannon_{soil,bioethanol,rice}.csv` | Writes `fig2_r1.pdf` / `.png` to the current directory |
| Figure 3 (q2-boots alpha stability) | `figure3_alpha_stability_v2_1.py` | `data/alpha_cv_summary_{soil,bioethanol,rice}.csv` | Unchanged from submission |
| Figure 4 (DA schematic) | `make_figure4_r1.py` | none (drawn programmatically) | |

The supplementary DA benchmark figures are generated separately from the
archived run summaries; see `../da_benchmark/README.md`.

## Input tables

`data/` contains the frozen summary tables the figures read:

- `shannon_*.csv` — per-sample Shannon diversity across rarefaction depths
  (phase-1 sensitivity analysis outputs; one file per dataset).
- `alpha_cv_summary_*.csv` — per-sample coefficient of variation of Shannon
  diversity over 100 q2-boots rarefaction iterations at the primary depth
  (inputs to Figure 3 and the manuscript's alpha-stability CV values:
  median CV Soil 0.17%, Bioethanol 0.45%, Rice 0.15%).

## Environment

Python with matplotlib 3.8.4, pandas, numpy, seaborn (the QIIME 2
amplicon 2025.10 environment described in the repository README was used;
fonts: Arial/Helvetica).

## Running elsewhere

`make_figure2_r1.py` and `figure3_alpha_stability_v2_1.py` embed the original
machine's data locations near the top of the file:

```python
DATA_DIR = '/Volumes/PS3000/benchmark_data/data_for_figure2'   # or ..._figure3
OUTPUT_DIR = '/Volumes/PS3000/benchmark_data/data_for_figure2/output'
```

To run from this repository, copy the script and point those two assignments
at this directory, e.g.:

```bash
cd figures
sed "s|^DATA_DIR = .*|DATA_DIR = 'data'|; s|^OUTPUT_DIR = .*|OUTPUT_DIR = 'output'|" \
  make_figure2_r1.py > /tmp/make_figure2_local.py
python /tmp/make_figure2_local.py
```

(`make_figure1_r1.py` and `make_figure4_r1.py` need no data and run as-is.)

## Regeneration verification

Both data-driven figures were re-generated from these archived scripts and
input tables and compared against the frozen figure PDFs (2026-08-19):

- **Figure 2**: regenerated PDF is content-identical to the frozen `fig2_r1.pdf`
  — identical `pdftotext` extraction and pixel-identical 150-dpi rasterization.
- **Figure 3**: identical text/number content and panel geometry; the only
  difference from the submitted `fig3.pdf` is a uniform sub-point (~0.2 pt)
  glyph-position drift attributable to matplotlib version differences since
  the original January render (PDF byte equality is not expected for
  regenerated matplotlib output).
