#!/usr/bin/env python3
"""
D-sim calibration step 1: extract the Rice Rhizoplane ASV count table.

Reads the denoised feature table (Rice_v2/02_denoised/table.qza) and the
Foundation-A sample ledger, keeps Rhizoplane samples that are present in the
feature table, drops ASVs that are all-zero within that subset, and writes a
plain features-x-samples TSV for the Dirichlet-multinomial fit (fit_dm.R).

Design reference: D_HARNESS_DESIGN.md v2 section 4 (Rice-calibrated parametric
simulation; ASV-level because the current pipeline has no taxonomy artifacts).

Environment contract (hard stop when missing):
    ANALYSIS_ROOT   e.g. .../benchmark_data
    REVISION_ROOT   e.g. ${ANALYSIS_ROOT}/revision_r1

Run with the qiime2 env python (needs biom-format + pandas):
    ANALYSIS_ROOT=... REVISION_ROOT=... \
        ~/miniforge3/envs/qiime2-amplicon-2025.10/bin/python extract_rice_rhizoplane.py
"""

import os
import sys
import zipfile
import tempfile
from datetime import datetime, timezone

import biom
import pandas as pd


def require_env(name: str) -> str:
    value = os.environ.get(name, "").strip()
    if not value:
        sys.stderr.write(f"ERROR: environment variable {name} is MISSING — stopping.\n")
        sys.exit(1)
    if not os.path.isdir(value):
        sys.stderr.write(f"ERROR: {name}={value} is not a directory — stopping.\n")
        sys.exit(1)
    return value


ANALYSIS_ROOT = require_env("ANALYSIS_ROOT")
REVISION_ROOT = require_env("REVISION_ROOT")

TABLE_QZA = os.path.join(ANALYSIS_ROOT, "Rice_v2", "02_denoised", "table.qza")
LEDGER = os.path.join(REVISION_ROOT, "sample_ledger", "output", "rice_sample_ledger.tsv")
OUT_DIR = os.path.join(REVISION_ROOT, "simulation", "harness", "calibrate", "output")
OUT_TSV = os.path.join(OUT_DIR, "rice_rhizoplane_asv_counts.tsv")
OUT_PROV = os.path.join(OUT_DIR, "extract_provenance.md")

for path in (TABLE_QZA, LEDGER):
    if not os.path.exists(path):
        sys.stderr.write(f"ERROR: required input MISSING: {path} — stopping.\n")
        sys.exit(1)
os.makedirs(OUT_DIR, exist_ok=True)


def load_biom_from_qza(qza_path: str) -> biom.Table:
    with zipfile.ZipFile(qza_path) as zf:
        names = [n for n in zf.namelist() if n.endswith("data/feature-table.biom")]
        if len(names) != 1:
            sys.stderr.write(f"ERROR: expected exactly 1 feature-table.biom in {qza_path}, "
                             f"found {len(names)} — stopping.\n")
            sys.exit(1)
        with tempfile.TemporaryDirectory() as td:
            zf.extract(names[0], td)
            return biom.load_table(os.path.join(td, names[0]))


ledger = pd.read_csv(LEDGER, sep="\t", dtype=str)
rhizo = ledger[(ledger["compartment"] == "Rhizoplane") & ledger["table_depth"].notna()
               & (ledger["table_depth"] != "")]
rhizo_ids = set(rhizo["sample-id"])

table = load_biom_from_qza(TABLE_QZA)
table_ids = set(table.ids("sample"))
keep = sorted(rhizo_ids & table_ids)
missing = sorted(rhizo_ids - table_ids)
if missing:
    sys.stderr.write(f"ERROR: {len(missing)} ledger Rhizoplane samples absent from table "
                     f"(e.g. {missing[:3]}) — stopping (ledger/table mismatch).\n")
    sys.exit(1)

sub = table.filter(keep, axis="sample", inplace=False)
sub = sub.remove_empty(axis="observation", inplace=False)

df = pd.DataFrame(sub.matrix_data.toarray(),
                  index=sub.ids("observation"), columns=sub.ids("sample"))
df = df.astype(int)
df.index.name = "asv_id"
df.to_csv(OUT_TSV, sep="\t")

total_reads = int(df.values.sum())
with open(OUT_PROV, "w") as fh:
    fh.write(
        "# Extraction provenance — rice_rhizoplane_asv_counts.tsv\n\n"
        f"- generated: {datetime.now(timezone.utc).isoformat(timespec='seconds')}\n"
        "- source table: `Rice_v2/02_denoised/table.qza` (relative to ANALYSIS_ROOT)\n"
        "- sample selection: `sample_ledger/output/rice_sample_ledger.tsv`, "
        "compartment == Rhizoplane AND present in feature table\n"
        f"- samples kept: {len(keep)}\n"
        f"- ASVs (non-empty within subset): {df.shape[0]}\n"
        f"- total reads: {total_reads}\n"
        "- script: `simulation/harness/calibrate/extract_rice_rhizoplane.py`\n"
    )

print(f"OK: {df.shape[0]} ASVs x {df.shape[1]} samples, {total_reads} reads -> {OUT_TSV}")
