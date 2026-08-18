#!/usr/bin/env python3
"""D-real input prep: extract the full Rice ASV count table + analysis metadata.

Reads the denoised feature table (Rice_v2/02_denoised/table.qza) and the
Foundation-A sample ledger, keeps every sample that is present in the feature
table (all compartments; n = 493 as of the 2026-08 ledger), drops ASVs that
are all-zero over that set, and writes:

    realdata_typeI/input/rice_asv_counts.tsv   (features x samples, integer)
    realdata_typeI/input/rice_meta.tsv         (sample-id, compartment, depth)
    realdata_typeI/input/extract_provenance.md

Design reference: D_HARNESS_DESIGN.md v2 section 5 (D-real; stratified split
uses `compartment` as the stratum, robustness cell restricts to Rhizoplane).

Environment contract (hard stop when missing):
    ANALYSIS_ROOT   e.g. .../benchmark_data
    REVISION_ROOT   e.g. ${ANALYSIS_ROOT}/revision_r1

Run with the qiime2 env python (needs biom-format + pandas):
    ANALYSIS_ROOT=... REVISION_ROOT=... \
        ~/miniforge3/envs/qiime2-amplicon-2025.10/bin/python extract_rice_full.py
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
OUT_DIR = os.path.join(REVISION_ROOT, "realdata_typeI", "input")
OUT_TSV = os.path.join(OUT_DIR, "rice_asv_counts.tsv")
OUT_META = os.path.join(OUT_DIR, "rice_meta.tsv")
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
in_table = ledger[ledger["table_depth"].notna() & (ledger["table_depth"] != "")]
ledger_ids = set(in_table["sample-id"])

table = load_biom_from_qza(TABLE_QZA)
table_ids = set(table.ids("sample"))
if ledger_ids != table_ids:
    only_ledger = sorted(ledger_ids - table_ids)
    only_table = sorted(table_ids - ledger_ids)
    sys.stderr.write(f"ERROR: ledger/table sample mismatch — {len(only_ledger)} only in "
                     f"ledger (e.g. {only_ledger[:3]}), {len(only_table)} only in table "
                     f"(e.g. {only_table[:3]}) — stopping.\n")
    sys.exit(1)

keep = sorted(ledger_ids)
sub = table.filter(keep, axis="sample", inplace=False)
sub = sub.remove_empty(axis="observation", inplace=False)

df = pd.DataFrame(sub.matrix_data.toarray(),
                  index=sub.ids("observation"), columns=sub.ids("sample"))
df = df.astype(int)
df.index.name = "asv_id"
df.to_csv(OUT_TSV, sep="\t")

meta = in_table[["sample-id", "compartment"]].copy()
meta["depth"] = df.sum(axis=0).reindex(meta["sample-id"]).values
lm = in_table.set_index("sample-id")["table_depth"].astype(float)
if not (meta.set_index("sample-id")["depth"] == lm.reindex(meta["sample-id"]).values).all():
    sys.stderr.write("ERROR: extracted depths do not match ledger table_depth — stopping.\n")
    sys.exit(1)
meta = meta.sort_values("sample-id")
meta.to_csv(OUT_META, sep="\t", index=False)

total_reads = int(df.values.sum())
comp_counts = meta["compartment"].value_counts().sort_index()
with open(OUT_PROV, "w") as fh:
    fh.write(
        "# Extraction provenance — rice_asv_counts.tsv / rice_meta.tsv\n\n"
        f"- generated: {datetime.now(timezone.utc).isoformat(timespec='seconds')}\n"
        "- source table: `Rice_v2/02_denoised/table.qza` (relative to ANALYSIS_ROOT)\n"
        "- sample selection: `sample_ledger/output/rice_sample_ledger.tsv`, "
        "table_depth non-null (= imported into the feature table)\n"
        f"- samples kept: {len(keep)}\n"
        f"- compartments: {'; '.join(f'{k} {v}' for k, v in comp_counts.items())}\n"
        f"- ASVs (non-empty over kept samples): {df.shape[0]}\n"
        f"- total reads: {total_reads}\n"
        "- depth check: column sums equal ledger `table_depth` for all samples\n"
        "- script: `realdata_typeI/harness/extract_rice_full.py`\n"
    )

print(f"OK: {df.shape[0]} ASVs x {df.shape[1]} samples, {total_reads} reads -> {OUT_TSV}")
