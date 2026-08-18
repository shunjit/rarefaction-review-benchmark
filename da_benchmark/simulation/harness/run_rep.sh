#!/bin/zsh
# One rep end-to-end: generate -> wilcoxon -> ancombc2 -> maaslin3 -> collect.
# usage: run_rep.sh <run_dir> <cell_id> <rep>
# Idempotent: exits 0 immediately when the completion marker exists.
# On any failure no marker is written, so the next orchestrate re-runs the rep.
set -euo pipefail

run_dir=$1; cell_id=$2; rep=$3
rep4=$(printf '%04d' "$rep")
marker="$run_dir/markers/${cell_id}_rep${rep4}.done"
[[ -f "$marker" ]] && exit 0

source "$run_dir/run.env"    # BASE_SEED, PSEUDO_SENS, REVISION_ROOT, HARNESS_R_DIR,
                             # ANCOMBC_RSCRIPT, MAASLIN3_RSCRIPT, WORKROOT

log="$run_dir/log/${cell_id}_rep${rep4}.log"
exec >>"$log" 2>&1
echo "=== $(date -u +%Y-%m-%dT%H:%M:%SZ) start ${cell_id} rep ${rep}"

row=$(awk -F'\t' -v c="$cell_id" 'NR>1 && $2==c' "$run_dir/cells.tsv")
[[ -n "$row" ]] || { echo "ERROR: cell_id $cell_id not found in cells.tsv"; exit 1; }
IFS=$'\t' read -r CELL_INDEX CELL_ID SCENARIO K N_PER_GROUP RHO D0 THETA FC DA_FRAC B REPS <<< "$row"
export CELL_INDEX CELL_ID SCENARIO K N_PER_GROUP RHO D0 THETA FC DA_FRAC B
export REP=$rep BASE_SEED PSEUDO_SENS REVISION_ROOT HARNESS_R_DIR

# Pin every numeric library to 1 thread: rep-level parallelism only. Without
# this, 8 workers x multithreaded BLAS oversubscribe the cores ~10x (measured
# 52 s/fit vs 5.5 s/fit for ancombc2 on 2026-08-11).
export OMP_NUM_THREADS=1 OPENBLAS_NUM_THREADS=1 MKL_NUM_THREADS=1 \
       VECLIB_MAXIMUM_THREADS=1 NUMEXPR_NUM_THREADS=1

export WORKDIR="$WORKROOT/${cell_id}_rep${rep4}"
rm -rf "$WORKDIR"; mkdir -p "$WORKDIR"

# All stages run on arm64-native R. The qiime2 env R is x86_64 (Rosetta) and
# collapses under concurrency on Apple Silicon (measured 2026-08-11: ancombc2
# 5.5 s/fit at j=2 -> 52 s/fit at j=8, while native maaslin3 stayed ~5 s).
# generate/wilcoxon need base R only, so they borrow the maaslin3 env.
# stdin is force-closed on every stage: children of xargs inherit its stdin
# (the job queue), and ANCOMBC 2.12 prompts interactively when pseudo_sens is
# off — an un-redirected read would swallow queue lines.
"$MAASLIN3_RSCRIPT" "$HARNESS_R_DIR/generate_rep.R"  < /dev/null
"$MAASLIN3_RSCRIPT" "$HARNESS_R_DIR/fit_wilcoxon.R"  < /dev/null
"$ANCOMBC_RSCRIPT"  "$HARNESS_R_DIR/fit_ancombc2.R"  < /dev/null
"$MAASLIN3_RSCRIPT" "$HARNESS_R_DIR/fit_maaslin3.R"  < /dev/null

outdir="$run_dir/results/$cell_id"
mkdir -p "$outdir"
all="$WORKDIR/all_results.tsv"
head -1 "$WORKDIR/results_wilcoxon.tsv" > "$all"
awk 'FNR > 1' "$WORKDIR"/results_wilcoxon.tsv "$WORKDIR"/results_ancombc2.tsv \
              "$WORKDIR"/results_maaslin3.tsv >> "$all"
gzip -c "$all" > "$outdir/rep_${rep4}.tsv.gz.tmp"
mv "$outdir/rep_${rep4}.tsv.gz.tmp" "$outdir/rep_${rep4}.tsv.gz"
gzip -c "$WORKDIR/truth.tsv" > "$outdir/rep_${rep4}_truth.tsv.gz.tmp"
mv "$outdir/rep_${rep4}_truth.tsv.gz.tmp" "$outdir/rep_${rep4}_truth.tsv.gz"
cp "$WORKDIR/geninfo.tsv" "$outdir/rep_${rep4}_geninfo.tsv"
cp "$WORKDIR/timings.tsv" "$outdir/rep_${rep4}_timings.tsv"

touch "$marker"
rm -rf "$WORKDIR"
echo "=== $(date -u +%Y-%m-%dT%H:%M:%SZ) done ${cell_id} rep ${rep}"
