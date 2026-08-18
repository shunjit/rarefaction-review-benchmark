#!/bin/zsh
# One D-real rep end-to-end: generate (subsample + paired assignments) ->
# 3 methods x 2 assignments -> collect. The two assignments (depth-ordered
# split and restricted permutation) share ONE set of count tables; only
# meta.tsv differs, so each fit runs twice on symlinked tables.
# usage: run_rep_real.sh <run_dir> <cell_id> <rep>
# Idempotent: exits 0 immediately when the completion marker exists.
# On any failure no marker is written, so the next orchestrate re-runs the rep.
set -euo pipefail

run_dir=$1; cell_id=$2; rep=$3
rep4=$(printf '%04d' "$rep")
marker="$run_dir/markers/${cell_id}_rep${rep4}.done"
[[ -f "$marker" ]] && exit 0

source "$run_dir/run.env"    # BASE_SEED, PSEUDO_SENS, REVISION_ROOT, HARNESS_R_DIR,
                             # REAL_R_DIR, ANCOMBC_RSCRIPT, MAASLIN3_RSCRIPT, WORKROOT

log="$run_dir/log/${cell_id}_rep${rep4}.log"
exec >>"$log" 2>&1
echo "=== $(date -u +%Y-%m-%dT%H:%M:%SZ) start ${cell_id} rep ${rep}"

row=$(awk -F'\t' -v c="$cell_id" 'NR>1 && $2==c' "$run_dir/cells.tsv")
[[ -n "$row" ]] || { echo "ERROR: cell_id $cell_id not found in cells.tsv"; exit 1; }
IFS=$'\t' read -r CELL_INDEX CELL_ID DATASET STRAT_COL STRAT_MODE FILTER_VAL SUBSAMPLE_FRAC B REPS <<< "$row"
export CELL_INDEX CELL_ID DATASET STRAT_COL STRAT_MODE FILTER_VAL SUBSAMPLE_FRAC B
export REP=$rep BASE_SEED PSEUDO_SENS REVISION_ROOT HARNESS_R_DIR

# Pin every numeric library to 1 thread: rep-level parallelism only (same
# rationale and measurement as the D-sim harness, 2026-08-11).
export OMP_NUM_THREADS=1 OPENBLAS_NUM_THREADS=1 MKL_NUM_THREADS=1 \
       VECLIB_MAXIMUM_THREADS=1 NUMEXPR_NUM_THREADS=1

WORKDIR_ROOT="$WORKROOT/${cell_id}_rep${rep4}"
rm -rf "$WORKDIR_ROOT"; mkdir -p "$WORKDIR_ROOT"

# stdin force-closed on every stage (ANCOMBC 2.12 prompts when pseudo_sens is
# off; children of xargs must never swallow the job queue).
t_gen0=$(date +%s)
WORKDIR="$WORKDIR_ROOT" "$MAASLIN3_RSCRIPT" "$REAL_R_DIR/generate_real_rep.R" < /dev/null
echo "generate: $(( $(date +%s) - t_gen0 )) s"

# Two assignments over the same tables: depth-ordered split vs restricted
# permutation. SCENARIO is stamped into every result row by the fit scripts.
for assign in depth perm; do
  ad="$WORKDIR_ROOT/$assign"
  mkdir -p "$ad"
  ln -s "$WORKDIR_ROOT/raw.tsv"      "$ad/raw.tsv"
  ln -s "$WORKDIR_ROOT/rarefied.tsv" "$ad/rarefied.tsv"
  ln -s "$WORKDIR_ROOT/boots"        "$ad/boots"
  cp "$WORKDIR_ROOT/meta_${assign}.tsv" "$ad/meta.tsv"
  scenario=$([[ "$assign" == depth ]] && echo "real_depth" || echo "real_perm")

  WORKDIR="$ad" SCENARIO="$scenario" "$MAASLIN3_RSCRIPT" "$HARNESS_R_DIR/fit_wilcoxon.R"  < /dev/null
  WORKDIR="$ad" SCENARIO="$scenario" "$ANCOMBC_RSCRIPT"  "$HARNESS_R_DIR/fit_ancombc2.R"  < /dev/null
  WORKDIR="$ad" SCENARIO="$scenario" "$MAASLIN3_RSCRIPT" "$HARNESS_R_DIR/fit_maaslin3.R"  < /dev/null
done

outdir="$run_dir/results/$cell_id"
mkdir -p "$outdir"
all="$WORKDIR_ROOT/all_results.tsv"
head -1 "$WORKDIR_ROOT/depth/results_wilcoxon.tsv" > "$all"
awk 'FNR > 1' "$WORKDIR_ROOT"/depth/results_wilcoxon.tsv "$WORKDIR_ROOT"/depth/results_ancombc2.tsv \
              "$WORKDIR_ROOT"/depth/results_maaslin3.tsv \
              "$WORKDIR_ROOT"/perm/results_wilcoxon.tsv "$WORKDIR_ROOT"/perm/results_ancombc2.tsv \
              "$WORKDIR_ROOT"/perm/results_maaslin3.tsv >> "$all"
gzip -c "$all" > "$outdir/rep_${rep4}.tsv.gz.tmp"
mv "$outdir/rep_${rep4}.tsv.gz.tmp" "$outdir/rep_${rep4}.tsv.gz"
gzip -c "$WORKDIR_ROOT/truth.tsv" > "$outdir/rep_${rep4}_truth.tsv.gz.tmp"
mv "$outdir/rep_${rep4}_truth.tsv.gz.tmp" "$outdir/rep_${rep4}_truth.tsv.gz"
cp "$WORKDIR_ROOT/geninfo.tsv" "$outdir/rep_${rep4}_geninfo.tsv"

# Merge the per-assignment timing files, tagging rows with the assignment so
# cost accounting can separate the pair (schema: + assign column vs D-sim).
{ head -1 "$WORKDIR_ROOT/depth/timings.tsv" | awk '{print $0 "\tassign"}'
  awk 'FNR > 1 {print $0 "\tdepth"}' "$WORKDIR_ROOT/depth/timings.tsv"
  awk 'FNR > 1 {print $0 "\tperm"}'  "$WORKDIR_ROOT/perm/timings.tsv"
} > "$outdir/rep_${rep4}_timings.tsv"

touch "$marker"
rm -rf "$WORKDIR_ROOT"
echo "=== $(date -u +%Y-%m-%dT%H:%M:%SZ) done ${cell_id} rep ${rep}"
