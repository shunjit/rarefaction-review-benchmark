#!/bin/zsh
# D-real orchestrator: build the (cell x rep) queue, skip completed markers,
# run under xargs -P + caffeinate. Resumable at rep granularity.
# Mirrors simulation/harness/orchestrate.sh; differences are the cells schema
# (9 columns, REPS = $9), run_rep_real.sh, a one-time tsv->rds input
# conversion, and a larger default timeout (real-scale fits).
#
# usage: REVISION_ROOT=... ./orchestrate_real.sh <run_id> <cells_tsv> [jobs]
#   env overrides: BASE_SEED (default 20260813), PSEUDO_SENS (default TRUE),
#                  TIMEOUT_SEC (default 14400), DRY=1 (print queue size and exit)
set -euo pipefail

[[ -n "${REVISION_ROOT:-}" && -d "${REVISION_ROOT:-}" ]] || { echo "ERROR: REVISION_ROOT MISSING - stopping." >&2; exit 1; }
run_id=$1; cells_src=$2; jobs=${3:-8}
[[ -f "$cells_src" ]] || { echo "ERROR: cells file $cells_src MISSING" >&2; exit 1; }

module="$REVISION_ROOT/realdata_typeI"
harness="$module/harness"
sim_r="$REVISION_ROOT/simulation/harness/R"    # fit scripts + lib_common live here
# Hot run state on LOCAL APFS; canonical copy on PS3000 is written once by a
# single-writer rsync afterwards (night_launch_real.sh / sync_runs_real.sh).
# Same rationale as D-sim (exFAT visibility incident, 2026-08-11).
runs_root="${RUNS_ROOT:-$HOME/dreal_runs}"
run_dir="$runs_root/$run_id"
mkdir -p "$run_dir"/{results,log,markers,sessionInfo}

# cells.tsv is frozen per run: first launch copies it, later launches must match
if [[ -f "$run_dir/cells.tsv" ]]; then
  cmp -s "$cells_src" "$run_dir/cells.tsv" || { echo "ERROR: cells.tsv differs from the frozen copy in $run_dir - refusing to resume with different parameters." >&2; exit 1; }
else
  cp "$cells_src" "$run_dir/cells.tsv"
fi

ancombc_rscript=${ANCOMBC_RSCRIPT:-$HOME/miniforge3/envs/ancombc_native/bin/Rscript}
maaslin3_rscript=${MAASLIN3_RSCRIPT:-$HOME/miniforge3/envs/maaslin3_env/bin/Rscript}
[[ -x "$ancombc_rscript" ]] || { echo "ERROR: ANCOMBC_RSCRIPT $ancombc_rscript not executable" >&2; exit 1; }
[[ -x "$maaslin3_rscript" ]] || { echo "ERROR: MAASLIN3_RSCRIPT $maaslin3_rscript not executable" >&2; exit 1; }

# One-time input conversion: plain-tsv counts -> .rds for fast per-rep loads.
# Runs before the queue starts (single writer, no worker races).
for dataset in $(awk -F'\t' 'NR>1 {print $3}' "$cells_src" | sort -u); do
  tsv="$module/input/${dataset}_asv_counts.tsv"
  rds="$module/input/${dataset}_asv_counts.rds"
  [[ -f "$tsv" ]] || { echo "ERROR: input $tsv MISSING (run extract_${dataset}_full.py first)" >&2; exit 1; }
  if [[ ! -f "$rds" || "$tsv" -nt "$rds" ]]; then
    echo "converting $tsv -> $rds ..."
    "$maaslin3_rscript" -e '
      a <- commandArgs(TRUE)
      m <- as.matrix(read.delim(a[1], row.names = 1, check.names = FALSE))
      storage.mode(m) <- "integer"
      saveRDS(m, paste0(a[2], ".tmp"))
      file.rename(paste0(a[2], ".tmp"), a[2])
      cat(sprintf("rds written: %d taxa x %d samples\n", nrow(m), ncol(m)))
    ' "$tsv" "$rds" < /dev/null
  fi
done

if [[ ! -f "$run_dir/run.env" ]]; then
  cat > "$run_dir/run.env" <<EOF
BASE_SEED=${BASE_SEED:-20260813}
PSEUDO_SENS=${PSEUDO_SENS:-TRUE}
REVISION_ROOT="$REVISION_ROOT"
HARNESS_R_DIR="$sim_r"
REAL_R_DIR="$harness/R"
ANCOMBC_RSCRIPT="$ancombc_rscript"
MAASLIN3_RSCRIPT="$maaslin3_rscript"
WORKROOT="\$HOME/.cache/dreal_work/$run_id"
EOF
fi
source "$run_dir/run.env"
mkdir -p "$HOME/.cache/dreal_work/$run_id"

if [[ ! -f "$run_dir/config.yaml" ]]; then
  cat > "$run_dir/config.yaml" <<EOF
run_id: $run_id
created_utc: $(date -u +%Y-%m-%dT%H:%M:%SZ)
design: D_HARNESS_DESIGN.md v2 section 5 (D-real; depth-associated empirical stress test)
base_seed: $BASE_SEED
pseudo_sens: $PSEUDO_SENS
cells: cells.tsv (frozen copy in this directory; 9-column real schema)
scenarios: "real_depth = stratified median-depth split; real_perm = paired restricted permutation (same subsample, labels shuffled within stratum)"
seed_hierarchy: "dataset = base + cell_index*1e7 + rep*1e3; +901 single rarefaction; +1..B boots; +902 label permutation; +903 method"
prefilter: "prevalence >= 0.10 on the raw subsample table; intention-to-test set shared by all arms and both assignments"
environments:
  ancombc_rscript: $ancombc_rscript
  maaslin3_rscript: $maaslin3_rscript
host: $(hostname -s)
EOF
  (cd "$harness" && shasum -a 256 R/*.R run_rep_real.sh orchestrate_real.sh extract_rice_full.py) > "$run_dir/sessionInfo/script_checksums.sha256"
  (cd "$REVISION_ROOT/simulation/harness" && shasum -a 256 R/lib_common.R R/fit_wilcoxon.R R/fit_ancombc2.R R/fit_maaslin3.R) >> "$run_dir/sessionInfo/script_checksums.sha256"
  shasum -a 256 "$module"/input/*_asv_counts.tsv "$module"/input/*_meta.tsv > "$run_dir/sessionInfo/input_checksums.sha256"
  "$ancombc_rscript"  -e 'suppressMessages(library(ANCOMBC)); writeLines(capture.output(sessionInfo()))' > "$run_dir/sessionInfo/ancombc_env_R.txt" 2>/dev/null
  "$maaslin3_rscript" -e 'suppressMessages(library(maaslin3)); writeLines(capture.output(sessionInfo()))' > "$run_dir/sessionInfo/maaslin3_env_R.txt" 2>/dev/null
fi

queue=$(mktemp)
awk -F'\t' -v md="$run_dir/markers" 'NR>1 {
  for (r = 1; r <= $9; r++) {
    m = sprintf("%s/%s_rep%04d.done", md, $2, r)
    if ((getline _ < m) < 0) printf "%s %d\n", $2, r
    close(m)
  }
}' "$run_dir/cells.tsv" | sort -s -k2,2n > "$queue"
# round-robin by rep index: an interrupted night still leaves every cell with a
# usable (smaller) rep count instead of finishing some cells and starving others

total=$(awk -F'\t' 'NR>1 {s += $9} END {print s}' "$run_dir/cells.tsv")
todo=$(wc -l < "$queue" | tr -d ' ')
echo "run $run_id: $todo of $total reps to do (jobs=$jobs, timeout=${TIMEOUT_SEC:-14400}s, pseudo_sens=$PSEUDO_SENS)"
if [[ "${DRY:-0}" == "1" || "$todo" == "0" ]]; then rm -f "$queue"; exit 0; fi

start_ts=$(date +%s)
caffeinate -s -i xargs -P "$jobs" -n 2 \
  gtimeout "${TIMEOUT_SEC:-14400}" "$harness/run_rep_real.sh" "$run_dir" \
  < "$queue" || true
rm -f "$queue"
echo "wall: $(( $(date +%s) - start_ts )) s"

done_n=$(ls "$run_dir/markers" 2>/dev/null | wc -l | tr -d ' ')
echo "run $run_id finished: $done_n / $total reps complete"
