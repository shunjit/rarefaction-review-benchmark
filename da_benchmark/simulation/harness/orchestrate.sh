#!/bin/zsh
# D-sim orchestrator: build the (cell x rep) queue, skip completed markers,
# run under xargs -P + caffeinate. Resumable at rep granularity.
# (xargs -P replaces GNU parallel, which hangs in this execution environment;
#  one dependency fewer, identical semantics for an embarrassingly parallel queue.)
#
# usage: REVISION_ROOT=... ./orchestrate.sh <run_id> <cells_tsv> [jobs]
#   env overrides: BASE_SEED (default 20260811), PSEUDO_SENS (default TRUE),
#                  TIMEOUT_SEC (default 3600), DRY=1 (print queue size and exit)
set -euo pipefail

[[ -n "${REVISION_ROOT:-}" && -d "${REVISION_ROOT:-}" ]] || { echo "ERROR: REVISION_ROOT MISSING - stopping." >&2; exit 1; }
run_id=$1; cells_src=$2; jobs=${3:-8}
[[ -f "$cells_src" ]] || { echo "ERROR: cells file $cells_src MISSING" >&2; exit 1; }

sim="$REVISION_ROOT/simulation"
harness="$sim/harness"
# Hot run state lives on LOCAL APFS: the PS3000 exFAT (fskit) volume loses
# recently-created files under 8-way concurrent small-file writes (observed
# 2026-08-11: markers/results vanished, directory rolled back to launch time).
# The canonical copy on PS3000 is produced by a single-writer rsync afterwards
# (night_launch.sh / sync_runs.sh).
runs_root="${RUNS_ROOT:-$HOME/dsim_runs}"
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

if [[ ! -f "$run_dir/run.env" ]]; then
  cat > "$run_dir/run.env" <<EOF
BASE_SEED=${BASE_SEED:-20260811}
PSEUDO_SENS=${PSEUDO_SENS:-TRUE}
REVISION_ROOT="$REVISION_ROOT"
HARNESS_R_DIR="$harness/R"
ANCOMBC_RSCRIPT="$ancombc_rscript"
MAASLIN3_RSCRIPT="$maaslin3_rscript"
WORKROOT="\$HOME/.cache/dsim_work/$run_id"
EOF
fi
source "$run_dir/run.env"
mkdir -p "$HOME/.cache/dsim_work/$run_id"

if [[ ! -f "$run_dir/config.yaml" ]]; then
  cat > "$run_dir/config.yaml" <<EOF
run_id: $run_id
created_utc: $(date -u +%Y-%m-%dT%H:%M:%SZ)
design: D_HARNESS_DESIGN.md v2 (2026-08-11)
base_seed: $BASE_SEED
pseudo_sens: $PSEUDO_SENS
cells: cells.tsv (frozen copy in this directory)
seed_hierarchy: "dataset = base + cell_index*1e7 + rep*1e3; +901 single rarefaction; +1..B boots; +903 method"
prefilter: "prevalence >= 0.10 on the raw table; intention-to-test set shared by all arms"
environments:
  ancombc_rscript: $ancombc_rscript
  maaslin3_rscript: $maaslin3_rscript
host: $(hostname -s)
EOF
  (cd "$harness" && shasum -a 256 R/*.R run_rep.sh orchestrate.sh) > "$run_dir/sessionInfo/script_checksums.sha256"
  "$ancombc_rscript"  -e 'suppressMessages(library(ANCOMBC)); writeLines(capture.output(sessionInfo()))' > "$run_dir/sessionInfo/ancombc_env_R.txt" 2>/dev/null
  "$maaslin3_rscript" -e 'suppressMessages(library(maaslin3)); writeLines(capture.output(sessionInfo()))' > "$run_dir/sessionInfo/maaslin3_env_R.txt" 2>/dev/null
fi

queue=$(mktemp)
awk -F'\t' -v md="$run_dir/markers" 'NR>1 {
  for (r = 1; r <= $12; r++) {
    m = sprintf("%s/%s_rep%04d.done", md, $2, r)
    if ((getline _ < m) < 0) printf "%s %d\n", $2, r
    close(m)
  }
}' "$run_dir/cells.tsv" | sort -s -k2,2n > "$queue"
# round-robin by rep index: an interrupted night still leaves every cell with a
# usable (smaller) rep count instead of finishing some cells and starving others

total=$(awk -F'\t' 'NR>1 {s += $12} END {print s}' "$run_dir/cells.tsv")
todo=$(wc -l < "$queue" | tr -d ' ')
echo "run $run_id: $todo of $total reps to do (jobs=$jobs, timeout=${TIMEOUT_SEC:-3600}s, pseudo_sens=$PSEUDO_SENS)"
if [[ "${DRY:-0}" == "1" || "$todo" == "0" ]]; then rm -f "$queue"; exit 0; fi

start_ts=$(date +%s)
caffeinate -s -i xargs -P "$jobs" -n 2 \
  gtimeout "${TIMEOUT_SEC:-3600}" "$harness/run_rep.sh" "$run_dir" \
  < "$queue" || true
rm -f "$queue"
echo "wall: $(( $(date +%s) - start_ts )) s"

done_n=$(ls "$run_dir/markers" 2>/dev/null | wc -l | tr -d ' ')
echo "run $run_id finished: $done_n / $total reps complete"
