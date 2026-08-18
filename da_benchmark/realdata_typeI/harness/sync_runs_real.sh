#!/bin/zsh
# Manual one-way sync of local D-real run directories to the PS3000 canonical
# location (single writer). usage: REVISION_ROOT=... ./sync_runs_real.sh [run_id]
set -euo pipefail
[[ -n "${REVISION_ROOT:-}" && -d "${REVISION_ROOT:-}" ]] || { echo "ERROR: REVISION_ROOT MISSING" >&2; exit 1; }
runs_root="${RUNS_ROOT:-$HOME/dreal_runs}"
canon_root="$REVISION_ROOT/realdata_typeI/runs"
mkdir -p "$canon_root"
if [[ $# -ge 1 ]]; then
  dirs=("$runs_root/$1")
else
  dirs=("$runs_root"/*(N/))
fi
for d in "${dirs[@]}"; do
  [[ -d "$d" ]] || { echo "skip: $d not a directory"; continue; }
  rid=$(basename "$d")
  echo "sync $rid ..."
  rsync -a "$d/" "$canon_root/$rid/"
done
echo "done"
