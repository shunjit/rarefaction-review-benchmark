#!/bin/zsh
# One-way sync of local hot run directories to the canonical PS3000 location.
# Single sequential writer -> safe for the exFAT/fskit volume.
# usage: REVISION_ROOT=... ./sync_runs.sh [run_id ...]   (default: all local runs)
set -euo pipefail
[[ -n "${REVISION_ROOT:-}" && -d "$REVISION_ROOT" ]] || { echo "ERROR: REVISION_ROOT MISSING" >&2; exit 1; }
runs_root="${RUNS_ROOT:-$HOME/dsim_runs}"
canon_root="$REVISION_ROOT/simulation/runs"
if (( $# > 0 )); then ids=("$@"); else ids=($(ls "$runs_root" 2>/dev/null)); fi
for id in "${ids[@]}"; do
  [[ -d "$runs_root/$id" ]] || { echo "skip $id (no local dir)"; continue; }
  mkdir -p "$canon_root/$id"
  rsync -a "$runs_root/$id/" "$canon_root/$id/"
  echo "synced $id -> $canon_root/$id"
done
