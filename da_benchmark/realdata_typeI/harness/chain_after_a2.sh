#!/bin/zsh
# A-2 -> A-1 chain (2026-08-16): wait for dsim_a2_nsens_v1 (200 reps) to
# finish, aggregate it, sync the summary to the canonical location, then
# launch the A-1 night run (dreal_a1_x25_v1, R=100, B=0).
# Mirrors chain_after_dsim.sh. Start from a LOCAL CWD; absolute paths only.
# usage: REVISION_ROOT=... nohup ./chain_after_a2.sh &
set -euo pipefail
[[ -n "${REVISION_ROOT:-}" ]] || { echo "ERROR: REVISION_ROOT MISSING" >&2; exit 1; }
sim_run="$HOME/dsim_runs/dsim_a2_nsens_v1"
total=200
while :; do
  n=$(ls "$sim_run/markers" 2>/dev/null | wc -l | tr -d ' ')
  if [[ "$n" -ge "$total" ]] && ! pgrep -f "orchestrate.sh dsim_a2_nsens_v1" > /dev/null; then
    break
  fi
  sleep 300
done
echo "[chain] A-2 complete ($n/$total reps) at $(date)"

"$HOME/miniforge3/envs/maaslin3_env/bin/Rscript" \
  "$REVISION_ROOT/simulation/harness/R/aggregate_run.R" "$sim_run" \
  > "$sim_run/log/aggregate_$(date +%Y%m%d_%H%M%S).log" 2>&1
echo "[chain] A-2 aggregated"

REVISION_ROOT="$REVISION_ROOT" zsh "$REVISION_ROOT/simulation/harness/sync_runs.sh" dsim_a2_nsens_v1
echo "[chain] A-2 synced to canonical"

REVISION_ROOT="$REVISION_ROOT" zsh "$REVISION_ROOT/realdata_typeI/harness/night_launch_real.sh" \
  dreal_a1_x25_v1 "$REVISION_ROOT/realdata_typeI/harness/configs/cells_a1_x25_v1.tsv" 8
echo "[chain] A-1 night run launched at $(date)"
