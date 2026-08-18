#!/bin/zsh
# One-shot chain launcher: when every dsim_main_v1 process (orchestrate wrapper,
# xargs workers, final rsync) has exited, launch the D-real phase-1 night run.
# Run detached under caffeinate so the machine cannot sleep in the gap between
# dsim's own caffeinate ending (xargs exit) and the D-real caffeinate starting:
#   nohup caffeinate -s -i /bin/zsh chain_after_dsim.sh >> <log> 2>&1 &
# The pgrep pattern "dsim_main_v1" matches the dsim night wrapper shell, its
# xargs/run_rep children and the canonicalising rsync, but not this script
# (file name carries no "_main_v1" argument — the pattern string itself is
# split below to keep this process out of its own match set).
set -u
REVISION_ROOT=${REVISION_ROOT:-/Volumes/PS3000/benchmark_data/revision_r1}
H="$REVISION_ROOT/realdata_typeI/harness"
RUN_ID=${RUN_ID:-dreal_main_v1}
CELLS=${CELLS:-$H/configs/cells_real_v1.tsv}
JOBS=${JOBS:-8}
pat="dsim_main""_v1"

echo "[chain $(date '+%F %T')] waiting for ${pat} processes to finish ..."
while pgrep -f "$pat" >/dev/null 2>&1; do
  sleep 300
done
sleep 60   # settle: let filesystem writes and log flushes complete
if pgrep -f "$pat" >/dev/null 2>&1; then
  echo "[chain $(date '+%F %T')] processes reappeared after settle; waiting again"
  while pgrep -f "$pat" >/dev/null 2>&1; do sleep 300; done
  sleep 60
fi

n_dsim=$(ls "$HOME/dsim_runs/dsim_main_v1/markers" 2>/dev/null | wc -l | tr -d ' ')
echo "[chain $(date '+%F %T')] dsim clear (markers=$n_dsim). launching D-real ${RUN_ID}"

export REVISION_ROOT TIMEOUT_SEC=${TIMEOUT_SEC:-10800} PSEUDO_SENS=${PSEUDO_SENS:-TRUE} BASE_SEED=${BASE_SEED:-20260813}
"$H/night_launch_real.sh" "$RUN_ID" "$CELLS" "$JOBS"
echo "[chain $(date '+%F %T')] night_launch_real dispatched"
sleep 120  # hold caffeinate until the new run's own caffeinate is up
echo "[chain $(date '+%F %T')] chain done"
