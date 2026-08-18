#!/bin/zsh
# Detached overnight launch for D-real (mirrors simulation night_launch.sh).
# Hot run state on local APFS (RUNS_ROOT, default ~/dreal_runs); on completion
# the run directory is rsynced once to the canonical location
# PS3000 revision_r1/realdata_typeI/runs/<run_id>/.
# usage: REVISION_ROOT=... ./night_launch_real.sh <run_id> <cells_tsv> [jobs]
set -euo pipefail
[[ -n "${REVISION_ROOT:-}" ]] || { echo "ERROR: REVISION_ROOT MISSING" >&2; exit 1; }
run_id=$1; cells=$2; jobs=${3:-8}
here=$(cd "$(dirname "$0")" && pwd)
runs_root="${RUNS_ROOT:-$HOME/dreal_runs}"
run_dir="$runs_root/$run_id"
canon="$REVISION_ROOT/realdata_typeI/runs/$run_id"
mkdir -p "$run_dir/log"
out="$run_dir/log/night_$(date +%Y%m%d_%H%M%S).log"
nohup /bin/zsh -c "
  export REVISION_ROOT='$REVISION_ROOT' RUNS_ROOT='$runs_root'
  export TIMEOUT_SEC='${TIMEOUT_SEC:-14400}' PSEUDO_SENS='${PSEUDO_SENS:-TRUE}' BASE_SEED='${BASE_SEED:-20260813}'
  '$here/orchestrate_real.sh' '$run_id' '$cells' '$jobs'
  mkdir -p '$canon'
  rsync -a '$run_dir/' '$canon/' && echo 'SYNCED to PS3000 canonical'
" > "$out" 2>&1 &
pid=$!
disown
echo "launched D-real run $run_id (pid $pid, jobs $jobs)"
echo "local run dir: $run_dir"
echo "log:           $out"
echo "progress:      ls $run_dir/markers | wc -l"
