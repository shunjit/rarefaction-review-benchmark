#!/bin/bash
# =============================================================================
# A-3 (CG-19): Rice 共通検体集合 (n=486) 固定での 4 深度 alpha/beta 再計算
# run: a3_common486_v1
#
# フェーズ1 感度解析 (rarefaction-review-benchmark/script/06_sensitivity_analysis.sh
# Rice セクション) を、4 深度 (10000/15000/20000/22714) すべてで保持される共通検体
# 集合 (基盤 A 台帳 rice_sample_ledger.tsv の retained_d* 4 列 = True、n=486) に
# 入力 table を固定して再実行する。フェーズ1 との差分は「入力 table の検体集合のみ」
# (metadata・深度・コマンド・環境系列は同一)。目的 = R3-9 の「同一検体集合上の
# stability」を、脱落ゼロの固定集合で直接実証する (dropout 影響の分離)。
#
# 規約: PHASE3_PLAN §5 (ルート変数・MISSING 停止)、夜間ラン規約 (ローカル CWD・
# marker 再開・完了時 rsync 正本化)。
# 起動: cd ~/a3_runs/a3_common486_v1 && nohup bash script/run_a3_common486.sh \
#         > log/run.log 2>&1 &
# =============================================================================
set -euo pipefail

RAW_ROOT="${RAW_ROOT:-/Volumes/PS3000/benchmark_data}"
REVISION_ROOT="${REVISION_ROOT:-/Volumes/PS3000/benchmark_data/revision_r1}"
LOCAL_RUN="${LOCAL_RUN:-$HOME/a3_runs/a3_common486_v1}"
CONDA_SH="${CONDA_SH:-$HOME/miniforge3/etc/profile.d/conda.sh}"
QIIME_ENV="${QIIME_ENV:-qiime2-amplicon-2025.10}"
THREADS="${THREADS:-8}"
DEPTHS=(22714 20000 15000 10000)   # 最重深度を先頭に置きスモークを兼ねる
N_EXPECTED=486

SRC_TABLE="$RAW_ROOT/Rice_v2/02_denoised/table.qza"
SRC_TREE="$RAW_ROOT/Rice_v2/02_denoised/rooted-tree.qza"
SRC_META="$RAW_ROOT/00_github/Schloss_Rarefaction_mSphere_2024/data/rice/metadata_rice.tsv"
SRC_LEDGER="$REVISION_ROOT/sample_ledger/output/rice_sample_ledger.tsv"
DEST="$REVISION_ROOT/common_set/runs/a3_common486_v1"

log() { printf '[%s] %s\n' "$(date '+%Y-%m-%d %H:%M:%S')" "$*"; }
die() { log "ERROR: $*"; exit 1; }

log "=== A-3 a3_common486_v1 start (threads=$THREADS) ==="

# --- Phase 0: 入力チェックと環境 -------------------------------------------
for f in "$SRC_TABLE" "$SRC_TREE" "$SRC_META" "$SRC_LEDGER"; do
  [ -f "$f" ] || die "MISSING input: $f"
done
[ -f "$CONDA_SH" ] || die "MISSING conda.sh: $CONDA_SH"
# conda の activate フックは未定義変数を参照するため、activate の間だけ -u を外す
set +u
# shellcheck disable=SC1090
source "$CONDA_SH"
conda activate "$QIIME_ENV" || die "conda env $QIIME_ENV not found"
set -u
command -v qiime >/dev/null || die "qiime CLI not on PATH after activate"

cd "$LOCAL_RUN" || die "LOCAL_RUN missing: $LOCAL_RUN"
mkdir -p input output summary sessionInfo markers log

# 入力をローカルへ確定コピー (PS3000 読みは起動時のみ)
for f in "$SRC_TABLE" "$SRC_TREE" "$SRC_META" "$SRC_LEDGER"; do
  cp -f "$f" input/
done
( cd input && shasum -a 256 table.qza rooted-tree.qza metadata_rice.tsv \
    rice_sample_ledger.tsv > input_checksums.sha256 )
log "inputs copied and checksummed"

# --- Phase 1: 共通検体 486 の抽出と検証 ------------------------------------
# 台帳列: 8=retained_d10000 9=retained_d15000 10=retained_d20000 11=retained_d22714
awk -F'\t' 'NR==1{print "sample-id"; next}
  $8=="True" && $9=="True" && $10=="True" && $11=="True" {print $1}' \
  input/rice_sample_ledger.tsv > input/ids_common486.tsv
N_IDS=$(( $(wc -l < input/ids_common486.tsv) - 1 ))
[ "$N_IDS" -eq "$N_EXPECTED" ] || die "common set size $N_IDS != $N_EXPECTED"

# 全 486 検体の table_depth (列7) が最大深度以上であること (rarefy で脱落しない)
MIN_DEPTH=$(awk -F'\t' 'NR==FNR{if(FNR>1) keep[$1]=1; next}
  keep[$1]{print int($7)}' input/ids_common486.tsv input/rice_sample_ledger.tsv | sort -n | head -1)
[ "$MIN_DEPTH" -ge 22714 ] || die "min table_depth $MIN_DEPTH < 22714 in common set"
log "common set OK: n=$N_IDS, min table_depth=$MIN_DEPTH"

# --- Phase 2: table を 486 検体へ固定 ---------------------------------------
if [ ! -f markers/filter.done ]; then
  qiime feature-table filter-samples \
    --i-table input/table.qza \
    --m-metadata-file input/ids_common486.tsv \
    --o-filtered-table input/table_common486.qza
  python - <<'PY' || die "filtered table sample count != 486"
from qiime2 import Artifact
import biom, sys
t = Artifact.load("input/table_common486.qza").view(biom.Table)
n = len(t.ids("sample"))
print(f"filtered table samples = {n}", flush=True)
sys.exit(0 if n == 486 else 1)
PY
  touch markers/filter.done
fi
log "filtered table ready (n=486)"

# --- Phase 3: 深度ごとの core-metrics + PERMANOVA + alpha 群有意性 ----------
for D in "${DEPTHS[@]}"; do
  MARK="markers/depth_${D}.done"
  if [ -f "$MARK" ]; then log "depth $D: marker present, skip"; continue; fi
  T0=$SECONDS
  OUTDIR="output/depth_${D}"
  rm -rf "$OUTDIR"; mkdir -p "$OUTDIR"
  log "depth $D: core-metrics-phylogenetic start"

  qiime diversity core-metrics-phylogenetic \
    --i-table input/table_common486.qza \
    --i-phylogeny input/rooted-tree.qza \
    --p-sampling-depth "$D" \
    --m-metadata-file input/metadata_rice.tsv \
    --p-n-jobs-or-threads "$THREADS" \
    --output-dir "$OUTDIR/core_metrics"

  python - "$OUTDIR/core_metrics/rarefied_table.qza" <<'PY' || die "rarefied table != 486 samples"
from qiime2 import Artifact
import biom, sys
t = Artifact.load(sys.argv[1]).view(biom.Table)
n = len(t.ids("sample"))
print(f"rarefied table samples = {n}", flush=True)
sys.exit(0 if n == 486 else 1)
PY

  qiime diversity beta-group-significance \
    --i-distance-matrix "$OUTDIR/core_metrics/weighted_unifrac_distance_matrix.qza" \
    --m-metadata-file input/metadata_rice.tsv \
    --m-metadata-column compartment \
    --p-method permanova \
    --p-pairwise \
    --p-permutations 999 \
    --o-visualization "$OUTDIR/permanova_wunifrac_compartment.qzv"

  qiime diversity beta-group-significance \
    --i-distance-matrix "$OUTDIR/core_metrics/unweighted_unifrac_distance_matrix.qza" \
    --m-metadata-file input/metadata_rice.tsv \
    --m-metadata-column compartment \
    --p-method permanova \
    --p-pairwise \
    --p-permutations 999 \
    --o-visualization "$OUTDIR/permanova_uunifrac_compartment.qzv"

  qiime diversity alpha-group-significance \
    --i-alpha-diversity "$OUTDIR/core_metrics/shannon_vector.qza" \
    --m-metadata-file input/metadata_rice.tsv \
    --o-visualization "$OUTDIR/shannon_group_significance.qzv"

  qiime diversity alpha-group-significance \
    --i-alpha-diversity "$OUTDIR/core_metrics/faith_pd_vector.qza" \
    --m-metadata-file input/metadata_rice.tsv \
    --o-visualization "$OUTDIR/faith_pd_group_significance.qzv"

  touch "$MARK"
  log "depth $D: done in $((SECONDS - T0)) s"
done

# --- Phase 4: PERMANOVA/alpha 要約の抽出 ------------------------------------
python script/extract_a3_permanova.py output summary
log "summary extracted"

# --- Phase 5: 環境記録 -------------------------------------------------------
qiime info > sessionInfo/qiime_info.txt 2>&1
conda list > sessionInfo/conda_list.txt 2>&1
( cd script && shasum -a 256 run_a3_common486.sh extract_a3_permanova.py \
    > ../sessionInfo/script_checksums.sha256 )

# --- Phase 6: rsync 正本化 ---------------------------------------------------
mkdir -p "$DEST"
rsync -rt --exclude '.DS_Store' \
  script input/ids_common486.tsv input/input_checksums.sha256 \
  input/table_common486.qza output summary sessionInfo markers log "$DEST/" \
  || die "rsync to PS3000 failed"
log "rsynced to $DEST"
log "=== A-3 a3_common486_v1 COMPLETE ==="
