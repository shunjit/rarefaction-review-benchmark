#!/usr/bin/env bash
# run_with_base_dir.sh — run a frozen pipeline script with the data root overridden.
#
# The numbered pipeline scripts in this directory are archived exactly as they
# were executed for the manuscript (frozen; see README "Reproducibility").
# They therefore carry the absolute data root of the original analysis machine:
#
#   BASE_DIR="/Volumes/PS3000/benchmark_data"      (02_ ... 07_)
#   OUTPUT_BASE="/Volumes/PS3000/benchmark_data"   (01_download_data.sh)
#
# To keep those files byte-identical, this wrapper does NOT edit them in place.
# It mirrors the whole script/ directory into a temporary location (the
# pipeline scripts locate companion files, e.g. the Step 0 Python diagnostic,
# relative to their own path), replaces exactly one assignment line in the
# requested script, verifies that every other mirrored file is byte-identical
# to the frozen original, and then executes the patched copy.
#
# Usage:
#   BASE_DIR=/path/to/your/data_root bash script/run_with_base_dir.sh 04_confounding_diagnosis.sh
#
# Full pipeline (the frozen run_all.sh calls the numbered scripts directly and
# would bypass this wrapper, so loop over the steps instead):
#   for s in 01_download_data.sh 02_dada2_processing.sh 03_build_phylogeny.sh \
#            04_confounding_diagnosis.sh 05_alpha_rarefaction.sh \
#            06_sensitivity_analysis.sh 07_q2boots_analysis.sh; do
#     BASE_DIR=/path/to/your/data_root bash script/run_with_base_dir.sh "$s" || break
#   done
#
# Notes:
#   - BASE_DIR must not contain a double-quote character.
#   - Requires the QIIME 2 environment described in the repository README.

set -euo pipefail

HERE=$(cd "$(dirname "$0")" && pwd)
SELF=$(basename "$0")

if [[ $# -lt 1 ]]; then
  echo "usage: BASE_DIR=/path/to/data_root bash $0 <script-name> [args...]" >&2
  exit 2
fi
if [[ -z "${BASE_DIR:-}" ]]; then
  echo "ERROR: set the BASE_DIR environment variable to your data root" >&2
  exit 2
fi
case "$BASE_DIR" in
  *\"*) echo "ERROR: BASE_DIR must not contain double quotes" >&2; exit 2 ;;
esac

name=$1; shift
src="$HERE/$name"
if [[ ! -f "$src" ]]; then
  echo "ERROR: $src not found" >&2
  exit 2
fi

tmp=$(mktemp -d)
trap 'rm -rf "$tmp"' EXIT

# Mirror the frozen script directory (scripts resolve siblings via dirname "$0").
( cd "$HERE" && tar -cf - --exclude='._*' --exclude='.DS_Store' . ) | ( cd "$tmp" && tar -xf - )

# Replace exactly one literal assignment line in the requested script; fail otherwise.
if ! awk -v nb="$BASE_DIR" '
  $0 == "BASE_DIR=\"/Volumes/PS3000/benchmark_data\""    { print "BASE_DIR=\"" nb "\"";    c++; next }
  $0 == "OUTPUT_BASE=\"/Volumes/PS3000/benchmark_data\"" { print "OUTPUT_BASE=\"" nb "\""; c++; next }
  { print }
  END { if (c != 1) exit 3 }
' "$src" > "$tmp/$name.patched"; then
  echo "ERROR: expected exactly one BASE_DIR/OUTPUT_BASE assignment line in $name;" >&2
  echo "       refusing to run a copy that differs from the frozen script elsewhere." >&2
  exit 3
fi
mv "$tmp/$name.patched" "$tmp/$name"

# Integrity check: every mirrored file except the patched one and this wrapper
# must be byte-identical to the frozen original.
while IFS= read -r f; do
  rel=${f#"$tmp"/}
  [[ "$rel" == "$name" || "$rel" == "$SELF" ]] && continue
  if ! cmp -s "$f" "$HERE/$rel"; then
    echo "ERROR: mirrored copy of $rel differs from the frozen original; aborting" >&2
    exit 3
  fi
done < <(find "$tmp" -type f)

echo "[run_with_base_dir] $name with BASE_DIR=$BASE_DIR (frozen scripts untouched; one line patched in a temporary mirror)" >&2
exec bash "$tmp/$name" "$@"
