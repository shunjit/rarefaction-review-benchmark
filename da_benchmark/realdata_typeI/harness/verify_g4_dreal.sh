#!/bin/bash
# G4 reproducibility check (D_HARNESS_DESIGN.md v2, section 14).
#
# Phase 2 (dreal_boots_v1, B=10, R=25) shares BASE_SEED=20260813 with phase 1
# (dreal_main_v1, B=0, R=100) and recomputes the raw/rarefied rows for reps
# 1-25. If the harness is deterministic, those recomputed rows must be
# byte-identical to phase 1 — which licenses pairing the phase-2 boots rows
# with the phase-1 raw/rarefied distributions.
#
# Per (cell, rep) pair, compares md5 of:
#   1. result rows with condition raw/rarefied (header kept, sorted to
#      neutralise write order)
#   2. the *_truth.tsv.gz file (identical subsample => identical tested set)
#   3. geninfo minus column 25 ("B", the only intended difference: 10 vs 0)
#
# usage: ./verify_g4_dreal.sh <main_run_dir> <boots_run_dir> [n_reps]
set -u
MAIN=${1:?main run dir}
BOOTS=${2:?boots run dir}
NREP=${3:-25}
PASS=0; FAIL=0
for cell in rice_strat rice_rhizo; do
  for i in $(seq 1 "$NREP"); do
    rep=$(printf "rep_%04d" "$i")
    m="$MAIN/results/$cell"; b="$BOOTS/results/$cell"
    h1=$(gzcat "$m/$rep.tsv.gz" | awk -F'\t' 'NR==1 || $5=="raw" || $5=="rarefied"' | sort | md5 -q)
    h2=$(gzcat "$b/$rep.tsv.gz" | awk -F'\t' 'NR==1 || $5=="raw" || $5=="rarefied"' | sort | md5 -q)
    t1=$(gzcat "$m/${rep}_truth.tsv.gz" | md5 -q)
    t2=$(gzcat "$b/${rep}_truth.tsv.gz" | md5 -q)
    g1=$(cut -f1-24,26-28 "$m/${rep}_geninfo.tsv" | md5 -q)
    g2=$(cut -f1-24,26-28 "$b/${rep}_geninfo.tsv" | md5 -q)
    if [ "$h1" = "$h2" ] && [ "$t1" = "$t2" ] && [ "$g1" = "$g2" ]; then
      PASS=$((PASS+1))
    else
      FAIL=$((FAIL+1))
      echo "MISMATCH $cell $rep: results=$([ "$h1" = "$h2" ] && echo ok || echo NG)" \
           "truth=$([ "$t1" = "$t2" ] && echo ok || echo NG)" \
           "geninfo=$([ "$g1" = "$g2" ] && echo ok || echo NG)"
    fi
  done
done
echo "G4 VERIFY: PASS=$PASS FAIL=$FAIL (of $((2 * NREP)) cell-rep pairs)"
[ "$FAIL" -eq 0 ]
