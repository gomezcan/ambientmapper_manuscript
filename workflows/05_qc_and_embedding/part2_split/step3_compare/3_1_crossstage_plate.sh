#!/bin/bash
# 3_1_crossstage_plate.sh -- run the plate cross-stage compare (Pre vs wd vs nd) per genome.
# Lightweight (base R + optional ggplot2) -- PDFs need ggplot2; without it only the TSVs are written.
# Resolves the single frozen v7 metadata per sample by glob, so it tracks the frozen config tag
# automatically (At carries .mclst_40, B73 does not). Skips a genome whose Pre object is not built yet.
#
#   bash 0_scripts/part2_split/step3_compare/3_1_crossstage_plate.sh
set -uo pipefail
PROJECT_ROOT="${PROJECT_ROOT:-.}"
BASE="${PROJECT_ROOT}/6_socrates"
cd "$BASE"
SCRIPTS="${SCRIPTS:-${BASE}/0_scripts}"; S2="SM2v2_plate/step2_cluster"; OUT="SM2v2_plate/step3_compare"

# conda activate ambientmapper-manuscript   (environment.yml at the repo root; needed for the ggplot2 PDFs)

run_one() {   # genome  pre_prefix  wd_prefix  nd_prefix
  local g="$1" pre="$2" wd="$3" nd="$4"
  local pf wf nf
  pf=$(ls "$S2/${pre}.mQCv6.updated_metadata_v7."*.txt 2>/dev/null | head -1)
  wf=$(ls "$S2/${wd}.mQCv6.updated_metadata_v7."*.txt  2>/dev/null | head -1)
  nf=$(ls "$S2/${nd}.mQCv6.updated_metadata_v7."*.txt  2>/dev/null | head -1)
  if [ -z "$pf" ]; then echo "SKIP $g: Pre v7 metadata not found (clustering not run yet)"; return; fi
  [ -z "$wf" ] && wf=NA; [ -z "$nf" ] && nf=NA
  echo "== $g ==  Pre=$(basename "$pf")  wd=$([ "$wf" = NA ] && echo NA || basename "$wf")  nd=$([ "$nf" = NA ] && echo NA || basename "$nf")"
  Rscript "${SCRIPTS}/part2_split/step3_compare/3_1_crossstage_plate.R" "$g" "$OUT" "$pf" "$wf" "$nf"
}

run_one TAIR10 SM2_At_TAIR10  Clean.SM2v2wd_At_TAIR10   Clean.SM2v2_At_TAIR10
run_one B73v5  SM2_B73_B73v5  Clean.SM2v2wd_B73_B73v5   Clean.SM2v2_B73_B73v5
echo "done. Outputs -> $OUT (TSVs) + $OUT/plots (PDFs if ggplot2 available)."
