#!/bin/bash
########## BATCH Lines for Resource Request ##########
#SBATCH --time=2:00:00
#SBATCH --nodes=1
#SBATCH --ntasks-per-node=1
#SBATCH --cpus-per-task=4
#SBATCH --mem=64G
#SBATCH --job-name=Step5_5_metacell_rZ
#SBATCH --partition=standard
#SBATCH --output=_logs/Step5_5_metacell_rZ_%A_%a.log
#SBATCH --array=0-5
#
# 5_5_metacell_rZ_plate.sh -- Step 5.5 (plate arm): annotate SEACells meta-cells with the
# 4_3i reciprocal-z idea (engine: 0_scripts/common/4_3p_metacell_rZ.R).
#
#   0: SM2_At_TAIR10            Pre    3: SM2_B73_B73v5            Pre
#   1: Clean.SM2v2wd_At_TAIR10  wd     4: Clean.SM2v2wd_B73_B73v5  wd
#   2: Clean.SM2v2_At_TAIR10    nd     5: Clean.SM2v2_B73_B73v5    nd
#
# Light enough to run locally (sparse perkb is 10-70 MB; the 4_3f-smoothed matrices, which
# this deliberately does NOT use, are 152 MB / 3.7 GB).
#
# MODE=pooled   all seeds' meta-cells stacked as columns (At N=70, Zi cap 8.25)   [default]
# MODE=perseed  one run per seed (At N=14, Zi cap 3.47); the across-seed spread is the
#               uncertainty, and it is the only seed-aware quantity that is NOT
#               pseudo-replicated
# MODE=both     both of the above
#
# POOLED SEEDS ARE PSEUDO-REPLICATES. The 5 At seeds re-partition the SAME 1061 cells, so
# 70 meta-cells are five overlapping views, not 70 independent units. Pooling is legitimate
# for raising the Zi cap and stabilising the row statistics; it is NOT a sample size. Never
# hand the pooled N to a null model (same trap as the per-cell tests).
# ---------------------------------------------------------------------------

set -euo pipefail

PROJECT_ROOT="${PROJECT_ROOT:-.}"
BASE="${PROJECT_ROOT}/6_socrates"
S3="${BASE}/SM2v2_plate/step3_compare"
S5="${BASE}/SM2v2_plate/step5_metacell"
MK="${BASE}/_data/markers"
ENGINE="${SCRIPTS:-${BASE}/0_scripts}/common/4_3p_metacell_rZ.R"
OUT="${S5}/rZ_annotation"
MODE="${MODE:-pooled}"
TOPN="${TOPN:-6}"

declare -a NAMES=( SM2_At_TAIR10 Clean.SM2v2wd_At_TAIR10 Clean.SM2v2_At_TAIR10 \
                   SM2_B73_B73v5 Clean.SM2v2wd_B73_B73v5 Clean.SM2v2_B73_B73v5 )
declare -a STAGES=( Pre wd nd Pre wd nd )
declare -a PANELS=( At At At maize maize maize )
# the two Step-5.3 pass-2 deliverables (count + bandwidth + depth matched)
RUN_At="seacells_cps50_F14235_N14"      # 5 seeds -> seed<k>/ subdirs
RUN_maize="seacells_cps50_F44658_N286"  # 1 seed  -> flat

mkdir -p "$OUT"

run_one () {                 # $1 index  $2 tag  $3 comma-separated c2s list
  local i="$1" tag="$2" c2s="$3"
  local name="${NAMES[$i]}" stage="${STAGES[$i]}" panel="${PANELS[$i]}"
  local mat="${S3}/${name}.plate.perkb.genes.sparse.rds"     # RAW perkb -- never .smoothed
  local bed="${MK}/markers.${panel}.informative_top15.bed"
  [[ -s "$mat" ]] || { echo "ERROR: missing $mat"; return 3; }
  [[ -s "$bed" ]] || { echo "ERROR: missing $bed"; return 3; }
  echo "--- [${i}] ${name} (${stage}) ${tag}"
  Rscript "$ENGINE" "${OUT}/${tag}" "$name" "$mat" "$c2s" "$bed" "$stage" "$TOPN"
}

collect_c2s () {             # $1 index -> echoes comma-separated assignment paths
  local i="$1" name="${NAMES[$i]}" run f
  if [[ "${PANELS[$i]}" == "At" ]]; then run="$RUN_At"; else run="$RUN_maize"; fi
  local acc=""
  for d in "${S5}/${run}"/seed*/ "${S5}/${run}"/; do
    f="${d}${name}.seacells.cell_to_seacell.tsv"
    [[ -s "$f" ]] && acc="${acc:+${acc},}${f}"
  done
  echo "$acc"
}

INDICES=( "${@:-0 1 2 3 4 5}" )
if [[ -n "${SLURM_ARRAY_TASK_ID:-}" ]]; then INDICES=( "$SLURM_ARRAY_TASK_ID" ); fi
read -r -a INDICES <<< "${INDICES[*]}"

for i in "${INDICES[@]}"; do
  all="$(collect_c2s "$i")"
  [[ -n "$all" ]] || { echo "ERROR: no cell_to_seacell.tsv found for ${NAMES[$i]}"; exit 4; }
  if [[ "$MODE" == "pooled" || "$MODE" == "both" ]]; then
    run_one "$i" "pooled" "$all"
  fi
  if [[ "$MODE" == "perseed" || "$MODE" == "both" ]]; then
    IFS=',' read -r -a arr <<< "$all"
    for f in "${arr[@]}"; do
      sd="$(sed -E 's|.*/seed([0-9]+)/.*|\1|' <<< "$f")"; [[ "$sd" == "$f" ]] && sd=1
      run_one "$i" "seed${sd}" "$f"
    done
  fi
done
echo "$(date): done -> ${OUT}"
