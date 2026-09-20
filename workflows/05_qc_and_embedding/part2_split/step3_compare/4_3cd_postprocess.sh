#!/bin/bash
########## BATCH Lines for Resource Request ##########
#SBATCH --time=0:30:00
#SBATCH --nodes=1
#SBATCH --ntasks-per-node=1
#SBATCH --cpus-per-task=1
#SBATCH --mem=8G
#SBATCH --job-name=Step4_diff_SM2v2_indep
#SBATCH --partition=standard
#SBATCH --output=_logs/Step4_diff_SM2v2_indep_%A.log
#
# 4_3cd_postprocess.sh  --  Part-2 Step-4 post-processing (light, base R).
# ---------------------------------------------------------------------------
# Runs, per genome (headline Pre vs Post-wd):
#   4_3c_annotation_diff.R        cross-stage contamination confirmation + coherence delta
#   4_3d_marker_informativeness.R empirical marker informativeness (canonical vs PlantscRNAdb), per stage
#
# Prereq: 4_3_eval_pergenome.sh finished (annotation dirs + marker_zscore.tsv present) and Step-3
#         (3_1) cross-stage outputs exist in SM2v2_indep/step3_compare/.
# Usage:  sbatch 0_scripts/part2_split/step3_compare/4_3cd_postprocess.sh
# ---------------------------------------------------------------------------

# conda activate ambientmapper-manuscript   (environment.yml at the repo root)
set -euo pipefail

PROJECT_ROOT="${PROJECT_ROOT:-.}"
BASE="${PROJECT_ROOT}/6_socrates"
SCR="${SCRIPTS:-${BASE}/0_scripts}/part2_split/step3_compare"
STEP3="${BASE}/SM2v2_indep/step3_compare"
MARK="${BASE}/_data/markers"
Z_THR=1.5

# genome : original canonical panel (for source tagging in 4_3d)
declare -a GEN=( B73v5 TAIR10 )
declare -A CANON=( [B73v5]="${MARK}/markers.maize.Marand2025.bed" [TAIR10]="${MARK}/markers.At.Ecker2025.bed" )

for g in "${GEN[@]}"; do
  pre="SM2_${g}"; post="Clean.SM2v2wd_${g}"
  preann="${STEP3}/4_3_annotation_${pre}/${pre}.cluster_annotation.tsv"
  postann="${STEP3}/4_3_annotation_${post}/${post}.cluster_annotation.tsv"

  echo "=== [$g] 4_3c cross-stage annotation diff ==="
  Rscript "${SCR}/4_3c_annotation_diff.R" "$g" "$preann" "$postann" "$STEP3" "$STEP3"

  echo "=== [$g] 4_3d marker informativeness (per stage) ==="
  for st in "$pre" "$post"; do
    zt="${STEP3}/4_3_annotation_${st}/${st}.marker_zscore.tsv"
    at="${STEP3}/4_3_annotation_${st}/${st}.cluster_annotation.tsv"
    Rscript "${SCR}/4_3d_marker_informativeness.R" "$zt" "${CANON[$g]}" "${STEP3}/${st}" "$Z_THR" "$at"
  done
done
echo "$(date): Step-4 post-processing DONE"
