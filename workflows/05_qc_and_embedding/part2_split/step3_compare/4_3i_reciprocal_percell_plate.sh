#!/bin/bash
########## BATCH Lines for Resource Request ##########
#SBATCH --time=2:00:00
#SBATCH --nodes=1
#SBATCH --ntasks-per-node=1
#SBATCH --cpus-per-task=4
#SBATCH --mem=64G
#SBATCH --job-name=Step4_rZpercell_SM2v2_plate
#SBATCH --partition=standard
#SBATCH --output=_logs/Step4_rZpercell_SM2v2_plate_%A_%a.log
#SBATCH --array=0-5
#
# 4_3i_reciprocal_percell_plate.sh -- Step-4 (plate arm): PER-CELL reciprocal z-score.
# ---------------------------------------------------------------------------
# Per-cell bidirectional marker z on the SMOOTHED PERKB gene-activity matrix:
#   Zi = gene z ACROSS CELLS, Zj = cell z ACROSS GENES ; both large-N so the Euclidean
#   combine sqrt(max0(Zi)^2+max0(Zj)^2) is balanced (no 5-cluster 1.79 cap). Also emits the
#   AND variant sqrt(max0(Zi)*max0(Zj)) for comparison. rZ needs perkb (Zj is length-sensitive)
#   -> NORM defaults to perkb here (the other 4_3* drivers default raw).
#
#   0: SM2_At_TAIR10            Pre    3: SM2_B73_B73v5            Pre
#   1: Clean.SM2v2wd_At_TAIR10  wd     4: Clean.SM2v2wd_B73_B73v5  wd
#   2: Clean.SM2v2_At_TAIR10    nd     5: Clean.SM2v2_B73_B73v5    nd
#   GENE_LENGTH_NORM=perkb sbatch --array=0 ...  # At Pre  ;  --array=3 ... # B73 Pre
#
# Prereq: SMOOTHED perkb matrix present (4_3f_smooth_plate.sh with GENE_LENGTH_NORM=perkb) +
#         Step-2 v7 metadata. "all markers" = the augmented per-genome panel (rZ ranks them itself).
# ---------------------------------------------------------------------------

# conda activate ambientmapper-manuscript   (environment.yml at the repo root)
set -euo pipefail

PROJECT_ROOT="${PROJECT_ROOT:-.}"
BASE="${PROJECT_ROOT}/6_socrates"
SCRIPTS="${SCRIPTS:-${BASE}/0_scripts}"
MARK="${BASE}/_data/markers"
STEP2="${BASE}/SM2v2_plate/step2_cluster"
OUT="${BASE}/SM2v2_plate/step3_compare"
s4_3i="${SCRIPTS}/common/4_3i_reciprocal_zscore_percell.R"
TOPN=6

NORM="${GENE_LENGTH_NORM:-perkb}"      # perkb (default) | raw ; rZ REQUIRES perkb to be meaningful
case "$NORM" in raw) VAR="" ;; perkb) VAR=".perkb" ;; *) echo "ERROR: GENE_LENGTH_NORM must be raw|perkb"; exit 1 ;; esac

declare -a NAMES=( SM2_At_TAIR10 Clean.SM2v2wd_At_TAIR10 Clean.SM2v2_At_TAIR10 SM2_B73_B73v5 Clean.SM2v2wd_B73_B73v5 Clean.SM2v2_B73_B73v5 )
declare -a STAGES=( Pre wd nd Pre wd nd )
i="${SLURM_ARRAY_TASK_ID:?set via --array}"
name="${NAMES[$i]:?bad array index $i}"; stage="${STAGES[$i]}"
species="${name##*_}"

case "$species" in
  B73v5)  panel="${MARK}/markers.maize.Marand2025_PlantscRNAdb4.bed" ;;
  TAIR10) panel="${MARK}/markers.At.Ecker2025_PlantscRNAdb4.bed" ;;
  *) echo "ERROR: unknown species '$species' from '$name'"; exit 1 ;;
esac

sm="${OUT}/${name}.plate${VAR}.genes.smoothed.rds"
meta="$(ls ${STEP2}/${name}.mQCv6.updated_metadata_v7.*.txt 2>/dev/null | head -1)"
adir="${OUT}/4_3i_reciprocal_percell_${name}${VAR}"

[[ -f "$sm" ]]                     || { echo "ERROR: missing smoothed matrix $sm (run 4_3f_smooth_plate.sh GENE_LENGTH_NORM=$NORM first)"; exit 2; }
[[ -n "${meta:-}" && -f "$meta" ]] || { echo "ERROR: no v7 metadata for $name in $STEP2"; exit 2; }
[[ -f "$panel" ]]                  || { echo "ERROR: missing panel $panel"; exit 2; }

echo "$(date): Step-4 per-cell rZ | $name ($stage norm=$NORM) | matrix=$(basename "$sm") panel=$(basename "$panel")"
Rscript "$s4_3i" "$adir" "$name" "$sm" "$meta" "$panel" "$stage" LouvainClusters "$TOPN"
echo "$(date): DONE $name -> $adir"
