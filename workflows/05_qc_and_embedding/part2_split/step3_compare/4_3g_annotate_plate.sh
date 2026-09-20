#!/bin/bash
########## BATCH Lines for Resource Request ##########
#SBATCH --time=1:00:00
#SBATCH --nodes=1
#SBATCH --ntasks-per-node=1
#SBATCH --cpus-per-task=4
#SBATCH --mem=48G
#SBATCH --job-name=Step4_annot_SM2v2_plate
#SBATCH --partition=standard
#SBATCH --output=_logs/Step4_annot_SM2v2_plate_%A_%a.log
#SBATCH --array=0-2
#
# 4_3g_annotate_plate.sh -- Step-4 STEP 3 (plate arm): per-cell annotation ON the smoothed activity.
# ---------------------------------------------------------------------------
# Reads *.plate.genes.smoothed.rds (from 4_3f_smooth_plate.sh) + the frozen v7 clusters, runs the
# Marand-style enrichment classifier (4_3g_cell_annotation.R). Re-run freely to tune z_thresh/ratio
# (env vars below) -- it does NOT re-smooth.
#
#   0: SM2_At_TAIR10   Pre    3: SM2_B73_B73v5           Pre
#   1: Clean.SM2v2wd_At_TAIR10 wd   4: Clean.SM2v2wd_B73_B73v5 wd
#   2: Clean.SM2v2_At_TAIR10   nd   5: Clean.SM2v2_B73_B73v5   nd
#   sbatch --array=0-2 ...  # At (decisive) ;  sbatch --array=3-5 ...  # B73
#   Z_THRESH=2 RATIO=1.5 sbatch --array=0-2 ...   # override the call gates
#
# ANNOTATION PANEL = informative_top15 (the mean-z call uses the informative panel + min_markers=3).
# Prereq: 4_3f_smooth_plate.sh done (smoothed matrices present).
# ---------------------------------------------------------------------------

# conda activate ambientmapper-manuscript   (environment.yml at the repo root)
set -euo pipefail

PROJECT_ROOT="${PROJECT_ROOT:-.}"
BASE="${PROJECT_ROOT}/6_socrates"
SCRIPTS="${SCRIPTS:-${BASE}/0_scripts}"
STEP2="${BASE}/SM2v2_plate/step2_cluster"
MARK="${BASE}/_data/markers"
OUT="${BASE}/SM2v2_plate/step3_compare"
s4_3g="${SCRIPTS}/common/4_3g_cell_annotation.R"
Z_THRESH="${Z_THRESH:-2}"; RATIO="${RATIO:-1.5}"; MIN_MARKERS="${MIN_MARKERS:-3}"
NORM="${GENE_LENGTH_NORM:-raw}"        # raw (default) | perkb -> consume the length-normalized smoothed matrix
case "$NORM" in raw) VAR="" ;; perkb) VAR=".perkb" ;; *) echo "ERROR: GENE_LENGTH_NORM must be raw|perkb"; exit 1 ;; esac

declare -a NAMES=( SM2_At_TAIR10 Clean.SM2v2wd_At_TAIR10 Clean.SM2v2_At_TAIR10 SM2_B73_B73v5 Clean.SM2v2wd_B73_B73v5 Clean.SM2v2_B73_B73v5 )
declare -a STAGES=( Pre wd nd Pre wd nd )
i="${SLURM_ARRAY_TASK_ID:?set via --array}"
name="${NAMES[$i]:?bad array index $i}"; stage="${STAGES[$i]}"; species="${name##*_}"

case "$species" in
  B73v5)  panel="${MARK}/markers.maize.informative_top15.bed" ;;
  TAIR10) panel="${MARK}/markers.At.informative_top15.bed" ;;
  *) echo "ERROR: unknown species '$species'"; exit 1 ;;
esac

sm="${OUT}/${name}.plate${VAR}.genes.smoothed.rds"
meta="$(ls ${STEP2}/${name}.mQCv6.updated_metadata_v7.*.txt 2>/dev/null | head -1)"
adir="${OUT}/4_3g_annotation_plate_${name}${VAR}"
[[ -f "$sm" ]]                     || { echo "ERROR: missing smoothed matrix $sm (run 4_3f_smooth_plate.sh first)"; exit 2; }
[[ -n "${meta:-}" && -f "$meta" ]] || { echo "ERROR: no v7 metadata for $name"; exit 2; }
[[ -f "$panel" ]]                  || { echo "ERROR: missing panel $panel"; exit 2; }

echo "$(date): Step-4 ANNOTATE | $name ($stage norm=$NORM) | smoothed=$(basename "$sm") panel=$(basename "$panel") z_thresh=$Z_THRESH ratio=$RATIO"
Rscript "$s4_3g" "$adir" "$name" "$sm" "$meta" "$panel" "$stage" "$Z_THRESH" "$RATIO" "$MIN_MARKERS"
echo "$(date): DONE $name -> $adir/${name}.cluster_majority.tsv (+ .cell_annotation.tsv)"
