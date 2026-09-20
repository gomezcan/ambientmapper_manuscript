#!/bin/bash
########## BATCH Lines for Resource Request ##########
#SBATCH --time=8:00:00
#SBATCH --nodes=1
#SBATCH --ntasks-per-node=1
#SBATCH --cpus-per-task=8
#SBATCH --mem=96G
#SBATCH --job-name=P1indep_3_2_grid
#SBATCH --partition=standard
#SBATCH --output=_logs/P1indep_3_2_grid_%A_%a.log
#SBATCH --array=0-1
#
# 3_2_coembed_gridscan.sh -- STAGE 2: Fig-1/S3 grid scan on the co-embedding.
# ---------------------------------------------------------------------------
#   Task 0: SM2v2_coembed_Pre      Task 1: SM2v2_coembed_Post_wd
#
# Runs 3_2_coembed_gridscan.R (== the combined Fig-1/S3 engine) on the joint object:
# pcs[20,25,30,40,50] x k[15,20,30] x min_dist[.05,.15,.3], scores on genome_mixing,
# min.c=50. genome_mixing here = plate-of-origin -> the SAME quantity as combined Fig S3.
# THE VERIFICATION: the Fig-S3 text claims (i) mixing stable across the grid,
# (ii) mixing INCREASES with knn_preservation, (iii) QC correlations low. Read those
# three off UMAP_grid_scan.minc_50.metrics.tsv here. Requires Stage 1 done.
# ---------------------------------------------------------------------------
set -euo pipefail
# conda activate ambientmapper-manuscript   (environment.yml at the repo root)

PROJECT_ROOT="${PROJECT_ROOT:-.}"
BASE="${PROJECT_ROOT}/6_socrates"
cd "$BASE"
OUT="${BASE}/SM2v2_indep/coembed"
SCR="${SCRIPTS:-${BASE}/0_scripts}/part1_indep"
MINC="${MINC:-50}"

i="${SLURM_ARRAY_TASK_ID:-${1:?set via --array or pass 0/1}}"
case "$i" in
  0) STAGE="Pre" ;;
  1) STAGE="Post_wd" ;;
  *) echo "bad array index $i"; exit 2 ;;
esac
PREFIX="${OUT}/SM2v2_coembed_${STAGE}"
[[ -s "${PREFIX}.coembed.soc.rds" ]] || { echo "MISSING Stage-1 output: ${PREFIX}.coembed.soc.rds -- run 1_1 first"; exit 2; }

GRID_OUT="${OUT}/${STAGE}_grid"
mkdir -p "${GRID_OUT}/plots"
echo "$(date): [$STAGE] grid scan (min.c=${MINC}) -> ${GRID_OUT}"
Rscript "${SCR}/3_2_coembed_gridscan.R" \
        "${PREFIX}.coembed.soc.rds" "${PREFIX}.coembed.meta.tsv" "$GRID_OUT" 1 "$MINC"
echo "$(date): DONE [$STAGE]"
echo " -> ${GRID_OUT}/UMAP_grid_scan.minc_${MINC}.metrics.tsv (genome_mixing = plate-of-origin)"
