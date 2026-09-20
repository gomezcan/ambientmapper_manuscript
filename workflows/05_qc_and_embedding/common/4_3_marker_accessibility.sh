#!/usr/bin/env bash
#SBATCH --time=02:00:00
#SBATCH --nodes=1
#SBATCH --ntasks-per-node=1
#SBATCH --cpus-per-task=4
#SBATCH --mem=80G
#SBATCH --job-name=marker_accessibility
#SBATCH --partition=standard
#SBATCH --output=_logs/%x-%j.log
#
# Cluster-annotation phase, Task 4_3: pseudobulk marker-accessibility annotation.
# Emits the figure-feeding TABLES (+ exploratory plots) per stage; the final figure is built in
# analysis/ from these tables. Adapted from the maize_282 reference annotation code (dependency-light).
#
# Array: 0 = PRE  (SM2/step5_markers/SM2.genes.sparse.rds            -> SM2/step5_markers/)
#        1 = POST (SM2v2_clean/step5_markers/Clean.SM2v2.genes.sparse.rds -> SM2v2_clean/step5_markers/)
#   sbatch --array=0-1 0_scripts/common/4_3_marker_accessibility.sh
#
# Entry objects = the minc_50 clustering (canonical, non-bal). Override the clustering tag with TAG=...
#   TAG=pcs_20.k_near_30.min_dis_0.3.minc_50 sbatch --array=0-1 0_scripts/common/4_3_marker_accessibility.sh

set -euo pipefail
PROJECT_ROOT="${PROJECT_ROOT:-.}"
BASE="${PROJECT_ROOT}/6_socrates"
cd "$BASE"
# conda activate ambientmapper-manuscript   (environment.yml at the repo root)

MARKERS="_data/markers/markers.SM2v2.bed"
TAG="${TAG:-pcs_20.k_near_30.min_dis_0.3.minc_50}"
mkdir -p _logs

case "${SLURM_ARRAY_TASK_ID:-0}" in
  0) OUTDIR="SM2/step5_markers";            PREFIX="SM2"
     MAT="SM2/step5_markers/SM2.genes.sparse.rds"
     META="SM2/step4_cluster/SM2.updated_metadata_v7.${TAG}.txt"
     STAGE="PreClean" ;;
  1) OUTDIR="SM2v2_clean/step5_markers";    PREFIX="Clean.SM2v2"
     MAT="SM2v2_clean/step5_markers/Clean.SM2v2.genes.sparse.rds"
     META="SM2v2_clean/step4_cluster/Clean.SM2v2.updated_metadata_v7.${TAG}.txt"
     STAGE="PostClean" ;;
  *) echo "ERROR: array index must be 0 (PRE) or 1 (POST), got ${SLURM_ARRAY_TASK_ID:-unset}"; exit 1 ;;
esac

for f in "$MAT" "$META" "$MARKERS"; do
  [[ -s "$f" ]] || { echo "ERROR: missing input: $f"; exit 1; }
done

echo " - 4_3 marker accessibility: STAGE=${STAGE}  OUTDIR=${OUTDIR}"
Rscript "${SCRIPTS:-${BASE}/0_scripts}/common/4_3_marker_accessibility.R" "$OUTDIR" "$PREFIX" "$MAT" "$META" "$MARKERS" "$STAGE"
echo " - done -> ${OUTDIR}"
