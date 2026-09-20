#!/usr/bin/env bash
#SBATCH --time=03:00:00
#SBATCH --nodes=1
#SBATCH --ntasks-per-node=1
#SBATCH --cpus-per-task=8
#SBATCH --mem=100G
#SBATCH --job-name=marker_smoothing
#SBATCH --partition=standard
#SBATCH --output=_logs/%x-%j.log
#
# Cluster-annotation phase, Task 4_3b: Markov-graph smoothing of MARKER gene activity for the
# marker-on-UMAP figure panels (companion to 4_3, which does cluster-level tables / no smoothing).
# Output (smoothed marker x cell RDS + coords) is consumed by the figure scripts.
#
# Array: 0 = PRE  (SM2)            1 = POST (SM2v2_clean)
#   sbatch --array=0-1 0_scripts/common/4_3b_marker_smoothing.sh
#
# Tunables (env): TAG (clustering tag), K (knn, default 25), STEP (diffusion hops, default 3).
#   K=25 STEP=3 sbatch --array=0-1 0_scripts/common/4_3b_marker_smoothing.sh

set -euo pipefail
PROJECT_ROOT="${PROJECT_ROOT:-.}"
BASE="${PROJECT_ROOT}/6_socrates"
cd "$BASE"
# conda activate ambientmapper-manuscript   (environment.yml at the repo root)

MARKERS="_data/markers/markers.SM2v2.bed"
TAG="${TAG:-pcs_20.k_near_30.min_dis_0.3.minc_50}"
K="${K:-25}"; STEP="${STEP:-3}"
mkdir -p _logs

case "${SLURM_ARRAY_TASK_ID:-0}" in
  0) OUTDIR="SM2/step5_markers";         PREFIX="SM2"
     MAT="SM2/step5_markers/SM2.genes.sparse.rds"
     META="SM2/step4_cluster/SM2.updated_metadata_v7.${TAG}.txt"
     RD="SM2/step4_cluster/SM2.reduced_dimensions_v7.${TAG}.txt"
     STAGE="PreClean" ;;
  1) OUTDIR="SM2v2_clean/step5_markers"; PREFIX="Clean.SM2v2"
     MAT="SM2v2_clean/step5_markers/Clean.SM2v2.genes.sparse.rds"
     META="SM2v2_clean/step4_cluster/Clean.SM2v2.updated_metadata_v7.${TAG}.txt"
     RD="SM2v2_clean/step4_cluster/Clean.SM2v2.reduced_dimensions_v7.${TAG}.txt"
     STAGE="PostClean" ;;
  *) echo "ERROR: array index must be 0 (PRE) or 1 (POST), got ${SLURM_ARRAY_TASK_ID:-unset}"; exit 1 ;;
esac

for f in "$MAT" "$META" "$RD" "$MARKERS"; do
  [[ -s "$f" ]] || { echo "ERROR: missing input: $f"; exit 1; }
done

echo " - 4_3b marker smoothing: STAGE=${STAGE}  OUTDIR=${OUTDIR}  k=${K} step=${STEP}"
Rscript "${SCRIPTS:-${BASE}/0_scripts}/common/4_3b_marker_smoothing.R" "$OUTDIR" "$PREFIX" "$MAT" "$META" "$RD" "$MARKERS" "$STAGE" "$K" "$STEP"
echo " - done -> ${OUTDIR}"
