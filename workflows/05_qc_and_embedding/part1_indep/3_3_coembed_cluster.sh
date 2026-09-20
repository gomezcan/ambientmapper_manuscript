#!/bin/bash
########## BATCH Lines for Resource Request ##########
#SBATCH --time=4:00:00
#SBATCH --nodes=1
#SBATCH --ntasks-per-node=1
#SBATCH --cpus-per-task=8
#SBATCH --mem=80G
#SBATCH --job-name=P1indep_3_3_cluster
#SBATCH --partition=standard
#SBATCH --output=_logs/P1indep_3_3_cluster_%A.log
#
# 3_3_coembed_cluster.sh -- STAGE 3: final Leiden clustering of the co-embedding -> Fig 1F/G/H.
# ---------------------------------------------------------------------------
# Runs 3_3_coembed_cluster.R (== 3_0_0b) on the co-embed object at the config CHOSEN FROM THE
# 3_2 grid (Fig S3). Emits <prefix>.updated_metadata_v7.<tag>.txt with umap1/umap2 +
# LouvainClusters + Genome -- the exact columns the Fig 1 script needs for panels F (by Genome),
# G (by Cluster), H (composition).
#
# SET THE CONFIG from 3_2's UMAP_grid_scan.minc_50.best.tsv (pick within the stable region):
#   PCS=.. KNN=.. MD=.. RES=.. sbatch 3_3_coembed_cluster.sh
# STAGE defaults to Pre (Fig 1F = the pre-clean problem). min.c=50.
#
# meta input: the ENRICHED .coembed.meta_full.tsv, built here by Stage 3.1b (see below).
#   The bare 3_1 .coembed.meta.tsv does NOT work: it lacks `nSites`, which callClusters needs
#   in its cluster-size filter, and total/pTSS/FRiP, which Fig 1F/G/H gate on. 3.1b runs
#   automatically and is skipped if its output is already newer than the 3_1 meta.
# Requires Stage 1 (and a config from Stage 2).
# ---------------------------------------------------------------------------
set -euo pipefail
# conda activate ambientmapper-manuscript   (environment.yml at the repo root)

PROJECT_ROOT="${PROJECT_ROOT:-.}"
BASE="${PROJECT_ROOT}/6_socrates"
cd "$BASE"
OUT="${BASE}/SM2v2_indep/coembed"
SCR="${SCRIPTS:-${BASE}/0_scripts}/part1_indep"

STAGE="${STAGE:-Pre}"
PCS="${PCS:?set PCS from 3_2 best.tsv}"
KNN="${KNN:?set KNN from 3_2 best.tsv}"
MD="${MD:?set MD (min_dist) from 3_2 best.tsv}"
RES="${RES:-0.5}"
MINC="${MINC:-50}"
MCLST="${MCLST:-50}"

PREFIX="${OUT}/SM2v2_coembed_${STAGE}"
[[ -s "${PREFIX}.coembed.soc.rds" ]] || { echo "MISSING Stage-1 output: ${PREFIX}.coembed.soc.rds"; exit 2; }
CL_OUT="${OUT}/${STAGE}_cluster"
mkdir -p "${CL_OUT}/plots"

# ---- Stage 3.1b: enrich the co-embed meta to the full v6 schema -------------------
# REQUIRED, not optional. The 3_1 meta lacks `nSites`, and Socrates::callClusters dies on that
# in its cluster-size filter AFTER cleanData/SVD/UMAP/graph all succeed. The same
# gap also hides total/pTSS/FRiP, which the Fig 1 script gates panels F/G/H on. 3_1b grafts the v6
# columns on (prefer B73v5, else TAIR10) into a NEW file; the 3_1 output is left untouched.
QC="${BASE}/SM2v2_indep/step0_qc"
case "$STAGE" in
  Pre)     V6_B="${QC}/SM2_B73v5.minDepth200.updated_metadata_v6.txt"
           V6_A="${QC}/SM2_TAIR10.minDepth200.updated_metadata_v6.txt" ;;
  Post_wd) V6_B="${QC}/Clean.SM2v2wd_B73v5.minDepth200.updated_metadata_v6.txt"
           V6_A="${QC}/Clean.SM2v2wd_TAIR10.minDepth200.updated_metadata_v6.txt" ;;
  *) echo "unknown STAGE '$STAGE' (expected Pre or Post_wd)"; exit 2 ;;
esac
META_FULL="${PREFIX}.coembed.meta_full.tsv"
if [[ -s "$META_FULL" && "$META_FULL" -nt "${PREFIX}.coembed.meta.tsv" ]]; then
  echo "$(date): [$STAGE] 3.1b meta already current -> $(basename "$META_FULL")"
else
  echo "$(date): [$STAGE] 3.1b enrich co-embed meta to the full v6 schema"
  Rscript "${SCR}/3_1b_coembed_meta_enrich.R" \
          "${PREFIX}.coembed.meta.tsv" "$V6_B" "$V6_A" "$META_FULL"
fi

echo "$(date): [$STAGE] cluster co-embed | pcs=${PCS} k=${KNN} min_dist=${MD} min_c=${MINC} res=${RES}"
Rscript "${SCR}/3_3_coembed_cluster.R" \
        "${PREFIX}.coembed.soc.rds" "$META_FULL" \
        "$CL_OUT" "SM2v2_coembed_${STAGE}" \
        "$PCS" "$KNN" "$MD" "$MINC" "$RES" 1 "coembed" "$MCLST"
echo "$(date): DONE [$STAGE] -> ${CL_OUT}/SM2v2_coembed_${STAGE}.updated_metadata_v7.*.txt  (umap1/umap2 + LouvainClusters + Genome)"
