#!/bin/bash
########## BATCH Lines for Resource Request ##########
#SBATCH --time=6:00:00
#SBATCH --nodes=1
#SBATCH --ntasks-per-node=1
#SBATCH --cpus-per-task=8
#SBATCH --mem=80G
#SBATCH --job-name=P1indep_2_3_cluster
#SBATCH --partition=standard
#SBATCH --output=_logs/P1indep_2_3_cluster_%A_%a.log
#SBATCH --array=0-1
#
# 2_3_cluster.sh -- STAGE 2c: per-genome Leiden clustering at the FROZEN config.
# ---------------------------------------------------------------------------
#   Task 0: SM2_B73v5   pcs=20 k=20 min_dist=0.05 min_c=50 res=0.5   (12 clusters)
#   Task 1: SM2_TAIR10  pcs=8  k=30 min_dist=0.3  min_c=50 res=1     ( 8 clusters)
#
# Config frozen after the 2_1 gridscan + 2_2 resolution scan (plus config-search sweeps, not shipped).
# Produces SM2v2_indep/step2_cluster/<prefix>.full.SocObj_v7.<tag>.rds -- the cleaned-tile +
# v7-cell objects that Stage 3 (3_1_coembed_build.R) stacks into the co-embedding.
# Runs 2_3_cluster.R (== 3_0_0b_cluster_pergenome.R).
# ---------------------------------------------------------------------------
set -euo pipefail
# conda activate ambientmapper-manuscript   (environment.yml at the repo root)

PROJECT_ROOT="${PROJECT_ROOT:-.}"
BASE="${PROJECT_ROOT}/6_socrates"
cd "$BASE"
SCR="${SCRIPTS:-${BASE}/0_scripts}/part1_indep"
QC="${BASE}/SM2v2_indep/step0_qc"
OUT="${BASE}/SM2v2_indep/step2_cluster"
mkdir -p "$OUT"

# frozen per-genome config: pcs k_near min_dist min_c resolution
PREFIXES=(SM2_B73v5 SM2_TAIR10)
GENOMES=(B73v5 TAIR10)
PCS=(20 8); KNN=(20 30); MD=(0.05 0.3); RES=(0.5 1); MINC=(50 50)

i="${SLURM_ARRAY_TASK_ID:-${1:?set via --array or pass 0/1}}"
PREFIX="${PREFIXES[$i]}"; GENOME="${GENOMES[$i]}"
SOC="${QC}/${PREFIX}.raw.soc.rds"
META="${QC}/${PREFIX}.minDepth200.updated_metadata_v6.txt"
for f in "$SOC" "$META"; do [[ -s "$f" ]] || { echo "MISSING: $f -- run Stage 1 (1_0_qc_run.sh) first"; exit 2; }; done

echo "$(date): [$PREFIX] cluster | pcs=${PCS[$i]} k=${KNN[$i]} min_dist=${MD[$i]} min_c=${MINC[$i]} res=${RES[$i]}"
Rscript "${SCR}/2_3_cluster.R" \
        "$SOC" "$META" "$OUT" "$PREFIX" \
        "${PCS[$i]}" "${KNN[$i]}" "${MD[$i]}" "${MINC[$i]}" "${RES[$i]}" 1 "$GENOME" 50
echo "$(date): DONE [$PREFIX] -> ${OUT}/${PREFIX}.full.SocObj_v7.*.rds"
