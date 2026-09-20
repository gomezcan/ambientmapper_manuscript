#!/bin/bash
########## BATCH Lines for Resource Request ##########
#SBATCH --time=10:00:00
#SBATCH --nodes=1
#SBATCH --ntasks-per-node=1
#SBATCH --cpus-per-task=10
#SBATCH --mem=120G
#SBATCH --job-name=Step1b_resscan
#SBATCH --partition=standard
#SBATCH --output=_logs/Step1b_resscan_%A_%a.log
#SBATCH --array=0-1
#
# Step 1b (Part 2) -- pcs x resolution stability scan on the PRE objects.
# 1a could not principledly pick pcs, so 1b ARBITRATES pcs jointly with resolution by
# bootstrap cluster stability (the non-circular objective). Fixed k_near + min_dist
# (min_dist is UMAP-viz-only; Leiden clusters on the SVD graph). Marker-free.
#
# Grid: pcs {20,30,50} x res {0.3,0.5,0.8,1.0,1.5}. k_near=30 (matches combined-SM2 3_0_0
# default; 1a did not separate k). min_dist=0.05 (1a; viz only).
# min.c per genome (MUST match the intended clustering floor): B73v5 NA(data-driven ~250),
# TAIR10 50 (the data-driven floor annihilates the low-coverage minor genome).
#
#   sbatch 0_scripts/part1_indep/2_2_resolution_scan.sh              # both genomes (array 0-1)
#   sbatch --array=1 0_scripts/part1_indep/2_2_resolution_scan.sh    # TAIR10 only (fast)
#
# B73v5 (~25k cells) is the heavy task: 3 pcs x 5 res x (1 + n_boot) Leiden clusterings
# (cleanData+tfidf shared; reduceDims+scan per pcs). TAIR10 (~3k) finishes fast.

set -euo pipefail
PROJECT_ROOT="${PROJECT_ROOT:-.}"
BASE="${PROJECT_ROOT}/6_socrates"
SCRIPTS="${SCRIPTS:-${BASE}/0_scripts}"
cd "${BASE}"
# conda activate ambientmapper-manuscript   (environment.yml at the repo root)

# task 0 = B73v5, task 1 = TAIR10
PREFIXES=(SM2_B73v5 SM2_TAIR10)
GENOMES=(B73v5 TAIR10)
MINCS=(NA 50)

# frozen embedding axes (1b arbitrates pcs)
PCS_LIST="20,30,50"
K=30
MD=0.05
NBOOT=10
BOOTFRAC=0.8

i="${SLURM_ARRAY_TASK_ID}"
PREFIX="${PREFIXES[$i]}"; GENOME="${GENOMES[$i]}"; MINC="${MINCS[$i]}"

SOC="SM2v2_indep/step0_qc/${PREFIX}.raw.soc.rds"
META="SM2v2_indep/step0_qc/${PREFIX}.minDepth200.updated_metadata_v6.txt"
OUTDIR="SM2v2_indep/step1_cluster_opt/${PREFIX}"
mkdir -p "${OUTDIR}/plots"

[[ -s "$SOC"  ]] || { echo "ERROR: missing SocObj: $SOC";       exit 1; }
[[ -s "$META" ]] || { echo "ERROR: missing v6 metadata: $META"; exit 1; }

echo " - Step 1b pcs x res scan: ${PREFIX} (${GENOME}) | pcs={${PCS_LIST}} k=${K} min_dist=${MD} | min.c=${MINC}"
Rscript "${SCRIPTS}/part1_indep/2_2_resolution_scan.R" \
  "$SOC" "$META" "$OUTDIR" "$PREFIX" "$PCS_LIST" "$K" "$MD" "$MINC" 1 "$NBOOT" "$BOOTFRAC"

echo " - done. Review ${OUTDIR}/${PREFIX}.resolution_scan.tsv + plots/${PREFIX}.resolution_scan.pdf ,"
echo "   then finalize SM2v2_indep/${GENOME}.cluster_config.tsv from ${PREFIX}.cluster_config.draft.tsv."
