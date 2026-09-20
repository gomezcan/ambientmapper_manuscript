#!/bin/bash
########## BATCH Lines for Resource Request ##########
#SBATCH --time=6:00:00
#SBATCH --nodes=1
#SBATCH --ntasks-per-node=1
#SBATCH --cpus-per-task=10
#SBATCH --mem=90G
#SBATCH --job-name=Step1b_resscan_plate
#SBATCH --partition=standard
#SBATCH --output=_logs/Step1b_resscan_plate_%A_%a.log
#SBATCH --array=0-1
#
# 3_0_1_scan_plate.sh  --  Step 1b (Part 2, PLATE arm): per-genome pcs x resolution scan.
# ---------------------------------------------------------------------------
# WHY THIS EXISTS (separate from the indep 3_0_1_resolution_scan.sh)
#   The plate cell sets are a DIFFERENT data regime from the indep objects, so the frozen
#   indep configs do not transfer. Proven: the borrowed indep TAIR10 config
#   (res=1, tuned on 2937 mixed cells) shatters the ~800 pure At plate cells into 94 clusters.
#   Each (arm x genome) regime needs its OWN config, optimized on its Pre object and frozen
#   across its 3 stages (plan decision B). This scan produces the plate configs.
#
# KEY DIFFERENCE from the indep scan: a LOWER resolution grid for At.
#   res=1 already gives 94 clusters on ~800 cells, so the sensible range is BELOW 0.3. The
#   R script now takes an optional res_grid (arg 12); At sweeps {0.05..0.5}, B73 {0.3..1.0}.
#
# Optimize on PRE (the superset), freeze, then apply the IDENTICAL config to Pre/nd/wd in
# Step 2 (re-run 3_0_eval_plate_v5v6.sh with the frozen params replacing the borrowed ones).
# Optimize on v6 metadata (v6 is canonical -- the qc_check-dropped cells were shown to be genuinely
# low quality: half depth/sites, FRiP 0.31 vs 0.57).
#
#   task 0 = At  (SM2_At_TAIR10)  : pcs {8,15}  k=30 md=0.3  min.c=50  res {0.05,0.1,0.2,0.3,0.5}
#   task 1 = B73 (SM2_B73_B73v5)  : pcs {20}    k=20 md=0.05 min.c=50  res {0.3,0.5,0.8,1.0}
#
#   sbatch --array=0 0_scripts/part2_split/step1_cluster_opt/3_0_1_scan_plate.sh   # At -- fast + decisive
#   sbatch --array=1 0_scripts/part2_split/step1_cluster_opt/3_0_1_scan_plate.sh   # B73 -- verification
#   sbatch          0_scripts/part2_split/step1_cluster_opt/3_0_1_scan_plate.sh   # both
#
# PREREQ: 1_QC_scifiATAC_SM2v2_plate.sh done (raw.soc.rds + v6 for both Pre configs).
# OUTPUT: SM2v2_plate/step1_cluster_opt/<prefix>/<prefix>.resolution_scan.{tsv,pdf} + .cluster_config.draft.tsv
# ---------------------------------------------------------------------------

set -euo pipefail
PROJECT_ROOT="${PROJECT_ROOT:-.}"
BASE="${PROJECT_ROOT}/6_socrates"
SCRIPTS="${SCRIPTS:-${BASE}/0_scripts}"
cd "${BASE}"
# conda activate ambientmapper-manuscript   (environment.yml at the repo root)

# per-task params (0 = At/TAIR10, 1 = B73/B73v5)
PREFIXES=(SM2_At_TAIR10          SM2_B73_B73v5)
GENOMES=(TAIR10                  B73v5)
PCS_LISTS=("8,15"                "20")
KS=(30                          20)
MDS=(0.3                        0.05)
MINCS=(50                       50)
RESGRIDS=("0.05,0.1,0.2,0.3,0.5" "0.3,0.5,0.8,1.0")
NBOOT=10
BOOTFRAC=0.8

i="${SLURM_ARRAY_TASK_ID}"
PREFIX="${PREFIXES[$i]}"; GENOME="${GENOMES[$i]}"
PCS_LIST="${PCS_LISTS[$i]}"; K="${KS[$i]}"; MD="${MDS[$i]}"; MINC="${MINCS[$i]}"; RESGRID="${RESGRIDS[$i]}"

SOC="SM2v2_plate/step0_qc/${PREFIX}.raw.soc.rds"
META="SM2v2_plate/step0_qc/${PREFIX}.minDepth200.updated_metadata_v6.txt"
OUTDIR="SM2v2_plate/step1_cluster_opt/${PREFIX}"
mkdir -p "${OUTDIR}/plots"

[[ -s "$SOC"  ]] || { echo "ERROR: missing SocObj: $SOC";       exit 1; }
[[ -s "$META" ]] || { echo "ERROR: missing v6 metadata: $META"; exit 1; }

echo "$(date): Step-1b plate res-scan ${PREFIX} (${GENOME}) | pcs={${PCS_LIST}} k=${K} md=${MD} min.c=${MINC} | res={${RESGRID}}"
Rscript "${SCRIPTS}/part2_split/step1_cluster_opt/3_0_1_resolution_scan.R" \
  "$SOC" "$META" "$OUTDIR" "$PREFIX" "$PCS_LIST" "$K" "$MD" "$MINC" 1 "$NBOOT" "$BOOTFRAC" "$RESGRID"

echo "$(date): done ${PREFIX}. Review ${OUTDIR}/${PREFIX}.resolution_scan.tsv + plots/${PREFIX}.resolution_scan.pdf ,"
echo "  then freeze SM2v2_plate/${GENOME}.plate.cluster_config.tsv from ${PREFIX}.cluster_config.draft.tsv and re-run Step 2 with the frozen params."
