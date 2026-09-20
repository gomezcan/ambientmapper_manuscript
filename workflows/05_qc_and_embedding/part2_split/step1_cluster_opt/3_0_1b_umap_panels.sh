#!/bin/bash
########## BATCH Lines for Resource Request ##########
#SBATCH --time=1:30:00
#SBATCH --nodes=1
#SBATCH --ntasks-per-node=1
#SBATCH --cpus-per-task=10
#SBATCH --mem=90G
#SBATCH --job-name=Step1b_umap_plate
#SBATCH --partition=standard
#SBATCH --output=_logs/Step1b_umap_plate_%A_%a.log
#SBATCH --array=0-1
#
# 3_0_1b_umap_panels.sh  --  Step 1b (Part 2, PLATE arm): UMAP visualization companion
#                            to the 3_0_1 res scan (which plots metric curves, no UMAP).
# ---------------------------------------------------------------------------
# WHY: the plan requires an eyeball call for At ("continuum vs real islands; pick bins by
#   eye"), impossible from stability curves alone -- and a single pcs is not a decision, so
#   pcs is SWEPT (like 2_0_2b / 3_0_4). This renders the SAME embedding Step 2 would freeze,
#   across a pcs x resolution GRID, colored by cluster AND by QC (depth/nSites/FRiP/pTSS/pOrg/dif).
#   If At's split tracks a depth/sites gradient across pcs -> continuum (root gradient), not biology.
#
#   task 0 = At  (SM2_At_TAIR10)  : pcs {5,8,10,15,20,30}  k=30 md=0.3  min.c=50  res {0.1,0.2,0.3,0.5,0.8,1.0}
#   task 1 = B73 (SM2_B73_B73v5)  : pcs {10,15,20,25,30}   k=20 md=0.05 min.c=50  res {0.3,0.5,0.8,1.0}
#
#   sbatch --array=0 0_scripts/part2_split/step1_cluster_opt/3_0_1b_umap_panels.sh   # At -- decisive
#   sbatch --array=1 0_scripts/part2_split/step1_cluster_opt/3_0_1b_umap_panels.sh   # B73 -- confirm 5 blobs
#   sbatch          0_scripts/part2_split/step1_cluster_opt/3_0_1b_umap_panels.sh   # both
#
# reduceDims recomputed per pcs (truncated SVD at that pcs, cor.max=0.6 -- may drop e.g.
# B73 pcs=20 -> 19; the R reports pcs_act and annotates facets "pcs=20 (r19)").
#
# PREREQ: 1_QC_scifiATAC_SM2v2_plate.sh done (raw.soc.rds + v6 for both Pre configs).
# OUTPUT: SM2v2_plate/step1_cluster_opt/<prefix>/plots/<prefix>.umap_panels.pdf  (pcs x res grid)
#         + <prefix>.umap_grid_summary.tsv  (n_clusters + n_unlabeled per pcs x res)
# ---------------------------------------------------------------------------

set -euo pipefail
PROJECT_ROOT="${PROJECT_ROOT:-.}"
BASE="${PROJECT_ROOT}/6_socrates"
SCRIPTS="${SCRIPTS:-${BASE}/0_scripts}"
cd "${BASE}"
# conda activate ambientmapper-manuscript   (environment.yml at the repo root)

# per-task params (0 = At/TAIR10, 1 = B73/B73v5) -- match 3_0_1_scan_plate.sh
PREFIXES=(SM2_At_TAIR10             SM2_B73_B73v5)
GENOMES=(TAIR10                     B73v5)
PCS_LISTS=("5,8,10,15,20,30"        "10,15,20,25,30")
KS=(30                             20)
MDS=(0.3                           0.05)
MINCS=(50                          50)
RESGRIDS=("0.1,0.2,0.3,0.5,0.8,1.0" "0.3,0.5,0.8,1.0")

i="${SLURM_ARRAY_TASK_ID}"
PREFIX="${PREFIXES[$i]}"; GENOME="${GENOMES[$i]}"
PCS="${PCS_LISTS[$i]}"; K="${KS[$i]}"; MD="${MDS[$i]}"; MINC="${MINCS[$i]}"; RESGRID="${RESGRIDS[$i]}"

SOC="SM2v2_plate/step0_qc/${PREFIX}.raw.soc.rds"
META="SM2v2_plate/step0_qc/${PREFIX}.minDepth200.updated_metadata_v6.txt"
OUTDIR="SM2v2_plate/step1_cluster_opt/${PREFIX}"
mkdir -p "${OUTDIR}/plots"

[[ -s "$SOC"  ]] || { echo "ERROR: missing SocObj: $SOC";       exit 1; }
[[ -s "$META" ]] || { echo "ERROR: missing v6 metadata: $META"; exit 1; }

echo "$(date): Step-1b UMAP grid ${PREFIX} (${GENOME}) | pcs={${PCS}} k=${K} md=${MD} min.c=${MINC} | res={${RESGRID}}"
Rscript "${SCRIPTS}/part2_split/step1_cluster_opt/3_0_1b_umap_panels.R" \
  "$SOC" "$META" "$OUTDIR" "$PREFIX" "$PCS" "$K" "$MD" "$MINC" "$RESGRID" 1

echo "$(date): done ${PREFIX}. Eyeball ${OUTDIR}/plots/${PREFIX}.umap_panels.pdf"
echo "  page 1 = pcs x res cluster grid; later pages = QC overlays (continuum check)."
echo "  decision table: ${OUTDIR}/${PREFIX}.umap_grid_summary.tsv (n_clusters + n_unlabeled per pcs x res)."
