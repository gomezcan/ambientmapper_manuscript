#!/bin/bash
########## BATCH Lines for Resource Request ##########
#SBATCH --time=1:30:00
#SBATCH --nodes=1
#SBATCH --ntasks-per-node=1
#SBATCH --cpus-per-task=10
#SBATCH --mem=90G
#SBATCH --job-name=Step1c_kgrid_plate
#SBATCH --partition=standard
#SBATCH --output=_logs/Step1c_kgrid_plate_%A_%a.log
#SBATCH --array=0-1
#
# 3_0_1c_umap_kgrid.sh  --  Step 1b (Part 2, PLATE arm): SECOND grid, sweeping k_near and
#                           lowering m.clst, per pcs, at every resolution. Extends 3_0_1b.
# ---------------------------------------------------------------------------
# WHY: for a SMALL cell set (plate At ~1090) the fixed k_near=30 / m.clst=50 are large
#   relative to n and can THEMSELVES force the merge to 2 -- so "At = 2" may be a k/m.clst
#   artifact, not biology. This sweeps k_near and drops m.clst=40 to test that directly.
#   (Grid: pcs {5,8,10} x k_near {10,15,20,25,30} x a 6-value res grid, m.clst=40.)
#
#   task 0 = At  (SM2_At_TAIR10)  : pcs {5,8,10}   k_near {10,15,20,25,30}  md=0.3  min.c=50  m.clst=40  res {0.1,0.2,0.3,0.5,0.8,1.0}
#   task 1 = B73 (SM2_B73_B73v5)  : pcs {10,15,20} k_near {15,20,30}        md=0.05 min.c=50  m.clst=40  res {0.3,0.5,0.8,1.0}
#
#   sbatch --array=0 0_scripts/part2_split/step1_cluster_opt/3_0_1c_umap_kgrid.sh   # At
#   sbatch --array=1 0_scripts/part2_split/step1_cluster_opt/3_0_1c_umap_kgrid.sh   # B73 -- optional analog
#
# k_near enters BOTH projectUMAP (layout) and callClusters (graph), so each (pcs,k) is its
#   own embedding (faithful to Step 2). reduceDims is k-independent -> once per pcs.
#
# PREREQ: 1_QC_scifiATAC_SM2v2_plate.sh done (raw.soc.rds + v6 for both Pre configs).
# OUTPUT: SM2v2_plate/step1_cluster_opt/<prefix>/plots/<prefix>.umap_panels.kgrid.pdf
#         + <prefix>.umap_grid_summary.kgrid.tsv  (n_clusters + n_unlabeled per pcs x k x res)
#   (distinct ".kgrid" tag -> does NOT overwrite the 3_0_1b pcs x res outputs)
# ---------------------------------------------------------------------------

set -euo pipefail
PROJECT_ROOT="${PROJECT_ROOT:-.}"
BASE="${PROJECT_ROOT}/6_socrates"
SCRIPTS="${SCRIPTS:-${BASE}/0_scripts}"
cd "${BASE}"
# conda activate ambientmapper-manuscript   (environment.yml at the repo root)

# per-task params (0 = At/TAIR10, 1 = B73/B73v5)
PREFIXES=(SM2_At_TAIR10             SM2_B73_B73v5)
GENOMES=(TAIR10                     B73v5)
PCS_LISTS=("5,8,10"                 "10,15,20")
K_LISTS=("10,15,20,25,30"           "15,20,30")
MDS=(0.3                            0.05)
MINCS=(50                           50)
RESGRIDS=("0.1,0.2,0.3,0.5,0.8,1.0" "0.3,0.5,0.8,1.0")
MCLST=40

i="${SLURM_ARRAY_TASK_ID}"
PREFIX="${PREFIXES[$i]}"; GENOME="${GENOMES[$i]}"
PCS="${PCS_LISTS[$i]}"; K="${K_LISTS[$i]}"; MD="${MDS[$i]}"; MINC="${MINCS[$i]}"; RESGRID="${RESGRIDS[$i]}"

SOC="SM2v2_plate/step0_qc/${PREFIX}.raw.soc.rds"
META="SM2v2_plate/step0_qc/${PREFIX}.minDepth200.updated_metadata_v6.txt"
OUTDIR="SM2v2_plate/step1_cluster_opt/${PREFIX}"
mkdir -p "${OUTDIR}/plots"

[[ -s "$SOC"  ]] || { echo "ERROR: missing SocObj: $SOC";       exit 1; }
[[ -s "$META" ]] || { echo "ERROR: missing v6 metadata: $META"; exit 1; }

echo "$(date): Step-1c kgrid ${PREFIX} (${GENOME}) | pcs={${PCS}} k={${K}} md=${MD} min.c=${MINC} m.clst=${MCLST} | res={${RESGRID}}"
Rscript "${SCRIPTS}/part2_split/step1_cluster_opt/3_0_1c_umap_kgrid.R" \
  "$SOC" "$META" "$OUTDIR" "$PREFIX" "$PCS" "$K" "$MD" "$MINC" "$RESGRID" "$MCLST" 1

echo "$(date): done ${PREFIX}. Eyeball ${OUTDIR}/plots/${PREFIX}.umap_panels.kgrid.pdf"
echo "  page 1 = n_clusters heatmap (k x res x pcs); page 2 = n_unlabeled; then per-pcs UMAP grids + depth overlays."
echo "  decision table: ${OUTDIR}/${PREFIX}.umap_grid_summary.kgrid.tsv"