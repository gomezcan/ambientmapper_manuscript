#!/bin/bash
########## BATCH Lines for Resource Request ##########
#SBATCH --time=2:00:00
#SBATCH --nodes=1
#SBATCH --ntasks-per-node=1
#SBATCH --cpus-per-task=6
#SBATCH --mem=48G
#SBATCH --job-name=Step1a_gridscan
#SBATCH --partition=standard
#SBATCH --output=_logs/Step1a_gridscan_%A_%a.log
#SBATCH --array=0-1
#
# Step 1a (Part 2) -- per-genome UMAP parameter scan on the PRE (superset) objects.
# Optimize on Pre; the chosen (pcs,k,min_dist) is frozen and applied identically to
# Pre + Post-wd in Step 2 (plan decision B: one frozen configuration per genome).
#
# NOTE on min.c: this is the cleanData "features/cell" knob applied at CLUSTERING time,
# NOT the minDepth=200 QC depth filter already baked into the v6 metadata. Per-genome:
#   B73v5  (major, high-cov) -> NA = Socrates data-driven floor (~250)
#   TAIR10 (minor, low-cov)  -> 50  (the data-driven floor annihilates the low-coverage minor genome)
# Both are trivially changeable here; 1b (3_0_1) reports cluster metrics to finalize.
#
#   sbatch 0_scripts/part1_indep/2_1_gridscan_pergenome.sh              # both genomes (array 0-1)
#   sbatch --array=0 0_scripts/part1_indep/2_1_gridscan_pergenome.sh    # B73v5 only

set -euo pipefail
PROJECT_ROOT="${PROJECT_ROOT:-.}"
BASE="${PROJECT_ROOT}/6_socrates"
SCRIPTS="${SCRIPTS:-${BASE}/0_scripts}"
cd "${BASE}"
# conda activate ambientmapper-manuscript   (environment.yml at the repo root)

# task 0 = B73v5 (major genome), task 1 = TAIR10 (minor genome)
PREFIXES=(SM2_B73v5 SM2_TAIR10)
GENOMES=(B73v5 TAIR10)
MINCS=(NA 50)

i="${SLURM_ARRAY_TASK_ID}"
PREFIX="${PREFIXES[$i]}"; GENOME="${GENOMES[$i]}"; MINC="${MINCS[$i]}"

SOC="SM2v2_indep/step0_qc/${PREFIX}.raw.soc.rds"
META="SM2v2_indep/step0_qc/${PREFIX}.minDepth200.updated_metadata_v6.txt"
OUTDIR="SM2v2_indep/step1_cluster_opt/${PREFIX}"
mkdir -p "${OUTDIR}/plots"

[[ -s "$SOC"  ]] || { echo "ERROR: missing SocObj: $SOC";   exit 1; }
[[ -s "$META" ]] || { echo "ERROR: missing v6 metadata: $META"; exit 1; }

echo " - Step 1a grid scan: ${PREFIX} (${GENOME}) | min.c=${MINC} | out=${OUTDIR}"
Rscript "${SCRIPTS}/part1_indep/2_1_gridscan_pergenome.R" "$SOC" "$META" "$OUTDIR" 1 "$MINC"

echo " - done. Inspect ${OUTDIR}/UMAP_grid_scan*.best.tsv for (pcs,k_near,min_dist);"
echo "   2_2_resolution_scan.sh reads that best.tsv automatically for Step 1b."
