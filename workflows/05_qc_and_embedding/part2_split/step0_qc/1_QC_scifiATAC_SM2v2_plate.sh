#!/bin/bash
########## BATCH Lines for Resource Request ##########
#SBATCH --time=3:00:00
#SBATCH --nodes=1
#SBATCH --ntasks-per-node=1
#SBATCH --cpus-per-task=10
#SBATCH --mem=40G
#SBATCH --job-name=Step1_SocQC_SM2v2_plate
#SBATCH --partition=standard
#SBATCH --output=_logs/Step1_SocQC_SM2v2_plate_%A_%a.log
#SBATCH --array=0-3
#
# 1_QC_scifiATAC_SM2v2_plate.sh
# ---------------------------------------------------------------------------
# PLATE-SPLIT ("normal user") Socrates QC for SM2v2 — the deployment-scenario
# counterpart to 1_QC_scifiATAC_SM2v2_indep.sh.
#
# THE DESIGN (and how it differs from SM2v2_indep — read this before comparing)
#   SM2v2_indep maps EVERY barcode to BOTH references. That is AmbientMapper's
#   substrate, not a user's analysis. Here each plate goes to its EXPECTED
#   genome only, which is what a scifi-ATAC user actually runs:
#       At  plate -> TAIR10          B73 plate -> B73v5
#   Pre vs Post-AM on that design measures what cleaning buys the normal
#   downstream workflow.
#
#   *** THIS DOES NOT REPLACE SM2v2_indep. ***
#   Cross-plate barcodes are excluded a priori here, so the dominant cleaning
#   effect (maize object sheds 8,775 -> 506 At-plate barcodes) is invisible BY
#   CONSTRUCTION, and the Step-4 cl5 contamination call cannot be reproduced --
#   cl5 *is* that At-plate-in-maize population. Steps 1-4 stay on SM2v2_indep.
#   What this arm measures instead is within-plate, read-level cleaning.
#
# WHY THE isCell REFIT MATTERS (the point of the arm, not a side effect)
#   In SM2v2_indep the TAIR10 object is 80% B73-plate barcodes, so its isCellv2
#   / v6 QC model is fit largely on the WRONG plate. Re-fitting on At-plate-only
#   data is the correct model for an At sample -- and it means the mixed-fit
#   cell counts (At 597 Pre / 770 Post-wd) will NOT carry over unchanged.
#   Expect them to move; that is the corrected number, not a bug.
#
# Array task layout (Pre + Post-nd + Post-wd x {At->TAIR10, B73->B73v5}):
#   0: SM2_At_TAIR10             PreClean   TAIR10
#   1: SM2_B73_B73v5             PreClean   B73v5
#   2: Clean.SM2v2wd_At_TAIR10   PostClean  TAIR10  (wd = design-guided; HEADLINE)
#   3: Clean.SM2v2wd_B73_B73v5   PostClean  B73v5   (wd)
#   4: Clean.SM2v2_At_TAIR10     PostClean  TAIR10  (nd = design-free; sensitivity)
#   5: Clean.SM2v2_B73_B73v5     PostClean  B73v5   (nd)
# Default --array=0-3 runs Pre + Post-wd (the headline pair, matching the
# Step-3 plan's choice). The nd sensitivity arm:
#   sbatch --array=4,5 0_scripts/part2_split/step0_qc/1_QC_scifiATAC_SM2v2_plate.sh
#
# Prereq: sbatch 0_scripts/part2_split/step0_qc/0_09_split_beds_by_plate.sh
# Usage:  sbatch 0_scripts/part2_split/step0_qc/1_QC_scifiATAC_SM2v2_plate.sh [OUTDIR]
#         (default OUTDIR: SM2v2_plate). Submit from the 6_socrates dir.
# ---------------------------------------------------------------------------

# conda activate ambientmapper-manuscript   (environment.yml at the repo root)
# module load macs2   (site-specific; Socrates' .preRunChecks() needs macs2 on PATH even with precomputed peaks)
set -euo pipefail

OUTDIR="${1:-SM2v2_plate}"

# Paths are anchored on PROJECT_ROOT ($BASH_SOURCE self-location does not survive sbatch).
PROJECT_ROOT="${PROJECT_ROOT:-.}"
BASE="${PROJECT_ROOT}/6_socrates"
SCRIPTS="${SCRIPTS:-${BASE}/0_scripts}"
INPUTDIR="${BASE}/_data/_BED_files"
GENOME_INFO="${BASE}/_data/_GenomeInfo"

step_1_1="${SCRIPTS}/common/1_1_QC_scifiATAC_data.R"
step_1_2="${SCRIPTS}/common/1_2_filter_lowQC_cells_scifiATAC_data.R"
step_1_3="${SCRIPTS}/common/1_3_metaQC_scifiATAC_data.R"

# SM2 canonical depth (matches the Fig-1 lineage / the indep arm). NOT 500:
# the At object is low-coverage and 500 over-filters it. Keeping DEPTH=200
# identical to SM2v2_indep is what makes the two arms comparable.
DEPTH=200

# Per-genome annotation, chr_sizes, MACS genome size, organelle scaffolds.
# Identical registry to the indep arm -- same references, same contig naming
# (the SM2 B73v5 mapping names organelles Mt:AGPv4/Pt:AGPv4).
declare -A SAMPLE_ANN SAMPLE_CHR SAMPLE_GSIZE SAMPLE_ORG
SAMPLE_ANN[B73v5]="${GENOME_INFO}/B73v5.gff3"   # original NAM5 gff3; the *.gtf variants have no gene_id -> TxDb fails
SAMPLE_CHR[B73v5]="${GENOME_INFO}/SM2v2_B73v5.chrs.size.txt"
SAMPLE_GSIZE[B73v5]="2.1e9"
SAMPLE_ORG[B73v5]="Mt:AGPv4,Pt:AGPv4"

SAMPLE_ANN[TAIR10]="${GENOME_INFO}/TAIR10.gff3"          # Ensembl TAIR10 r60 (chroms 1-5/Mt/Pt)
SAMPLE_CHR[TAIR10]="${GENOME_INFO}/SM2v2_TAIR10.chrs.size.txt"
SAMPLE_GSIZE[TAIR10]="0.12e9"
SAMPLE_ORG[TAIR10]="Mt,Pt"

mkdir -p "${BASE}/${OUTDIR}/step0_qc/plots"
cd "${BASE}/${OUTDIR}/step0_qc"

runQC() {
  local name="$1"
  local species
  species="${name#Clean.}"        # strip optional PostClean prefix
  species="${species##*_}"        # trailing token -> B73v5 or TAIR10

  local ann="${SAMPLE_ANN[$species]:?No annotation registered for species=$species}"
  local chr="${SAMPLE_CHR[$species]:?No chr_sizes registered for species=$species}"
  local gsize="${SAMPLE_GSIZE[$species]:?No genomesize registered for species=$species}"
  local org="${SAMPLE_ORG[$species]:?No org_scaffolds registered for species=$species}"
  local bed="${INPUTDIR}/${name}_scifiATAC.mq10.BC.rmdup.mm.tn5.bed.gz"

  [[ -f "$bed" ]] || { echo "ERROR: missing BED $bed  (run 0_09_split_beds_by_plate.sh first)"; exit 2; }
  [[ -f "$ann" ]] || { echo "ERROR: missing annotation $ann"; exit 2; }
  [[ -f "$chr" ]] || { echo "ERROR: missing chr_sizes $chr  (run 0_08_setup_SM2v2_indep.sh first)"; exit 2; }

  echo " - 1. Socrates build + isCellv2 for $name  (species=$species, gsize=$gsize, org=$org)"
  Rscript "$step_1_1" "$bed" "$name" "$ann" "$chr" "$gsize" "$org"

  echo " - 2. QC values for $name"
  Rscript "$step_1_2" "${name}.raw.soc.rds" "$name"

  echo " - 3. metaQC (depth=${DEPTH}) for $name"
  Rscript "$step_1_3" "${name}.updated_metadata.txt" "$name" "${DEPTH}"
}

samples_list=(
  SM2_At_TAIR10             # 0  PreClean   TAIR10
  SM2_B73_B73v5             # 1  PreClean   B73v5
  Clean.SM2v2wd_At_TAIR10   # 2  Post-wd    TAIR10  (headline)
  Clean.SM2v2wd_B73_B73v5   # 3  Post-wd    B73v5   (headline)
  Clean.SM2v2_At_TAIR10     # 4  Post-nd    TAIR10  (sensitivity)
  Clean.SM2v2_B73_B73v5     # 5  Post-nd    B73v5   (sensitivity)
)

sample="${samples_list[$SLURM_ARRAY_TASK_ID]:-}"
[[ -n "$sample" ]] || { echo "ERROR: no sample for array index $SLURM_ARRAY_TASK_ID"; exit 1; }

echo "$(date): SM2v2-plate QC on $sample  (OUTDIR=${OUTDIR}, DEPTH=${DEPTH})"
runQC "$sample"
echo "$(date): DONE $sample"
