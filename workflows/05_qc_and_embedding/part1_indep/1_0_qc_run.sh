#!/bin/bash
########## BATCH Lines for Resource Request ##########
#SBATCH --time=3:00:00
#SBATCH --nodes=1
#SBATCH --ntasks-per-node=1
#SBATCH --cpus-per-task=10
#SBATCH --mem=40G
#SBATCH --job-name=Step1_SocQC_SM2v2_indep
#SBATCH --partition=standard
#SBATCH --output=_logs/Step1_SocQC_SM2v2_indep_%A_%a.log
#SBATCH --array=0-3
#
# 1_0_qc_run.sh
# ---------------------------------------------------------------------------
# SPLIT-genome (independent-mapping) Socrates QC for SM2v2 — the per-genome
# counterpart to the combined ZmATcombined chunked run. This is the substrate
# for the biological-impact (PreClean vs PostClean) analysis: each genome mapped
# to itself, standard single-genome Socrates, run before and after cleaning.
#
# Array task layout (Pre + Post-nd + Post-wd x {B73v5, TAIR10}):
#   0: SM2_B73v5             PreClean   B73v5
#   1: SM2_TAIR10            PreClean   TAIR10
#   2: Clean.SM2v2_B73v5     PostClean  B73v5   (nd = design-free decontam)
#   3: Clean.SM2v2_TAIR10    PostClean  TAIR10  (nd)
#   4: Clean.SM2v2wd_B73v5   PostClean  B73v5   (wd = design-guided; sensitivity)
#   5: Clean.SM2v2wd_TAIR10  PostClean  TAIR10  (wd)
# The default #SBATCH --array=0-3 runs Pre+Post-nd; submit the
# wd fast-follow explicitly:  sbatch --array=4,5 0_scripts/part1_indep/1_0_qc_run.sh
#
# Prereq: bash 0_scripts/part1_indep/0_2_setup_indep.sh   (chr_sizes + symlinks)
# Usage:  sbatch 0_scripts/part1_indep/1_0_qc_run.sh [OUTDIR]   (default: SM2v2_indep)
#         Submit from the 6_socrates dir so _logs/ resolves.
# ---------------------------------------------------------------------------

# conda activate ambientmapper-manuscript   (environment.yml at the repo root)
# module load macs2   (site-specific; Socrates' .preRunChecks() needs macs2 on PATH even with precomputed peaks)
set -euo pipefail

OUTDIR="${1:-SM2v2_indep}"

# Paths are anchored on PROJECT_ROOT ($BASH_SOURCE self-location does not survive sbatch).
PROJECT_ROOT="${PROJECT_ROOT:-.}"
BASE="${PROJECT_ROOT}/6_socrates"
SCRIPTS="${SCRIPTS:-${BASE}/0_scripts}"
INPUTDIR="${BASE}/_data/_BED_files"
GENOME_INFO="${BASE}/_data/_GenomeInfo"

step_1_1="${SCRIPTS}/part1_indep/1_1_qc_build.R"
step_1_2="${SCRIPTS}/part1_indep/1_2_qc_filter.R"
step_1_3="${SCRIPTS}/part1_indep/1_3_qc_metaqc.R"

# SM2 canonical depth (matches the Fig-1 lineage). NOT the B73Mo17 default of 500:
# the TAIR10/At object is low-coverage and 500 over-filters it (the At-minDepth500 collapse).
DEPTH=200

# Per-genome annotation, chr_sizes, MACS genome size, and organelle scaffolds.
# chr_sizes are SM2-specific (the SM2 B73v5 mapping names organelles Mt:AGPv4/Pt:AGPv4).
declare -A SAMPLE_ANN SAMPLE_CHR SAMPLE_GSIZE SAMPLE_ORG
SAMPLE_ANN[B73v5]="${GENOME_INFO}/B73v5.gff3"   # original NAM5 gff3 (proper ID/Parent); the *.gtf variants have GFF3 keys in GTF syntax with NO gene_id -> TxDb fails
SAMPLE_CHR[B73v5]="${GENOME_INFO}/SM2v2_B73v5.chrs.size.txt"
SAMPLE_GSIZE[B73v5]="2.1e9"
SAMPLE_ORG[B73v5]="Mt:AGPv4,Pt:AGPv4"

SAMPLE_ANN[TAIR10]="${GENOME_INFO}/TAIR10.gff3"          # Ensembl TAIR10 r60 (chroms 1-5/Mt/Pt)
SAMPLE_CHR[TAIR10]="${GENOME_INFO}/SM2v2_TAIR10.chrs.size.txt"
SAMPLE_GSIZE[TAIR10]="0.12e9"
SAMPLE_ORG[TAIR10]="Mt,Pt"

mkdir -p "${BASE}/${OUTDIR}/step0_qc/plots"
cd "${BASE}/${OUTDIR}/step0_qc"        # QC objects (raw.soc/metadata/minDepth) live under step0_qc/

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

  [[ -f "$bed" ]] || { echo "ERROR: missing BED $bed"; exit 2; }
  [[ -f "$ann" ]] || { echo "ERROR: missing annotation $ann"; exit 2; }
  [[ -f "$chr" ]] || { echo "ERROR: missing chr_sizes $chr  (run 0_2_setup_indep.sh first)"; exit 2; }

  echo " - 1. Socrates build + isCellv2 for $name  (species=$species, gsize=$gsize, org=$org)"
  Rscript "$step_1_1" "$bed" "$name" "$ann" "$chr" "$gsize" "$org"

  echo " - 2. QC values for $name"
  Rscript "$step_1_2" "${name}.raw.soc.rds" "$name"

  echo " - 3. metaQC (depth=${DEPTH}) for $name"
  Rscript "$step_1_3" "${name}.updated_metadata.txt" "$name" "${DEPTH}"
}

samples_list=(
  SM2_B73v5             # 0  PreClean   B73v5
  SM2_TAIR10            # 1  PreClean   TAIR10
  Clean.SM2v2_B73v5     # 2  Post-nd    B73v5
  Clean.SM2v2_TAIR10    # 3  Post-nd    TAIR10
  Clean.SM2v2wd_B73v5   # 4  Post-wd    B73v5  (sensitivity)
  Clean.SM2v2wd_TAIR10  # 5  Post-wd    TAIR10
)

sample="${samples_list[$SLURM_ARRAY_TASK_ID]:-}"
[[ -n "$sample" ]] || { echo "ERROR: no sample for array index $SLURM_ARRAY_TASK_ID"; exit 1; }

echo "$(date): SM2v2-indep QC on $sample  (OUTDIR=${OUTDIR}, DEPTH=${DEPTH})"
runQC "$sample"
echo "$(date): DONE $sample"
