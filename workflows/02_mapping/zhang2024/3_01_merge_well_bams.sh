#!/usr/bin/env bash
#SBATCH --time=4:00:00
#SBATCH --nodes=1
#SBATCH --cpus-per-task=8
#SBATCH --mem=16G
#SBATCH --job-name=merge_wells
#SBATCH --partition=standard
#SBATCH --output=_logs/merge_%j.log

# =============================================================================
# 3_01_merge_well_bams.sh
#
# Merge per-well BAMs (from step2) into replicate-level BAMs per genome.
# Output goes to 3_Mapping/ambientmapper_input/ for AmbientMapper.
#
# Usage (from ${PROJECT_ROOT}/3_Mapping; LISTS = <repo>/workflows/02_mapping/zhang2024):
#   sbatch 3_01_merge_well_bams.sh B73Mo17_rep1        ${LISTS}/Genome_list_scifi_B73_Mo17
#   sbatch 3_01_merge_well_bams.sh B73Mo17_rep2        ${LISTS}/Genome_list_scifi_B73_Mo17
#   sbatch 3_01_merge_well_bams.sh multiGenotypes_rep1  ${LISTS}/Genome_list_MultipleGenome
#
# Input:  3_Mapping/{DATASET_REP}/{WELL}/{WELL}_{GENOME}_scifiATAC.mq10.BC.rmdup.mm.bam
# Output: 3_Mapping/ambientmapper_input/{DATASET_REP}_{GENOME}_scifiATAC.mq10.BC.rmdup.mm.bam
# =============================================================================

set -euo pipefail
: "${PROJECT_ROOT:?set PROJECT_ROOT to the analysis tree that holds 3_Mapping/}"

########## Environment #########
# conda activate <env from environment.yml>   (samtools)

########## Arguments #########
DATASET_REP="${1:?Usage: sbatch 3_01_merge_well_bams.sh <DATASET_REP> <GENOME_LIST>}"
GENOME_LIST="${2:?Missing GENOME_LIST argument}"

########## Paths #########
BASEDIR="${PROJECT_ROOT}"
MAPPING_DIR="${BASEDIR}/3_Mapping"
WELL_DIR="${MAPPING_DIR}/${DATASET_REP}"
OUTDIR="${MAPPING_DIR}/ambientmapper_input"
THREADS=8
MAPQ_MIN=10

cd "${MAPPING_DIR}"
mkdir -p "${OUTDIR}" _logs

########## Merge per genome #########
while IFS= read -r GENOME; do
  [[ -z "$GENOME" ]] && continue

  OUT_BAM="${OUTDIR}/${DATASET_REP}_${GENOME}_scifiATAC.mq${MAPQ_MIN}.BC.rmdup.mm.bam"

  if [[ -f "$OUT_BAM" ]]; then
    echo "[$(date)] Skipping ${GENOME}: ${OUT_BAM} already exists"
    continue
  fi

  # Collect per-well BAMs for this genome
  BAM_PATTERN="*_${GENOME}_scifiATAC.mq${MAPQ_MIN}.BC.rmdup.mm.bam"
  BAM_LIST=$(mktemp)
  find "${WELL_DIR}" -name "$BAM_PATTERN" -type f | sort > "$BAM_LIST"
  N_BAMS=$(wc -l < "$BAM_LIST")

  if [[ "$N_BAMS" -eq 0 ]]; then
    echo "[WARN] No BAMs found for genome=${GENOME} in ${WELL_DIR}"
    rm "$BAM_LIST"
    continue
  fi

  echo "[$(date)] Merging ${N_BAMS} well BAMs for genome=${GENOME} -> ${OUT_BAM}"

  # Use a file list to avoid E2BIG with 96 BAM paths
  samtools merge -@ "${THREADS}" -b "$BAM_LIST" "$OUT_BAM"
  samtools index -@ "${THREADS}" "$OUT_BAM"

  rm "$BAM_LIST"
  echo "[$(date)] Done: ${OUT_BAM} ($(du -h "$OUT_BAM" | cut -f1))"

done < "${GENOME_LIST}"

echo "[$(date)] All genomes merged for ${DATASET_REP}."
