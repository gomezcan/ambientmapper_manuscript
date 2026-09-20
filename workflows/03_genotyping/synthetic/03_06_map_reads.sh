#!/bin/bash
#SBATCH --job-name=bwa_synth
#SBATCH --partition=standard
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=8
#SBATCH --mem=16G
#SBATCH --time=02:00:00
#SBATCH --output=_logs/03_06_map_reads_%A_%a.log
#SBATCH --array=0-44

# =============================================================================
# 03_06_map_reads.sh — Phase 4: BWA mapping for synthetic benchmark
#
# Maps each of the 15 titration datasets against 3 reference genomes.
# SLURM array 0-44: task = dataset_index * 3 + genome_index
#
# Requires: Phase 3 output (synthetic/barcoded/{dataset}/all_reads_{R1,R2}.fq.gz)
# Output:   synthetic/mapping/{dataset}/syn_{genome}.sorted.bam
# =============================================================================

PROJECT_ROOT="${PROJECT_ROOT:-$(pwd)}"   # root of the data tree (3_Mapping/, 5_AmbientDetection/, ...)
cd "${PROJECT_ROOT}/5_AmbientDetection"

# conda activate <env from environment.yml>
module load Bioinformatics bwa/0.7.17-mil4ns7

# --- Config ---
# BWA indexes of the NAM assemblies (Index_Zm_<genome>_bwa). BWA_INDEX_ROOT holds the Zea/ tree.
BWA_INDEX_ROOT="${BWA_INDEX_ROOT:?set BWA_INDEX_ROOT to the directory holding Zea/NAN_Indexes/}"
BWA_INDEX_DIR="${BWA_INDEX_ROOT}/Zea/NAN_Indexes"
BARCODED_DIR="synthetic/barcoded"
MAP_DIR="synthetic/mapping"
THREADS=8

# Dataset list (must match Phase 3 output)
DATASETS=(
    alpha_000
    alpha_002_Il14H
    alpha_005_Il14H
    alpha_010_Il14H
    alpha_020_Il14H
    alpha_030_Il14H
    alpha_040_Il14H
    alpha_050_Il14H
    alpha_002_Ki11
    alpha_005_Ki11
    alpha_010_Ki11
    alpha_020_Ki11
    alpha_030_Ki11
    alpha_040_Ki11
    alpha_050_Ki11
)

# Genome list (BWA index names)
GENOMES=(B73v5 Il14H Ki11)

N_DATASETS=${#DATASETS[@]}
N_GENOMES=${#GENOMES[@]}

# --- Decode array task ID ---
TASK_ID=${SLURM_ARRAY_TASK_ID}
DATASET_IDX=$((TASK_ID / N_GENOMES))
GENOME_IDX=$((TASK_ID % N_GENOMES))

DATASET="${DATASETS[$DATASET_IDX]}"
GENOME="${GENOMES[$GENOME_IDX]}"

if [ -z "${DATASET}" ] || [ -z "${GENOME}" ]; then
    echo "Error: Invalid task ID ${TASK_ID} (dataset_idx=${DATASET_IDX}, genome_idx=${GENOME_IDX})"
    exit 1
fi

# --- Paths ---
R1="${BARCODED_DIR}/${DATASET}/all_reads_R1.fq.gz"
R2="${BARCODED_DIR}/${DATASET}/all_reads_R2.fq.gz"
REF="${BWA_INDEX_DIR}/Index_Zm_${GENOME}_bwa"
OUTDIR="${MAP_DIR}/${DATASET}"
BAM="${OUTDIR}/syn_${GENOME}.sorted.bam"

mkdir -p "${OUTDIR}"

echo "============================================"
echo "[$(date)] Phase 4: BWA mapping"
echo "  Task:    ${TASK_ID} (dataset ${DATASET_IDX}/${N_DATASETS}, genome ${GENOME_IDX}/${N_GENOMES})"
echo "  Dataset: ${DATASET}"
echo "  Genome:  ${GENOME}"
echo "  R1:      ${R1}"
echo "  Ref:     ${REF}"
echo "  Output:  ${BAM}"
echo "============================================"

# --- Verify inputs ---
if [ ! -f "${R1}" ]; then
    echo "Error: FASTQ not found: ${R1}"
    exit 1
fi

# --- BWA mapping ---
echo "[$(date)] Running bwa mem..."

bwa mem -M -t "${THREADS}" "${REF}" "${R1}" "${R2}" \
  | samtools sort -@ 4 -m 2G -o "${BAM}" -

if [ $? -ne 0 ]; then
    echo "Error: BWA/samtools failed"
    exit 1
fi

echo "[$(date)] Indexing BAM..."
samtools index "${BAM}"

# --- Stats ---
TOTAL=$(samtools view -c "${BAM}")
MAPPED=$(samtools view -c -F 4 "${BAM}")
echo ""
echo "=== Mapping stats ==="
echo "  Total alignments:  ${TOTAL}"
echo "  Mapped:            ${MAPPED}"
echo "  BAM size:          $(ls -lh "${BAM}" | awk '{print $5}')"

echo ""
echo "============================================"
echo "[$(date)] Done: ${DATASET} x ${GENOME}"
echo "============================================"
