#!/bin/bash
########## BATCH Lines for Resource Request ##########
#SBATCH --time=3:00:00
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=1
#SBATCH --mem=20G
#SBATCH --job-name=CorrectChunk_ATAC
#SBATCH --partition=standard
#SBATCH --output=_logs/process_chunk_atac_%A_%a.log
#SBATCH --array=1-120
#
# Legacy chain, step 4 (SM2 only): per chunk (60 chunk pairs = 120 files), correct the Tn5 well barcode
# (each 5 bp half within one mismatch, N counts as a mismatch, closest well of config/96well_Tn5_bc_layout.txt),
# assign the read to a sample (SM2_B73 / SM2_MUDR / SM2_At) by the plate design config/PlateDesign_SM2_ATAC.legacy_demux.txt
# and rewrite the corrected barcode into the read name (2_1_1_fastq_sample_assign_and_fixTn5.ATAC.py).
# Usage (from ${PROJECT_ROOT}/2_CleanReads): sbatch 2_1_1_fastq_sample_assign_and_fixTn5.ATAC.sh SM2_ATAC

: "${PROJECT_ROOT:?set PROJECT_ROOT to the analysis tree that holds 2_CleanReads/}"
REPO_ROOT="${REPO_ROOT:-$(git rev-parse --show-toplevel 2>/dev/null || pwd)}"
# conda activate <env from environment.yml>   (python 3.8+)
cd "${PROJECT_ROOT}/2_CleanReads"
mkdir -p _logs

########## Define Variables #########
# Read SAMPLE from command-line argument.
SAMPLE="$1"
if [ -z "$SAMPLE" ]; then
    echo "Usage: $0 SAMPLE"
    exit 1
fi

# Read all FASTQ files for the sample into an array.
readarray -t SAMPLES < <(find "${SAMPLE}_chunks" -name "*fastq.gz" -exec basename {} \; | sort -u)
if [ ${#SAMPLES[@]} -eq 0 ]; then
    echo "Error: No FASTQ result files found in ${SAMPLE}_chunks/"
    exit 1
fi

# Compute the index for the current array task (subtracting 1 because SLURM_ARRAY_TASK_ID is 1-indexed)
INDEX=$((SLURM_ARRAY_TASK_ID - 1))

# Check that the index is within the bounds of the SAMPLES array.
TOTAL_SAMPLES=${#SAMPLES[@]}
if [ $INDEX -ge $TOTAL_SAMPLES ]; then
    echo "Error: Array index $INDEX exceeds total available files ($TOTAL_SAMPLES)"
    exit 1
fi

# Get the file name corresponding to this job's array task.
FILE="${SAMPLES[$INDEX]}"

# Format the chunk number with leading zeros for logging consistency.
CHUNK_NUM=$(printf "%03d" "$SLURM_ARRAY_TASK_ID")

# Barcode layout and plate design (sample -> well ranges).
LAYOUT="${REPO_ROOT}/config/96well_Tn5_bc_layout.txt"
DESIGN="${DESIGN:-${REPO_ROOT}/config/PlateDesign_SM2_ATAC.legacy_demux.txt}"
ASSIGN_PY="${REPO_ROOT}/workflows/01_preprocessing/legacy_sm2_root1/2_1_1_fastq_sample_assign_and_fixTn5.ATAC.py"

# Define input (split) and output (corrected) directories.
SPLIT_OUTPUT_DIR="${SAMPLE}_chunks"
CORRECTED_OUTPUT_DIR="${SAMPLE}_chunks_Corrected"

# Create the output directory if it doesn't already exist.
mkdir -p "${CORRECTED_OUTPUT_DIR}"

echo "[$(date)] Starting processing for chunk ${CHUNK_NUM} with file: ${FILE}"

########## Function to process a single sample file ##########
process_sample() {
    local file="$1"
    echo "[$(date)] Starting barcode correction for file: $file (chunk ${CHUNK_NUM})"

    python "${ASSIGN_PY}" "${LAYOUT}" "${SPLIT_OUTPUT_DIR}/${file}" "${DESIGN}" --output "${CORRECTED_OUTPUT_DIR}"

    if [ $? -ne 0 ]; then
        echo "[$(date)] Error processing file: $file"
    else
        echo "[$(date)] Processing completed for file: $file"
    fi
}

# Process the single file for this array task.
process_sample "${FILE}"

echo "[$(date)] Processing of chunk ${CHUNK_NUM} completed."
