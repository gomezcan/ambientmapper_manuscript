#!/bin/bash
#SBATCH --job-name=syn_pipeline
#SBATCH --partition=standard
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=16
#SBATCH --mem=32G
#SBATCH --time=08:00:00
#SBATCH --output=_logs/03_08_run_pipeline_%A_%a.log
#SBATCH --array=0-14

# =============================================================================
# 03_08_run_pipeline.sh — Phase 6: AmbientMapper filter+chunks+assign
#
# For each dataset: copies QCMapping to filtered_QCFiles/ format, runs
# ambientmapper chunks and assign. Skips extract (done in Phase 5) and
# filter (synthetic barcodes don't need frequency filtering).
#
# Requires: Phase 5 output (synthetic/{dataset}/qc/{genome}_QCMapping.txt)
# Output:   synthetic/{dataset}/cell_map_ref_chunks/*_filtered.tsv.gz
# =============================================================================

PROJECT_ROOT="${PROJECT_ROOT:-$(pwd)}"   # root of the data tree (3_Mapping/, 5_AmbientDetection/, ...)
cd "${PROJECT_ROOT}/5_AmbientDetection"

# conda activate <env from environment.yml>

export OMP_NUM_THREADS=1
export MKL_NUM_THREADS=1
export OPENBLAS_NUM_THREADS=1
export BLIS_NUM_THREADS=1
export NUMEXPR_NUM_THREADS=1
export GOTO_NUM_THREADS=1

# --- Dataset list ---
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

DATASET="${DATASETS[$SLURM_ARRAY_TASK_ID]}"
GENOMES=(B73 Il14H Ki11)
WORKDIR="${PROJECT_ROOT}/5_AmbientDetection/synthetic"
SAMPLE="${DATASET}"
THREADS=16

if [ -z "${DATASET}" ]; then
    echo "Error: Invalid task ID ${SLURM_ARRAY_TASK_ID}"
    exit 1
fi

echo "============================================"
echo "[$(date)] Phase 6: AmbientMapper pipeline (chunks + assign)"
echo "  Task:    ${SLURM_ARRAY_TASK_ID}"
echo "  Dataset: ${DATASET}"
echo "  Workdir: ${WORKDIR}"
echo "  Sample:  ${SAMPLE}"
echo "============================================"

# --- Step 1: Copy QCMapping to filtered_QCFiles/ ---
# Phase 5 wrote to: synthetic/{dataset}/qc/{genome}_QCMapping.txt (with header)
# ambientmapper assign reads: filtered_QCFiles/filtered_{genome}_QCMapping.txt (with header)
# We just rename/copy to match the expected naming.

FILT_DIR="${WORKDIR}/${SAMPLE}/filtered_QCFiles"
mkdir -p "${FILT_DIR}"

echo "[$(date)] Copying QCMapping files to filtered_QCFiles/..."
for G in "${GENOMES[@]}"; do
    SRC="${WORKDIR}/${SAMPLE}/qc/${G}_QCMapping.txt"
    DST="${FILT_DIR}/filtered_${G}_QCMapping.txt"
    if [ ! -f "${SRC}" ]; then
        echo "Error: QCMapping not found: ${SRC}"
        exit 1
    fi
    cp "${SRC}" "${DST}"
    echo "  ${SRC} -> ${DST} ($(wc -l < "${DST}") lines)"
done

# --- Step 2: Create config JSON for this dataset ---
CONFIG="${WORKDIR}/${SAMPLE}/config.json"
BAM_DIR="${WORKDIR}/../synthetic/mapping/${DATASET}"

cat > "${CONFIG}" <<CONFIGEOF
{
    "sample": "${SAMPLE}",
    "workdir": "${WORKDIR}",
    "min_barcode_freq": 3,
    "chunk_size_cells": 100,
    "genomes": {
        "B73": "${BAM_DIR}/syn_B73v5.sorted.bam",
        "Il14H": "${BAM_DIR}/syn_Il14H.sorted.bam",
        "Ki11": "${BAM_DIR}/syn_Ki11.sorted.bam"
    }
}
CONFIGEOF

echo "[$(date)] Config: ${CONFIG}"

# --- Step 3: Run chunks ---
echo ""
echo "[$(date)] Running ambientmapper chunks..."
ambientmapper chunks \
    --config "${CONFIG}" \
    --chunk-size-cells 100

if [ $? -ne 0 ]; then
    echo "Error: ambientmapper chunks failed"
    exit 1
fi

CHUNK_DIR="${WORKDIR}/${SAMPLE}/cell_map_ref_chunks"
N_CHUNKS=$(ls "${CHUNK_DIR}"/*_cell_map_ref_chunk_*.txt 2>/dev/null | wc -l)
echo "  Created ${N_CHUNKS} chunk files in ${CHUNK_DIR}/"

# --- Step 4: Run assign ---
echo ""
echo "[$(date)] Running ambientmapper assign..."
ambientmapper assign \
    --config "${CONFIG}" \
    --threads "${THREADS}" \
    --alpha 0.05 \
    --k 10 \
    --mapq-min 10 \
    --xa-max 2 \
    --chunksize 500000 \
    --batch-size 6

if [ $? -ne 0 ]; then
    echo "Error: ambientmapper assign failed"
    exit 1
fi

# --- Summary ---
echo ""
echo "=== Phase 6 summary ==="
N_FILTERED=$(ls "${CHUNK_DIR}"/*_filtered.tsv.gz 2>/dev/null | wc -l)
echo "  Filtered chunks: ${N_FILTERED}"
echo "  Chunk dir: ${CHUNK_DIR}/"

echo ""
echo "============================================"
echo "[$(date)] Phase 6 complete for ${DATASET}"
echo "============================================"
