#!/bin/bash
#SBATCH --job-name=syn_reassign
#SBATCH --partition=standard
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=16
#SBATCH --mem=32G
#SBATCH --time=06:00:00
#SBATCH --output=_logs/03_09_rerun_assign_%A_%a.log
#SBATCH --array=0-29

# =============================================================================
# 03_09_rerun_assign.sh — Full re-run (extract→filter→chunks→assign)
#                                    on Track B + Track B-disc
#
# Re-runs the full ambientmapper pipeline (steps 1-4) so that:
#   1) extract produces frag_loc column (needed by friend rescue)
#   2) assign generates the rescued flag via _friend_rescue
#
# Array 0-14:  Track B (synthetic/)       — 3,831 peaks, 3 genomes
# Array 15-29: Track B-disc (synthetic_disc/) — 843 peaks, 3 genomes
#
# Cleanup: removes old qc/, filtered_QCFiles/, cell_map_ref_chunks/,
#          ExplorationReadLevel/, raw_cell_map_ref_chunks/, _sentinels/,
#          and genotyping_runs/ from previous runs.
#          Keeps: BAMs (mapping/), FASTQs (barcoded/), orthologs/, reads/
#
# Before submitting:
#   pip install .   (from ambientmapper repo)
# =============================================================================

set -euo pipefail

REPO_ROOT="$(git rev-parse --show-toplevel 2>/dev/null || pwd)"   # this repository (sourced helpers, sibling scripts)
PROJECT_ROOT="${PROJECT_ROOT:-$(pwd)}"   # root of the data tree (3_Mapping/, 5_AmbientDetection/, ...)

cd "${PROJECT_ROOT}/5_AmbientDetection"

# conda activate <env from environment.yml>

export OMP_NUM_THREADS=1
export MKL_NUM_THREADS=1
export OPENBLAS_NUM_THREADS=1
export BLIS_NUM_THREADS=1
export NUMEXPR_NUM_THREADS=1
export GOTO_NUM_THREADS=1

DATASETS=(
    alpha_000
    alpha_002_Il14H alpha_005_Il14H alpha_010_Il14H alpha_020_Il14H
    alpha_030_Il14H alpha_040_Il14H alpha_050_Il14H
    alpha_002_Ki11  alpha_005_Ki11  alpha_010_Ki11  alpha_020_Ki11
    alpha_030_Ki11  alpha_040_Ki11  alpha_050_Ki11
)

GENOMES=(B73 Il14H Ki11)
GENOME_KEYS=(B73v5 Il14H Ki11)  # BAM naming: syn_{key}.sorted.bam
TASK=${SLURM_ARRAY_TASK_ID}

if [ "${TASK}" -lt 15 ]; then
    TRACK="trackB"
    BASEDIR="synthetic"
    DS_IDX=${TASK}
else
    TRACK="trackB_disc"
    BASEDIR="synthetic_disc"
    DS_IDX=$((TASK - 15))
fi

DATASET="${DATASETS[$DS_IDX]}"
WORKDIR_ABS="${PROJECT_ROOT}/5_AmbientDetection/${BASEDIR}"
DS_DIR="${BASEDIR}/${DATASET}"
MAPPING_DIR="${BASEDIR}/mapping/${DATASET}"
THREADS=16

echo "============================================"
echo "[$(date)] Re-run pipeline: ${TRACK} / ${DATASET}"
echo "  Task:       ${TASK}"
echo "  Basedir:    ${BASEDIR}"
echo "  Dataset:    ${DATASET}"
echo "  Mapping:    ${MAPPING_DIR}"
echo "============================================"

# =====================================================================
# STEP 0: Clean old outputs
# =====================================================================
echo "[$(date)] Cleaning old outputs..."

# Directories to remove entirely (will be regenerated)
for d in qc filtered_QCFiles cell_map_ref_chunks raw_cell_map_ref_chunks \
         ExplorationReadLevel _sentinels final genotyping_runs; do
    if [ -d "${DS_DIR}/${d}" ]; then
        echo "  rm -rf ${DS_DIR}/${d}/"
        rm -rf "${DS_DIR}/${d}"
    fi
done

# Remove old config (will be regenerated)
rm -f "${DS_DIR}/config.json"

echo "[$(date)] Cleanup done"

# =====================================================================
# STEP 1: Create config JSON
# =====================================================================
CONFIG="${DS_DIR}/config.json"
BAM_DIR="${WORKDIR_ABS}/../${MAPPING_DIR}"

# Verify BAMs exist
for GIDX in 0 1 2; do
    GK="${GENOME_KEYS[$GIDX]}"
    BAM="${BAM_DIR}/syn_${GK}.sorted.bam"
    if [ ! -f "${BAM}" ]; then
        echo "ERROR: BAM not found: ${BAM}"
        exit 1
    fi
done

mkdir -p "${DS_DIR}"
cat > "${CONFIG}" <<CONFIGEOF
{
    "sample": "${DATASET}",
    "workdir": "${WORKDIR_ABS}",
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

# =====================================================================
# STEP 2: Custom extract (produces frag_loc from BAMs)
# =====================================================================
# Uses 03_07_extract_metrics.py instead of ambientmapper extract because
# synthetic read names use {barcode}|{genome}|{art_id} format (no CB tag).
echo ""
echo "[$(date)] Step 1: extract (custom, with frag_loc)..."
QC_DIR="${WORKDIR_ABS}/${DATASET}/qc"
mkdir -p "${QC_DIR}"

python3 "${REPO_ROOT}/workflows/03_genotyping/synthetic/03_07_extract_metrics.py" \
    --bam-dir "${MAPPING_DIR}" \
    --outdir "${QC_DIR}" \
    --genome-map B73v5:B73 Il14H:Il14H Ki11:Ki11

echo "[$(date)] Done: extract"

# =====================================================================
# STEP 3: Copy QCMapping to filtered_QCFiles/ + Filter via chunks
# =====================================================================
# Synthetic barcodes don't need frequency filtering — just copy and chunk.
echo ""
echo "[$(date)] Step 2: copy QCMapping to filtered_QCFiles/..."
FILT_DIR="${WORKDIR_ABS}/${DATASET}/filtered_QCFiles"
mkdir -p "${FILT_DIR}"
for G in "${GENOMES[@]}"; do
    cp "${QC_DIR}/${G}_QCMapping.txt" "${FILT_DIR}/filtered_${G}_QCMapping.txt"
    echo "  ${G}: $(wc -l < "${FILT_DIR}/filtered_${G}_QCMapping.txt") lines"
done
echo "[$(date)] Done: copy"

# =====================================================================
# STEP 4: Chunks
# =====================================================================
echo ""
echo "[$(date)] Step 3: chunks..."
ambientmapper chunks --config "${CONFIG}" --chunk-size-cells 100
echo "[$(date)] Done: chunks"

# =====================================================================
# STEP 5: Assign (with friend rescue + rescued flag)
# =====================================================================
echo ""
echo "[$(date)] Step 4: assign..."
ambientmapper assign \
    --config "${CONFIG}" \
    --threads "${THREADS}" \
    --alpha 0.05 \
    --k 10 \
    --mapq-min 10 \
    --xa-max 2 \
    --chunksize 500000 \
    --batch-size 6

echo "[$(date)] Done: assign"

# =====================================================================
# VERIFY
# =====================================================================
echo ""
echo "=== Verification ==="
CHUNK_DIR="${WORKDIR_ABS}/${DATASET}/cell_map_ref_chunks"
N_FILTERED=$(ls "${CHUNK_DIR}"/*_filtered.tsv.gz 2>/dev/null | wc -l)
echo "  Filtered chunks: ${N_FILTERED}"

# Check for rescued flag
python3 -c "
import pandas as pd, glob
files = sorted(glob.glob('${CHUNK_DIR}/*_filtered.tsv.gz'))
if files:
    df = pd.read_csv(files[0], sep='\t', compression='gzip', nrows=1000)
    print('  Columns:', list(df.columns))
    vc = df['assigned_class'].value_counts()
    print('  assigned_class:', dict(vc))
    has_frag = 'frag_loc' in df.columns
    has_rescued = 'rescued' in vc.index
    print('  frag_loc: ' + ('PRESENT' if has_frag else 'MISSING'))
    print('  rescued flag: ' + ('PRESENT' if has_rescued else 'not in sample (may appear in other chunks)'))
else:
    print('  ERROR: no filtered files found')
"

echo ""
echo "============================================"
echo "[$(date)] Done: ${TRACK} / ${DATASET}"
echo "============================================"
