#!/usr/bin/env bash
#SBATCH --job-name=pipeline_newds
#SBATCH --partition=standard
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=16
#SBATCH --mem=128G
#SBATCH --time=24:00:00
#SBATCH --array=0-2
#SBATCH --output=_logs/04_01_pipeline_%A_%a.log
#
# 04_01_pipeline_B73Mo17_Multi.sh
#
# Full AmbientMapper pipeline (extract → filter ��� chunks → assign) for:
#   Task 0: B73Mo17_rep1  (2 genomes: B73 + Mo17)
#   Task 1: B73Mo17_rep2  (2 genomes: B73 + Mo17)
#   Task 2: multiGenotypes_rep1  (7 genomes: B73, B97, Ky21, M162W, Mo18W, Oh7B, Tzi8)
#
# Assign uses the α=0.05 / k=10 reconciled bundle (same as Root1 Phase 1-4).
# Friend rescue is ON (default). Genotyping is NOT run here — it will be
# configured separately based on the synthetic disc analysis.
#
# Resume-safe: each step checks for existing output before re-running.
#
# =============================================================================

set -euo pipefail

PROJECT_ROOT="${PROJECT_ROOT:-$(pwd)}"   # root of the data tree (3_Mapping/, 5_AmbientDetection/, ...)

cd "${PROJECT_ROOT}/5_AmbientDetection"

# conda activate <env from environment.yml>

export OMP_NUM_THREADS=1
export MKL_NUM_THREADS=1
export OPENBLAS_NUM_THREADS=1
export BLIS_NUM_THREADS=1
export NUMEXPR_NUM_THREADS=1
export GOTO_NUM_THREADS=1

# -----------------------------------------------------------------------------
# Dataset selection by array task ID
# -----------------------------------------------------------------------------
CONFIGS=(
    "configs/B73Mo17_rep1.ambientmapper.json"
    "configs/B73Mo17_rep2.ambientmapper.json"
    "configs/multiGenotypes_rep1.ambientmapper.json"
)
NAMES=(
    "B73Mo17_rep1"
    "B73Mo17_rep2"
    "multiGenotypes_rep1"
)

IDX=${SLURM_ARRAY_TASK_ID}
CONFIG=${CONFIGS[$IDX]}
SAMPLE=${NAMES[$IDX]}
THREADS=16

# Pin ambientmapper version for reproducibility
# AMBIENTMAPPER_REPO: optional local checkout of the ambientmapper tool, used only to log its commit
AM_VERSION=$(pip show ambientmapper 2>/dev/null | awk '/^Version/ {print $2}')
AM_COMMIT=$(git -C "${AMBIENTMAPPER_REPO:-.}" \
    rev-parse --short HEAD 2>/dev/null || echo "unknown")

echo "============================================================"
echo "[$(date)] Pipeline: extract → filter → chunks → assign"
echo "  task          = ${IDX} / ${SAMPLE}"
echo "  config        = ${CONFIG}"
echo "  ambientmapper : v${AM_VERSION} @ ${AM_COMMIT}"
echo "============================================================"

if [[ ! -f "${CONFIG}" ]]; then
    echo "ERROR: config not found: ${CONFIG}" >&2
    exit 1
fi

# -----------------------------------------------------------------------------
# Step 1/4: extract — QCMapping from BAMs
# -----------------------------------------------------------------------------
echo ""
echo "[$(date)] Step 1/4: extract"
ambientmapper extract \
    --config "${CONFIG}" \
    --threads "${THREADS}"

# -----------------------------------------------------------------------------
# Step 2/4: filter — barcode frequency filter
# -----------------------------------------------------------------------------
echo ""
echo "[$(date)] Step 2/4: filter (min_barcode_freq=5)"
ambientmapper filter \
    --config "${CONFIG}" \
    --threads 2 \
    --min-barcode-freq 5

# -----------------------------------------------------------------------------
# Step 3/4: chunks — partition barcodes into chunks
# -----------------------------------------------------------------------------
echo ""
echo "[$(date)] Step 3/4: chunks (chunk_size_cells=100)"
ambientmapper chunks \
    --config "${CONFIG}" \
    --chunk-size-cells 100

# -----------------------------------------------------------------------------
# Step 4/4: assign — score reads at α=0.05 (friend rescue ON)
# -----------------------------------------------------------------------------
echo ""
echo "[$(date)] Step 4/4: assign (α=0.05, k=10, friend rescue ON)"
# score-batch-size/workers: B73Mo17 (2 genomes) can use defaults (200/4).
# multiGenotypes (7 genomes, task 2) is heavier but 128G should handle 200/4.
# If multiGenotypes OOMs, reduce to --score-batch-size 100 --score-workers 2.
ambientmapper assign \
    --config "${CONFIG}" \
    --threads "${THREADS}" \
    --alpha 0.05 \
    --k 10 \
    --mapq-min 10 \
    --xa-max 2 \
    --chunksize 500000 \
    --batch-size 6 \
    --no-prepare

# -----------------------------------------------------------------------------
# Summary
# -----------------------------------------------------------------------------
echo ""
echo "============================================================"
echo "[$(date)] Pipeline COMPLETE for ${SAMPLE}"

WORKDIR=${PROJECT_ROOT}/5_AmbientDetection/${SAMPLE}
if [[ -d "${WORKDIR}/cell_map_ref_chunks" ]]; then
    N_CHUNKS=$(find "${WORKDIR}/cell_map_ref_chunks" -maxdepth 1 -name '*_cell_map_ref_chunk_*.txt' | wc -l)
    N_FILTERED=$(find "${WORKDIR}/cell_map_ref_chunks" -maxdepth 1 -name '*_filtered.tsv.gz' | wc -l)
    echo "  chunks dir: ${N_CHUNKS} chunk .txt files, ${N_FILTERED} filtered files"
fi

echo ""
echo "Next: configure genotyping based on synthetic disc analysis,"
echo "      then run genotyping as a separate step."
echo "============================================================"
