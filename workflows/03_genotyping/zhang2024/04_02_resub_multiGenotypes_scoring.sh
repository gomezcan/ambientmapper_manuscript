#!/usr/bin/env bash
#SBATCH --job-name=multi_rescore
#SBATCH --partition=standard
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=16
#SBATCH --mem=128G
#SBATCH --time=7-00:00:00
#SBATCH --output=_logs/04_02_multi_rescore_%j.log
#
# 04_02_resub_multiGenotypes_scoring.sh
#
# Resubmit multiGenotypes_rep1 assign/scoring with tuned Parquet params.
# Original 04_01 task 2 used score_batch_size=200 / score_workers=4
# (new defaults), giving ~7,700s per batch of 200 chunks = ~32 days total.
#
# Fix: batch_size=50 (better DuckDB row-group skipping on 200 GB Parquet)
#      + score_workers=6 (uses more of the 16 CPUs; 6*17.5 GB = 105 GB < 128 GB)
#
# Resume-safe: skips chunks with existing *_filtered.tsv.gz > 64 bytes.
#
# If resubmitting over a cancelled run, delete the partially written
# *_filtered.tsv.gz / *_raw.tsv.gz files of the cancelled batch first.
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

CONFIG="configs/multiGenotypes_rep1.ambientmapper.json"
SAMPLE="multiGenotypes_rep1"
THREADS=16

# AMBIENTMAPPER_REPO: optional local checkout of the ambientmapper tool, used only to log its commit
AM_VERSION=$(pip show ambientmapper 2>/dev/null | awk '/^Version/ {print $2}')
AM_COMMIT=$(git -C "${AMBIENTMAPPER_REPO:-.}" \
    rev-parse --short HEAD 2>/dev/null || echo "unknown")

echo "============================================================"
echo "[$(date)] Rescore: multiGenotypes_rep1 (tuned Parquet params)"
echo "  config        = ${CONFIG}"
echo "  ambientmapper : v${AM_VERSION} @ ${AM_COMMIT}"
echo "  score_batch_size = 50  (was 200)"
echo "  score_workers    = 6   (was 4)"
echo "============================================================"

# Count existing filtered files (resume baseline)
N_EXIST=$(find ${SAMPLE}/cell_map_ref_chunks -maxdepth 1 -name '*_filtered.tsv.gz' | wc -l)
N_TOTAL=$(find ${SAMPLE}/cell_map_ref_chunks -maxdepth 1 -name '*_cell_map_ref_chunk_*.txt' | wc -l)
echo "[$(date)] Resume: ${N_EXIST}/${N_TOTAL} chunks already done"

# Extract, filter, chunks: all done — skip directly to assign
ambientmapper assign \
    --config "${CONFIG}" \
    --threads "${THREADS}" \
    --alpha 0.05 \
    --k 10 \
    --mapq-min 10 \
    --xa-max 2 \
    --chunksize 500000 \
    --batch-size 6 \
    --no-prepare \
    --score-batch-size 50 \
    --score-workers 6

echo ""
echo "============================================================"
echo "[$(date)] Rescore COMPLETE for ${SAMPLE}"
N_FILTERED=$(find ${SAMPLE}/cell_map_ref_chunks -maxdepth 1 -name '*_filtered.tsv.gz' | wc -l)
echo "  chunks scored: ${N_FILTERED}/${N_TOTAL}"
echo "============================================================"
