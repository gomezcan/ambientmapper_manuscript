#!/usr/bin/env bash
#SBATCH --job-name=multi_score_arr
#SBATCH --partition=standard
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=16
#SBATCH --mem=128G
#SBATCH --time=2-00:00:00
#SBATCH --array=0-7
#SBATCH --output=_logs/04_03_multi_score_%A_%a.log
#
# 04_03_scoring_array_multiGenotypes.sh
#
# SLURM array scoring for multiGenotypes_rep1 using --score-chunk-range.
# 8 tasks each process ~9.1K chunks in parallel (73K total).
# Expected wall: ~50h per task at 20s/chunk.
#
# Resume-safe: ambientmapper skips chunks with existing *_filtered.tsv.gz > 64 bytes.
#
# Requires an ambientmapper build with --score-chunk-range. If resubmitting over
# a cancelled run, delete the partially written *_filtered.tsv.gz / *_raw.tsv.gz
# files of the cancelled batch first.
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
N_TASKS=8

# AMBIENTMAPPER_REPO: optional local checkout of the ambientmapper tool, used only to log its commit
AM_VERSION=$(pip show ambientmapper 2>/dev/null | awk '/^Version/ {print $2}')
AM_COMMIT=$(git -C "${AMBIENTMAPPER_REPO:-.}" \
    rev-parse --short HEAD 2>/dev/null || echo "unknown")

# Count total chunks (sorted glob order matches ambientmapper's sorted() call)
TOTAL=$(find ${SAMPLE}/cell_map_ref_chunks -maxdepth 1 -name '*_cell_map_ref_chunk_*.txt' | wc -l)
PER_TASK=$(( (TOTAL + N_TASKS - 1) / N_TASKS ))
START=$(( SLURM_ARRAY_TASK_ID * PER_TASK ))
END=$(( START + PER_TASK ))
if [ ${END} -gt ${TOTAL} ]; then END=${TOTAL}; fi
if [ ${START} -ge ${TOTAL} ]; then
    echo "[$(date)] Task ${SLURM_ARRAY_TASK_ID}: START=${START} >= TOTAL=${TOTAL}, nothing to do"
    exit 0
fi

N_EXIST=$(find ${SAMPLE}/cell_map_ref_chunks -maxdepth 1 -name '*_filtered.tsv.gz' | wc -l)

echo "============================================================"
echo "[$(date)] Score array: multiGenotypes_rep1"
echo "  task          = ${SLURM_ARRAY_TASK_ID} / ${N_TASKS}"
echo "  config        = ${CONFIG}"
echo "  ambientmapper : v${AM_VERSION} @ ${AM_COMMIT}"
echo "  total chunks  = ${TOTAL}"
echo "  chunk range   = ${START}:${END} ($(( END - START )) chunks)"
echo "  resume base   = ${N_EXIST}/${TOTAL} already done"
echo "  batch_size=50, workers=6"
echo "============================================================"

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
    --score-workers 6 \
    --score-chunk-range ${START}:${END}

echo ""
echo "============================================================"
echo "[$(date)] Score array task ${SLURM_ARRAY_TASK_ID} COMPLETE"
N_FILTERED=$(find ${SAMPLE}/cell_map_ref_chunks -maxdepth 1 -name '*_filtered.tsv.gz' | wc -l)
echo "  total filtered: ${N_FILTERED}/${TOTAL}"
echo "============================================================"
