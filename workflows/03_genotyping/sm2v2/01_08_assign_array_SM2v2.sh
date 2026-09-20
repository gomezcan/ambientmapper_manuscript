#!/usr/bin/env bash
#SBATCH --job-name=sm2v2_assign
#SBATCH --partition=standard
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=16
#SBATCH --mem=64G
#SBATCH --time=6:00:00
#SBATCH --array=0-15
#SBATCH --output=_logs/01_08_sm2v2_assign_%A_%a.log
#
# 01_08_assign_array_SM2v2.sh — SM2v2 step 5/7: assign + score (SLURM array).
#
# Splits the chunk list into N_TASKS disjoint slices via --score-chunk-range.
# Each task scores its slice independently and writes:
#   SM2v2/cell_map_ref_chunks/SM2v2_cell_map_ref_chunk_<N>_filtered.tsv.gz
#
# 2-genome dataset (B73 + At), ~37.8K chunks → ~2,362/task at 16 tasks.
# Resume-safe: ambientmapper skips chunks with existing *_filtered.tsv.gz > 64 bytes.
#
# alpha=0.05 / k=10 / mapq-min=10 / xa-max=2: matches sub1k C0 reconciled
# bundle. (Genotyping-side mq=20/xa=0 is a SEPARATE filter applied later.)
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

CONFIG="configs/SM2v2.ambientmapper.json"
SAMPLE="SM2v2"
THREADS=16
N_TASKS=16

# AMBIENTMAPPER_REPO: optional local checkout of the ambientmapper tool, used only to log its commit
AM_VERSION=$(pip show ambientmapper 2>/dev/null | awk '/^Version/ {print $2}')
AM_COMMIT=$(git -C "${AMBIENTMAPPER_REPO:-.}" \
    rev-parse --short HEAD 2>/dev/null || echo "unknown")

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
echo "[$(date)] Score array: SM2v2"
echo "  task          = ${SLURM_ARRAY_TASK_ID} / ${N_TASKS}"
echo "  config        = ${CONFIG}"
echo "  ambientmapper : v${AM_VERSION} @ ${AM_COMMIT}"
echo "  total chunks  = ${TOTAL}"
echo "  chunk range   = ${START}:${END} ($(( END - START )) chunks)"
echo "  resume base   = ${N_EXIST}/${TOTAL} already done"
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
    --score-chunk-range ${START}:${END}

echo ""
echo "============================================================"
echo "[$(date)] Score array task ${SLURM_ARRAY_TASK_ID} COMPLETE"
N_FILTERED=$(find ${SAMPLE}/cell_map_ref_chunks -maxdepth 1 -name '*_filtered.tsv.gz' | wc -l)
echo "  total filtered: ${N_FILTERED}/${TOTAL}"
echo "============================================================"
