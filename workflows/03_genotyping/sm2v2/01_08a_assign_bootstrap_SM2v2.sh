#!/usr/bin/env bash
#SBATCH --job-name=sm2v2_assign_boot
#SBATCH --partition=standard
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=16
#SBATCH --mem=64G
#SBATCH --time=4:00:00
#SBATCH --output=_logs/01_08a_sm2v2_assign_boot_%A.log
#
# 01_08a_assign_bootstrap_SM2v2.sh — SM2v2 step 5a/7: edges + ECDF bootstrap.
#
# Runs `ambientmapper assign` WITHOUT --score-chunk-range to learn:
#   SM2v2/ExplorationReadLevel/global_edges.npz
#   SM2v2/ExplorationReadLevel/global_ecdf.npz
#
# Required because --score-chunk-range implies --only-score, which implies
# --skip-edges (assign_streaming requires the edges/ECDF models to exist).
#
# In addition to learning the models, this single-task run will also start
# scoring chunks (resume-safe). Whatever it gets done before walltime is fine
# — the 16-task array (01_08) resume-scores the remainder.
#
# Once the edges + ECDF .npz files exist, this bootstrap is unnecessary and
# 01_08 can be re-submitted any number of times.
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

# AMBIENTMAPPER_REPO: optional local checkout of the ambientmapper tool, used only to log its commit
AM_VERSION=$(pip show ambientmapper 2>/dev/null | awk '/^Version/ {print $2}')
AM_COMMIT=$(git -C "${AMBIENTMAPPER_REPO:-.}" \
    rev-parse --short HEAD 2>/dev/null || echo "unknown")

TOTAL=$(find ${SAMPLE}/cell_map_ref_chunks -maxdepth 1 -name '*_cell_map_ref_chunk_*.txt' | wc -l)
N_EXIST=$(find ${SAMPLE}/cell_map_ref_chunks -maxdepth 1 -name '*_filtered.tsv.gz' | wc -l)

echo "============================================================"
echo "[$(date)] SM2v2 step 5a/7: assign bootstrap (edges + ECDFs)"
echo "  config        = ${CONFIG}"
echo "  ambientmapper : v${AM_VERSION} @ ${AM_COMMIT}"
echo "  total chunks  = ${TOTAL}"
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
    --no-prepare

echo ""
echo "============================================================"
echo "[$(date)] bootstrap DONE (or hit walltime — whichever)"
N_FILTERED=$(find ${SAMPLE}/cell_map_ref_chunks -maxdepth 1 -name '*_filtered.tsv.gz' | wc -l)
echo "  total filtered: ${N_FILTERED}/${TOTAL}"
echo "  edges npz: $(ls ${SAMPLE}/ExplorationReadLevel/global_edges.npz 2>/dev/null && echo OK || echo MISSING)"
echo "  ecdf  npz: $(ls ${SAMPLE}/ExplorationReadLevel/global_ecdf.npz  2>/dev/null && echo OK || echo MISSING)"
echo "============================================================"
echo "Next: sbatch workflows/03_genotyping/sm2v2/01_08_assign_array_SM2v2.sh  (16-task array, resumes)"
