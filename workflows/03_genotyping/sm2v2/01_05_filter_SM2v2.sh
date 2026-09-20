#!/usr/bin/env bash
#SBATCH --job-name=sm2v2_filter
#SBATCH --partition=standard
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=16
#SBATCH --mem=40G
#SBATCH --time=2:00:00
#SBATCH --output=_logs/01_05_sm2v2_filter_%A.log
#
# 01_05_filter_SM2v2.sh — SM2v2 step 2/7: filter QCMapping (--min-barcode-freq).
#
# Reads SM2v2/qc/{B73,At}_QCMapping.txt and produces:
#   SM2v2/filtered_QCFiles/filtered_B73_QCMapping.txt
#   SM2v2/filtered_QCFiles/filtered_At_QCMapping.txt
#
# min_barcode_freq is read from the JSON config (= 3 for SM2v2).
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
THREADS=16

# AMBIENTMAPPER_REPO: optional local checkout of the ambientmapper tool, used only to log its commit
AM_VERSION=$(pip show ambientmapper 2>/dev/null | awk '/^Version/ {print $2}')
AM_COMMIT=$(git -C "${AMBIENTMAPPER_REPO:-.}" \
    rev-parse --short HEAD 2>/dev/null || echo "unknown")

echo "============================================================"
echo "[$(date)] SM2v2 step 2/7: filter"
echo "  config        = ${CONFIG}"
echo "  ambientmapper : v${AM_VERSION} @ ${AM_COMMIT}"
echo "============================================================"

[[ -f "${CONFIG}" ]] || { echo "ERROR: config not found: ${CONFIG}" >&2; exit 1; }

ambientmapper filter \
    --config "${CONFIG}" \
    --threads "${THREADS}"

echo ""
echo "[$(date)] filter DONE"
ls -lh SM2v2/filtered_QCFiles/ 2>/dev/null || true
