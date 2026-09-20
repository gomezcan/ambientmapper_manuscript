#!/usr/bin/env bash
#SBATCH --job-name=sm2v2_extract
#SBATCH --partition=standard
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=16
#SBATCH --mem=40G
#SBATCH --time=4:00:00
#SBATCH --output=_logs/01_04_sm2v2_extract_%A.log
#
# 01_04_extract_SM2v2.sh — SM2v2 step 1/7: extract QCMapping from input BAMs.
#
# Reads the two BAMs in configs/SM2v2.ambientmapper.json and produces:
#   SM2v2/qc/B73_QCMapping.txt
#   SM2v2/qc/At_QCMapping.txt
#
# Single task. Resume-safe (skips if outputs already exist).
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
echo "[$(date)] SM2v2 step 1/7: extract"
echo "  config        = ${CONFIG}"
echo "  ambientmapper : v${AM_VERSION} @ ${AM_COMMIT}"
echo "============================================================"

[[ -f "${CONFIG}" ]] || { echo "ERROR: config not found: ${CONFIG}" >&2; exit 1; }

ambientmapper extract \
    --config "${CONFIG}" \
    --threads "${THREADS}"

echo ""
echo "[$(date)] extract DONE"
ls -lh SM2v2/qc/ 2>/dev/null || true
