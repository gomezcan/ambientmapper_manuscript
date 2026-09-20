#!/usr/bin/env bash
#SBATCH --job-name=sm2v2_prepare
#SBATCH --partition=standard
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=8
#SBATCH --mem=64G
#SBATCH --time=2:00:00
#SBATCH --output=_logs/01_07_sm2v2_prepare_%A.log
#
# 01_07_prepare_SM2v2.sh — SM2v2 step 4/7: TSV → Parquet conversion.
#
# Converts SM2v2/filtered_QCFiles/filtered_*_QCMapping.txt → .parquet
# (BC-sorted with row groups). The assign step's DuckDB scoring path
# auto-detects and prefers Parquet (~100-700x faster than TSV scan).
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

echo "============================================================"
echo "[$(date)] SM2v2 step 4/7: prepare (Parquet)"
echo "  config        = ${CONFIG}"
echo "============================================================"

[[ -f "${CONFIG}" ]] || { echo "ERROR: config not found: ${CONFIG}" >&2; exit 1; }

ambientmapper prepare \
    --config "${CONFIG}" \
    --duckdb-threads 8

echo ""
echo "[$(date)] prepare DONE"
ls -lh SM2v2/filtered_QCFiles/*.parquet 2>/dev/null || true
