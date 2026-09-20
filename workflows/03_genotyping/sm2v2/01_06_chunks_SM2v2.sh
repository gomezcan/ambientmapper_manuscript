#!/usr/bin/env bash
#SBATCH --job-name=sm2v2_chunks
#SBATCH --partition=standard
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=4
#SBATCH --mem=16G
#SBATCH --time=1:00:00
#SBATCH --output=_logs/01_06_sm2v2_chunks_%A.log
#
# 01_06_chunks_SM2v2.sh — SM2v2 step 3/7: barcode chunking.
#
# Splits the union barcode set into chunks of size chunk_size_cells (= 50)
# and writes:
#   SM2v2/cell_map_ref_chunks/SM2v2_cell_map_ref_chunk_<N>.txt
#
# These barcode lists drive the assign step (one chunk per task slice).
# =============================================================================

set -euo pipefail

PROJECT_ROOT="${PROJECT_ROOT:-$(pwd)}"   # root of the data tree (3_Mapping/, 5_AmbientDetection/, ...)

cd "${PROJECT_ROOT}/5_AmbientDetection"

# conda activate <env from environment.yml>

CONFIG="configs/SM2v2.ambientmapper.json"

echo "============================================================"
echo "[$(date)] SM2v2 step 3/7: chunks"
echo "  config        = ${CONFIG}"
echo "============================================================"

[[ -f "${CONFIG}" ]] || { echo "ERROR: config not found: ${CONFIG}" >&2; exit 1; }

ambientmapper chunks --config "${CONFIG}"

echo ""
echo "[$(date)] chunks DONE"
N_CHUNK=$(find SM2v2/cell_map_ref_chunks -maxdepth 1 -name '*_cell_map_ref_chunk_*.txt' 2>/dev/null | wc -l)
echo "  ${N_CHUNK} chunk .txt files in SM2v2/cell_map_ref_chunks/"
