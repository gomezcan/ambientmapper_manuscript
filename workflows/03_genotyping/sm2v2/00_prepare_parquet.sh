#!/usr/bin/env bash
#SBATCH --job-name=prepare_pq
#SBATCH --partition=standard
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=8
#SBATCH --mem=64G
#SBATCH --time=12:00:00
#SBATCH --array=0-3
#SBATCH --output=_logs/00_prepare_pq_%A_%a.log
#
# 00_prepare_parquet.sh
#
# Convert filtered QCMapping TSVs to BC-sorted Parquet files.
# One-time preparation that speeds up the assign/score step ~100-700x
# by enabling DuckDB row-group predicate pushdown.
#
# DuckDB memory is capped at 16 GB (SET memory_limit) with external
# sort for files larger than that. 64 GB SLURM allocation gives
# headroom for the OS + DuckDB overhead.
#
# Task array:
#   0: Root1_rep1      (26 genomes, 226 GB)   ~30-60 min
#   1: B73Mo17_rep1    (2 genomes, 175 GB)    ~20-30 min
#   2: B73Mo17_rep2    (2 genomes, 88 GB)     ~10-15 min
#   3: multiGenotypes  (7 genomes, 1.2 TB)    ~2-4 h
#
# Resume-safe: skips genomes whose .parquet is already up to date.
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

CONFIGS=(
    "configs/Root1_rep1.ambientmapper.json"
    "configs/B73Mo17_rep1.ambientmapper.json"
    "configs/B73Mo17_rep2.ambientmapper.json"
    "configs/multiGenotypes_rep1.ambientmapper.json"
)
NAMES=(
    "Root1_rep1"
    "B73Mo17_rep1"
    "B73Mo17_rep2"
    "multiGenotypes_rep1"
)

IDX=${SLURM_ARRAY_TASK_ID}
CONFIG=${CONFIGS[$IDX]}
SAMPLE=${NAMES[$IDX]}

# AMBIENTMAPPER_REPO: optional local checkout of the ambientmapper tool, used only to log its commit
AM_VERSION=$(pip show ambientmapper 2>/dev/null | awk '/^Version/ {print $2}')
AM_COMMIT=$(git -C "${AMBIENTMAPPER_REPO:-.}" \
    rev-parse --short HEAD 2>/dev/null || echo "unknown")

echo "============================================================"
echo "[$(date)] Prepare: TSV → Parquet"
echo "  task          = ${IDX} / ${SAMPLE}"
echo "  config        = ${CONFIG}"
echo "  ambientmapper : v${AM_VERSION} @ ${AM_COMMIT}"
echo "============================================================"

if [[ ! -f "${CONFIG}" ]]; then
    echo "ERROR: config not found: ${CONFIG}" >&2
    exit 1
fi

ambientmapper prepare \
    --config "${CONFIG}" \
    --duckdb-threads 8

echo ""
echo "============================================================"
echo "[$(date)] Prepare COMPLETE for ${SAMPLE}"

# Report Parquet file sizes
WORKDIR=$(python3 -c "import json; c=json.load(open('${CONFIG}')); print(c['workdir'])")
FDIR="${WORKDIR}/${SAMPLE}/filtered_QCFiles"
echo "  Parquet files:"
ls -lh "${FDIR}"/filtered_*_QCMapping.parquet 2>/dev/null || echo "  (none)"
echo ""
echo "  Total:"
du -sh "${FDIR}"/filtered_*_QCMapping.parquet 2>/dev/null | tail -1 || echo "  (none)"
echo "============================================================"
